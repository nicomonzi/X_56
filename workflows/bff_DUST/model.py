"""Render a single X-56 DUST BFF experiment from the repository's proven decks."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import re
import shutil

import h5py
import numpy as np
from scipy.signal import butter

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
OPEN = REPO / 'workflows/bff_open_loop'
COUPLED = REPO / 'workflows/coupled_dust'
FEM = REPO / 'assets/fem/mbdyn_modal_60.fem'
PID_NAMES = ('ALT_PID', 'VZ_PID', 'PITCH_PID', 'Q_PID', 'ROLL_PID', 'P_PID', 'VY_PID')
TAGS = ('BFR', 'WF1R', 'WF2R', 'WF3R', 'WF4R', 'BFL', 'WF1L', 'WF2L', 'WF3L', 'WF4L')
PARENTS = (881004, 881008, 881011, 881014, 881017, 880004, 880008, 880011, 880014, 880017)


def replace(pattern, replacement, text, expected=1):
    result, count = re.subn(pattern, lambda m: replacement(m) if callable(replacement) else replacement,
                            text, flags=re.M)
    if count != expected:
        raise ValueError(f'Expected {expected} matches, got {count}: {pattern}')
    return result


def constant(text, name, value):
    return replace(rf'^set:\s*const\s+real\s+{name}\s*=\s*[^;]+;',
                   f'set: const real {name} = {value};', text)


def without_comments(text):
    return re.sub(r'(?m)#.*$', '', text)


def spectrum():
    text = FEM.read_text()
    def matrix(group):
        body = re.search(rf'\*\* RECORD GROUP {group},[^\n]*\n(.*?)(?=\*\* RECORD GROUP|\Z)', text, re.S).group(1)
        values = ' '.join(l for l in body.splitlines() if not l.lstrip().startswith('**'))
        a = np.fromstring(values.replace('D', 'E'), sep=' ')
        if a.size != 3600:
            raise ValueError(f'FEM60 group {group}: expected 3600 matrix entries')
        return a.reshape(60, 60, order='F')
    mass, stiffness = matrix(9), matrix(10)
    freq = np.sqrt(np.maximum(0, np.diag(stiffness)/np.diag(mass)))/(2*math.pi)
    if np.max(freq[:6]) > .01 or np.any(freq[6:] <= .01):
        raise ValueError('Unexpected rigid/elastic mode numbering')
    return freq, hashlib.sha256(text.encode()).hexdigest()


def load_config(path):
    c = json.loads(path.read_text())
    positive = ('velocity_mps', 'dt_s', 'final_time_s', 'actuator_tau_s', 'output_dt_s',
                'rho_lbf_s2_in4', 'sas_off_duration_s', 'rap_frequency_hz', 'hinge_endpoint_inset_in',
                'wake_length_chords', 'wake_reference_chord_in')
    for name in positive:
        if not math.isfinite(c[name]) or c[name] <= 0:
            raise ValueError(f'{name} must be positive and finite')
    if not isinstance(c['elastic_modes'], int) or not 1 <= c['elastic_modes'] <= 25:
        raise ValueError('Choose 1..25 elastic modes')
    if c['nelem_chord'] != 30:
        raise ValueError('This workflow uses the selected chord-30 mesh')
    for name in ('n_wake_panels', 'n_wake_particles'):
        if type(c[name]) is not int or c[name] < 1:
            raise ValueError(f'{name} must be a positive integer')
    if not math.isfinite(c['wake_reference_te_x_in']) or c['wake_reference_te_x_in'] < 221:
        raise ValueError('Wake reference must lie at or downstream of the selected mesh trailing edge (x~220.725 in)')
    if set(c['controller_gain_scale']) != set(PID_NAMES):
        raise ValueError('All seven non-yaw controllers must be configured')
    for v in c['controller_gain_scale'].values():
        if not math.isfinite(v) or v <= 0:
            raise ValueError('Controller gain scales must be finite and positive')
    for name in ('final_time_s', 'output_dt_s', 'sas_off_start_s', 'sas_off_duration_s'):
        if not math.isclose(c[name]/c['dt_s'], round(c[name]/c['dt_s']), abs_tol=1.e-8):
            raise ValueError(f'{name} must be a multiple of dt_s')
    off, end = c['sas_off_start_s'], c['sas_off_start_s']+c['sas_off_duration_s']
    if not 1 < off < off+c['rap_delay_s']+.742/c['rap_frequency_hz']+.1 < end < c['final_time_s']:
        raise ValueError('Invalid settling / RAP / observation / recovery sequence')
    return c


def build(output, c, smoke=False):
    output.mkdir(parents=True, exist_ok=False)
    inc = output/'INCLUDE'; inc.mkdir()
    dustdir = output/'dust'; dustdir.mkdir()
    c = dict(c)
    if smoke:
        # Exercise all phases in a short technical test, not a flutter measurement.
        c.update(final_time_s=.04, sas_off_start_s=.008, sas_off_duration_s=.024,
                 rap_delay_s=.002, rap_frequency_hz=.742/.008, output_dt_s=.004)
    freq, digest = spectrum()
    modes = list(range(7, 7+c['elastic_modes']))
    samples = 1/(c['dt_s']*freq[modes[-1]-1])
    if samples < 12:
        raise ValueError(f'Only {samples:.2f} samples/period at highest mode; reduce dt')
    c['mode_ids'] = modes
    c['mode_frequencies_hz'] = [float(freq[i-1]) for i in modes]
    c['smoke'] = smoke
    c['rap_duration_s'] = .742/c['rap_frequency_hz']
    c['sas_on_start_s'] = c['sas_off_start_s']+c['sas_off_duration_s']
    (output/'case.json').write_text(json.dumps(c, indent=2)+'\n')

    # Snapshot just the required non-aerodynamic structural includes.
    for name in ('node.mbd', 'reference.mbd', 'driver.mbd', 'control_surfaces.mbd'):
        shutil.copy2(OPEN/'INCLUDE'/name, inc/name)
    # MBDyn's external structural force accepts dynamic interface nodes.
    node_text = (inc/'node.mbd').read_text()
    node_text = re.sub(r'(structural:\s*\d+,)\s*static,', r'\1 dynamic,', node_text)
    (inc/'node.mbd').write_text(node_text)
    (inc/'mbdyn_modal_60.fem').symlink_to(FEM)
    const = without_comments((OPEN/'INCLUDE/setconst.mbd').read_text())
    const = re.sub(r'^set: const (?:real|integer) NASTRAN_DLM_ROM_\w+\s*=[^;]+;\s*', '', const, flags=re.M)
    for name in ('YAW_PID', 'YAW_DRIVE', 'R_PID', 'R_DRIVE'):
        const = replace(rf'^set: const integer {name} =[^;]+;', '', const)
    filters = {}
    for prefix, cutoff in c['filter_cutoffs_hz'].items():
        order = 2 if prefix == 'LP' else 1
        b, a = butter(order, cutoff, fs=1/c['dt_s'])
        for i in range(1, order+1):
            const = constant(const, f'{prefix}_A{i}', f'{-a[i]:.17g}')
        for i in range(order+1):
            const = constant(const, f'{prefix}_B{i}', f'{b[i]:.17g}')
        filters[prefix] = {'cutoff_hz': cutoff, 'b': b.tolist(), 'a': a.tolist()}
    dt, tau = c['dt_s'], c['actuator_tau_s']
    aa, bb = (2*tau-dt)/(2*tau+dt), dt/(2*tau+dt)
    for name, value in {'ACT_A1': aa, 'ACT_B0': bb, 'ACT_B1': bb,
                        'SAS_MODAL_DAMPING': c['sas_modal_damping']}.items():
        const = constant(const, name, f'{value:.17g}')
    filters['ACT'] = {'tau_s': tau, 'b': [bb, bb], 'a': [1, -aa]}
    (inc/'setconst.mbd').write_text('# Filter coefficients regenerated for dt; gains are DUST tuning seeds.\n'+const)

    modal = (OPEN/'INCLUDE/modaljoint.mbd').read_text()
    tail = modal[modal.index('        origin position,'):]
    init = ',\n'.join(f'            mode, {i}, '+(f'Q_INIT_{i}' if i <= 12 else '0.')+', 0.' for i in modes)
    (inc/'modaljoint.mbd').write_text(
        'joint: MODAL_JOINT, modal, BASE_NODE,\n'
        f'    {len(modes)}, list, '+', '.join(map(str, modes))+',\n'
        '    initial value,\n'+init+',\n'
        f'    8527, single factor damping, {c["structural_damping_ratio"]},\n'
        f'    "{inc / "mbdyn_modal_60.fem"}",\n'+tail)

    main = without_comments((OPEN/'main_bff_open_loop.mbd').read_text())
    for pattern, value in [(r'aerodynamic elements: 58;', ''), (r'structural nodes: 56;', 'structural nodes: 76;'),
                           (r'loadable elements: 9;', 'loadable elements: 7;'), (r'joints: 12;', 'joints: 32;'),
                           (r'^\s*air properties;', ''), (r'^\s*default aerodynamic output:[^;]+;', ''),
                           (r'orientation constraint, inactive, inactive, inactive, null;',
                            'orientation constraint, inactive, inactive, active, null;'),
                           (r'output results: netcdf, no text;', 'output results: netcdf, sync, no text;')]:
        main = replace(pattern, value, main)
    for pid, drive in [('YAW_PID', 'YAW_DRIVE'), ('R_PID', 'R_DRIVE')]:
        main = replace(rf'\s*user defined: {pid}, pid,[^;]+;', '', main)
        main = replace(rf'\s*drive caller: {drive},[^;]+;', '', main)
    main = replace(r'drive caller: DIR_RAW_DRIVE, array, 3,[^;]+;',
                   'drive caller: DIR_RAW_DRIVE, reference, VY_DRIVE;', main)
    main = replace(r'\s*force: NASTRAN_DLM_ROM_FORCE,[^;]+;', '', main)
    main = replace(r'\s*air properties: RHO_AIR,[^;]+;', '', main)
    main = replace(r'include: "./INCLUDE/aerobody.mbd";',
                   'include: "./INCLUDE/coupling_hinge_joints.mbd";\n'
                   '    include: "./INCLUDE/external_dust_force.mbd";', main)
    main = main.replace('end: nodes;', '    include: "./INCLUDE/coupling_hinge_nodes.mbd";\nend: nodes;')
    for pid, scale in c['controller_gain_scale'].items():
        main = replace(rf'user defined: {pid}, pid,[^;]+;',
                       lambda m: re.sub(r'(K[pid],)\s*([^,]+),',
                                        lambda g: f'{g[1]} ({g[2]})*{scale:.12g},', m[0]), main)
    for name, value in {'V_INF': c['velocity_mps'], 'RHO_AIR': c['rho_lbf_s2_in4'],
                        'TIME_STEP': dt, 'FINAL_TIME': c['final_time_s'],
                        'SAS_OFF_START': c['sas_off_start_s'], 'SAS_OFF_DURATION': c['sas_off_duration_s'],
                        'BFF_RAP_TARGET_FREQUENCY': c['rap_frequency_hz'],
                        'BFF_RAP_START': f'SAS_OFF_START + {c["rap_delay_s"]}',
                        'BFF_RAP_AMPLITUDE': f'{c["rap_amplitude_deg"]}*deg2rad'}.items():
        main = constant(main, name, value)
    # The filtered RAP has a tail. Allow >= 8 actuator time constants before identification.
    main = constant(main, 'IDENTIFICATION_START', f'BFF_RAP_END + {max(.05, 8*tau)}')
    main = main.replace('./INCLUDE/', str(inc)+'/')
    (output/'main.mbd').write_text('# BFF DUST: IPS, yaw constrained, seven non-yaw PIDs, no MBDyn aero.\n'+main)

    nodes = np.loadtxt(COUPLED/'model/dust/coupling_nodes.in')
    # Coincident endpoints on adjacent independently moving flaps cause NN ties.
    # Move only duplicate auxiliary sampling points inside their own hinge line.
    # These are force/kinematic probes, not aerodynamic surface vertices.
    original = nodes.copy()
    duplicate = np.sum(np.all(original[:,None,:] == original[None,:,:], axis=2), axis=1)>1
    for i in np.flatnonzero(duplicate):
        j = i+1 if (i-39)%2 == 0 else i-1
        direction = original[j]-original[i]
        nodes[i] += c['hinge_endpoint_inset_in']*direction/np.linalg.norm(direction)
    if len(np.unique(nodes, axis=0)) != 59:
        raise ValueError('Ambiguous coupling coordinates remain')
    np.savetxt(output/'coupling_nodes.in', nodes, fmt='%.12f')
    hn, hj = [], []
    for i, xyz in enumerate(nodes[39:]):
        label, parent = 882001+i, PARENTS[i//2]
        pos = ', '.join(f'{x:.12f}' for x in xyz)
        hn.append(f'structural: {label}, dynamic, reference, GROUND, {pos},\n'
                  '    reference, GROUND, eye, reference, GROUND, null, reference, GROUND, null;\n')
        hj.append(f'joint: {9601+i}, total joint,\n'
                  f'    {parent}, position, reference, GROUND, {pos},\n'
                  '        position orientation, reference, node, eye, rotation orientation, reference, node, eye,\n'
                  f'    {label}, position, reference, GROUND, {pos},\n'
                  '        position orientation, reference, other node, eye, rotation orientation, reference, other node, eye,\n'
                  '    position constraint, 1, 1, 1, null, orientation constraint, 1, 1, 1, null;\n')
    (inc/'coupling_hinge_nodes.mbd').write_text('\n'.join(hn))
    (inc/'coupling_hinge_joints.mbd').write_text('\n'.join(hj))
    shutil.copy2(COUPLED/'model/mbdyn/INCLUDE/external_dust_force.mbd', inc)
    mesh = (COUPLED/'meshes/FINE/parametric_mesh.in').read_text()
    mesh = mesh.replace('model/dust/coupling_nodes.in', str(output/'coupling_nodes.in'))
    mesh = mesh.replace('model/dust/airfoilsection/', str(COUPLED/'model/dust/airfoilsection')+'/')
    (output/'mesh.in').write_text(mesh)
    shutil.copy2(COUPLED/'model/dust/References.in', output)
    (output/'dust_pre.in').write_text('comp_name = X56\ngeo_file = mesh.in\nref_tag = centerbody\nfile_name = geo_input.h5\n')
    dust = (COUPLED/'model/dust/dust_production.in').read_text()
    speed = c['velocity_mps']/.0254
    for key, value in {'basename': 'dust/case', 'geometry_file': 'geo_input.h5', 'reference_file': 'References.in',
                       'precice_config': str(output/'precice.xml'), 'tend': c['final_time_s'], 'dt': dt,
                       'dt_out': c['output_dt_s'], 'rho_inf': c['rho_lbf_s2_in4'],
                       'u_inf': f'(/{speed:.12f},0.,0./)'}.items():
        dust = replace(rf'^{key}\s*=.*$', f'{key} = {value}', dust)
    # Fixed global downstream deletion plane: 15 reference chords behind the
    # most downstream TE (~220.725 in, rounded to 221 in). Not a particle age.
    downstream = c['wake_reference_te_x_in'] + c['wake_length_chords']*c['wake_reference_chord_in']
    dust = replace(r'^particles_box_max\s*=.*$', f'particles_box_max = (/{downstream},1500.,1500./)', dust)
    # FMM boxes enclose the deletion plane, but must not round that plane outward.
    dust = replace(r'^n_box\s*=.*$', f'n_box = (/{math.ceil((downstream+1000)/500)},6,6/)', dust)
    dust = replace(r'^n_wake_panels\s*=.*$', f'n_wake_panels = {c["n_wake_panels"]}', dust)
    # Capacity is separate from spatial deletion; do not allocate a full-run wake.
    dust = replace(r'^n_wake_particles\s*=.*$', f'n_wake_particles = {c["n_wake_particles"]}', dust)
    dust = re.sub(r'(?m)^! 9.5 s.*\n', '', dust)
    (output/'dust.in').write_text(dust)
    xml = (COUPLED/'model/precice-config-v3.xml').read_text()
    xml = xml.replace('max-time value="9.50"', f'max-time value="{c["final_time_s"]}"')
    xml = xml.replace('exchange-directory="."', f'exchange-directory="{output / "exchange"}"')
    (output/'precice.xml').write_text(xml)
    # Local snapshots prevent later changes to the old workflow altering a run.
    for name in ('coupling_adapter.py', 'mbdyn_interface.py'):
        shutil.copy2(COUPLED/'adapters'/name, output/name)
    invariants = {'controllers': re.findall(r'user defined: (\w+), pid,', main),
                  'no_mbdyn_aero': not re.search(r'aerodynamic elements:|air properties:|NASTRAN_DLM_ROM|aerobody.mbd', main),
                  'yaw_locked': 'orientation constraint, inactive, inactive, active, null;' in main,
                  'hold_surfaces': (inc/'control_surfaces.mbd').read_text().count('sample and hold'),
                  'unique_coupling_vertices': int(len(np.unique(nodes, axis=0)))}
    if set(invariants['controllers']) != set(PID_NAMES) or not invariants['no_mbdyn_aero'] or not invariants['yaw_locked'] or invariants['hold_surfaces'] != 10:
        raise ValueError(f'Invalid model: {invariants}')
    report = {'config': c, 'invariants': invariants, 'filters': filters, 'fem_sha256': digest,
              'samples_per_highest_period': samples, 'coupling_probe_adjustments_in': (nodes-original).tolist(),
              'source_sha256': {str(p.relative_to(REPO)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in [OPEN/'main_bff_open_loop.mbd', OPEN/'INCLUDE/setconst.mbd', COUPLED/'meshes/FINE/parametric_mesh.in']}}
    (output/'model_audit.json').write_text(json.dumps(report, indent=2)+'\n')
    return c


def check_geometry(output):
    with h5py.File(output/'geo_input.h5') as h:
        g = h['Components/Comp001']
        shape = [len(g['Geometry/ee']), len(g['Geometry/rr']), len(g['Geometry/CouplingNodes']), int(g['Hinges/n_hinges'][()])]
        np.testing.assert_allclose(g['Geometry/CouplingNodes'][:], np.loadtxt(output/'coupling_nodes.in'), atol=1.e-10, rtol=0)
    if shape != [7320, 7503, 59, 10]:
        raise ValueError(f'Unexpected mesh: {shape}')
    return shape
