"""Acceptance checks and body-frame diagnostics for one DUST BFF run."""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path
import sys

import h5py
import numpy as np
from scipy.io import netcdf_file
from scipy.spatial.transform import Rotation

SURFACES = {'BFL':1004,'WF1L':1008,'WF2L':1011,'WF3L':1014,'WF4L':1017,
            'BFR':2004,'WF1R':2008,'WF2R':2011,'WF3R':2014,'WF4R':2017}


def read_state(path):
    with netcdf_file(path/'case.nc', 'r', mmap=False) as d:
        def get(key):
            return np.asarray(d.variables[key].data, dtype=float).copy()
        result = {'time':get('time'), 'x':get('node.struct.990000.X'),
                  'v':get('node.struct.990000.XP'), 'phi':get('node.struct.990000.Phi'),
                  'omega':get('node.struct.990000.Omega'), 'a':get('elem.joint.5.a'),
                  'ap':get('elem.joint.5.aPrime'), 'yaw_constraint':get('elem.joint.1.Phi')[:,2],
                  'left':get('node.struct.990020.X'), 'right':get('node.struct.991020.X')}
        result['surfaces'] = {name:get(f'elem.joint.{label}.Phi')[:,1] for name,label in SURFACES.items()}
    return result


def release_check(state, c):
    t = state['time']; off = c['sas_off_start_s']
    mask = (t >= off-1) & (t < off)
    if mask.sum() < max(3, round(.9/c['dt_s'])):
        return {'passed':False, 'reason':'less than 0.9 s of pre-release history'}
    v, w = state['v'][mask], state['omega'][mask]
    checks = {'vz_inps':abs(float(v[-1,2])), 'pitch_rate_radps':abs(float(w[-1,1])),
              'roll_rate_radps':abs(float(w[-1,0])), 'std_vz_inps':float(v[:,2].std()),
              'std_pitch_rate_radps':float(w[:,1].std())}
    failed = [k for k,v in checks.items() if not math.isfinite(v) or v > c['release_limits'][k]]
    return {'passed':not failed, 'values':checks, 'limits':c['release_limits'], 'failed':failed}


def analyse(path):
    c = json.loads((path/'case.json').read_text())
    s = read_state(path); t=s['time']
    finite = all(np.isfinite(a).all() for k,a in s.items() if k != 'surfaces') and all(np.isfinite(a).all() for a in s['surfaces'].values())
    complete = bool(len(t)>1 and t[-1] >= c['final_time_s']-c['dt_s']*1.e-5)
    count_ok = s['a'].shape[1] == c['elastic_modes']
    yaw = float(np.max(np.abs(s['yaw_constraint'])))
    inv = Rotation.from_rotvec(s['phi']).inv()
    left = inv.apply(s['left']-s['x']); right = inv.apply(s['right']-s['x'])
    symmetric = .5*(left[:,2]+right[:,2]); symmetric -= symmetric[0]
    pitch = Rotation.from_rotvec(s['phi']).as_euler('xyz')[:,1]
    air_body = inv.apply(np.array([c['velocity_mps']/.0254, 0., 0.])-s['v'])
    alpha = np.arctan2(air_body[:,2], air_body[:,0])
    bff_end=c['sas_off_start_s']+c['rap_delay_s']+c['rap_duration_s']
    id_start=bff_end+max(.05,8*c['actuator_tau_s'])
    id_end=c['sas_on_start_s']-.05
    use=(t >= id_start)&(t <= id_end)
    # Smoke hold check excludes WF4 RAP and checks the other eight held surfaces.
    held=(t >= c['sas_off_start_s']+2*c['dt_s'])&(t < c['sas_on_start_s']-c['dt_s'])
    hold_changes = {name:float(np.ptp(value[held])) for name,value in s['surfaces'].items()
                    if held.sum() and (not c['smoke'] or not name.startswith('WF4'))}
    if not c['smoke']:
        hold_changes = {name:float(np.ptp(value[use])) for name,value in s['surfaces'].items()} if use.sum() else {}
    hold_ok = bool(hold_changes) and max(hold_changes.values()) < math.radians(.001)
    with (path/'coupled_response.csv').open() as f:
        coupled = list(csv.DictReader(f))
    loads_finite = bool(coupled) and all(math.isfinite(float(row[key])) for row in coupled for key in row if key not in ('state',))
    max_iterations = max((int(row['coupling_iterations']) for row in coupled), default=0)
    h5_finite=True; hinge_max={}; frame_count=0
    for frame in sorted((path/'dust').glob('case_res_*.h5')):
        frame_count+=1
        with h5py.File(frame) as h:
            for key in ('Geometry/rr','Geometry/rr_virtual','Solution/dF','Solution/dMom'):
                h5_finite &= bool(np.isfinite(h['Components/Comp001/'+key][:]).all())
            for i in range(1,11):
                a=h[f'Components/Comp001/Hinges/Hinge_{i:02d}/theta'][:]
                h5_finite &= bool(np.isfinite(a).all())
                hinge_max[str(i)]=max(hinge_max.get(str(i),0.),float(np.max(np.abs(a))))
    release = {'passed':False,'reason':'compressed technical smoke; no equilibrium validation'} if c['smoke'] else release_check(s,c)
    technical = complete and finite and count_ok and loads_finite and h5_finite and yaw < 1.e-5 and hold_ok and 0<max_iterations<20 and frame_count>0
    report={'technical_pass':bool(technical),'complete':complete,'final_time_s':float(t[-1]),
            'states':len(t), 'modes':s['a'].shape[1], 'finite':bool(finite and loads_finite and h5_finite),
            'yaw_constraint_max_rad':yaw,'max_implicit_iterations':max_iterations,
            'hold_changes_rad':hold_changes, 'surface_hold_pass':hold_ok,'release':release,
            'dust_frames':frame_count,'dust_hinge_max_abs_theta':hinge_max,
            'dust_theta_note':'theta is unused for coupled hinges; zero is not a motion failure (mod_hinges.f90).',
            'flutter_validated':False,'load_timing_validated':False,
            'identification_valid':False,'classification':'technical smoke' if c['smoke'] else 'single-speed BFF candidate'}
    columns = {'time_s':t,'pitch_rad':pitch,'alpha_rad':alpha,'pitch_rate_radps':s['omega'][:,1],
               'vz_inps':s['v'][:,2], 'symmetric_tip_body_in':symmetric,
               **{f'q_{mode}':s['a'][:,i] for i,mode in enumerate(c['mode_ids'])},
               **{f'{name}_rad':value for name,value in s['surfaces'].items()}}
    with (path/'response.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(columns);writer.writerows(zip(*columns.values()))
    if technical and not c['smoke'] and release['passed'] and use.sum()>100:
        sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'bff_open_loop'))
        from modal_identification import identify_multimodal, WINDOWS_S
        sig={'symmetric_tip':symmetric[use], 'swb1':s['a'][use,0], 'swb1_velocity':s['ap'][use,0],
             'q':s['omega'][use,1], 'alpha':alpha[use], 'pitch':pitch[use]}
        _,clusters,candidate,sp=identify_multimodal(t[use],sig,{w:w<=t[use][-1]-t[use][0] for w in WINDOWS_S})
        report.update(candidate=candidate,short_period=sp,clusters=clusters,
                      identification_valid=bool(candidate.get('accepted',False)))
    (path/'analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
