#!/usr/bin/env python3
"""One-velocity BFF DUST: prepare/check/smoke/run/analyse, with no sweep."""
from __future__ import annotations
import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time

import numpy as np
from model import ROOT, COUPLED, build, check_geometry, load_config
from analyse import analyse, read_state, release_check


def binaries():
    defaults={'MBDYN_BIN':'/usr/local/mbdyn/bin/mbdyn',
              'DUST_BIN':str(Path.home()/'dust-patched/build-user/bin/dust'),
              'DUST_PRE_BIN':str(Path.home()/'dust-patched/build-user/bin/dust_pre'),
              'DUST_POST_BIN':str(Path.home()/'dust-patched/build-user/bin/dust_post'),
              'MBDYN_PYTHON_PATH':'/usr/local/mbdyn/libexec/mbpy'}
    local=ROOT/'machine.env'
    if local.exists():
        for line in local.read_text().splitlines():
            if line.strip() and not line.lstrip().startswith('#'):
                key,value=line.split('=',1); defaults[key.strip()]=value.strip()
    defaults={key:os.environ.get(key,value) for key,value in defaults.items()}
    for key,value in defaults.items():
        if key.endswith('_BIN') and not os.access(value,os.X_OK):
            raise RuntimeError(f'Executable unavailable: {key}={value}')
    for directory in defaults['MBDYN_PYTHON_PATH'].split(os.pathsep):
        sys.path.insert(0,directory)
    return defaults


def stop(p):
    if p is not None and p.poll() is None:
        os.killpg(p.pid,signal.SIGTERM)
        try:p.wait(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(p.pid,signal.SIGKILL);p.wait()


def execute(cmd,path,log,env,timeout=120):
    with (path/log).open('w') as f:
        r=subprocess.run(cmd,cwd=path,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=timeout)
    text=(path/log).read_text(errors='replace')
    if r.returncode or 'ERROR in' in text:
        raise RuntimeError(f'Command failed: see {path/log}')


def worker(path,threads):
    """Run in a separate process, giving the parent a timeout for initialization."""
    bins=binaries();c=json.loads((path/'case.json').read_text())
    env=os.environ.copy()
    env.update(OMP_NUM_THREADS=str(threads),OPENBLAS_NUM_THREADS='1',OMP_PLACES='cores',OMP_PROC_BIND='close')
    env.pop('DISPLAY',None);env.pop('WAYLAND_DISPLAY',None)
    sys.path.insert(0,str(path))
    from mbdyn_interface import MBDynInterface
    from coupling_adapter import CouplingAdapter
    class CheckedAdapter(CouplingAdapter):
        def _write_kinematics(self):
            for key in ('Position','Rotation','Velocity','AngularVelocity'):
                if not np.isfinite(self.mbdyn.data[key]).all():
                    raise RuntimeError(f'Non-finite {key}')
            super()._write_kinematics()
        def _read_loads(self):
            super()._read_loads()
            for key in ('Force','Moment'):
                if not np.isfinite(self.mbdyn.data[key]).all():
                    raise RuntimeError(f'Non-finite {key}')
        def _write_diagnostics(self,t,iterations):
            super()._write_diagnostics(t,iterations)
            if iterations>=20:
                raise RuntimeError('Implicit iteration limit reached; run not accepted')
            (path/'progress.json').write_text(json.dumps({'time_s':t,'iterations':iterations})+'\n')
            if not c['smoke'] and not getattr(self,'release_checked',False) and t>=c['sas_off_start_s']-2*c['dt_s']:
                check=release_check(read_state(path),c)
                (path/'release_check.json').write_text(json.dumps(check,indent=2)+'\n')
                self.release_checked=True
                if not check['passed']:
                    raise RuntimeError('DUST settling failed release gate; retune non-yaw controllers. See release_check.json')
    processes=[]
    with tempfile.TemporaryDirectory(prefix='bff_dust_socket_') as temp:
        sock=Path(temp)/'mbdyn.sock';env['MBSOCK']=str(sock)
        with (path/'dust.log').open('w') as df,(path/'mbdyn.log').open('w') as mf:
            try:
                processes.append(subprocess.Popen([bins['DUST_BIN'],'dust.in'],cwd=path,env=env,stdout=df,stderr=subprocess.STDOUT,start_new_session=True))
                processes.append(subprocess.Popen([bins['MBDYN_BIN'],'-f','main.mbd','-o','case'],cwd=path,env=env,stdout=mf,stderr=subprocess.STDOUT,start_new_session=True))
                (path/'solver_pids.json').write_text(json.dumps([p.pid for p in processes]))
                deadline=time.monotonic()+90
                while not sock.exists():
                    if any(p.poll() is not None for p in processes):raise RuntimeError('Solver exited before handshake; see logs')
                    if time.monotonic()>deadline:raise RuntimeError('MBDyn socket timeout')
                    time.sleep(.1)
                structure=MBDynInterface(str(sock),str(path/'coupling_nodes.in'))
                def state(t):
                    if t<c['sas_off_start_s']:return 'SETTLING'
                    if t<c['sas_off_start_s']+c['rap_delay_s']:return 'HOLD'
                    if t<c['sas_off_start_s']+c['rap_delay_s']+c['rap_duration_s']:return 'RAP'
                    if t<c['sas_on_start_s']:return 'OPEN_LOOP'
                    return 'RECOVERY'
                adapter=CheckedAdapter(structure,str(path/'precice.xml'),str(path/'coupled_response.csv'),state)
                adapter.run(c['dt_s'])
                codes=[p.wait(timeout=60) for p in processes]
                (path/'solver_exit_codes.json').write_text(json.dumps(codes))
                if codes!=[0,0]:raise RuntimeError(f'Solver return codes {codes}')
            finally:
                for p in reversed(processes):stop(p)


def run_worker(path,threads,stall_seconds):
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS=str(threads),PYTHONUNBUFFERED='1')
    p=subprocess.Popen([sys.executable,str(ROOT/'run_case.py'),'--worker',str(path),'--threads',str(threads)],
                       cwd=path,env=env,start_new_session=True)
    try:
        last_progress=time.monotonic();stamp=None
        while p.poll() is None:
            file=path/'progress.json'
            new=file.stat().st_mtime_ns if file.exists() else None
            if new!=stamp:stamp=new;last_progress=time.monotonic()
            if time.monotonic()-last_progress>stall_seconds:
                raise RuntimeError(f'No coupling progress for {stall_seconds} seconds; stopped test')
            time.sleep(1)
        if p.returncode:raise RuntimeError(f'Coupling failed ({p.returncode}); inspect {path}')
    finally:
        stop(p)
        # Worker solvers are separate process groups. Stop only recorded, still-live
        # solver PIDs whose command line points to the expected binaries.
        pidfile=path/'solver_pids.json'
        if pidfile.exists():
            for pid in json.loads(pidfile.read_text()):
                try:
                    cmd=Path(f'/proc/{pid}/cmdline').read_bytes().split(b'\0')[0].decode()
                    cwd=Path(f'/proc/{pid}/cwd').resolve()
                    if Path(cmd).name in ('dust','mbdyn') and cwd==path:
                        os.killpg(pid,signal.SIGTERM)
                except (FileNotFoundError,ProcessLookupError):pass


def check(path,bins,env):
    execute([bins['DUST_PRE_BIN'],'dust_pre.in'],path,'dust_pre.log',env)
    shape=check_geometry(path)
    validator='/usr/local/bin/precice-config-validate'
    execute([validator,str(path/'precice.xml'),'MBDyn','1'],path,'precice_validate.log',env)
    with tempfile.TemporaryDirectory(prefix='bff_dust_parse_') as temp:
        socket=Path(temp)/'parse.sock';env=env.copy();env['MBSOCK']=str(socket)
        with (path/'parse.log').open('w') as f:
            p=subprocess.Popen([bins['MBDYN_BIN'],'-f','main.mbd','-o','parse'],cwd=path,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
            try:
                deadline=time.monotonic()+60
                while not socket.exists() and p.poll() is None and time.monotonic()<deadline:time.sleep(.1)
                if not socket.exists():raise RuntimeError('MBDyn parse did not reach socket; see parse.log')
            finally:stop(p)
    result={'check_pass':True,'geometry':shape,'scope':'FEM25, input and XML validation; no time step'}
    (path/'check.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


def postprocess(path,bins,env,c):
    frames=len(list((path/'dust').glob('case_res_*.h5')))
    if not frames:return
    (path/'paraview').mkdir(exist_ok=True)
    text=f'''basename = paraview/x56
data_basename = dust/case
analysis = {{
 type = viz
 name = bff
 start_res = 1
 end_res = {frames}
 step_res = 1
 format = vtk
 wake = T
 separate_wake = T
 variable = pressure
 component = all
}}
'''
    (path/'dust_post.in').write_text(text)
    execute([bins['DUST_POST_BIN'],'dust_post.in'],path,'dust_post.log',env,timeout=600)
    execute([sys.executable,str(COUPLED/'tools/fix_vtu_xml.py'),str(path/'paraview')],path,'vtk_fix.log',env)
    for suffix in ('','_wpan','_wpart'):
        lines=['<?xml version="1.0"?>','<VTKFile type="Collection" version="0.1" byte_order="LittleEndian"><Collection>']
        for i in range(1,frames+1):
            lines.append(f'<DataSet timestep="{i*c["output_dt_s"]:.8f}" file="x56_bff{suffix}-{i:04d}.vtu"/>')
        lines.append('</Collection></VTKFile>')
        (path/'paraview'/f'x56{suffix}.pvd').write_text('\n'.join(lines))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    for name in ('prepare','check','smoke','run'):mode.add_argument('--'+name,action='store_true')
    mode.add_argument('--analyse',type=Path)
    mode.add_argument('--worker',type=Path,help=argparse.SUPPRESS)
    parser.add_argument('--config',type=Path,default=ROOT/'case.json')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--threads',type=int,default=12)
    parser.add_argument('--stall-seconds',type=float,default=900)
    args=parser.parse_args()
    if args.threads<1:parser.error('threads must be positive')
    if args.worker:worker(args.worker.resolve(),args.threads);return
    if args.analyse:print(json.dumps(analyse(args.analyse.resolve()),indent=2));return
    config=load_config(args.config)
    kind='smoke' if args.smoke else 'check' if args.check else 'run' if args.run else 'prepared'
    output=(args.output or ROOT/'runs'/f'{kind}_{datetime.now():%Y%m%d_%H%M%S}').resolve()
    c=build(output,config,smoke=args.smoke)
    print(f'Case: {output}\nV={c["velocity_mps"]} m/s = {c["velocity_mps"]/.0254:.6f} in/s; '
          f'{c["elastic_modes"]} elastic modes; dt={c["dt_s"]} s',flush=True)
    if args.prepare:return
    bins=binaries();env=os.environ.copy();env.update(OMP_NUM_THREADS=str(args.threads),OPENBLAS_NUM_THREADS='1')
    env.pop('DISPLAY',None);env.pop('WAYLAND_DISPLAY',None)
    check(output,bins,env)
    if args.check:return
    run_worker(output,args.threads,args.stall_seconds)
    report=analyse(output)
    postprocess(output,bins,env,c)
    print(json.dumps(report,indent=2))
    if not report['technical_pass']:raise RuntimeError('Acceptance failed; see analysis.json')
    if args.run and not report['identification_valid']:raise RuntimeError('Run completed but BFF identification not accepted; see analysis.json')


if __name__=='__main__':
    try:main()
    except KeyboardInterrupt:raise SystemExit(130)
