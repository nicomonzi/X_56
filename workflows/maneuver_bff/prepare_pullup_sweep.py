#!/usr/bin/env python3
"""Prepare in TESI; execute only on explicit user request; results on Desktop.

28 additional cases: 18 causal controls, 8 lower-speed paired trajectories,
2 fine-step 1g trajectories. Existing primary/timestep results are read only.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import json
from pathlib import Path
import re
import signal
import subprocess
import sys

from campaign import StudyCase, render_case, source_fingerprints
import run_pullup_causal_controls as execution

ROOT=Path(__file__).resolve().parent
PACKAGE=ROOT/'sweep_pullup_v3'
DESKTOP=Path('/mnt/c/Users/Utente/Desktop/BFF_PULLUP_V2')
RESULTS=DESKTOP/'pullup_bff_sweep_v3'


def write_same(path,text):
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists() and path.read_text()!=text:
        raise RuntimeError(f'Refusing changed-file overwrite: {path}')
    if not path.exists():path.write_text(text)


def number(text,name):
    return float(re.search(rf'set: const real {name} = ([^;]+);',text)[1])


def tag(value):return f'{value:.3f}'.replace('.','p')


def prepare(package=PACKAGE,results=RESULTS):
    rows=[];designs=[]
    def add(name,group,text,info):
        input_path=package/'inputs'/group/(name+'.mbd')
        text='# PULLUP_SWEEP_V3 '+json.dumps(info,sort_keys=True)+'\n'+text
        write_same(input_path,text)
        row=dict(name=name,group=group,input_path=str(input_path),
                 prefix=str(results/group/'cases'/name),
                 input_sha256=hashlib.sha256(text.encode()).hexdigest(),
                 required_final_s=number(text,'FINAL_TIME'),
                 reference_release_s=number(text,'SAS_OFF_START'),**info)
        rows.append(row);return row
    for v in (66.25,66.75,67.25):
        for n in (1.3,1.6):
            matches=list((DESKTOP/'primary/cases').glob(f'*V_{v:07.3f}'.replace('.','p')+f'*nnom_0{n:04.2f}'.replace('.','p')+'_shadow.mbd'))
            if len(matches)!=1:raise RuntimeError(('Missing/ambiguous original',v,n,matches))
            source=matches[0];original=source.read_text();key=f'V_{tag(v)}_n_{tag(n)}'
            group='causal/'+key;control_rows=[]
            original_meta=json.loads(re.search(r'(?m)^# STIFFNESS_STUDY_METADATA (.+)$',original)[1])
            if original.count(execution.GATE)!=2:raise RuntimeError('Unexpected SAS gate layout')
            for name,maneuver,release_sas in (('sham_release',False,True),
                                             ('pullup_sas_continuous',True,False),
                                             ('sham_sas_continuous',False,False)):
                text=original
                if not maneuver:
                    for label in ('PULLUP_PITCH_ANGLE','DIVE_PITCH_RATE_COMMAND'):
                        text=execution.replace_constant(text,label,'0.')
                    text=execution.replace_constant(text,'DIVE_NOMINAL_LOAD_CLASS','1.')
                if not release_sas:text=text.replace(execution.GATE,'1.')
                text=text.replace('closest next, SAS_OFF_START - STUDY_OUTPUT_PREHISTORY,','closest next, 7.,')
                meta=dict(original_meta,campaign='causal_controls',excited=False,
                          nominal_load_factor=n if maneuver else 1.,
                          pitch_amplitude_deg=original_meta['pitch_amplitude_deg'] if maneuver else 0.,
                          pitch_rate_command_deg_s=original_meta['pitch_rate_command_deg_s'] if maneuver else 0.,
                          modification=name)
                text=re.sub(r'(?m)^# STIFFNESS_STUDY_METADATA .+$','# STIFFNESS_STUDY_METADATA '+json.dumps(meta,sort_keys=True),text)
                text=re.sub(r'(?m)^# MANEUVER_METADATA .+\n','',text)
                info=dict(velocity_mps=v,nominal_load_factor=n,control=name,
                          maneuver=maneuver,sas_release=release_sas,dt=.01,excited=False,
                          source=str(source),source_sha256=hashlib.sha256(original.encode()).hexdigest())
                row=add(key+'__'+name,group,text,info)
                control_rows.append(dict(row,name=name))
            design=dict(velocity_mps=v,nominal_load_factor=n,
                        reused_pullup_release=str(source.with_suffix('.nc')),cases=control_rows)
            write_same(package/'controls'/key/'design.json',json.dumps(design,indent=2)+'\n')
            designs.append(str(package/'controls'/key/'design.json'))
    for group,velocities,loads,dt in (
        ('paired_extension',(65.,65.75),(1.,1.6),.01),
        ('timestep_extension',(66.75,),(1.,),.005)):
        for v in velocities:
            for n in loads:
                for excited in (False,True):
                    case=StudyCase('primary',v,n,dt,1.,excited)
                    text=render_case(case)
                    metadata=case.metadata();metadata['campaign']=group
                    text=re.sub(r'(?m)^# STIFFNESS_STUDY_METADATA .+$','# STIFFNESS_STUDY_METADATA '+json.dumps(metadata,sort_keys=True),text)
                    text=text.replace('closest next, SAS_OFF_START - STUDY_OUTPUT_PREHISTORY,','closest next, 7.,')
                    name=f'V_{tag(v)}_n_{tag(n)}_dt_{dt:.4f}_'+('excited' if excited else 'shadow')
                    # A decimal point in a solver prefix would be mistaken for a suffix.
                    name=name.replace('.','p')
                    add(name,group,text,dict(velocity_mps=v,nominal_load_factor=n,
                                          control='paired',maneuver=n>1,sas_release=True,
                                          dt=dt,excited=excited))
    manifest=dict(version=3,case_count=len(rows),input_root=str(package),result_root=str(results),
                  reused_campaigns=[str(DESKTOP/'primary'),str(DESKTOP/'timestep')],
                  controls=designs,source_sha256=source_fingerprints(),cases=rows)
    write_same(package/'manifest.json',json.dumps(manifest,indent=2)+'\n')
    import io
    buf=io.StringIO();fields=['name','group','velocity_mps','nominal_load_factor','dt','excited','input_path','prefix']
    writer=csv.DictWriter(buf,fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerows(rows)
    write_same(package/'manifest.csv',buf.getvalue())
    print(f'[prepared only] {len(rows)} inputs in {package}/inputs',flush=True)
    print(f'[results when YOU execute] {results}',flush=True)
    return manifest


def execute_case(row):
    if execution.STOP.is_set():return dict(name=row['name'],complete=False,cancelled_before_start=True)
    text=Path(row['input_path']).read_text()
    if hashlib.sha256(text.encode()).hexdigest()!=row['input_sha256']:
        raise RuntimeError(f"Input modified after preparation: {row['input_path']}")
    prefix=Path(row['prefix']);write_same(prefix.with_suffix('.mbd'),text)
    # Existing successful outputs can be reused only with exactly matching inputs.
    return execution.run(row)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--package',type=Path,default=PACKAGE)
    p.add_argument('--results',type=Path,default=RESULTS)
    p.add_argument('--execute',action='store_true')
    p.add_argument('--jobs',type=int,default=2)
    p.add_argument('--stage',choices=('all','paired','causal','timestep'),default='all')
    p.add_argument('--analyse',action='store_true',help='Run postprocessing after a successful explicit execution, or on existing results')
    args=p.parse_args()
    if args.jobs<1:p.error('--jobs must be positive')
    manifest=prepare(args.package.resolve(),args.results.resolve())
    if args.execute:
        if manifest['source_sha256']!=source_fingerprints():raise RuntimeError('Model dependencies changed')
        signal.signal(signal.SIGINT,execution.stop_children)
        signal.signal(signal.SIGTERM,execution.stop_children)
        chosen=[r for r in manifest['cases'] if args.stage=='all' or
                (args.stage=='causal' and r['group'].startswith('causal/')) or
                (args.stage=='paired' and r['group']=='paired_extension') or
                (args.stage=='timestep' and r['group']=='timestep_extension')]
        # Establish the growth bracket before spending time on causal controls.
        chosen.sort(key=lambda r:(0 if r['group']=='paired_extension' else 1 if r['group']=='timestep_extension' else 2,r['name']))
        outcomes=[]
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures=[pool.submit(execute_case,row) for row in chosen]
            for i,future in enumerate(as_completed(futures),1):
                try:result=future.result()
                except Exception as error:
                    execution.stop_children(signal.SIGTERM,None)
                    raise RuntimeError('Sweep stopped after an execution error') from error
                outcomes.append(result);print(f'[{i}/{len(chosen)}] '+json.dumps(result),flush=True)
        if execution.STOP.is_set() or not all(r.get('complete') for r in outcomes):raise SystemExit(1)
    if args.analyse:
        subprocess.run([sys.executable,str(ROOT/'analyse_pullup_bff.py'),
                        '--sweep-root',str(args.results),'--controls',str(args.package/'controls'),
                        '--output',str(args.results/'analysis')],check=True)


if __name__=='__main__':main()
