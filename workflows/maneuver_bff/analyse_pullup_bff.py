#!/usr/bin/env python3
"""Reassess BFF excitation/growth using raw trajectories, not short-window Hilbert slopes.

Variable projection fits polynomial drift plus an exponentially modulated sinusoid.
Separate q/qdot-constrained fits and cycle amplitudes check growth independently.
No solver calls and no changes to previous analysis outputs.
"""
from __future__ import annotations
import argparse
import csv
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/pullup_bff_matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from scipy.optimize import least_squares
from scipy.signal import find_peaks

from analyse_sweep import discover, constant
from analyse_time_domain_pairs import read_symmetric_tip
from compare_paired_response import growth_metrics

ROOT = Path(__file__).resolve().parent
DATA = Path('/mnt/c/Users/Utente/Desktop/BFF_PULLUP_V2')


def read(path):
    with Dataset(path) as nc:
        out = {key: np.asarray(nc[name][:], dtype=float) for key,name in (
            ('t','time'),('q','elem.joint.5.a'),('qd','elem.joint.5.aPrime'),
            ('omega','node.struct.990000.Omega'))}
    out['t'] = out['t'].squeeze()
    _, out['tip'] = read_symmetric_tip(path)
    return out


def fit(t, y, degree=2, yd=None, pole=None):
    t=np.asarray(t); y=np.asarray(y); x=t-t.mean()
    scale=max(float(np.std(y)),1e-15)
    dscale=max(float(np.std(yd)),1e-15) if yd is not None else None
    def design(p):
        sig,fr=p; w=2*np.pi*fr; e=np.exp(sig*x); c=np.cos(w*x); s=np.sin(w*x)
        D=np.column_stack([x**i for i in range(degree+1)]+[e*c,e*s])
        if yd is None:
            A=D/scale; b=y/scale
        else:
            P=np.column_stack([np.zeros_like(x) if i==0 else i*x**(i-1) for i in range(degree+1)]+
                              [e*(sig*c-w*s),e*(sig*s+w*c)])
            A=np.vstack([D/scale,P/dscale]);b=np.r_[y/scale,yd/dscale]
        coeff=np.linalg.lstsq(A,b,rcond=None)[0]
        return D,coeff,A@coeff-b
    if pole is None:
        candidates=[least_squares(lambda p:design(p)[2],[.3,f],bounds=([-3.,1.],[3.,3.5]),
                                  xtol=1e-10,ftol=1e-10,gtol=1e-10)
                    for f in (1.4,1.9,2.5,3.1)]
        best=min(candidates,key=lambda r:r.cost);pole=best.x
    D,c,res=design(pole); trend=D[:,:degree+1]@c[:degree+1]; osc=D[:,-2:]@c[-2:]
    drift_only=np.column_stack([x**i for i in range(degree+1)])
    drift_only=drift_only@np.linalg.lstsq(drift_only,y,rcond=None)[0]
    sse=float(np.sum((y-trend-osc)**2)); sst=float(np.sum((y-y.mean())**2))
    detrended_sse=float(np.sum((y-drift_only)**2))
    return dict(sigma=float(pole[0]),frequency=float(pole[1]),amplitude=float(np.hypot(*c[-2:])),
                center_time=float(t.mean()),phase=float(np.arctan2(-c[-1],c[-2])),
                r2=1-sse/max(sst,1e-30),oscillatory_variance_explained=1-sse/max(detrended_sse,1e-30),
                normalized_residual=float(np.sqrt(np.mean(res**2))),trend=trend,osc=osc)


def scalar(result):
    return {k:v for k,v in result.items() if np.isscalar(v)}


def cycles(t,y,fr):
    # Peak-to-adjacent-trough amplitudes remove slow offset without filtering.
    dt=float(np.median(np.diff(t))); distance=max(1,int(.7/fr/dt))
    maxima=find_peaks(y,distance=distance)[0]; minima=find_peaks(-y,distance=distance)[0]
    times=[]; amps=[]
    for i in maxima:
        before=minima[minima<i];after=minima[minima>i]
        if len(before) and len(after):
            lo,hi=before[-1],after[0]
            trough=np.interp(t[i],[t[lo],t[hi]],[y[lo],y[hi]])
            times.append(float(t[i]));amps.append(float((y[i]-trough)/2))
    sig=float(np.polyfit(times,np.log(amps),1)[0]) if len(times)>=2 and min(amps)>0 else None
    return dict(times=times,amplitudes=amps,sigma=sig)


def synthetic_check():
    t=np.arange(0,1.54,.02);rows=[]
    for phase in (0.,.7,1.4,2.1):
        y=np.exp(.4*t)*np.sin(2*np.pi*1.9*t+phase)
        old,_,_=growth_metrics(t,y,0.,float(t[-1]),2.0575)
        direct=fit(t,y)
        assert abs(direct['sigma']-.4)<1e-7 and abs(direct['frequency']-1.9)<1e-7
        rows.append(dict(true_sigma=.4,phase=phase,old_hilbert_sigma=old['sigma_per_s'],
                         direct_sigma=direct['sigma']))
    return rows


def analyze_pairs(output):
    rows=[];variants=[];cache={}
    for campaign in ('primary','timestep'):
        groups=discover(DATA/campaign/'cases')
        for key,pair in sorted(groups.items()):
            if set(pair)!={False,True}:continue
            s=read(pair[False][0]);e=read(pair[True][0]);mbd=pair[False][1]
            if not np.allclose(s['t'],e['t'],rtol=0,atol=1e-9):raise ValueError('Unmatched pair')
            rel=constant(mbd,'SAS_OFF_START');t=s['t']-rel
            dq=e['q'][:,0]-s['q'][:,0];dqd=e['qd'][:,0]-s['qd'][:,0]
            dtip=e['tip']-s['tip'];dw=e['omega'][:,1]-s['omega'][:,1]
            start=.1+.742/constant(mbd,'BFF_RAP_TARGET_FREQUENCY');end=2.
            mask=(t>=start)&(t<=end)
            qfit=fit(t[mask],dq[mask]);joint=fit(t[mask],dq[mask],yd=dqd[mask])
            tipfit=fit(t[mask],dtip[mask]);pitchfit=fit(t[mask],dw[mask])
            pole=(qfit['sigma'],qfit['frequency'])
            shadow=fit(t[mask],s['q'][mask,0],pole=pole)
            shadowtip=fit(t[mask],s['tip'][mask],pole=pole)
            peak=cycles(t[mask],dq[mask],qfit['frequency'])
            old,_,_=growth_metrics(s['t'],dq,rel+start,rel+end,constant(mbd,'BFF_RAP_TARGET_FREQUENCY'))
            row=dict(campaign=campaign,velocity=key[1],nominal_n=key[2],dt=key[3],
                     saved_dt=float(np.median(np.diff(t))),reference_release=rel,
                     source_shadow=str(pair[False][0]),source_excited=str(pair[True][0]),
                     direct=scalar(qfit),joint_q_qdot=scalar(joint),tip=scalar(tipfit),pitch_rate=scalar(pitchfit),
                     peak_to_trough=peak,old_hilbert_sigma=old['sigma_per_s'],
                     shadow_q=scalar(shadow),shadow_tip=scalar(shadowtip))
            sensitivity=[]
            for degree in (1,2,3):
                for trim in (0.,.1,.2):
                    m=(t>=start+trim)&(t<=end-trim)
                    f=fit(t[m],dq[m],degree=degree)
                    sensitivity.append(dict(degree=degree,trim_each_end=trim,**scalar(f)))
            row['window_trend_sensitivity']=sensitivity
            rows.append(row)
            if campaign=='primary':cache[(key[1],key[2])]=(t,dq,s,mask,qfit,shadowtip)
            print(campaign,key[1:3], 'sigma',round(qfit['sigma'],4),'Hz',round(qfit['frequency'],4),flush=True)
    for row in rows:
        baseline=next((r for r in rows if r['campaign']=='primary' and r['velocity']==row['velocity'] and r['nominal_n']==1.),None)
        row['delta_sigma_vs_primary_1g']=row['direct']['sigma']-baseline['direct']['sigma'] if baseline else None
        row['delta_frequency_vs_primary_1g']=row['direct']['frequency']-baseline['direct']['frequency'] if baseline else None
        if baseline:
            row['shadow_tip_amplitude_ratio_vs_primary_1g']=row['shadow_tip']['amplitude']/baseline['shadow_tip']['amplitude']
    fig,axes=plt.subplots(2,3,figsize=(13,7))
    for j,v in enumerate((66.25,66.75,67.25)):
        for n in (1.,1.3,1.6):
            t,dq,s,mask,f,st=cache[(v,n)];idx=t>=.45
            axes[0,j].plot(t[idx],dq[idx],label=f'nnom={n:g}')
            # Show raw shadow after subtracting only its fitted slow drift.
            axes[1,j].plot(t[mask],1000*(s['tip'][mask]-st['trend']),label=f'nnom={n:g}')
        axes[0,j].set_title(f'V={v:g} m/s');axes[0,j].set_ylabel('Delta q7: excited - shadow')
        axes[1,j].set_ylabel('Shadow symmetric tip [mm]\nslow polynomial drift removed')
        for ax in axes[:,j]:ax.set_xlabel('Time after SAS release [s]');ax.grid(alpha=.25);ax.legend()
    fig.suptitle('Existing simulations: incremental growth and oscillation without added rap')
    fig.tight_layout();fig.savefig(output/'existing_trajectories.png',dpi=170);plt.close(fig)
    return rows


def causal_controls(root,rows,output):
    design_path=root/'design.json'
    if not design_path.exists():return dict(status='not_prepared')
    design=json.loads(design_path.read_text());missing=[];cases={}
    reference=next(r for r in rows if r['campaign']=='primary' and r['velocity']==67.25 and r['nominal_n']==1.6)
    pole=(reference['direct']['sigma'],reference['direct']['frequency'])
    definitions=[dict(name='pullup_release',prefix=str(Path(design['reused_pullup_release']).with_suffix('')),
                      reference_release_s=reference['reference_release'],required_final_s=reference['reference_release']+2.05)] + design['cases']
    result=[]
    for row in definitions:
        path=Path(row['prefix']).with_suffix('.nc')
        if not path.exists():missing.append(row['name']);continue
        # Avoid reading an actively written NetCDF.
        if row['name']!='pullup_release':
            status=Path(row['prefix']).with_suffix('.status.json')
            if not status.exists() or not json.loads(status.read_text()).get('complete'):
                missing.append(row['name']);continue
        d=read(path);t=d['t']-row['reference_release_s'];cases[row['name']]=(t,d)
        metrics={}
        for label,a,b in [('pre_release',-1.8,-.1),('post_release',.47,2.)]:
            m=(t>=a)&(t<=b)
            # Open-loop pole is only a projection basis for continuous-SAS data,
            # not an assertion that the closed-loop pole is the same.
            f=fit(t[m],d['tip'][m],pole=pole if label=='post_release' else (0.,pole[1]))
            metrics[label]=scalar(f)
            metrics[label]['amplitude_mm']=1000*f['amplitude']
        result.append(dict(name=row['name'],path=str(path),metrics=metrics))
    if len(cases)==4:
        fig,axes=plt.subplots(2,2,figsize=(12,7))
        for j,names in enumerate((('sham_release','pullup_release'),('sham_sas_continuous','pullup_sas_continuous'))):
            for name in names:
                t,d=cases[name];m=(t>=-2.)&(t<=2.)
                axes[0,j].plot(t[m],d['q'][m,0],label=name)
                m=(t>=.47)&(t<=2.);f=fit(t[m],d['tip'][m],pole=pole)
                axes[1,j].plot(t[m],1000*(d['tip'][m]-f['trend']),label=name)
            axes[0,j].set_title(('SAS released','SAS continuous')[j])
            axes[0,j].set_ylabel('Raw q7 (includes mean wing bending)')
            axes[1,j].set_ylabel('Tip oscillation [mm], slow drift removed')
            for ax in axes[:,j]:ax.axvline(0,color='gray',ls=':');ax.grid(alpha=.25);ax.legend(fontsize=8);ax.set_xlabel('Time relative to reference release [s]')
        fig.tight_layout();fig.savefig(output/'causal_controls.png',dpi=170);plt.close(fig)
    return dict(status='complete' if not missing else 'pending',missing=missing,cases=result,
                warning='Continuous-SAS amplitudes are projections onto the open-loop BFF signature, not an identified closed-loop flutter pole.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'results/pullup_bff')
    p.add_argument('--controls',type=Path,default=Path('/home/nicomonzi/ZENO/BFF_PULLUP_CAUSAL_READY'))
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    synth=synthetic_check();rows=analyze_pairs(args.output)
    controls=causal_controls(args.controls,rows,args.output)
    result=dict(method='Unfiltered polynomial-plus-damped-sinusoid variable projection; independent derivative, tip, pitch and extrema checks.',
                synthetic=synth,paired=rows,causal_controls=controls,
                limitations=['Finite-window effective pole, not a stationary eigenvalue certificate.',
                             'DLM calibration and constant structural stiffness retained.',
                             'Shadow contains maneuver history and possible SAS-release transients; causal controls distinguish their combined effects.',
                             'Window/trend spreads are robustness diagnostics, not statistical confidence intervals.'])
    (args.output/'analysis.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    flat=[dict(campaign=r['campaign'],V=r['velocity'],n=r['nominal_n'],dt=r['dt'],
               sigma=r['direct']['sigma'],frequency=r['direct']['frequency'],
               sigma_joint=r['joint_q_qdot']['sigma'],sigma_tip=r['tip']['sigma'],
               sigma_pitch_rate=r['pitch_rate']['sigma'],old_sigma=r['old_hilbert_sigma'],
               delta_sigma_vs_1g=r['delta_sigma_vs_primary_1g'],
               shadow_tip_amplitude_mm=1000*r['shadow_tip']['amplitude'],
               shadow_oscillation_explained=r['shadow_tip']['oscillatory_variance_explained']) for r in rows]
    with (args.output/'paired_summary.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(flat)
    print(json.dumps(dict(output=str(args.output),controls=controls['status'],synthetic=synth),indent=2))


if __name__=='__main__':main()
