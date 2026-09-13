#!/usr/bin/env python3
"""Exact real-slice invariant attribution, with explicit interference.

One Hessian per invariant suffices: block multihomogeneity determines all
radial derivatives and their continuation to any positive radial point.
The attribution freezes propagators; candidate physics must rebuild them.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
import verify_p54_common_renormalization as common
import verify_p54_goldstone_ir as ir

RF=Path(__file__).resolve().parents[1]
OUT=RF/'output/p54_invariant_decomposition'
LOOP2=32*np.pi**2
# Total real degrees in (Phi54,Sigma126,H10,S), including conjugates.
DEGREES={
 'mu2':(2,0,0,0),'c':(3,0,0,0),'a':(4,0,0,0),'b':(4,0,0,0),
 'nu2':(0,2,0,0),'lambda0':(0,4,0,0),'lambda2':(0,4,0,0),
 'lambda4':(0,4,0,0),'lambda4p':(0,4,0,0),
 'alpha':(2,2,0,0),'beta':(2,2,0,0),'xi02':(0,0,2,0),
 'xi1':(0,0,4,0),'xi2':(0,0,4,0),'xi3':(1,0,2,0),
 'gamma1':(0,2,2,0),'gamma2':(0,2,2,0),
 'eta0':(2,0,2,0),'eta1':(0,3,1,0),'eta2':(2,0,2,0),'eta3':(0,2,2,0),
 'mus2':(0,0,0,2),'chi1':(0,0,0,4),'chi2':(0,2,0,2),
 'chi3':(2,0,0,2),'chi4':(1,2,0,1),'chi5':(0,0,2,2),
 'chi6':(0,0,2,1),'chi7':(1,0,2,1)}
SLICES=(slice(0,54),slice(54,306),slice(306,326),slice(326,328))
RADIAL_SECTORS=(0,1,3)


class InvariantJets:
    def __init__(self,p1,parameters,r0,cache_dir,check=None):
        if set(parameters)!=set(DEGREES):
            raise ValueError('Invariant catalogue no longer matches the action')
        if any(np.iscomplexobj(v) and np.imag(v)!=0 for v in parameters.values()):
            raise ValueError('This verifier is explicitly restricted to real coupling directions')
        self.p1=p1;self.names=list(DEGREES);self.r0=np.asarray(r0)
        self.units=[];self.hits=0;self.evaluated=0
        x=p1.vacuum_vector(*r0)
        for name in self.names:
            p={k:0. for k in parameters};p[name]=1.
            store=common.scalar.ActionJets(p1,p,cache_dir)
            if name in ('xi1','xi2'):
                h=np.zeros((328,328)) # H^4 has identically zero Hessian at H=0.
            else:
                h=store.hessian(x,'unit_invariant_'+name)
                p1.jax.clear_caches()
            self.hits+=store.hits;self.evaluated+=store.evaluated
            clean=h.copy()
            for a,sa in enumerate(SLICES):
                for b,sb in enumerate(SLICES):
                    e=np.array(DEGREES[name]);e[a]-=1;e[b]-=1
                    if min(e)<0 or e[2]!=0:
                        if check:check(name+f'_forbidden_block_{a}{b}',h[sa,sb])
                        clean[sa,sb]=0.
            self.units.append(clean)
        self.units=np.asarray(self.units)

    def component(self,index,r):
        """Unit-coefficient ambient H,T,U at arbitrary strictly positive r."""
        r=np.asarray(r)
        if np.any(r<=0):raise ValueError('Positive radial coordinates required')
        h=np.zeros((328,328));t=np.zeros((3,328,328));u=np.zeros((3,3,328,328))
        for a,sa in enumerate(SLICES):
            for b,sb in enumerate(SLICES):
                e=np.array(DEGREES[self.names[index]]);e[a]-=1;e[b]-=1
                if min(e)<0 or e[2]!=0:continue
                er=e[list(RADIAL_SECTORS)]
                block=self.units[index,sa,sb]*np.prod((r/self.r0)**er)
                h[sa,sb]=block
                for j in range(3):
                    t[j,sa,sb]=er[j]/r[j]*block
                    for k in range(3):
                        u[j,k,sa,sb]=er[j]*(er[k]-(j==k))/(r[j]*r[k])*block
        return h,t,u

    def assemble(self,parameters,r):
        h=np.zeros((328,328));t=np.zeros((3,328,328));u=np.zeros((3,3,328,328))
        for i,name in enumerate(self.names):
            hi,ti,ui=self.component(i,r)
            h+=parameters[name]*hi;t+=parameters[name]*ti;u+=parameters[name]*ui
        return h,t,u


def weights(lam,mu2):
    m=np.array(lam);m[:38]=0.
    mask=np.ones((328,328),bool);mask[:38,:38]=False
    slope=np.zeros((328,328));log=np.zeros_like(slope)
    xx,yy=np.broadcast_arrays(m[:,None],m[None,:])
    slope[mask]=common.scalar.scalar.bubble_slope(xx[mask],yy[mask])
    # Degeneracies are handled by the analytic equal-mass branch.
    for i in range(328):
        for j in range(i,328):
            if mask[i,j]:log[i,j]=log[j,i]=ir.log_bubble(m[i],m[j],0.,mu2)
    a0=np.zeros(328);a0[38:]=m[38:]*(np.log(m[38:]/mu2)-1)
    return a0,log,slope


def run(cache_dir):
    oldpath=RF/'output/p54_full_doublet_cw.json';irpath=RF/'output/p54_goldstone_ir.json'
    old=json.loads(oldpath.read_text());oldir=json.loads(irpath.read_text())
    p1=common.yuk.module('inv_p1',RF/'code/verify_p54_p1_hessian_spectrum.py')
    cw=common.yuk.module('inv_cw',RF/'code/verify_p54_p2_bosonic_cw.py')
    p=old['tree_parameters'];r=np.array([old['vacuum'][k] for k in ('omega','sigma','vs')])
    mu=old['scheme']['mu_over_omega'];mu2=mu*mu
    rad=np.column_stack([p1.vacuum_vector(*e) for e in np.eye(3)]);metric=rad.T@rad
    checks=[]
    def check(name,a,b=0.,tol=3e-8):
        err=common.yuk.error(np.asarray(a),np.asarray(b))
        checks.append(dict(name=name,residual=err,tolerance=tol,pass_=bool(err<tol)))
    inv=InvariantJets(p1,p,r,cache_dir,check)
    h,t,u=inv.assemble(p,r);store=common.scalar.ActionJets(p1,p,cache_dir);x=rad@r
    oldh=store.hessian(x,'full_base')
    check('sum_of_invariants_reconstructs_full_action_Hessian',h,oldh)
    for j in range(3):
        tj,uj=store.jets(x,oldh,rad[:,j],f'independent_radial_{j}')
        check(f'homogeneity_T_{j}_equals_independent_action_jets',t[j],tj)
        check(f'homogeneity_U_{j}{j}_equals_independent_action_jets',u[j,j],uj)
        for k in range(j+1,3):
            step=.02
            hp=store.hessian(x+step*(rad[:,j]+rad[:,k]),f'independent_pair_{j}{k}')
            pair=(hp-oldh-step*(t[j]+t[k])-.5*step**2*(u[j,j]+u[k,k]))/step**2
            check(f'homogeneity_U_{j}{k}_equals_independent_action_jets',u[j,k],pair)
    # An independent new positive radial point tests the continuation,
    # not just the derivatives used in old caches.
    probe=np.array([.91,.21,.31]);hp,_,_=inv.assemble(p,probe)
    check('nonproportional_radial_continuation_independent_full_action',hp,
          store.hessian(rad@probe,'independent_nonproportional_r'))
    lam,q=np.linalg.eigh(h)
    check('38_tree_zero_modes',lam[:38]);check('290_tree_positive_modes',float(lam[38]<=0))
    a0,log,slope=weights(lam,mu2)
    n=len(inv.names);te=np.zeros((n,3,328,328));linear=np.zeros((n,3,3));treeparts=[]
    for a,name in enumerate(inv.names):
        ha,ta,ua=inv.component(a,r);ha*=p[name];ta*=p[name];ua*=p[name]
        te[a]=np.array([q.T@v@q for v in ta])
        tad=np.einsum('i,rii->r',a0,te[a])/LOOP2
        dp=-np.linalg.solve(common.radial_mass_map(r),tad)
        for j in range(3):
            for k in range(3):linear[a,j,k]=np.dot(a0,np.diag(q.T@ua[j,k]@q))/LOOP2
        linear[a]+=rad.T@common.mass_operator(dp)@rad
        treeparts.append((rad.T@ha@rad).tolist())
    # Tensor axes [invariant_a,invariant_b,radial_r,radial_s].
    def pairs(weight):
        flat=te.reshape(n*3,-1)
        result=(flat*weight.reshape(1,-1))@flat.T/LOOP2
        return result.reshape(n,3,n,3).transpose(0,2,1,3)
    masspair=pairs(log);kinpair=pairs(slope)
    # Each cross block may be nonsymmetric in r,s; only its paired sum
    # with reversed invariant labels is a physical symmetric matrix.
    hs=linear.sum(0)+masspair.sum((0,1));dz=kinpair.sum((0,1))
    # Authoritative old scalar+CT matrix is also reconstructed independently
    # from the full trace functional, avoiding fragile report key aliases.
    full_tad=np.einsum('i,rii->r',a0,te.sum(0))/LOOP2
    dp=-np.linalg.solve(common.radial_mass_map(r),full_tad)
    direct=np.zeros((3,3))+rad.T@common.mass_operator(dp)@rad
    for j in range(3):
        for k in range(3):direct[j,k]+=cw.hard_trace_hessian(lam,q,t[j],t[k],u[j,k],290,mu2,1.5,1.)
    check('linear_plus_all_interferences_equals_full_scalar_fixed_VEV_Hessian',hs,direct)
    check('full_mass_interference_sum_equals_P_IR1',hs,
          oldir['uniform_scalar_parameter_repair']['scalar_fixed_matrix'])
    check('full_kinetic_interference_sum_equals_P_IR1',dz,oldir['hard_scalar_kinetic']['matrix'])
    check('kinetic_radial_sum_PSD',max(0.,-eigh(dz,metric,eigvals_only=True).min()))
    joint=kinpair.transpose(0,2,1,3).reshape(n*3,n*3)
    check('kinetic_full_joint_invariant_radial_Gram_PSD',max(0.,-np.linalg.eigvalsh(joint).min()))
    check('ordered_mass_pair_transpose_symmetry',masspair,masspair.transpose(1,0,3,2))
    check('ordered_kinetic_pair_transpose_symmetry',kinpair,kinpair.transpose(1,0,3,2))
    # Signed symmetric row allocation is an accounting convention, not an
    # independently measurable per-coupling contribution.
    massalloc=linear+.5*(masspair.sum(1)+masspair.sum(0))
    kinalloc=.5*(kinpair.sum(1)+kinpair.sum(0))
    check('mass_signed_allocations_sum',massalloc.sum(0),hs)
    check('kinetic_signed_allocations_sum',kinalloc.sum(0),dz)
    zg,zv=eigh(dz,metric);witness=zv[:,-1]
    # The sigma-axis is separately fixed by the original failure, not by
    # an optimized allocation definition.
    sigma=np.array([0.,1/np.sqrt(2),0.])
    rows=[]
    for a,name in enumerate(inv.names):
        rows.append(dict(invariant=name,coefficient=p[name],degree=DEGREES[name],
            tree_radial_matrix=treeparts[a],linear_seagull_plus_CT=linear[a].tolist(),
            signed_mass_allocation=massalloc[a].tolist(),signed_kinetic_allocation=kinalloc[a].tolist(),
            sigma_mass_allocation=float(sigma@massalloc[a]@sigma),
            leading_kinetic_allocation=float(witness@kinalloc[a]@witness),
            leading_kinetic_self=float(witness@kinpair[a,a]@witness)))
    # On the fixed-weight vertex functional, Euler identity supplies a
    # separate finite-difference check of quadratic interference.
    rng=np.random.default_rng(20260913);v=rng.normal(size=n);eps=1e-5
    def frozen(coeff):
        return np.einsum('a,ars->rs',coeff,linear)+np.einsum('a,b,abrs->rs',coeff,coeff,masspair)
    derivative=np.einsum('a,ars->rs',v,linear+masspair.sum(1)+masspair.sum(0))
    check('frozen_vertex_derivative_not_allocation_finite_difference',
          (frozen(1+eps*v)-frozen(1-eps*v))/(2*eps),derivative)
    for trial in range(2):
        kap=rng.uniform(.2,1.7,n)
        _,kt,ku=inv.assemble({key:p[key]*kap[a] for a,key in enumerate(inv.names)},r)
        kte=np.array([q.T@v@q for v in kt])
        ktad=np.einsum('i,rii->r',a0,kte)/LOOP2
        independent=-np.diag(ktad/r)
        for j in range(3):
            for k in range(3):
                independent[j,k]+=cw.hard_trace_hessian(lam,q,kt[j],kt[k],ku[j,k],290,mu2,1.5,1.)
        check(f'random_vertex_rescaling_independent_Frechet_{trial}',frozen(kap),independent)
        independent_z=np.einsum('rij,sij,ij->rs',kte,kte,slope,optimize=True)/LOOP2
        check(f'random_vertex_rescaling_independent_kinetic_{trial}',
              np.einsum('a,b,abrs->rs',kap,kap,kinpair),independent_z)
    sources=[Path(__file__),Path(p1.__file__),Path(cw.__file__),Path(common.__file__),
             Path(common.scalar.__file__),Path(common.scalar.scalar.__file__),
             Path(common.yuk.__file__),Path(ir.__file__),oldpath,irpath]
    return dict(schema='p54-invariant-attribution-v1',date='2026-09-14',
        scope='fixed mass basis/weights, real coupling slice, radial scalar hard one-loop; not pole matching',
        invariant_order=inv.names,reference_parameters=p,reference_radii=r.tolist(),mu=mu,
        radial_metric=metric.tolist(),mass_scalar_plus_CT=hs.tolist(),hard_scalar_kinetic=dz.tolist(),
        hard_scalar_kinetic_eigenvalues=zg.tolist(),kinetic_witness=witness.tolist(),
        contributions=rows,mass_bubble_pairs=masspair.tolist(),kinetic_pairs=kinpair.tolist(),
        exclusions=['xi1 and xi2 have zero radial Hessian jets at H=0, not zero general interactions',
                    'mass parameters affect propagators despite zero cubic/quartic jets',
                    'imaginary coupling directions are outside this real parameter audit',
                    'gauge/ghost/fermion momentum, Wilson/box, lower matching and finite seesaw incomplete'],
        selected_physical_benchmark=None,physical_fit_enabled=False,default_parameters_changed=False,
        cache=dict(unit_hits=inv.hits,new_unit_Hessians=inv.evaluated,new_full_Hessians=store.evaluated),
        source_sha256={str(p.relative_to(RF)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        checks=checks,summary=dict(passed=sum(c['pass_'] for c in checks),total=len(checks),all_pass=all(c['pass_'] for c in checks)))


def markdown(r):
    lines=['# P54 invariant mass / kinetic attribution','',
      'Fixed propagators; signed allocations include interference and are not observables.',
      f"Checks: {r['summary']['passed']}/{r['summary']['total']}.",'',
      '| Invariant | Sigma-axis scalar mass + CT allocation | Leading kinetic allocation | Kinetic self term |',
      '|---|---:|---:|---:|']
    for row in sorted(r['contributions'],key=lambda a:abs(a['leading_kinetic_allocation']),reverse=True):
        lines.append(f"| {row['invariant']} | {row['sigma_mass_allocation']:.8g} | {row['leading_kinetic_allocation']:.8g} | {row['leading_kinetic_self']:.8g} |")
    lines+=['',f"Kinetic generalized eigenvalues: {r['hard_scalar_kinetic_eigenvalues']}.",'',
      'All ordered interference matrices, full 3x3 allocations, actual normalizations, and source hashes are retained in JSON.',
      'A new parameter card must rebuild its spectrum, weights, tadpoles and light-doublet constraint. No physical fit is enabled.']
    return '\n'.join(lines)+'\n'


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--cache-dir',type=Path,required=True)
    args=ap.parse_args();result=run(args.cache_dir)
    OUT.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    OUT.with_suffix('.md').write_text(markdown(result))
    print(markdown(result));print('failed',[c for c in result['checks'] if not c['pass_']])
    if not result['summary']['all_pass']:raise SystemExit(1)
