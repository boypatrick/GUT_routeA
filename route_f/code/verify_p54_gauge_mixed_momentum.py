#!/usr/bin/env python3
"""Same-action radial background gauge/scalar-vector momentum kernel.

DR/MSbar, F^a=partial.a^a+xi*g*(T^a X).eta, Euclidean s=p^2.
All tensor numerators are reduced BEFORE epsilon -> 0. No coupling scan.
This computes the radial one-loop functional, not an all-orders vacuum,
Nielsen-closed physical pole or complete finite matching certificate.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import sys
from functools import lru_cache
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.linalg import eigh
from numpy.polynomial.legendre import leggauss
import verify_p54_common_renormalization as common
import verify_p54_goldstone_ir as ir
import verify_p54_radial_fermion_momentum as fermion

RF=Path(__file__).resolve().parents[1]
OUT=RF/'output/p54_gauge_mixed_momentum'
LOOP=16*np.pi**2
mp.mp.dps=45


def aa(x,mu2):
    x,mu2=mp.mpf(x),mp.mpf(mu2)
    return x*(mp.log(x/mu2)-1) if x else mp.mpf(0)


@lru_cache(maxsize=8192)
def bb(a,b,s,mu2):
    """Finite Euclidean bubble, with +1/epsilon UV residue; no 1/(16pi²)."""
    a,b,s,mu2=map(mp.mpf,(a,b,s,mu2))
    if min(a,b,s)<0 or mu2<=0:raise ValueError('Euclidean nonnegative masses and momentum required')
    if not a and not b:
        if not s:raise ValueError('Massless zero-momentum bubble is not a local finite kernel')
        return 2-mp.log(s/mu2)
    if not a or not b:
        m=max(a,b)
        if not s:return 1-mp.log(m/mu2)
        return 2-mp.log(m/mu2)-(1+m/s)*mp.log1p(s/m)
    if not s:
        return -mp.log(a/mu2) if a==b else (aa(b,mu2)-aa(a,mu2))/(a-b)
    return -mp.quad(lambda x:mp.log(((1-x)*a+x*b+x*(1-x)*s)/mu2),[0,.5,1])


def dot_numerator(a,b,s,mu2):
    """4*int (k.q)^2/[(k²+a)(q²+b)], q=k+p, finite part."""
    a,b,s=map(mp.mpf,(a,b,s));c=a+b+s
    value=-(3*a+b+s)*aa(a,mu2)-(a+3*b+s)*aa(b,mu2)
    if c:value+=c*c*bb(a,b,s,mu2)
    return value


def pk_numerator(a,b,s,mu2):
    """int (p.k)^2/[(k²+a)(q²+b)], finite part."""
    a,b,s=map(mp.mpf,(a,b,s))
    if not s:return mp.mpf(0)
    c=a-b-s
    return (c*aa(a,mu2)+(b-a+3*s)*aa(b,mu2)+c*c*bb(a,b,s,mu2))/4


@lru_cache(maxsize=8192)
def vector_bubble(a,b,s,mu2,xi=0.):
    """Finite int tr[D_a(k)D_b(k+p)]; UV pole 3+xi².

    D=delta/(k²+a)+kk/a[1/(k²+a)-1/(k²+xi*a)].
    The -2 rational term comes from (d-2)B, NOT a fitted CW offset.
    """
    a,b,s,xi=map(mp.mpf,(a,b,s,xi))
    if min(a,b)<=0 or xi<0:raise ValueError('Active radial vector vertices require positive masses')
    if not s:
        return 3*bb(a,b,0,mu2)-2+(xi*xi*bb(xi*a,xi*b,0,mu2) if xi else 0)
    c,d=xi*a,xi*b
    value=2*bb(a,b,s,mu2)-2
    if xi:value+=xi*(bb(a,d,s,mu2)+bb(c,b,s,mu2))
    value+=(dot_numerator(a,b,s,mu2)-dot_numerator(a,d,s,mu2)
            -dot_numerator(c,b,s,mu2)+dot_numerator(c,d,s,mu2))/(4*a*b)
    return value


def vector_seagull(a,mu2,xi=0.):
    a,xi=map(mp.mpf,(a,xi))
    return 3*aa(a,mu2)+2*a+(xi*aa(xi*a,mu2) if xi else 0)


@lru_cache(maxsize=8192)
def mixed_bubble(a,b,s,mu2,xi=0.):
    """Finite int p.D_a(k).p/[(k+p)²+b]; UV pole s*(3+xi)/4."""
    a,b,s,xi=map(mp.mpf,(a,b,s,xi))
    if a<=0 or b<0 or s<0 or xi<0:raise ValueError('Invalid mixed kernel arguments')
    if not s:return mp.mpf(0)
    return s*bb(a,b,s,mu2)+(pk_numerator(a,b,s,mu2)-pk_numerator(xi*a,b,s,mu2))/a


def mixed_slope(a,b,mu2,xi=0.):
    """Finite s derivative at zero. Retains epsilon from 1/d."""
    a,b,xi=map(mp.mpf,(a,b,xi))
    value=mp.mpf(3)/4*bb(a,b,0,mu2)-mp.mpf(1)/8
    if xi:value+=xi*bb(xi*a,b,0,mu2)/4+xi/8
    return value


@lru_cache(maxsize=128)
def vector_slope(a,b):
    """UV-finite Landau derivative, independent of tensor reduction / mu.

    Angularly expanding tr(P_T(k)P_T(k+p)) in d=4 gives this convergent
    rational integral. The result is symmetric in a,b by integration by parts.
    """
    a,b=map(mp.mpf,(a,b))
    return -mp.quad(lambda t:3*b*t/((t+a)*(t+b)**3)+mp.mpf(3)/(4*(t+a)*(t+b)),
                    [0,min(a,b),max(a,b),mp.inf])


def independent_massless_mixed_difference(a,s,xi):
    """Exact S^3 angular integration, then UV-convergent radial quadrature.

    <1/q²>=1/max(t,s), <z²/q²>=(t+s)/(4 max(t,s)²).
    These shell-theorem identities avoid the q=0 angular endpoint cusp.
    """
    a,s,xi=map(mp.mpf,(a,s,xi))
    def integrand(t):
        h=max(t,s);zero=1/h;second=(t+s)/(4*h*h)
        value=(t*(zero-second)-mp.mpf(3)/4)/(t+a)
        if xi:value+=xi*(t*second-mp.mpf(1)/4)/(t+xi*a)
        return value
    return mp.quad(integrand,[0,s,max(a,s),mp.inf])


def independent_four_dimensional_difference(a,b,s,kind,xi=0.,n=240):
    """UV-convergent momentum difference, direct radial/angular integral.

    Returns 16pi² times the integral. For VV subtract s=0; for SV
    subtract the zero-momentum slope from I(s)/s. Rational MS constants
    cancel in these differences, and are independently checked in d dims.
    """
    nodes,weights=leggauss(n)
    scale=max(a,b,s,.01)
    # Resolve the massless q=0 integrable cusp (t=s,z=-1) explicitly.
    split=[0.,s/(s+scale),1.] if b==0 else [0.,1.]
    y=np.concatenate([(hi+lo)/2+(hi-lo)*nodes/2 for lo,hi in zip(split,split[1:])])
    wy=np.concatenate([(hi-lo)*weights/2 for lo,hi in zip(split,split[1:])])
    t=scale*y/(1-y);dt=scale/(1-y)**2
    theta=np.pi*np.arange(1,n+1)/(n+1)
    z=np.cos(theta);wz=2*np.sin(theta)**2/(n+1)
    tt=t[:,None];zz=z[None,:];q2=tt+s+2*np.sqrt(tt*s)*zz
    if kind=='SV':
        first=(1-zz*zz)/(tt+a)
        if xi:first+=xi*zz*zz/(tt+xi*a)
        integrand=first*(1/(q2+b)-1/(tt+b))
    elif kind=='VV':
        cos2=(tt+np.sqrt(tt*s)*zz)**2/(tt*q2)
        # TT, TL, LT, LL projector contractions at d=4.
        integrand=(2+cos2)/((tt+a)*(q2+b))
        if xi:
            integrand+=xi*(1-cos2)/((tt+a)*(q2+xi*b))
            integrand+=xi*(1-cos2)/((tt+xi*a)*(q2+b))
            integrand+=xi*xi*cos2/((tt+xi*a)*(q2+xi*b))
        integrand-=3/((tt+a)*(tt+b))
        if xi:integrand-=xi*xi/((tt+xi*a)*(tt+xi*b))
    else:raise ValueError(kind)
    return float(np.dot(wy*dt*t,integrand@wz))


def groups(values,zero_count):
    values=np.array(values);values[:zero_count]=0
    result=[];masses=[];start=0
    while start<len(values):
        stop=start+1
        while stop<len(values) and abs(values[stop]-values[start])<2e-11*max(1.,abs(values[start])):stop+=1
        result.append(np.arange(start,stop));masses.append(float(values[start:stop].mean()));start=stop
    return result,np.array(masses)


def run(cache_dir):
    oldpath=RF/'output/p54_full_doublet_cw.json';renpath=RF/'output/p54_common_renormalization.json'
    old=json.loads(oldpath.read_text());ren=json.loads(renpath.read_text())
    irpath=RF/'output/p54_goldstone_ir.json';old_ir=json.loads(irpath.read_text())
    p1=common.yuk.module('gm_p1',RF/'code/verify_p54_p1_hessian_spectrum.py')
    cw=common.yuk.module('gm_cw',RF/'code/verify_p54_p2_bosonic_cw.py')
    pars=old['tree_parameters'];r=np.array([old['vacuum'][k] for k in ('omega','sigma','vs')])
    mu=old['scheme']['mu_over_omega'];mu2=mu*mu;g=mu
    rad=np.column_stack([p1.vacuum_vector(*e) for e in np.eye(3)]);metric=rad.T@rad;x=rad@r
    store=common.scalar.ActionJets(p1,pars,cache_dir);h=store.hessian(x,'reference')
    t=[];u=np.zeros((3,3,328,328))
    for i in range(3):
        ti,ui=store.jets(x,h,rad[:,i],f'radial_{i}');t.append(ti);u[i,i]=ui
    t=np.array(t)
    for i in range(3):
        for j in range(i+1,3):
            step=.02;hp=store.hessian(x+step*(rad[:,i]+rad[:,j]),f'pair_{i}{j}')
            u[i,j]=(hp-h-step*(t[i]+t[j])-.5*step**2*(u[i,i]+u[j,j]))/step**2;u[j,i]=u[i,j]
    gens=cw.generators();v=cw.field_orbit(p1,x,gens)
    d=np.array([cw.field_orbit(p1,q,gens) for q in rad.T])
    mv=g*g*v.T@v;vl,vu=np.linalg.eigh(mv);vl[:12]=0.
    bg=np.array([g*g*(dr.T@v+v.T@dr) for dr in d])
    cg=np.array([[g*g*(dr.T@ds+ds.T@dr) for ds in d] for dr in d])
    bg=np.array([vu.T@z@vu for z in bg]);cg=np.array([[vu.T@z@vu for z in row] for row in cg])
    vg,vm=groups(vl,12)
    checks=[]
    def check(name,a,b=0.,tol=3e-8):
        error=common.yuk.error(np.asarray(a),np.asarray(b))
        checks.append(dict(name=name,residual=float(error),tolerance=tol,passed=bool(error<tol)))
    check('zero_new_action_Hessian_evaluations',store.evaluated)
    check('radial_vector_vertices_annihilate_unbroken_subspace',bg[:,:12,:])
    check('radial_gauge_mass_second_vertices_annihilate_unbroken_subspace',cg[:,:,:12,:])
    rotated_d=np.einsum('ria,ab->rib',d,vu)
    check('mixed_vertices_annihilate_unbroken_subspace',rotated_d[:,:,:12])
    check('tree_gauge_orbit_Ward_identity',h@v)
    nv=len(vm);vw=np.zeros((3,3,nv,nv));vc=np.zeros((3,3,nv))
    for i,gi in enumerate(vg):
        vc[:,:,i]=np.einsum('rsii->rs',cg[:,:,gi[:,None],gi])
        for j,gj in enumerate(vg):
            block=bg[:,gi[:,None],gj]
            vw[:,:,i,j]=np.einsum('rij,sij->rs',block,block)

    def vector(xi,s):
        bubble=np.zeros((3,3));seagull=np.zeros((3,3));ghost=np.zeros((3,3));tad=np.zeros(3)
        for i,a in enumerate(vm):
            if not a:continue
            av=float(vector_seagull(a,mu2,xi));seagull+=vc[:,:,i]*av/(2*LOOP)
            trb=np.trace(bg[:,vg[i][:,None],vg[i]],axis1=1,axis2=2)
            tad+=trb*av/(2*LOOP)
            if xi:
                ghost-=vc[:,:,i]*xi*float(aa(xi*a,mu2))/LOOP
                tad-=trb*xi*float(aa(xi*a,mu2))/LOOP
            for j,b in enumerate(vm):
                if not b or np.linalg.norm(vw[:,:,i,j])<1e-25:continue
                bubble-=vw[:,:,i,j]*float(vector_bubble(a,b,s,mu2,xi))/(2*LOOP)
                if xi:ghost+=vw[:,:,i,j]*xi*xi*float(bb(xi*a,xi*b,s,mu2))/LOOP
        return dict(seagull=seagull,bubble=bubble,ghost=ghost,tadpole=tad,total=seagull+bubble+ghost)

    gauge_mass=g*g*v@v.T
    gauge_first=np.array([g*g*(dr@v.T+v@dr.T) for dr in d])
    gauge_second=np.array([[g*g*(dr@ds.T+ds@dr.T) for ds in d] for dr in d])
    def prepare_scalar(xi):
        hs=h+xi*gauge_mass;ts=t+xi*gauge_first;us=u+xi*gauge_second
        lam,q=np.linalg.eigh(hs);soft=38 if xi==0 else 5
        check(f'xi_{xi}_scalar_soft_count',lam[:soft])
        check(f'xi_{xi}_positive_complement',float(lam[soft]<=0))
        lam[:soft]=0.;sg,sm=groups(lam,soft)
        te=np.array([q.T@z@q for z in ts]);nw=len(sm)
        sw=np.zeros((3,3,nw,nw));sc=np.zeros((3,3));st=np.zeros(3)
        a0=np.array([float(aa(z,mu2)) for z in lam])
        st=np.einsum('i,rii->r',a0,te)/(2*LOOP)
        for i in range(3):
            for j in range(3):sc[i,j]=np.dot(a0,np.diag(q.T@us[i,j]@q))/(2*LOOP)
        for k,gk in enumerate(sg):
            for l,gl in enumerate(sg):
                block=te[:,gk[:,None],gl]
                sw[:,:,k,l]=np.einsum('rij,sij->rs',block,block)/(2*LOOP)
        # Vertex contractions use both the scalar and vector eigenbases.
        ds=np.einsum('ik,rka->ria',q.T,rotated_d)
        mixed=np.zeros((3,3,nv,nw))
        for k,gk in enumerate(vg):
            for l,gl in enumerate(sg):
                block=ds[:,gl[:,None],gk]
                mixed[:,:,k,l]=np.einsum('ria,sia->rs',block,block)
        return dict(masses=sm,weights=sw,seagull=sc,tadpole=st,mixed=mixed,
                    scalar_eigenvalues=lam,scalar_eigenvectors=q,first_vertices=te)

    def scalar(prep,s):
        sm,sw=prep['masses'],prep['weights'];kernel=prep['seagull'].copy()
        for i,a in enumerate(sm):
            for j,b in enumerate(sm):
                if np.linalg.norm(sw[:,:,i,j])<1e-25:continue
                if s==0 and a==0 and b==0:continue # only for explicit hard CW check
                # The same high-precision master also resolves tiny finite-xi
                # Goldstone masses without a fixed-node endpoint error.
                kernel-=sw[:,:,i,j]*float(bb(a,b,s,mu2))
        return kernel

    def mixed(prep,xi,s):
        result=np.zeros((3,3))
        for i,a in enumerate(vm):
            if not a:continue
            for j,b in enumerate(prep['masses']):
                weight=prep['mixed'][:,:,i,j]
                if np.linalg.norm(weight)<1e-24:continue
                result-=4*g*g*weight*float(mixed_bubble(a,b,s,mu2,xi))/LOOP
        return result

    # Independent master checks: Feynman gauge, s=0, dimensional constants
    # and direct UV-convergent four-dimensional angular integrations.
    for i,(a,b,s) in enumerate(((.2,.31,.07),(.013,.14,.08),(.03,0.,.015))):
        for xi in (0.,.1,1.):
            got=float(mixed_bubble(a,b,s,mu2,xi)/s-mixed_slope(a,b,mu2,xi))
            direct=(float(independent_massless_mixed_difference(a,s,xi)) if not b else
                    independent_four_dimensional_difference(a,b,s,'SV',xi,320))
            check(f'independent_mixed_integral_{i}_{xi}',got,direct,tol=3e-6)
            if not b:
                lo=independent_four_dimensional_difference(a,b,s,'SV',xi,320)
                hi=independent_four_dimensional_difference(a,b,s,'SV',xi,640)
                check(f'massless_angular_endpoint_convergence_{xi}',float(abs(hi-direct)>abs(lo-direct)+1e-10))
                check(f'massless_angular_endpoint_integral_{xi}',hi,direct,tol=3e-6)
            small=mp.mpf('1e-8')*min(a,b if b else a)
            check(f'mixed_small_momentum_dimensional_slope_{i}_{xi}',
                  float(mixed_bubble(a,b,small,mu2,xi)/small),float(mixed_slope(a,b,mu2,xi)),tol=2e-7)
        check(f'Feynman_gauge_mixed_master_{i}',float(mixed_bubble(a,b,s,mu2,1)),float(s*bb(a,b,s,mu2)))
        if b:
            for xi in (0.,.1,1.):
                got=float(vector_bubble(a,b,s,mu2,xi)-vector_bubble(a,b,0.,mu2,xi))
                direct=independent_four_dimensional_difference(a,b,s,'VV',xi,320)
                check(f'independent_vector_integral_{i}_{xi}',got,direct,tol=3e-6)
                check(f'vector_bubble_exchange_symmetry_{i}_{xi}',float(vector_bubble(a,b,s,mu2,xi)),float(vector_bubble(b,a,s,mu2,xi)))
            check(f'Feynman_gauge_vector_master_{i}',float(vector_bubble(a,b,s,mu2,1)),float(4*bb(a,b,s,mu2)-2))
    a=mp.mpf('.17')
    for eps in (mp.mpf('1e-6'),mp.mpf('2e-7')):
        bare_b=mp.exp(mp.euler*eps)*mp.gamma(eps)*(a/mu2)**(-eps)
        bare_a=mp.exp(mp.euler*eps)*mp.gamma(-1+eps)*a*(a/mu2)**(-eps)
        check('dimensional_vector_seagull_rational_'+str(eps),float((3-2*eps)*bare_a+3*a/eps),float(3*aa(a,mu2)+2*a),tol=2e-6)
        check('dimensional_vector_bubble_rational_'+str(eps),float((3-2*eps)*bare_b-3/eps),float(3*bb(a,a,0,mu2)-2),tol=4e-6)
        check('dimensional_mixed_slope_rational_'+str(eps),float((1-1/(4-2*eps))*bare_b-mp.mpf(3)/(4*eps)),float(mixed_slope(a,a,mu2)),tol=2e-6)

    landau=prepare_scalar(0.);v0=vector(0.,0.)
    check('zero_momentum_vector_CW_same_MS_constant',v0['total'],ren['radial_bosonic']['vector_loop_Hessian'])
    check('zero_momentum_vector_tadpole_same_scheme',v0['tadpole'],ren['broken_tadpoles']['vector'])
    check('zero_momentum_scalar_hard_CW_same_scheme',scalar(landau,0.),ren['radial_bosonic']['scalar_loop_Hessian'])
    check('zero_momentum_scalar_tadpole_same_scheme',landau['tadpole'],ren['broken_tadpoles']['scalar'])
    ct=-np.diag((landau['tadpole']+v0['tadpole'])/r)
    tree=rad.T@h@rad
    # Verify the original fixed-VEV invariant mass CT, not a new subtraction.
    check('original_common_mass_CT_unchanged',ct,ren['radial_bosonic']['CT_Hessian'])
    check('original_tree_radial_Hessian_unchanged',tree,ren['radial_bosonic']['tree_Hessian'])
    print('model groups',dict(vector_masses=vm.tolist(),scalar_group_count=len(landau['masses'])),flush=True)
    momenta=[]
    for s in (.001,.01,.05,.1,.5,1.):
        vs=vector(0.,s);sv=mixed(landau,0.,s);sc=scalar(landau,s)
        kernel=s*metric+tree+ct+sc+vs['total']+sv
        previous=s*metric+tree+ct+sc+v0['total']
        old_subset=next(row['radial_Euclidean_subset'] for row in old_ir['scalar_momentum_probes'] if row['pE2']==s)
        check(f'previous_scalar_plus_static_vector_recovered_{s}',previous,old_subset)
        momenta.append(dict(pE2=s,scalar=sc.tolist(),vector_seagull=vs['seagull'].tolist(),
            vector_bubble=vs['bubble'].tolist(),mixed=sv.tolist(),ghost=vs['ghost'].tolist(),
            gauge_momentum_change=(vs['total']-v0['total']+sv).tolist(),
            complete_radial_bosonic_kernel=kernel.tolist(),
            generalized_bosonic_eigenvalues=eigh(kernel,metric,eigvals_only=True).tolist(),
            scalar_plus_static_vector_eigenvalues=eigh(previous,metric,eigvals_only=True).tolist(),
            fermion_addition='sigma-sigma entry: fixed_vev_majorana_kernel(M_i,sigma,mu,s); masses not fitted',
            is_Nielsen_certified_physical_pole=False))
        check(f'Landau_total_kernel_symmetric_{s}',kernel,kernel.T)
        print('s',s,'bosonic eigs',momenta[-1]['generalized_bosonic_eigenvalues'],flush=True)
    # Finite xi is a regulator/implementation check at nonzero s, NOT a
    # gauge-independent pole comparison. Hold the original finite CT fixed;
    # record, rather than cancel, the finite-xi tadpole change.
    xi_rows=[];s_test=.05
    baseline=np.asarray(momenta[2]['complete_radial_bosonic_kernel'])
    for xi in (.1,.01,.001,.0001,.00001,.000001):
        sp=prepare_scalar(xi);vv=vector(xi,s_test);sv=mixed(sp,xi,s_test)
        sc=scalar(sp,s_test);kernel=s_test*metric+tree+ct+sc+vv['total']+sv
        xi_rows.append(dict(xi=xi,pE2=s_test,scalar=sc.tolist(),vector=vv['total'].tolist(),
            ghost=vv['ghost'].tolist(),mixed=sv.tolist(),
            common_CT_held_fixed=True,finite_xi_tadpole_residual=(sp['tadpole']+vv['tadpole']+np.diag(ct)*r).tolist(),
            kernel=kernel.tolist(),landau_difference_norm=float(np.linalg.norm(kernel-baseline)),
            difference_over_xi_logxi=float(np.linalg.norm(kernel-baseline)/(xi*abs(np.log(xi)))),
            ghost_norm=float(np.linalg.norm(vv['ghost'])),
            full_fixed_bare_Nielsen_identity_checked=False))
        print('xi limit',xi,xi_rows[-1]['landau_difference_norm'],flush=True)
    check('finite_xi_to_Landau_contraction',float(any(b['landau_difference_norm']>=a['landau_difference_norm'] for a,b in zip(xi_rows,xi_rows[1:]))))
    check('direct_scalar_ghost_Landau_limit_contraction',float(any(b['ghost_norm']>=a['ghost_norm'] for a,b in zip(xi_rows,xi_rows[1:]))))
    check('finite_xi_endpoint_near_Landau',xi_rows[-1]['landau_difference_norm'],tol=1.2e-4)
    # Local gauge derivative: vector bubble finite, mixed wave UV pole
    # removed minimally. No arbitrary finite wave-function subtraction.
    dzv=np.zeros((3,3));dzsv=np.zeros((3,3))
    for i,a in enumerate(vm):
        if not a:continue
        for j,b in enumerate(vm):
            if b and np.linalg.norm(vw[:,:,i,j])>1e-25:
                step=mp.mpf('1e-8')*min(a,b)
                zero=vector_bubble(a,b,0,mu2)
                slope1=(vector_bubble(a,b,step,mu2)-zero)/step
                slope2=(vector_bubble(a,b,step/2,mu2)-zero)/(step/2)
                check(f'vector_local_slope_convergence_{i}_{j}',float(slope1),float(slope2),tol=2e-7)
                slope=vector_slope(a,b)
                check(f'vector_local_slope_independent_integral_{i}_{j}',float(2*slope2-slope1),float(slope),tol=2e-7)
                check(f'vector_local_slope_exchange_symmetry_{i}_{j}',float(slope),float(vector_slope(b,a)))
                dzv-=vw[:,:,i,j]*float(slope)/(2*LOOP)
        for j,b in enumerate(landau['masses']):
            dzsv-=4*g*g*landau['mixed'][:,:,i,j]*float(mixed_slope(a,b,mu2))/LOOP
    gauge_gram=g*g*np.einsum('ria,sia->rs',d,d)
    wave_pole=-3*gauge_gram/LOOP
    hard_scalar_z=np.asarray(old_ir['hard_scalar_kinetic']['matrix'])
    gauge_completed_hard_z=hard_scalar_z+dzv+dzsv
    combined_z_eigenvalues=eigh(gauge_completed_hard_z,metric,eigvals_only=True)
    sources=[Path(__file__),Path(p1.__file__),Path(cw.__file__),Path(common.__file__),
        Path(common.scalar.__file__),Path(common.yuk.__file__),Path(ir.__file__),Path(fermion.__file__),oldpath,renpath,irpath]
    return dict(schema='p54-same-prescription-radial-gauge-momentum-v1',date='2026-09-17',
        runtime=dict(python=sys.version.split()[0],numpy=np.__version__,jax=p1.jax.__version__,
                     mpmath=mp.__version__,master_decimal_precision=mp.mp.dps),
        scheme=dict(subtraction='DR/MSbar',gauge_fixing='F=partial.a+xi*g*(T X).eta',
                    landau=True,external_fields='scalar backgrounds',tadpole='original fixed-VEV invariant mass CT',
                    dimensional_rational_terms_retained=True,mu=mu,g=g),
        vacuum=r.tolist(),radial_metric=metric.tolist(),vector_group_masses2=vm.tolist(),
        vector_vertices_grouped=vw.tolist(),mixed_vertices_grouped=landau['mixed'].tolist(),
        scalar_group_masses2=landau['masses'].tolist(),common_mass_CT=ct.tolist(),
        gauge_local_kinetic=dict(vector=dzv.tolist(),mixed=dzsv.tolist(),sum=(dzv+dzsv).tolist(),
            generalized_eigenvalues=eigh(dzv+dzsv,metric,eigvals_only=True).tolist(),
            loop_wave_UV_residue=wave_pole.tolist(),minimal_wave_CT_residue=(-wave_pole).tolist()),
        hard_scalar_plus_gauge_kinetic=dict(matrix=gauge_completed_hard_z.tolist(),
            generalized_eigenvalues=combined_z_eigenvalues.tolist(),
            omitted_scalar_soft_nonlocal_part=True,unfitted_fermion_part_not_added=True,
            small_insertion_control_passed=bool(max(abs(combined_z_eigenvalues))<.3),
            is_physical_residue=False),
        momenta=momenta,finite_xi_limit=xi_rows,
        complete_radial_Landau_bosonic_one_loop_kernel=all(c['passed'] for c in checks),
        radial_fermion_kernel_available_as_three_mass_function=True,
        finite_xi_limit_is_Nielsen_certificate=False,full_physical_stability_decided=False,
        physical_fit_enabled=False,default_parameters_changed=False,
        missing=['fixed-bare finite-xi Nielsen/BRST pole test and physical background/quantum identification',
                 'perturbative control, resummation/higher-loop remainder and fitted Majorana inputs',
                 'full Wilson/box, lower gauge/ghost/Yukawa matching and finite CHN seesaw'],
        cache=dict(hits=store.hits,new_Hessians=store.evaluated),
        source_sha256={str(p.relative_to(RF)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        checks=checks,summary=dict(passed=sum(c['passed'] for c in checks),total=len(checks),all_pass=all(c['passed'] for c in checks)))


def markdown(report):
    lines=['# Same-prescription radial gauge / mixed momentum kernel','',
       'Strict-Landau one-loop background kernel, not a gauge-independent pole or all-orders stability proof.',
       f"Checks: {report['summary']['passed']}/{report['summary']['total']}.",'',
       '| Euclidean p2 | Bosonic generalized eigenvalues including gauge momentum |',
       '|---:|---|']
    for row in report['momenta']:
        lines.append(f"| {row['pE2']} | "+', '.join(f'{v:.9g}' for v in row['generalized_bosonic_eigenvalues'])+' |')
    lines+=['','The gauge-only local kinetic insertion has generalized eigenvalues '+str(report['gauge_local_kinetic']['generalized_eigenvalues'])+'.',
      'Adding it to the previous hard scalar kinetic insertion gives '+str(report['hard_scalar_plus_gauge_kinetic']['generalized_eigenvalues'])+'. The soft scalar nonlocal part and unfitted fermions are not included in this local diagnostic; it is not a physical residue.',
      '', 'Finite-xi checks hold the original finite mass CT fixed and report residual tadpoles; convergence to Landau is not a Nielsen identity certificate.',
      'The MS dimensional rational terms +2a (vector seagull), -2 (vector bubble), and -1/8 (mixed slope) are retained and independently checked.',
      'The exact radial Majorana function can be added with three specified masses; no fitted values are invented.','',
      '## Still open','']+['- '+s for s in report['missing']]
    return '\n'.join(lines)+'\n'


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--cache-dir',type=Path,required=True)
    args=ap.parse_args();report=run(args.cache_dir)
    OUT.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    OUT.with_suffix('.md').write_text(markdown(report))
    print(markdown(report));print('failed',[c for c in report['checks'] if not c['passed']])
    if not report['summary']['all_pass']:raise SystemExit(1)
