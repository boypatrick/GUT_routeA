# Route G G2 — joint output mode and detector recoil

A conditional event-rate prediction from the unchanged G2 tree action and preparation. No empirical masses, new parameters, time-dependent background or absolute collision probability enter.

## Preparation and normalization

R=1, M5=0.5, alpha=0.25, MD=5, g=0.05, initial Phi mode n=1 and common initial COM momentum p_i=0.2. The exact finite detector packet has l=-4,...,4, width 1.2, center 0.7. Values use engineering units, not experimental measurements.

$$W_{mj}=|c_{m+j-n}|^2 K_{m+j-n,m},\quad K_{lm}=\sigma_{lm}v_{lm},\quad j-l=n-m.$$

$$P_A(m,j)=W_{mj}/K_A,\quad K_A=\sum_{m,j}W_{mj};\qquad P_C(m,j)=\mathbf1_{m\ne n}W_{mj}/K_C,\quad K_C=\sum_{m\ne n,j}W_{mj}.$$

A conditions on scattering; C conditions on mode-changing scattering. Both use the declared common-overlap dilute rate prescription. Other flux/overlap preparations can change these weights. Scattered elastic events m=n are included only in A; neither table includes the unscattered identity amplitude.

| Quantity | Value |
|---|---:|
| All open signed rows | 25 |
| Mode-changing rows | 16 |
| K_A (mass unit^-2) | 1.85661444734e-06 |
| K_C (mass unit^-2) | 1.63742847959e-06 |
| K_C/K_A | 0.881943196088 |

All open outputs are enumerated, not only the previous m=0 example. Signed j is a theoretical label; neutral detector mass alone identifies only r=|j|.

## All Scattering


### Signed joint P(m,j)

| m / j | -4 | -3 | -2 | -1 | 0 | 1 | 2 | 3 | 4 | 5 | sum |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| -2 | 0 | 0 | 0 | 0.00081339106 | 0.0070464615 | 0 | 0 | 0 | 0 | 0 | 0.0078598526 |
| -1 | 0 | 0 | 0.00094694418 | 0.011045122 | 0.058392772 | 0.12813932 | 0.072669892 | 0 | 0 | 0 | 0.27119405 |
| 0 | 0 | 0.00082438759 | 0.010298995 | 0.06109596 | 0.16801155 | 0.2099302 | 0.11843175 | 0.030464183 | 0.0036308502 | 0.00020141526 | 0.60288929 |
| 1 | 0.00010168495 | 0.0013697947 | 0.0089720786 | 0.028038794 | 0.041092099 | 0.028038794 | 0.0089720786 | 0.0013697947 | 0.00010168495 | 0 | 0.1180568 |
| sum | 0.00010168495 | 0.0021941823 | 0.020218018 | 0.10099327 | 0.27454288 | 0.36610832 | 0.20007372 | 0.031833978 | 0.0037325352 | 0.00020141526 | 1 |

### Mass-readout joint P(m,|j|)

| m / r | 0 | 1 | 2 | 3 | 4 | 5 | sum |
|---:|---:|---:|---:|---:|---:|---:|---:|
| -2 | 0.0070464615 | 0.00081339106 | 0 | 0 | 0 | 0 | 0.0078598526 |
| -1 | 0.058392772 | 0.13918445 | 0.073616837 | 0 | 0 | 0 | 0.27119405 |
| 0 | 0.16801155 | 0.27102616 | 0.12873074 | 0.031288571 | 0.0036308502 | 0.00020141526 | 0.60288929 |
| 1 | 0.041092099 | 0.056077588 | 0.017944157 | 0.0027395894 | 0.0002033699 | 0 | 0.1180568 |
| sum | 0.27454288 | 0.46710159 | 0.22029174 | 0.03402816 | 0.0038342201 | 0.00020141526 | 1 |

| Event moment | Value |
|---|---:|
| mean m | -0.168856956059 |
| mean j | 0.729276316469 |
| mean l | -0.43958063959 |
| mean delta_j | 1.16885695606 |
| m_and_signed_j, bits | 0.0951708986407 |
| m_and_abs_j, bits | 0.0393651383423 |
| loss_from_sign_coarse_graining, bits | 0.0558057602984 |

Mutual information quantifies correlation in this conditional event ensemble, not confidence in the theory. JSON provides all marginals, Bayes conditionals and covariance matrices.

| Initial l | Incident weight | Weight conditioned on this ensemble |
|---:|---:|---:|
| -4 | 0.001285380886 | 0.002686407779 |
| -3 | 0.01460860355 | 0.02976037369 |
| -2 | 0.08290718676 | 0.1284608112 |
| -1 | 0.2349536867 | 0.3241896659 |
| 0 | 0.3324902842 | 0.3236921927 |
| 1 | 0.2349536867 | 0.1464705419 |
| 2 | 0.08290718676 | 0.03943626163 |
| 3 | 0.01460860355 | 0.00500064496 |
| 4 | 0.001285380886 | 0.0003031002148 |

## Mode Conversion


### Signed joint P(m,j)

| m / j | -4 | -3 | -2 | -1 | 0 | 1 | 2 | 3 | 4 | 5 | sum |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| -2 | 0 | 0 | 0 | 0.00092227148 | 0.0079896999 | 0 | 0 | 0 | 0 | 0 | 0.0089119714 |
| -1 | 0 | 0 | 0.001073702 | 0.012523621 | 0.066209221 | 0.14529204 | 0.082397475 | 0 | 0 | 0 | 0.30749606 |
| 0 | 0 | 0.00093474 | 0.011677617 | 0.069274258 | 0.19050155 | 0.23803143 | 0.13428501 | 0.034542115 | 0.0041168754 | 0.00022837668 | 0.68359197 |
| 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| sum | 0 | 0.00093474 | 0.012751319 | 0.08272015 | 0.26470047 | 0.38332347 | 0.21668248 | 0.034542115 | 0.0041168754 | 0.00022837668 | 1 |

### Mass-readout joint P(m,|j|)

| m / r | 0 | 1 | 2 | 3 | 4 | 5 | sum |
|---:|---:|---:|---:|---:|---:|---:|---:|
| -2 | 0.0079896999 | 0.00092227148 | 0 | 0 | 0 | 0 | 0.0089119714 |
| -1 | 0.066209221 | 0.15781566 | 0.083471177 | 0 | 0 | 0 | 0.30749606 |
| 0 | 0.19050155 | 0.30730569 | 0.14596262 | 0.035476855 | 0.0041168754 | 0.00022837668 | 0.68359197 |
| 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| sum | 0.26470047 | 0.46604362 | 0.2294338 | 0.035476855 | 0.0041168754 | 0.00022837668 | 1 |

| Event moment | Value |
|---|---:|
| mean m | -0.325319999342 |
| mean j | 0.826897151318 |
| mean l | -0.498422848024 |
| mean delta_j | 1.32531999934 |
| m_and_signed_j, bits | 0.058136448339 |
| m_and_abs_j, bits | 0.0394287084725 |
| loss_from_sign_coarse_graining, bits | 0.0187077398665 |

Mutual information quantifies correlation in this conditional event ensemble, not confidence in the theory. JSON provides all marginals, Bayes conditionals and covariance matrices.

| Initial l | Incident weight | Weight conditioned on this ensemble |
|---:|---:|---:|
| -4 | 0.001285380886 | 0.002930713496 |
| -3 | 0.01460860355 | 0.03219093826 |
| -2 | 0.08290718676 | 0.1354834791 |
| -1 | 0.2349536867 | 0.3357935897 |
| 0 | 0.3324902842 | 0.3204289062 |
| 1 | 0.2349536867 | 0.1342850064 |
| 2 | 0.08290718676 | 0.03454211468 |
| 3 | 0.01460860355 | 0.004116875393 |
| 4 | 0.001285380886 | 0.0002283766847 |

## Conservation, event selection and the differential prediction

$$\langle j\rangle_E=\langle l\rangle_E+n-\langle m\rangle_E,\quad \mathrm{Var}_E(j)=\mathrm{Var}_E(l)+\mathrm{Var}_E(m)-2\mathrm{Cov}_E(l,m).$$

$$\mathrm{Cov}_E(m,j)=\mathrm{Cov}_E(m,l)-\mathrm{Var}_E(m),\quad j-l=n-m.$$

Both sides must use the same selected-event ensemble E. The incident mean l is zero, but selection by channel rates biases l. Raw m-versus-j covariance is therefore not the transfer law. The transfer j-l and m have perfect anticorrelation whenever their variance is nonzero.

$$dP_E(m,j,\Omega)=P_E(m,j)\frac{d\Omega}{4\pi},\quad\mathbf p_\Phi=p_{f,lm}\hat u,\quad\mathbf p_X=-p_{f,lm}\hat u.$$

Each JSON event records p_f, both outgoing energies, incoming sqrt(s_l), incident weight, channel rate and both ensemble fractions. The contact amplitude is isotropic. Mean three-momenta vanish, each covariance is p_f^2 times the unit matrix / 3, and the cross-covariance is its negative. The discrete and angular laws are normalized separately.

## Mass-readout ambiguity and ideal kinematic sign inference

Neutral X has mu_j=mu_-j. A mass measurement projects onto the full orthogonal degenerate subspace and sums probabilities; it does not implement a coherent sign-mixing readout. Thus measuring m and detector mass does not generally identify the sign. Additional ordinary COM momentum information can distinguish conversion branches.

$$S_f=\sqrt{m_m^2+p_f^2}+\sqrt{M_D^2+r^2/R^2+p_f^2},\quad L^2=R^2\left[\left(S_f-\sqrt{m_n^2+p_i^2}\right)^2-M_D^2-p_i^2\right].$$

$$d=m-n\ne0:\qquad j=\frac{L^2-d^2-r^2}{2d},\quad l=j+d.$$

The final identity follows from l=j+d and j^2=r^2. It reconstructs each conversion row at ideal precision. This is model-dependent inference using the known common incoming COM momentum and shared spectrum. It uses compact conservation, so reconstruction is not an independent experimental verification of that premise. The predicted outgoing momentum lines, conditional fractions and isotropy are testable. Finite resolution or unknown incoming momentum can restore ambiguity; elastic d=0 remains sign-degenerate except r=0.

### Concrete example: m=0 and |j|=1

| Incoming l | Outgoing j | Conditional sign weight given (m=0,r=1) | Outgoing p_f |
|---:|---:|---:|---:|
| -2 | -1 | 0.225424585912 | 1.36191317264 |
| 0 | 1 | 0.774575414088 | 1.02175457609 |

The lines differ by **0.340158596547** in the common mass unit. Sign weights agree when conditioned from A or C, because m=0 is already mode-changing. These are ideal lines, not detector-resolution or event-count forecasts.

### All ambiguous mass-readout pairs

| m | r | Possible j | p_f separation | Ideal-p_f result |
|---:|---:|---|---:|---|
| -1 | 1 | -1, 1 | 0.698653755781 | resolved for declared conversion |
| -1 | 2 | -2, 2 | 1.50881326735 | resolved for declared conversion |
| 0 | 1 | -1, 1 | 0.340158596547 | resolved for declared conversion |
| 0 | 2 | -2, 2 | 0.658957636517 | resolved for declared conversion |
| 0 | 3 | -3, 3 | 0.942365815039 | resolved for declared conversion |
| 1 | 1 | -1, 1 | 0 | elastic sign unresolved |
| 1 | 2 | -2, 2 | 0 | elastic sign unresolved |
| 1 | 3 | -3, 3 | 0 | elastic sign unresolved |
| 1 | 4 | -4, 4 | 0 | elastic sign unresolved |

Among the **5 nonelastic sign pairs**, the smallest momentum separation is **0.340158596547**. Elastic sign pairs have zero separation and are excluded from that minimum.

## Limits and next physical choice

Translation, arbitrary packet phases and dephasing preserve every inclusive joint rate here: each signed output fixes its incoming l. These observables test the spectrum/interaction but do not certify coherent localization or fix off-diagonal recoil coherence. A coherence-sensitive readout would be a separate experiment.

Next choose finite-resolution discrimination or an absolute encounter probability. Discrimination needs a declared detector response and incoming momentum spread. Absolute probabilities additionally need ordinary-space collision packets or luminosity/overlap and duration. Neither is silently supplied. A powered time-dependent background remains unused.

## Verification

**357/357 checks passed.** This validates the implementation, not the empirical truth of this model.

Checks cover the frozen G2 totals/hash, joint normalization and Bayes conditionals, coarse readout, event conservation, on-shell/back-to-back recoil, isotropic angular moments, selected mean/covariance identities, full-joint large-gauge relabelling, phase/translation/dephasing invariance, ideal model-based reconstruction and information loss under |j| coarse graining.

Reproduce: python3 -B route_g/code/verify_g2_recoil_joint.py. JSON keeps every event, full-precision results, checks and source hashes. This imports frozen G2 functions without executing or rewriting their main.
