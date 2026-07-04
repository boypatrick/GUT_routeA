# DYN-0: Dynamics Conventions Addendum and Inventory

`13/13` checks passed.

- Conventions addendum recorded: superpotential parameters, PS
  singlet-vev registry (mirrored from HEAD), Aulakh-Girdhar branch
  convention, adopted decisions (root code/ + whitelisted output
  lanes; M_S benchmark 3 TeV, DYN-2 fit window [1, 10] TeV).
- Vacuum-cubic special points verified exactly: SU(5) x=1/2 (xi=-5),
  SU(5) x=-1 (xi=10), G_LR x=0 (xi=3), flipped SU(5) x=1/3 (xi=-2/3).
- Inventory: 38 history assets across the DYN items, all present in HEAD (read-only contract).
- Rerun-from-history gates: the four audit-4a1 scripts, extracted from
  HEAD into a temporary mini-tree, regenerate ledgers matching HEAD up
  to timestamps/sha stamps/source manifests.

Boundary: no dynamics results are claimed.  The heavy spectrum is
still a placeholder and the scalar-Hessian Goldstone eigenvector audit
is still pending; both are DYN-1 targets.
