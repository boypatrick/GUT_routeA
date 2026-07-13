# DYN-5 model-validity gate

- arithmetic checks: **9/9 pass**
- scientific status: **invalid_pending_rederivation**
- the historical one-loop promotion is not valid for the displayed quadratic action

## Correct tree matching

- `delta Z_tree = |zeta|/3 = 0.0434773011`
- `Z_N = 1.0434773011`
- unretuned `|zeta|/Z_N = 0.1249973556` (`-4.167%` shift)
- tree / historical putative-loop-per-efold ratio = `52.637890`
- full-6x6 versus Schur Weinberg residual = `7.147e-16`

## Fail-closed blockers

- no cubic or gauge messenger interaction is specified, so the one-loop anomalous dimension does not exist in the displayed model
- the displayed post-B-L R charges allow X L H_u
- a full branch-local action and an operator-complete selection rule are required

## Checks

- [PASS] K_tr^dagger K_tr = I/3: residual=1.923e-16
- [PASS] tree Kahler correction is delta Z_N = |zeta|/3: deltaZ_tree=0.043477301084
- [PASS] tree correction is 16 pi^2/3 times the historical putative loop coefficient: ratio=52.637890139
- [PASS] NN block of the full inverse equals the Schur-complement inverse: max relative residual=4.487e-16
- [PASS] full and Schur-complement Weinberg tensors agree at zero momentum: max relative residual=7.147e-16
- [PASS] consistent canonical normalization leaves Y M^-1 Y^T invariant: relative residual=4.063e-16
- [PASS] displayed action has no cubic messenger tensor: W is quadratic; the claimed |lambda|^2/(16 pi^2) anomalous dimension is undefined
- [PASS] the displayed R charges allow the dangerous X L H_u operator: R(XLH_u)=2=R(W)
- [PASS] one additive Abelian charge cannot allow XN, XX, LNH_u while forbidding XLH_u: q(XLH)=2 q_X=w follows algebraically (also modulo n)
