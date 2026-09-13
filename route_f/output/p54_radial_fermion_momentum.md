# P54 radial Majorana momentum subset

Checks: **56/56**.

The three Majorana masses are input parameters, not fitted values. All five numerical cards below are synthetic regression tests. No physical stability conclusion follows.

Exact fixed-VEV radial kernel: `Pi_F(s)=-sum[(M_i/sigma)^2 (4 M_i^2+s) L(M_i^2,M_i^2;s)]/(16 pi^2)`.

Only the sigma-sigma coordinate is nonzero. Its canonical generalized component has an extra factor 1/2 because G_sigma,sigma=2.

Zero-momentum universal three-Majorana cap: `0.0875112373081`.

The actual gauge mass first/second vertices and mixed scalar-vector derivative Gram are stored in JSON. Their tensors are fixed by the existing canonical representation and coupling; they do not require a flavor fit.

Still missing: finite vector/scalar-vector momentum integrals with dimensional rational terms, explicit ghost/gauge-fixing limit checks, common wave-function counterterms, Nielsen tests, complete field-space kernel and physical-pole analysis.

```json
[
  {
    "synthetic_test_masses": [
      0.02,
      0.06,
      0.11
    ],
    "finite_sigma_kinetic": 0.018816204043282585,
    "generalized_radial_kinetic": 0.009408102021641293,
    "values": [
      {
        "pE2": 0.0,
        "sigma_coordinate_kernel": 0.0008500164353828568
      },
      {
        "pE2": 0.001,
        "sigma_coordinate_kernel": 0.0008687205038736651
      },
      {
        "pE2": 0.01,
        "sigma_coordinate_kernel": 0.0010292982406057052
      },
      {
        "pE2": 0.1,
        "sigma_coordinate_kernel": 0.0022471941630754955
      },
      {
        "pE2": 1.0,
        "sigma_coordinate_kernel": 0.004905768917119621
      }
    ]
  },
  {
    "synthetic_test_masses": [
      0.25,
      0.35,
      0.45
    ],
    "finite_sigma_kinetic": 0.019703421131180517,
    "generalized_radial_kinetic": 0.009851710565590258,
    "values": [
      {
        "pE2": 0.0,
        "sigma_coordinate_kernel": 0.06157173055975894
      },
      {
        "pE2": 0.001,
        "sigma_coordinate_kernel": 0.061591315343973815
      },
      {
        "pE2": 0.01,
        "sigma_coordinate_kernel": 0.06175697405407376
      },
      {
        "pE2": 0.1,
        "sigma_coordinate_kernel": 0.062427910186881215
      },
      {
        "pE2": 1.0,
        "sigma_coordinate_kernel": 0.00270495316587264
      }
    ]
  },
  {
    "synthetic_test_masses": [
      0.6,
      0.8,
      1.0
    ],
    "finite_sigma_kinetic": -1.1795293836356397,
    "generalized_radial_kinetic": -0.5897646918178199,
    "values": [
      {
        "pE2": 0.0,
        "sigma_coordinate_kernel": -2.303834911476482
      },
      {
        "pE2": 0.001,
        "sigma_coordinate_kernel": -2.305014559564093
      },
      {
        "pE2": 0.01,
        "sigma_coordinate_kernel": -2.315642062159423
      },
      {
        "pE2": 0.1,
        "sigma_coordinate_kernel": -2.4229602850477656
      },
      {
        "pE2": 1.0,
        "sigma_coordinate_kernel": -3.589675080214757
      }
    ]
  },
  {
    "synthetic_test_masses": [
      0.0,
      0.03,
      0.2
    ],
    "finite_sigma_kinetic": 0.024040747483887688,
    "generalized_radial_kinetic": 0.012020373741943844,
    "values": [
      {
        "pE2": 0.0,
        "sigma_coordinate_kernel": 0.0052463857877380485
      },
      {
        "pE2": 0.001,
        "sigma_coordinate_kernel": 0.005270350262884032
      },
      {
        "pE2": 0.01,
        "sigma_coordinate_kernel": 0.005480462979198756
      },
      {
        "pE2": 0.1,
        "sigma_coordinate_kernel": 0.007230595350036515
      },
      {
        "pE2": 1.0,
        "sigma_coordinate_kernel": 0.009851657589428253
      }
    ]
  },
  {
    "synthetic_test_masses": [
      0.4381555712788046,
      0.4381555712788046,
      0.4381555712788046
    ],
    "finite_sigma_kinetic": -0.037986204483764933,
    "generalized_radial_kinetic": -0.018993102241882467,
    "values": [
      {
        "pE2": 0.0,
        "sigma_coordinate_kernel": 0.08751123730813674
      },
      {
        "pE2": 0.001,
        "sigma_coordinate_kernel": 0.08747313242873107
      },
      {
        "pE2": 0.01,
        "sigma_coordinate_kernel": 0.08711954727358544
      },
      {
        "pE2": 0.1,
        "sigma_coordinate_kernel": 0.08256720878275418
      },
      {
        "pE2": 1.0,
        "sigma_coordinate_kernel": -0.040326362458434566
      }
    ]
  }
]
```
