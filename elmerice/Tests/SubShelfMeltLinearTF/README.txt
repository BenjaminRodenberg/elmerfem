Test for SubShelfMeltLinearTF solver
====================================

This test now checks branch-sensitive behaviour, not only the base melt formula.

Base formula:

  m0 = gammaT * (rhow * cp / (rhoi * Lf)) * (Tw - Tf)  [m/yr]

where:

  Tf = lambda1*S + lambda2 + cc*z

Common setup:
- 2D mesh (1000 m x 500 m)
- Uniform ocean temperature: temp_oce_post = -1.0 degC
- Uniform ocean salinity: sal_oce_post = 34.0 PSU
- Ice base elevation: Zb = -500 m

Scenario 1 (floating + water column scaling):
- GroundedMask = -1
- water column scaling = True
- bedrock = -510 m, scaling factor = 75 m
- Expected:
  Tf      = -2.2455 degC
  Tw - Tf = 1.2455 degC
  m0      = 52.19234573 m/yr
  wct     = Zb - bedrock = 10 m
  scaling = tanh(10 / (75/e)) = 0.3473593091
  m       = 18.12949715 m/yr

Scenario 2 (grounding-line masking):
- GroundedMask = 0 everywhere
- grounding line melt = False
- Expected melt = 0 everywhere

Scenario 3 (grounded masking):
- GroundedMask = 1 everywhere
- grounding line melt = True
- Expected melt = 0 everywhere
