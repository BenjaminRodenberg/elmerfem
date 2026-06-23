Test for SubShelfMeltLinearTF solver
=====================================

Verifies the linear thermal forcing basal melt parameterisation:

  m = gammaT * (rhow * cp / (rhoi * Lf)) * (Tw - Tf)  [m/yr]

where Tf = lambda1*S + lambda2 + cc*z is the freezing point temperature.

Setup:
- 2D mesh (1000 m x 500 m)
- All nodes floating (GroundedMask = -1)
- Uniform ocean temperature: temp_oce_post = -1.0 degC
- Uniform ocean salinity: sal_oce_post = 34.0 PSU
- Ice base elevation: Zb = -500 m
- No water column scaling
- Grounding line melt = False

Expected result:
  Tf  = -0.0573*34 + 0.0832 + 7.61e-4*(-500) = -2.2455 degC
  TF  = Tw - Tf = -1.0 - (-2.2455) = 1.2455 degC
  m   = 1e-4 * (1028*3974 / (917*335500)) * 1.2455 * 31557600 ≈ 52.19 m/yr
