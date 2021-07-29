# Solve_Generic_HJB_KFE
Julia code to solve generic HJBs &amp; KFEs

Current code: deterministic, one state, one control, infinite horizon, time/state separable

Extensions:     
endog := endogenous & exog := exogenous    
- Multiple states (n_s): 
  - (a,z): wealth `a` endog, income `z` is exog 2-state Poisson, [HJB_stateconstraint_implicit](https://benjaminmoll.com/wp-content/uploads/2020/06/HJB_stateconstraint_implicit.m)
  - (a,z): wealth `a` endog, income `z` exog diffusion, [HJB_diffusion_implicit](https://benjaminmoll.com/wp-content/uploads/2020/06/HJB_diffusion_implicit.m)
  - (a,z,t): wealth `a` endog, income `z` exog diffusion, age `t` exog, [Lifecycle.pdf](https://benjaminmoll.com/wp-content/uploads/2020/06/lifecycle.pdf), [lifecycle.m](https://benjaminmoll.com/wp-content/uploads/2020/06/lifecycle.m)
  - (k,z): capital `k` endog, tech `z` exog diffusion, [firm.pdf](https://benjaminmoll.com/wp-content/uploads/2020/06/firm.pdf), [firm.m](https://benjaminmoll.com/wp-content/uploads/2020/06/firm.m)
- Multiple controls (n_c):
  - (a,z): wealth `a` endog, income `z` exog 2-state Poisson, choices (c,ℓ), [labor_supply.pdf](https://benjaminmoll.com/wp-content/uploads/2020/06/labor_supply.pdf),  [HJB_labor_supply.m](https://benjaminmoll.com/wp-content/uploads/2020/06/HJB_labor_supply.m)
  - (a,h): wealth `a` endog, human capital `h` endog, choices (c,s), [human_capital.pdf](https://benjaminmoll.com/wp-content/uploads/2020/02/human_capital.pdf),  [human_capital.m](https://benjaminmoll.com/wp-content/uploads/2020/03/human_capital.m)
  - (a,b,z): liquid assets `a` endog, illiquid assets `b` endog, income `z` exog Poiss, choices (c,d), [two_asset_kinked.pdf](https://benjaminmoll.com/wp-content/uploads/2020/06/two_asset_kinked.pdf),  [two_asset_kinked.m](https://benjaminmoll.com/wp-content/uploads/2020/06/two_asset_kinked.m)
- Stochastic: see above
