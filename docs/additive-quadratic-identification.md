# Identification and regularization of the additive quadratic MTUEM

Companion to `Quadratic additive.md`, which derives the model. This note covers what the
model can and cannot identify from time-use and expenditure data, why an unregularized fit
produces negative values of leisure, and which restrictions to impose so the estimated
system keeps microeconomic meaning.

All numbers below come from the `maed` data used in `examples/maed_3eq_cholesky_aq.r`
(712 observations, `EcI > 0`), unless stated otherwise.

## 1. Summary

The additive quadratic MTUEM has two exact invariances. Both must be removed before
estimation, and only one of the two is the familiar utility-scale normalization.

| Invariance | Transformation leaving all allocations unchanged | Removed by |
|---|---|---|
| Utility scale | multiply every structural parameter by `c > 0` | fix one parameter, e.g. `eta_1 = 1` |
| Theta location | add `c` to `theta_w` and to every `theta_i` | fix one theta, or anchor it |

The second one has a consequence with no counterpart in the Cobb-Douglas versions of the
model: **the level of the value of leisure is not identified**. Choosing the theta
normalization chooses the reported VoL. Only `VoL - VTAW = w` is invariant.

The recommended regularization is three layers: two normalizations for identification,
concavity constraints for a well-posed maximum, and a check (not an imposition) that both
shadow prices come out positive.

## 2. The two invariances

**Scale.** Every optimality condition is homogeneous of degree zero in the parameter vector.
Multiplying `theta_w, beta_w, theta_i, beta_i, phi_j, eta_j` by `c` scales `S_theta` and
`S_phi` not at all, scales `B` and `H` by `1/c`, and scales `mu` and `lambda` by `c`. The
allocations, and the ratio `mu/lambda`, are unchanged.

**Theta location.** The time budget makes `T_w + sum_i T_i = tau - T_c` a constant. Adding
`c` to `theta_w` and to every `theta_i` therefore adds `c(tau - T_c)` to utility, a constant,
and leaves preferences over the feasible set untouched. Formally `S_theta` rises by `cB`, so
`mu` rises by exactly `c`, and

```
T_i* = theta_i/beta_i - mu/beta_i = (theta_i + c)/beta_i - (mu + c)/beta_i
```

is unchanged. `lambda` does not move at all, because the money budget has no analogous
structure.

Numerical confirmation (2 free activities, 3 free goods, real data):

```
multiply all 12 structural parameters by 3   -> max change in allocations  7e-13
add 2 to theta_w and to every theta_i        -> max change in allocations  5e-13
add 2 to theta_1 only            (control)   -> max change in allocations  30.1
add 2 to every phi_j             (control)   -> max change in allocations  2706
```

## 3. Rank evidence

Numerical rank of the Jacobian of the stacked mean function with respect to the structural
parameters, over all 712 observations.

Specification with 1 free activity and 1 free good (6 structural parameters):

| Fixed | Rank | Smallest normalized singular value |
|---|---|---|
| nothing | 4 / 6 | 3e-13 |
| `eta_1` | 4 / 5 | 4e-11 |
| `eta_1` + `theta_w` | 4 / 4 | 1.4e-2 |
| `eta_1` + `theta_1` | 4 / 4 | 1.8e-3 |
| `eta_1` + `beta_w` | 3 / 4 | 1.2e-10 |

Specification with 2 free activities and 3 free goods (12 structural parameters): 10/12 with
nothing fixed, 10/11 with `eta_1` fixed, 10/10 with `eta_1` and `theta_1` fixed.

Two readings matter. First, fixing `eta_1` alone is never enough. Second, the second
restriction has to fall on a theta: fixing `beta_w` instead leaves the rank deficient.

Estimation confirms it directly. With `eta_1` as the only fixed parameter, apollo reports
`Singular convergence (unfavorable)`, computes no Hessian, and returns `NA` for every
standard error, while reaching exactly the same maximum:

```
eta_1 only          theta_w = -6033.153   theta_1 = -4856.511   LL = -8907.87
eta_1 + theta_1     theta_w = -1176.645   theta_1 =        0    LL = -8907.87
                    difference -6033.153 - (-4856.511) = -1176.642
```

The identified object is the contrast `theta_w - theta_1`, reproduced to three decimals
under both parameterizations.

## 4. Consequence: the value of leisure has no identified level

`mu` moves one for one with the theta shift and `lambda` does not, so

```
VoL  = mu / lambda                shifts by c / lambda
VTAW = (theta_w - beta_w T_w*) / lambda    shifts by c / lambda
VoL - VTAW = w                    invariant
```

This is the structural difference from the Cobb-Douglas parameterizations. There the only
invariance is a scaling of utility, which multiplies `mu` and `lambda` equally and cancels in
`mu/lambda`, so VoL is identified. Here the invariance is additive and touches the theta side
only, so it does not cancel.

Varying the normalization while holding the fit fixed moves mean VoL from 1.52 to 6.52 times
the mean wage with an identical log-likelihood and identical allocations. Section 5 gives the
full table for the shipped specification. The ranking of VoL across individuals is nearly
preserved across normalizations (Spearman 0.999), so relative comparisons survive; levels do
not. Mean wage in the sample is 12.254.

No estimation strategy recovers the level, because no allocation data can. Corner conditions
do not help either: the work corner condition bounds `mu - theta_w`, which is already
identified. Pinning the level requires either an explicit external assumption or a utility
function in which the additive theta shift stops being neutral.

## 5. The bliss point, and why an unregularized fit gives negative VoL

### What a bliss point is

The marginal utility of a quadratic activity is a straight line:

```
dU_i/dT_i = theta_i - beta_i T_i
```

It starts at `theta_i` when `T_i = 0` and falls with slope `-beta_i`, crossing zero at

```
B_i = theta_i / beta_i
```

`B_i` is the **bliss point**: the hours of activity `i` this person would choose if time were
free. Past `B_i`, more of the activity makes them worse off. Completing the square shows the
same thing as a parabola peaking at `B_i`:

```
theta_i T_i - 0.5 beta_i T_i^2 = -0.5 beta_i (T_i - B_i)^2 + 0.5 beta_i B_i^2
```

So `(theta_i, beta_i)` and `(B_i, beta_i)` carry exactly the same information, and the second
pair is the readable one: **how much of it you want**, and **how sharply you mind not getting
it**.

### Why the bliss point controls the shadow price of time

At an interior optimum the first order condition equalizes marginal utilities across free
activities at the time shadow price, `theta_i - beta_i T_i* = mu`, which rearranges to

```
mu = beta_i (B_i - T_i*)      for every freely allocated activity i
```

The marginal value of an extra hour is *how far below your bliss point you ended up*, times
*how sharply you feel it*. So `mu > 0` is precisely the statement that every free activity got
less than its bliss point, that is, that time is scarce. That is the microeconomic content of
the sign condition.

### What the theta shift does in these coordinates

Adding a constant `c` to `theta_w` and to every `theta_i` moves each bliss point to
`B_i + c/beta_i` and moves `mu` to `mu + c`, while leaving every `T_i*` where it was. The data
therefore pins down the shape and the relative positions but not the absolute level of "how
much people want". One bliss point has to be set by assumption. That assumption is the
**anchor**:

```
theta_1 = beta_1 * A     i.e.     B_1 = A
```

read as: *if the reference free activity were costless, this person would choose `A` hours of
it and be satiated there.* Everything else follows, `mu = beta_1 (A - T_1*)` and
`VoL = mu / lambda`.

### Why `theta_1 = 0` produced negative VoL

`theta_1 = 0` means `B_1 = 0`: the reference activity satiates at zero hours. Every hour
actually spent on it is then past satiation with negative marginal utility, so `mu < 0` for
the whole sample. On the `maed` fit that produced

```
bliss point of Tf1 = theta_1/beta_1 = 0 hours,   mean T_1* = 45.98 hours
100% of the sample past the bliss point
mean mu     = -533.9   (negative for 100% of the sample)
mean lambda =   54.0   (positive for 100%; X_1* = 72.5 below its bliss point 126.6)
mean VoL    =  -13.76
```

The goods side was never the problem. Negative VoL was a direct consequence of the
normalization value, not a property of the data. This is what motivates anchoring the theta
location on a bliss point rather than on zero.

### How to choose `A`

Three considerations, in order.

1. **Lower limit, non-negotiable.** `A > max T_i*` is what makes `mu > 0` for every
   respondent. Under the shipped specification `max T_f* = 74.21` hours, so `A >= 75`.
2. **Upper reference points.** Natural candidates are the largest discretionary time budget
   in the sample (`168 - min(Tc) = 112.73` hours here) or the full week. The shipped
   specification uses `A = 168`, "satiated only if free time filled the entire week", the most
   permissive upper anchor and the one that needs no reference to the sample at hand.
3. **`A` is not estimable.** The log-likelihood is flat in it, and `mu`, hence VoL, is linear
   in it. Choosing `A` chooses the reported level.

Sensitivity under the shipped specification, holding the identified quantities fixed
(allocations identical to 2e-14 throughout). With more than one free activity, lowering `A`
shifts every other bliss point too, `B_i -> B_i + beta_1(A' - A)/beta_i`, so what has to stay
above the largest predicted free time is the sum `S_theta = A + B_i`. The smallest admissible
anchor here is `A = 60.04`:

| anchor `A` | mean VoL | VoL / mean wage | share `mu <= 0` |
|---|---|---|---|
| 61 (floor) | 18.63 | 1.52 | 0 |
| 90 | 35.25 | 2.88 | 0 |
| 112 | 47.85 | 3.90 | 0 |
| 168 (used) | 79.93 | 6.52 | 0 |

The smallest anchor compatible with time being scarce for everyone gives
`mean VoL = 18.63`, about `1.52 x` the mean wage. That is a floor, not an estimate: the
specification cannot produce a VoL below it without negative shadow prices for part of the
sample. Reporting that floor alongside the chosen `A` is the honest way to present the
result.

## 6. Recommended regularization

### Layer 1: identification (mandatory, exactly two restrictions)

1. **Utility scale.** `eta_1 = 1`.
2. **Theta location.** `theta_1 = beta_1 * A` with `A` a fixed constant, that is, anchor the
   bliss point of the reference free activity at `A` hours.

`theta_1` stops being a free parameter, so the count drops by one and the rank deficiency is
removed. `theta_w` stays free.

### Layer 2: concavity (check, and impose only if needed)

`beta_w > 0`, `beta_i > 0`, `eta_j > 0` are required for the stationary point of (12)-(14) to
be a maximum. `beta_1` is estimated as `exp(lbeta_1)` in the shipped example so it cannot go
wrong. `beta_w` is left free and comes out positive on its own:

```
beta_w free   LL = -7992.29   beta_w = +11.03 (rob. t = 4.00)
```

Concavity therefore holds without being imposed. Estimating `exp(lbeta_w)` instead gives an
identical fit, so use it if you want the guarantee on other data.

**This conclusion is specification dependent, and that is worth dwelling on.** Under the
earlier misspecified free set, with `Tf1` modelled alone rather than the complete aggregate,
`beta_w` estimated at `-2.304` (rob. t `-1.31`), which violates concavity, and imposing
`beta_w >= 0` drove it to the boundary at zero. Completing the free sets flipped its sign and
made it significant. A negative satiation parameter is therefore a useful
diagnostic: treat it as a signal that the free sets are wrong before concluding anything
about preferences.

### Layer 3: positive shadow prices (choose `A` for `mu`, check `lambda`)

`mu > 0` for every observation is controlled by the anchor, which is what section 5 covers.
`lambda > 0` and `X_j* < phi_j/eta_j` are not controlled by `A` at all. Verify them after
estimation rather than imposing them. Under the shipped specification `mu > 0` everywhere,
but `lambda < 0` for 7 of 712 whose predicted free expenditure runs past the aggregate goods
bliss point, and `Ei* < 0` for 81 of 712 (11.4%). Those are goods-side corner problems, not
anchor problems, and the goods side is where this specification is weakest.

## 6b. The free sets must be complete

This is the difference from the Cobb-Douglas parameterizations that is easiest to get wrong.

In `get_ti_thph`, `T_i* = (theta_i / Theta) * (tau - Tw - Tc)`, and `Theta` is its **own
estimated parameter** covering all free activities. An activity left out of `free_times` is
still in the model; its share is `Theta - sum(modelled theta_i)`. Leaving one out is
therefore free, and does only one job: keeping the error covariance non-singular.

Here there is no aggregate parameter. `S_theta = sum(theta_i/beta_i)` and `B = sum(1/beta_i)`
are literal sums over what is passed in, so an activity left out is **deleted from the
model** and the remaining ones must absorb its time. Modelling `Tf1` alone on `maed` gave

```
predicted Tf1 = tau - Tc - Tw*   mean 45.963
observed  Tf1                    mean 28.859    residual mean -17.104
observed  Tf1 + Tf2              mean 40.241    residual mean  -5.722
mean Tf2                                11.382
```

and `-17.104 - (-5.722) = -11.382`, exactly the mean of the omitted activity. The bias in the
Tf1 equation *is* the omission.

### How many equations, then?

The daily cycle here is five quantities under two adding-up constraints:

```
Tw + T1 + Ti + Tc = tau            E1 + Ei + Ec = w Tw
```

so only **three** residuals are statistically independent. Estimating five equations gives a
rank-3 covariance. Three equations is the same count as the Cobb-Douglas versions; what
changes is only *how* the omitted categories are handled.

**Disaggregate, which is what the shipped example does.** Declare every free category with
its own parameters so the aggregates are complete, then remove `Ti` and `Ei` from the
likelihood with `omitted_times` / `omitted_goods`. They keep driving `Tw*` through
`S_theta`, `B`, `S_phi`, `H`, and they are still predicted from their own parameters rather
than left over. Estimated equations: `Tw`, `T1`, `E1`, hence a 3x3 Cholesky with 6 elements.

**Aggregate, the simpler alternative.** Collapse the free sets in the data,
`Tf = Tf1 + Tf2` and `Ef = Ef1 + Ef2 + Ef3`. Also correctly specified, but `T_f*` is then
the time budget itself, and the data satisfies that budget up to rounding:

```
err_Tf = -err_Tw - (rounding),    cor(err_Tw, err_Tf) = -0.999914
rounding sd 0.095 h  against an err_Tw sd of 7.19 h
smallest eigenvalue of the 3x3 residual correlation matrix: 8.6e-05
```

That passes the `1e-6` guard, so the model would estimate while fitting 0.095 hours of
rounding noise. `omitted_times` handles it, and only two equations survive, `Tw` and `Ef`.
This route loses the free-activity split, so prefer the disaggregated one unless the
categories genuinely are one composite.

Note the general rule for the covariance: `M` estimated equations need `M(M+1)/2` Cholesky
elements. Three equations means six; two equations means three.

## 7. What to report

- **Report `theta_w - theta_1`, not `theta_w`.** The sign of `theta_w` alone is meaningless:
  across anchors it reads 123.89, 774.15 and 2649.60 with an identical fit. The contrast is
  what supports the statement that the baseline marginal utility of work sits below that of
  the reference free activity.
- **Report the anchor `A` explicitly** alongside any VoL or VTAW figure, and treat those as
  conditional on it.
- **Report `VoL - VTAW = w`** as the invariant, and prefer relative VoL comparisons across
  individuals or groups over levels.
- **Report the regularity checks**: shares of the sample with `mu <= 0`, `lambda <= 0`,
  `T_i* < 0`, `X_j* < 0`, and allocations above their bliss points.

## 8. The shipped specification

`examples/maed_3eq_cholesky_aq.r` implements the aggregated version:

```
daily cycle   Tw, T1, Ti (times)   E1, Ei (expenditures)
free sets     A^f = {T1, Ti},  G^f = {E1, Ei}             (complete)
omitted       Ti, Ei                                      (recovered from the budgets)
equations     Tw, T1, E1                                  (3x3 Cholesky, 6 elements)
eta_1  = 1                                                (utility scale)
theta_1 = beta_1 * 168                                    (theta location, full week)
beta_1, beta_i, eta_i = exp(.)                            (concavity)
beta_w  free                                              (comes out positive)
```

Free activities are parameterized as `(satiation, bliss point)` rather than
`(theta, beta)`, which is the same model but far better scaled for the optimizer.

Result: `LL = -7992.29`, 14 parameters, negative definite Hessian, favourable convergence,
every robust standard error finite and every parameter significant.

```
theta_w  3229.42 (rob. t 6.82)   beta_w  11.03 (rob. t 4.00)
beta_1     23.98   beta_i 334.63   B_i 21.35 h
Xb_1      131.71   eta_i   0.437   Xb_i 181.90
theta_w - theta_1 = -799.91    theta_w - theta_i = -3913.82    <- identified contrasts
sig_Tw 6.63   sig_T1 7.18   sig_E1 44.23
rho_Tw_T1 -0.787   rho_Tw_E1 0.545   rho_T1_E1 -0.544

fit of the daily cycle (means)      predicted   observed
                            Tw         38.454     38.310
                            T1         28.736     28.859
                            Ti         11.368     11.382
                            E1         81.207     80.681
                            Ei         66.458     66.559
both budgets close to 1e-14

mean VoL 79.93 = 6.52 x mean wage at A = 168; floor is 18.63 = 1.52 x
regularity: beta_w > 0 (concave); VoL <= 0 for 7 of 712 (all lambda < 0);
            Tw* out of range 1; T1*, Ti*, E1* >= 0 everywhere; Ei* < 0 for 81
```

Compare against deleting `Ti`: that specification predicted `Tw` at 32.60 against an observed
38.31 and gave a concavity-violating `beta_w = -2.30`. With the free sets complete, `Tw`
predicts at 38.45 and `beta_w` is `+11.03` with rob. t 4.00.

The anchor remains free: re-running at `A = 112` leaves the log-likelihood, the allocations,
all satiation parameters and the whole Cholesky factor unchanged. Only the theta level moves,
and mean VoL goes from 79.93 to 47.85. That is the non-identification of the VoL level,
visible in one comparison.

## 9. Open issues

**Corner solutions.** `T_w* < 0` for 1 of 712 observations, `lambda < 0` for 7, and
`Ei* < 0` for 81 of 712 (11.4%), the last being much the most common.
Those violate the interior-solution assumption and belong in the corner branch of the KKT
conditions, which the implementation does not handle. See regularity condition 3 in
`Quadratic additive.md`.

**Heterogeneous parameters.** `get_additive_quadratic_essentials` reduces its inputs with
`sum()`, which collapses across observations as well as across activities. Specifying any
theta, beta, phi or eta as a function of covariates therefore produces silently wrong
aggregates. The reduction has to run across activities only before covariates can be used.
