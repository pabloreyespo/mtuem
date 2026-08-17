## Microeconomic time-use model with Quadratic Additive Utility Function

$$
Max\ U(X,T) = \left(\theta_w T_w - \frac{1}{2}\beta_w T_w^2 \right) + \sum_{i\in A} \left( \theta_i T_i - \frac{1}{2}\beta_i T_i^2  \right) + \sum_{j\in G} \left( \phi_j X_j - \frac{1}{2}\eta_j X_j^2  \right) \tag{1}
$$
subject to:
$$
\omega T_w + I_f - \sum_{j\in G}P_j X_j \geq 0 \to (\lambda) \tag{2}
$$
$$\tau - T_w - \sum_{i\in A}T_i = 0 \to (\mu)\tag{3}$$
$$T_i - T_i^{min} \geq 0 \to (\kappa_i) \ \forall i\in A^c \tag{4}$$
$$X_j - X_j^{min} \geq 0 \to (\psi_j) \ \forall j\in G^c \tag{5} $$

Given the optimization problem (1)-(5) we can formulate the following Lagrangian:
$$
L(X,T;\lambda,\mu,\kappa_i,\psi_j)=\left(\theta_w T_w - \frac{1}{2}\beta_w T_w^2 \right) + \sum_{i\in A} \left( \theta_i T_i - \frac{1}{2}\beta_i T_i^2  \right) + \sum_{j\in G} \left( \phi_j X_j - \frac{1}{2}\eta_j X_j^2  \right) \\
+ \lambda \left(\omega T_w + I_f - \sum_{j\in G}P_j X_j \right)
+ \mu \left(\tau - T_w - \sum_{i\in A}T_i\right) \\
+ \sum_{i\in A^c} \kappa_i \left(T_i - T_i^{min}\right) 
+ \sum_{j\in G^c} \psi_j \left(X_j - X_j^{min}\right) 
\tag{6}
$$

The first order conditions are:
$$\frac{\partial L}{\partial T_w}= \theta_w - \beta_wT_w+\lambda\omega-\mu=0 \tag{7}$$
$$\frac{\partial L}{\partial T_i}= 
\begin{cases} 
\theta_i - \beta_iT_i-\mu=0  &,i \in A^f \\ 
\theta_i - \beta_iT_i-\mu+\kappa_i=0  &,i\in A^c
\end{cases} \tag{8}$$
$$\frac{\partial L}{\partial X_j}= 
\begin{cases} 
\phi_j - \eta_jX_j-P_j\lambda=0 &,j\in G^f \\
\phi_j - \eta_jX_j-P_j\lambda+\psi_j=0 &,j\in G^c
\end{cases} \tag{9} $$

The activity set $A$ is partitioned into freely allocated activities $A^f$ (interior, $\kappa_i=0$) and committed activities $A^c$ (at their lower bound, $T_i=T_i^{min}$, $\kappa_i\geq 0$). The goods set $G$ is partitioned analogously into $G^f$ (interior, $\psi_j=0$) and $G^c$ ($X_j=X_j^{min}$, $\psi_j\geq 0$). Committed consumption is therefore data rather than a choice, and it enters the problem only through the two aggregates
$$
T_c=\sum_{i\in A^c} T_i^{min}
\qquad\text{and}\qquad
E_c=\sum_{j\in G^c} P_j X_j^{min}-I_f
$$
so that $E_c$ is committed expenditure net of fixed income. Assuming non-satiation over the free sets, the monetary constraint (2) binds ($\lambda>0$), and the system below solves for the interior allocations.

Then, from equation (8) and (3) we get:
$$
\sum_{i\in A^f} \frac{\theta_i-\mu}{\beta_i}=\tau-T_c-T_w
$$
Given $S_\theta=\sum_{i\in A^f} \frac{\theta_i}{\beta_i}$ and $\Beta=\sum_{i\in A^f} \frac{1}{\beta_i}$:
$$
S_\theta-\Beta\mu=\tau-T_c-T_w \to
\mu= \frac{S_\theta - (\tau-T_c-T_w)}{\Beta} 
\tag{10}$$

From equations (9) and (2). Note that the monetary budget adds *expenditures* $P_jX_j$, not quantities, so each free-good FOC must be weighted by its own price before summing:
$$
\sum_{j\in G^f} P_j\left(\frac{\phi_j-P_j\lambda}{\eta_j}\right)=\omega T_w-E_c
$$
Given $S_\phi=\sum_{j\in G^f} \frac{P_j\phi_j}{\eta_j}$ and $\Eta=\sum_{j\in G^f}\frac{P_j^2}{\eta_j}$:
$$
S_\phi-\Eta\lambda=\omega T_w-E_c \to \lambda = \frac{S_\phi-(\omega T_w - E_c)}{\Eta}
\tag{11}
$$

Reemplazando (10) y (11) en (7):
$$\beta_wT_w - \theta_w -\omega \left[\frac{S_\phi-(\omega T_w - E_c)}{\Eta}\right] +\left[\frac{S_\theta - (\tau-T_c-T_w)}{\Beta} \right]=0 
\\ \ \\
\beta_w T_w + \frac{\omega^2}{\Eta}T_w + \frac{1}{\Beta}T_w= \theta_w + \omega \left[\frac{S_\phi+ E_c}{\Eta}\right] - \left[\frac{S_\theta - (\tau-T_c)}{\Beta} \right]
\\ \ \\
$$ 
$$
T_w^*=\frac{\theta_w + \omega \left[\frac{S_\phi+ E_c}{\Eta}\right] - \left[\frac{S_\theta - (\tau-T_c)}{\Beta} \right]}{\left(\beta_w+\frac{\omega^2}{\Eta} + \frac{1}{\Beta}\right)}
\tag{12}$$

Remplazando (10) en (8):
$$
T_i^*=\frac{\theta_i-\left[\frac{S_\theta - (\tau-T_c-T_w^*)}{\Beta}\right]}{\beta_i}
$$
$$
T_i^*=\frac{\theta_i}{\beta_i} - \frac{S_\theta - (\tau-T_c-T_w^*)}{\beta_i\Beta}
\tag{13}
$$
Reempazando (11) en (9):
$$
X_j^*=\frac{\phi_j-P_j\left[\frac{S_\phi-(\omega T_w^* - E_c)}{\Eta}\right]}{\eta_j}
$$
$$
X_j^*= \frac{\phi_j}{\eta_j}-\frac{P_j \left[S_\phi-(\omega T_w^* - E_c)\right]}{\eta_j \Eta}
\tag{14}
$$

Solution system corresponds to (12)-(14). 

From $\lambda$ and $\mu$ we can derivate the Value of Leisure and the Value of time assigned to work as:
$$
VoL = \frac{\mu}{\lambda}=\frac{\Eta}{\Beta}\cdot\frac{S_\theta - (\tau-T_c-T_w^*)}{S_\phi-(\omega T_w^* - E_c)}
\tag{15}
$$
$$
VTAW = \frac{\partial U/\partial T_w}{\lambda}= \frac{\theta_w - \beta_wT_w^*}{\lambda}
\tag{16}
$$
Dividing the work FOC (7) by $\lambda$ gives the equivalent and computationally simpler identity
$$
VoL = \omega + VTAW
\tag{17}
$$

### Regularity conditions

The solution (12)-(14) is valid only where the following hold, and they should be checked rather than assumed during estimation:

1. **Concavity.** $\beta_i>0\ \forall i\in A^f$ and $\eta_j>0\ \forall j\in G^f$, which also makes $\Beta$ and $\Eta$ strictly positive, as the divisions in (10)-(14) require. For work the condition is weaker than $\beta_w>0$: the stationary point is a constrained maximum if and only if the Hessian projected on the feasible set is negative definite, which holds exactly when the denominator of (12) is positive,
$$\beta_w+\frac{1}{\Beta}+\frac{\omega^2}{\Eta}>0 \tag{18}$$
so a moderately negative $\beta_w$ is still admissible. $\beta_w>0$ is the sufficient condition and the one worth imposing in estimation, because a negative satiation parameter has no clean reading and in practice signals an incomplete free set rather than a preference.
2. **Positive shadow prices, work excluded.** $\mu>0$ requires every free activity to be short of its bliss point and $\lambda>0$ requires every free good to be short of its own, since (8) and (9) can be written as
$$\mu=\beta_i\left(\frac{\theta_i}{\beta_i}-T_i^*\right)\ \forall i\in A^f
\qquad\text{and}\qquad
\lambda=\frac{\eta_j}{P_j}\left(\frac{\phi_j}{\eta_j}-X_j^*\right)\ \forall j\in G^f$$
so the conditions are $T_i^*<\theta_i/\beta_i$ and $X_j^*<\phi_j/\eta_j$. Work is subject to no such condition. The quantity $\theta_w-\beta_wT_w^*$ is free in sign and carries the sign of $VTAW$ in (16), which is the quadratic counterpart of the Cobb-Douglas property that $\theta_i$ and $\phi_j$ must be positive while $\theta_w$ "could be negative, positive or zero, an important analytical property indeed as it carries the sign of the marginal utility of labour" (Jara-Díaz et al., 2008, p. 949). Requiring $T_w^*<\theta_w/\beta_w$ would force $VTAW>0$; it is a restriction on the results, not a regularity condition, and must not be imposed.
3. **Correct partition.** $A^f$ and $A^c$ (and $G^f$, $G^c$) are treated above as given, but they are endogenous. A candidate solution is a genuine optimum only if the multipliers of the binding constraints are non-negative, $\kappa_i=\mu-\theta_i+\beta_iT_i^{min}\geq0\ \forall i\in A^c$ and $\psi_j=P_j\lambda-\phi_j+\eta_jX_j^{min}\geq0\ \forall j\in G^c$, and if the free allocations respect their own bounds, $T_i^*\geq T_i^{min}$ and $X_j^*\geq X_j^{min}$. Otherwise the partition must be revised.
4. **Non-empty free sets.** $A^f\neq\emptyset$ and $G^f\neq\emptyset$, otherwise $\Beta$ or $\Eta$ is zero and $\mu$ or $\lambda$ is undefined.

### Identification rules

The data are the allocations $(T_w^*, T_i^*, X_j^*)$ together with $\omega$, $P_j$, $\tau$, $T_c$ and $E_c$. Any transformation of $(\theta_w,\beta_w,\theta_i,\beta_i,\phi_j,\eta_j)$ that leaves the solution (12)-(14) unchanged cannot be recovered from those data, no matter how large the sample or how the errors are specified.

#### The invariance group

Three such transformations exist. All of them follow from one principle: **adding a multiple of a binding constraint to the objective changes utility by a constant on the feasible set, so it cannot move the optimum, but it moves the multiplier of that constraint one for one.** Since the quadratic family is closed under adding a linear function of $T$ and $X$, both constraints generate an admissible reparameterization.

| | Transformation | $U$ changes by | $\mu$ | $\lambda$ | Allocations |
|---|---|---|---|---|---|
| $I_1$ | $(\theta_\bullet,\beta_\bullet,\phi_\bullet,\eta_\bullet)\to a\,(\theta_\bullet,\beta_\bullet,\phi_\bullet,\eta_\bullet)$, $a>0$ | factor $a$ | $a\mu$ | $a\lambda$ | unchanged |
| $I_2$ | $\theta_w\to\theta_w+c$ and $\theta_i\to\theta_i+c\ \forall i\in A$ | $+c\,\tau$ | $\mu+c$ | $\lambda$ | unchanged |
| $I_3$ | $\phi_j\to\phi_j+P_jd\ \forall j\in G$ and $\theta_w\to\theta_w-\omega d$ | $+d\,I_f$ | $\mu$ | $\lambda+d$ | unchanged |

$I_2$ adds $c$ times the time budget (3), whose coefficients are $1$ for every activity including work, so the shift must be common to all of them and must include the committed ones (their $\theta_i$ enters only $\kappa_i$, which is invariant as well). $I_3$ adds $d$ times the money budget (2), whose coefficients are the prices for goods and $-\omega$ for work, so the shift is proportional to $P_j$ on the goods side and to $\omega$ on the work side.

$I_1$ is harmless: it is the usual cardinalization of utility and $VoL=\mu/\lambda$ is invariant under it. $I_2$ and $I_3$ are not: they move one shadow price and not the other, so

$$
VoL\to\frac{\mu+c}{\lambda+d},
\qquad
VTAW\to\frac{\mu+c}{\lambda+d}-\omega
\tag{19}
$$

**The level of the value of leisure is therefore not identified by allocation data.** This is the structural difference from the Cobb-Douglas precursors, where marginal utilities are $\theta_i U/T_i$ and no additive shift of the parameters reproduces a linear term, so the only invariance is $I_1$ and $\mu/\lambda$ survives it. Corner conditions do not help: $\kappa_i$ and $\psi_j$ are invariant under $I_2$ and $I_3$ as well, because they compare marginal utilities of time with each other and of goods with each other.

$I_3$ has one saving grace. It requires $\theta_w$ to move with the respondent's own wage and $\phi_j$ to move with the respondent's own price. If $\omega$ varies across the sample and $\theta_w$ is not specified as a function of $\omega$, then $I_3$ leaves the parameter space and the level of $\lambda$ is identified by the wage variation in (7). Nothing similar rescues $I_2$: the coefficient of every activity in the time budget is $1$ for every respondent, so there is no variation to exploit and $\mu$ has no identified level under any specification.

#### Evidence

Rank of the Jacobian of the stacked mean function $(T_w^*,T_1^*,E_1^*)$ over the 712 `maed` observations, with two free activities and two free goods, columns normalized:

| Specification | Fixed | $p$ | Rank | $s_{min}$ |
|---|---|---|---|---|
| $\theta_w$ constant | nothing | 10 | 8 | 1.0e-11 |
| $\theta_w$ constant | $\eta_1$ | 9 | 8 | 1.0e-11 |
| $\theta_w$ constant | $\eta_1,\theta_1$ | 8 | 8 | 2.67e-2 |
| $\theta_w$ constant | $\eta_1,\theta_w$ | 8 | 8 | 2.68e-2 |
| $\theta_w$ constant | $\eta_1,\beta_w$ | 8 | 7 | 1.0e-11 |
| $\theta_w$ constant | $\eta_1,\phi_1$ | 8 | 7 | 1.0e-11 |
| $\theta_w=\theta_{w0}+\theta_{w1}\omega$ | $\eta_1,\theta_1$ | 9 | 8 | 1.9e-9 |
| $\theta_w=\theta_{w0}+\theta_{w1}\omega$ | $\eta_1,\theta_1,\theta_{w1}$ | 8 | 8 | 2.67e-2 |
| $\theta_w=\theta_{w0}+\theta_{wz}z$, $z$ not the wage | $\eta_1,\theta_1$ | 9 | 9 | 2.75e-2 |

Two deficiencies with a constant $\theta_w$ ($I_1$ and $I_2$), three as soon as the wage enters $\theta_w$ ($I_3$ becomes admissible), and none of them is removed by fixing a $\beta$ or a $\phi$.

#### The rules

**R1. Utility scale.** Fix one strictly positive curvature parameter, for example $\eta_1=1$. This removes $I_1$.

**R2. Location of the time-side parameters.** Fix one $\theta$, and it must be a $\theta$. Write it as a bliss point,
$$
\theta_k=\beta_k A_k
\quad\Longleftrightarrow\quad
\frac{\theta_k}{\beta_k}=A_k
\tag{20}
$$
for one $k\in A^f\cup\{w\}$ and a constant $A_k$ chosen by the modeller. This removes $I_2$. Fixing $\beta_w$, $\phi_j$ or $\eta_j$ instead leaves the model unidentified, as the table shows.

**R3. Exclusion restrictions that keep $I_3$ out of the parameter space.** The wage must not enter $\theta_w$ and the price $P_j$ must not enter $\phi_j$. Any other covariate is admissible, and covariates are what the location restriction is stated relative to: $A_k$ pins the constant, not the slopes. If the wage genuinely belongs in $\theta_w$, a third restriction is needed and the level of $\lambda$, hence of every value of time, becomes a second assumption rather than an estimate.

**R4. Complete free sets.** $S_\theta$, $\Beta$, $S_\phi$ and $\Eta$ are literal sums over the categories supplied. An activity or good left out is deleted from the model rather than absorbed into an aggregate parameter, which biases every remaining equation. Categories whose equation is a budget identity should be kept in the aggregates and dropped from the likelihood, not dropped from the model.

With R1 to R4 the identified objects are: all $\beta$ and $\eta$ (in the scale set by R1), all differences $\theta_k-\theta_l$ among activities and work, all $\phi_j$, all allocations, and the identity $VoL-VTAW=\omega$. The levels of $\theta$, of $\mu$, of $VoL$ and of $VTAW$ are conditional on $A_k$.

### Choosing the anchor, and the sign of the values of time

R2 can be stated on any activity, and the choices are equivalent reparameterizations of each other: anchoring the reference free activity at $A_1$ and anchoring work at $A_w$ are related by $A_1=\left.A_1\right|_{A_w=0}+\left(\beta_w/\beta_1\right)A_w$, so each value of one corresponds to exactly one value of the other. They are not equally useful, because the anchor is exactly the quantity that decides the sign of the values of time. Anchoring work is the better coordinate, since it puts the assumption where the sign lives.

With $\theta_w=\beta_wA_w$, the work first order condition (7) gives, for $\lambda>0$,

$$
VTAW=\frac{\beta_w\left(A_w-T_w^*\right)}{\lambda},
\qquad
VoL=\omega+VTAW,
\qquad
\mu=\lambda\left(\omega+VTAW\right)
\tag{21}
$$

so that under $\beta_w>0$

$$
VTAW<0\iff T_w^*>A_w,
\qquad
\mu>0\iff VTAW>-\omega\iff A_w>T_w^*-\frac{\lambda\omega}{\beta_w}
\tag{22}
$$

$A_w$ reads as the bliss point of work, the hours a respondent would choose if work carried no monetary reward. It gives the model exactly the property the Cobb-Douglas precursors have through the sign of $\theta_w$: work is a marginal bad beyond $A_w$ hours and a marginal good below it, and whether an individual is above or below is decided by the data through $T_w^*$. The two regularity requirements bracket the anchor from both sides, giving the admissible window per respondent

$$
T_w^*-\frac{\lambda\omega}{\beta_w}\;<\;A_w\;<\;T_w^*
\tag{23}
$$

whose left inequality is time scarcity ($\mu>0$) and whose right inequality is a negative value of time assigned to work.

On the `maed` sample (712 observations, mean wage 12.254, $\hat\beta_w=11.03$, allocations and log likelihood identical at every $A_w$, $LL=-7992.29$):

| $A_w$ (h/week) | share $\mu>0$ | share $VTAW<0$ | share both, with $\lambda>0$ | mean $VoL$ | $VoL$ / mean $\omega$ | mean $VTAW$ |
|---|---|---|---|---|---|---|
| 0 | 0.725 | 0.989 | 0.723 | 2.79 | 0.23 | -9.47 |
| 20 | 0.955 | 0.978 | 0.942 | 8.06 | 0.66 | -4.20 |
| 33.3 | 0.990 | 0.767 | 0.760 | 11.57 | 0.94 | -0.69 |
| 41.1 | 0.996 | 0.407 | 0.396 | 13.62 | 1.11 | 1.37 |
| 58.1 | 0.999 | 0.017 | 0.007 | 18.09 | 1.48 | 5.84 |

The trade-off is monotone and unavoidable: raising the anchor buys time scarcity and spends the negative value of work. $A_w=0$, the classical reading that work has no intrinsic reward at the entry margin, delivers $VTAW<0$ for 98.9 per cent of the sample but leaves a quarter of it past the bliss point of its free activities, that is, with time not scarce. Three defensible criteria:

1. **External anchor, recommended.** Choose $A_w$ so that the mean $VoL$ matches the Cobb-Douglas model estimated on the same sample, where the level is identified by functional form. On `maed` that model gives $\theta_w=-0.052<0$, hence $VTAW<0$ for the whole sample, and a mean $VoL$ of 11.57, which is 0.94 times the mean wage. Matching it requires $A_w=33.3$ hours per week, and at that anchor the quadratic model has $\mu>0$ for 99.0 per cent of the sample and $VTAW<0$ for 76.7 per cent. The level then comes from the model that can identify it and the quadratic model contributes what it does better, namely the curvature and the heterogeneity. The price of this route is explicit: the level is inherited together with the Cobb-Douglas functional form, so it should be reported as such and not as an independent finding of the quadratic model.
2. **Maximum regularity.** Choose $A_w$ to maximize the share of respondents satisfying $\mu>0$, $\lambda>0$ and, if the sign is the object of interest, $VTAW<0$. On `maed` that is $A_w\approx20$ hours with 94.2 per cent joint compliance. This is a minimum-violation criterion, not an estimate: the log likelihood is flat in $A_w$.
3. **Behavioural anchor.** State $A_w$ directly as the weekly hours of work a respondent would choose for its own sake, with $A_w=0$ as the pure instrumental-labour limit, and report the sensitivity.

Whichever is used, the anchor is an assumption and must travel with the result. Do not attempt to estimate $A_w$: the likelihood is flat in it, and an optimizer allowed to move it will either stall or drift to a boundary. Do not impose $\theta_w>0$ or estimate $\theta_w$ as $\exp(\cdot)$, which forces the work bliss point to be positive and, combined with a large $\beta_w$, silently restricts the sign of $VTAW$.

### What to report

- The anchor $A_w$ (or the equivalent $A_1$) explicitly, next to every value of time, which is conditional on it.
- The contrasts $\theta_w-\theta_i$ rather than $\theta_w$ or $\theta_i$ alone, since only the contrasts are identified.
- $VoL-VTAW=\omega$ as the invariant, and comparisons of values of time across individuals or groups in preference to levels.
- The regularity shares: $\mu\leq0$, $\lambda\leq0$, $T_i^*<T_i^{min}$, $X_j^*<X_j^{min}$, and $T_w^*$ outside $[0,\tau-T_c]$.

The companion note `docs/additive-quadratic-identification.md` carries the empirical detail behind these rules.


