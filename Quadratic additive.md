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

1. **Concavity.** $\beta_w>0$, $\beta_i>0\ \forall i\in A^f$ and $\eta_j>0\ \forall j\in G^f$. This gives a negative definite Hessian (so the stationary point is a maximum) and also makes $\Beta$ and $\Eta$ strictly positive, which is what the divisions in (10)-(14) require.
2. **Positive marginal utilities.** The quadratic form has a bliss point, so the solution is economically meaningful only below it: $T_i^*<\theta_i/\beta_i$, $X_j^*<\phi_j/\eta_j$ and $T_w^*<\theta_w/\beta_w$. These are what guarantee $\mu>0$ and $\lambda>0$; outside that region (15) can turn negative or explode.
3. **Correct partition.** $A^f$ and $A^c$ (and $G^f$, $G^c$) are treated above as given, but they are endogenous. A candidate solution is a genuine optimum only if the multipliers of the binding constraints are non-negative, $\kappa_i=\mu-\theta_i+\beta_iT_i^{min}\geq0\ \forall i\in A^c$ and $\psi_j=P_j\lambda-\phi_j+\eta_jX_j^{min}\geq0\ \forall j\in G^c$, and if the free allocations respect their own bounds, $T_i^*\geq T_i^{min}$ and $X_j^*\geq X_j^{min}$. Otherwise the partition must be revised.
4. **Non-empty free sets.** $A^f\neq\emptyset$ and $G^f\neq\emptyset$, otherwise $\Beta$ or $\Eta$ is zero and $\mu$ or $\lambda$ is undefined.


