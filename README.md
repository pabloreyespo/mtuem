
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mtuem

<!-- badges: start -->
<!-- badges: end -->

The goal of mtuem is to facilitate the estimation of Microeconomic
Time-Use-Expenditure models in apollo (see:
<https://www.apollochoicemodelling.com/>)

## Installation

You can install the development version of mtuem like so: 1. Install
devtools

``` r
install.packages("devtools")
```

2.  Install the repo via devtools

``` r
devtools::install_github("https://github.com/pabloreyespo/mtuem")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(mtuem)
#> Cargando paquete requerido: apollo
#> 
#> 
#>              . ,,                                                            
#>             ,      ,,                                                        
#>  ,,,,,,    ,         ,,                                                      
#> ,     ,,  ,            ,,,,.                                                 
#> ,,     , ,,   ,,,,,,    ,,,                                 //  //           
#>   ,     ,,,.   ,,,,,.   ,,      ////                        //  //           
#> ,,     ,,,,,.           ,,     // //     //////    /////    //  //    /////  
#> ,,,        ,,           ,      //  //    /    //  //   //   //  //   //   // 
#>               ,,       ,      ////////   /    //  //   //   //  //   //   // 
#>                 ,,   ,,      //     //   /   ///  //   //   //  //   //   // 
#>                    ,         //      //  /////      ///      //  //    ///   
#>                                          //                                  
#>                                          //                                  
#> 
#> Apollo 0.3.5
#> http://www.ApolloChoiceModelling.com
#> See url for a detailed manual, examples and a user forum.
#> Sign up to the user forum to receive updates on new releases.
#> 
#> Please cite Apollo in all written material you produce:
#> Hess S, Palma D (2019). "Apollo: a flexible, powerful and customisable
#> freeware package for choice model estimation and application." Journal
#> of Choice Modelling, 32. doi.org/10.1016/j.jocm.2019.100170
#> 
#> The developers of Apollo acknowledge the substantial support provided by
#> the European Research Council (ERC) through the consolidator grant DECISIONS,
#> the proof of concept grant APOLLO, and the advanced grant SYNERGY.

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName       = "enut-mtuem-1eq",
  indivID         = "id_persona"
)


# ################################################################# #
#### LOAD DATA AND APPLY ANY TRANSFORMATIONS                     ####
# ################################################################# #

database = mtuem::enut.ii
database['Tc'] = rowSums(database[,c("Tc", "Tc_sleep", "Tc_meals")])
is_not_retist = database["Ec"] > 0
is_worker = database["is_worker"] == 1
can_afford_expenses =  (database["Ec"] / (database["w"]*(168 -  database["Tc"]))) < 1
database = database[is_not_retist & can_afford_expenses & is_worker,]

apollo_beta = c(Theta   = 1,
                Phi     = 1,
                thw = 0)

apollo_fixed = c("Theta")

apollo_inputs = apollo_validateInputs(
  database = database,
  apollo_beta = apollo_beta,
  apollo_fixed = apollo_fixed,
  apollo_control = apollo_control 
)
#> All checks on apollo_control completed.
#> All checks on database completed.

#tricking apollo
#mtuem_likelihood = mtuem::mtuem_likelihood

apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  ### Create list of probabilities P
  P = list()

  ### Define individual alternatives
  work_times = c("Tw")

  # must use this names
  work_elasticities = list(
    Theta = Theta,
    Phi   = Phi,
    thw = thw
  )

  ### Define settings for Jara-Diaz model
  mtuem_settings = list(
    work_times = work_times,
    work_elasticities = work_elasticities,
    Tc = Tc,
    Ec = Ec,
    w = w,
    tau = 168
  )

  ### Compute probabilities using Jara-Diaz model
  P[["model"]] = mtuem_likelihood(mtuem_settings, functionality)

  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

suppressWarnings({
  model = apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)
  apollo_modelOutput(model)
})
#> Preparing user-defined functions.
#> 
#> Testing likelihood function...
#> 
#> Pre-processing likelihood function...
#> INFORMATION: Apollo was not able to compute analytical gradients for your model.
#>   This could be because you are using model components for which
#>   analytical gradients are not yet implemented, or because you coded
#>   your own model functions. If however you only used apollo_mnl,
#>   apollo_fmnl, apollo_normalDensity, apollo_ol or apollo_op, then there
#>   could be another issue. You might want to ask for help in the Apollo
#>   forum (http://www.apollochoicemodelling.com/forum) on how to solve
#>   this issue. If you do, please post your code and data (if not
#>   confidential). 
#> Analytical gradients could not be calculated for all components,
#>   numerical gradients will be used.
#> 
#> Testing influence of parameters..
#> Starting main estimation
#> 
#> BGW is using FD derivatives for model Jacobian. (Caller did not provide derivatives.)
#> 
#> 
#> Iterates will be written to: 
#>  C:/Users/pablo/OneDrive/UdeC/Investigación/apollo_timeuse/mtuem/enut-mtuem-1eq_iterations.csv
#>     it    nf     F            RELDF    PRELDF    RELDX    MODEL stppar
#>      0     1 1.964498496e+04
#>      1     4 1.905531350e+04 3.002e-02 2.495e-02 3.27e-02   G   1.37e+00
#>      2     5 1.834376482e+04 3.734e-02 2.751e-02 1.95e-01   G   2.26e-02
#>      3     6 1.830390289e+04 2.173e-03 3.926e-03 6.50e-02   S   0.00e+00
#>      4     7 1.829727878e+04 3.619e-04 9.650e-04 3.92e-02   S   0.00e+00
#>      5     8 1.829121979e+04 3.311e-04 3.198e-04 3.81e-02   S   0.00e+00
#>      6     9 1.829073812e+04 2.633e-05 1.735e-05 9.61e-03   S   0.00e+00
#>      7    10 1.829055553e+04 9.983e-06 1.072e-05 8.83e-03   G   0.00e+00
#>      8    11 1.829054638e+04 4.998e-07 3.824e-07 3.87e-04   G   0.00e+00
#>      9    12 1.829054434e+04 1.117e-07 7.359e-08 6.26e-04   G   0.00e+00
#>     10    13 1.829054353e+04 4.434e-08 4.418e-08 6.57e-04   S   0.00e+00
#>     11    14 1.829054353e+04 9.281e-12 8.879e-12 4.62e-06   S   0.00e+00
#> 
#> ***** Relative function convergence *****
#> 
#> Estimated parameters with approximate standard errors from BHHH matrix:
#>          Estimate     BHHH se BHH t-ratio (0)
#> Theta     1.00000          NA              NA
#> Phi       0.66546     0.01962          33.909
#> thw      -0.07547     0.02840          -2.658
#> 
#> Final LL: -18290.5435
#> 
#> Calculating log-likelihood at equal shares (LL(0)) for applicable
#>   models...
#> Calculating log-likelihood at observed shares from estimation data
#>   (LL(c)) for applicable models...
#> Calculating LL of each model component...
#> VoL: 4.653763 
#> VTAW: -0.314996 
#> Computing covariance matrix using numerical methods (numDeriv).
#>  0%....25%....50%...100%
#> Negative definite Hessian with maximum eigenvalue: -563.860558
#> Computing score matrix...
#> 
#> Your model was estimated using the BGW algorithm. Please acknowledge
#>   this by citing Bunch et al. (1993) - doi.org/10.1145/151271.151279
#> 
#> Please acknowledge the use of Apollo by citing Hess & Palma (2019) -
#>   doi.org/10.1016/j.jocm.2019.100170
#> Model run by pablo using Apollo 0.3.5 on R 4.4.3 for Windows.
#> Please acknowledge the use of Apollo by citing Hess & Palma (2019)
#>   DOI 10.1016/j.jocm.2019.100170
#>   www.ApolloChoiceModelling.com
#> 
#> Model name                                  : enut-mtuem-1eq
#> Model description                           : No model description provided in apollo_control
#> Model run at                                : 2025-06-25 13:34:12.834823
#> Estimation method                           : bgw
#> Model diagnosis                             : Relative function convergence
#> Optimisation diagnosis                      : Maximum found
#>      hessian properties                     : Negative definite
#>      maximum eigenvalue                     : -563.860558
#>      reciprocal of condition number         : 0.0144147
#> Number of individuals                       : 5043
#> Number of rows in database                  : 5043
#> Number of modelled outcomes                 : 0
#> 
#> Number of cores used                        :  1 
#> Model without mixing
#> 
#> LL(start)                                   : 0
#> LL at equal shares, LL(0)                   : NA
#> LL at observed shares, LL(C)                : NA
#> LL(final)                                   : -18290.54
#> Rho-squared vs equal shares                  :  Not applicable 
#> Adj.Rho-squared vs equal shares              :  Not applicable 
#> Rho-squared vs observed shares               :  Not applicable 
#> Adj.Rho-squared vs observed shares           :  Not applicable 
#> AIC                                         :  36585.09 
#> BIC                                         :  -Inf 
#> 
#> Estimated parameters                        : 2
#> Time taken (hh:mm:ss)                       :  00:00:0.57 
#>      pre-estimation                         :  00:00:0.23 
#>      estimation                             :  00:00:0.27 
#>      post-estimation                        :  00:00:0.06 
#> Iterations                                  :  11  
#> 
#> Unconstrained optimisation.
#> 
#> Estimates:
#>          Estimate        s.e.   t.rat.(0)    Rob.s.e. Rob.t.rat.(0)
#> Theta     1.00000          NA          NA          NA            NA
#> Phi       0.66546     0.02216      30.026     0.02752        24.180
#> thw      -0.07547     0.03616      -2.087     0.04937        -1.529

# Compute values of time
pred <- apollo_prediction(model, apollo_probabilities, apollo_inputs)
#> Running predictions from model using parameter estimates...
#> Prediction at user provided parameters
#>                  Tw        T        X      VoL     VTAW
#> Aggregate 210940.40 173970.3 551679.3 23468.93 -1588.53
#> Average       41.83     34.5    109.4     4.65    -0.31
#> 
#> The output from apollo_prediction is a matrix containing the
#>   predictions at the estimated values.
Tw = pred$Tw
w  = database$w
Tc = database$Tc
Ec = database$Ec
cteVoL  <- mean((w*Tw - Ec) / (168-Tw-Tc))
cteVTAW <- mean((w*Tw - Ec) / (Tw))
delta <- apollo_deltaMethod(model, list(
  expression=c(
    VoL      = paste0(cteVoL,  "*Theta/Phi"),
    VTAW     = paste0(cteVTAW, "*thw/Phi"))
))
#> The expression VoL includes parameters that were fixed in estimation:
#>   Theta
#> These have been replaced by their fixed values, giving:
#>   3.09688960243761*1/Phi
#> 
#> Running Delta method computation for user-defined function using robust standard errors
#> 
#>  Expression   Value   s.e. t-ratio (0)
#>         VoL  4.6538 0.1925       24.18
#>        VTAW -0.3150 0.1934       -1.63
#> INFORMATION: The results of the Delta method calculations are returned invisibly as
#>   an output from this function. Calling the function via
#>   result=apollo_deltaMethod(...) will save this output in an object
#>   called result (or otherwise named object).

z = 1.96
delta[,"CI-95%"] <- delta[,"Value"] - 1.96*delta[,"s.e."]
delta[,"CI+95%"] <- delta[,"Value"] + 1.96*delta[,"s.e."]
delta
#>   Expression   Value   s.e. t-ratio (0)    CI-95%   CI+95%
#> 1        VoL  4.6538 0.1925       24.18  4.276500 5.031100
#> 2       VTAW -0.3150 0.1934       -1.63 -0.694064 0.064064
```
