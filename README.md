
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
can_afford_expenses =  (database["Ec"] / (database["w"]*(168 -  database["Tc"]))) < 1
database = database[is_not_retist & can_afford_expenses,]

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
#>  C:/Users/pablo/OneDrive/mtuem/enut-mtuem-1eq_iterations.csv
#>     it    nf     F            RELDF    PRELDF    RELDX    MODEL stppar
#>      0     1 1.956485669e+04
#>      1     4 1.896555286e+04 3.063e-02 2.535e-02 3.27e-02   G   1.32e+00
#>      2     5 1.825564290e+04 3.743e-02 2.748e-02 1.94e-01   G   2.16e-02
#>      3     6 1.821465491e+04 2.245e-03 4.373e-03 7.40e-02   S   0.00e+00
#>      4     8 1.820361475e+04 6.061e-04 1.520e-03 3.16e-02   S   1.58e-01
#>      5     9 1.819871195e+04 2.693e-04 2.570e-04 2.50e-02   S   0.00e+00
#>      6    10 1.819822035e+04 2.701e-05 1.642e-05 8.48e-03   S   0.00e+00
#>      7    11 1.819785709e+04 1.996e-05 2.030e-05 1.23e-02   G   0.00e+00
#>      8    12 1.819782760e+04 1.621e-06 1.143e-06 1.11e-03   G   0.00e+00
#>      9    13 1.819781989e+04 4.234e-07 2.757e-07 1.13e-03   G   0.00e+00
#>     10    14 1.819781669e+04 1.760e-07 1.758e-07 1.27e-03   S   0.00e+00
#>     11    15 1.819781669e+04 1.773e-11 1.678e-11 7.70e-06   S   0.00e+00
#> 
#> ***** Relative function convergence *****
#> 
#> Estimated parameters with approximate standard errors from BHHH matrix:
#>          Estimate     BHHH se BHH t-ratio (0)
#> Theta     1.00000          NA              NA
#> Phi       0.67534     0.01952          34.599
#> thw      -0.08793     0.02815          -3.124
#> 
#> Final LL: -18197.8167
#> 
#> Calculating log-likelihood at equal shares (LL(0)) for applicable
#>   models...
#> Calculating log-likelihood at observed shares from estimation data
#>   (LL(c)) for applicable models...
#> Calculating LL of each model component...
#> VoL: 4.741416 
#> VTAW: -0.3744373 
#> Computing covariance matrix using numerical methods (numDeriv).
#>  0%....25%....50%...100%
#> Negative definite Hessian with maximum eigenvalue: -575.920528
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
#> Model run at                                : 2025-06-18 22:59:33.076758
#> Estimation method                           : bgw
#> Model diagnosis                             : Relative function convergence
#> Optimisation diagnosis                      : Maximum found
#>      hessian properties                     : Negative definite
#>      maximum eigenvalue                     : -575.920528
#>      reciprocal of condition number         : 0.0144156
#> Number of individuals                       : 5037
#> Number of rows in database                  : 5037
#> Number of modelled outcomes                 : 0
#> 
#> Number of cores used                        :  1 
#> Model without mixing
#> 
#> LL(start)                                   : 0
#> LL at equal shares, LL(0)                   : NA
#> LL at observed shares, LL(C)                : NA
#> LL(final)                                   : -18197.82
#> Rho-squared vs equal shares                  :  Not applicable 
#> Adj.Rho-squared vs equal shares              :  Not applicable 
#> Rho-squared vs observed shares               :  Not applicable 
#> Adj.Rho-squared vs observed shares           :  Not applicable 
#> AIC                                         :  36399.63 
#> BIC                                         :  -Inf 
#> 
#> Estimated parameters                        : 2
#> Time taken (hh:mm:ss)                       :  00:00:0.82 
#>      pre-estimation                         :  00:00:0.43 
#>      estimation                             :  00:00:0.31 
#>      post-estimation                        :  00:00:0.08 
#> Iterations                                  :  11  
#> 
#> Unconstrained optimisation.
#> 
#> Estimates:
#>          Estimate        s.e.   t.rat.(0)    Rob.s.e. Rob.t.rat.(0)
#> Theta     1.00000          NA          NA          NA            NA
#> Phi       0.67534     0.02197      30.745     0.02726        24.774
#> thw      -0.08793     0.03576      -2.459     0.04894        -1.796

# Compute values of time
pred <- apollo_prediction(model, apollo_probabilities, apollo_inputs)
#> Running predictions from model using parameter estimates...
#> Prediction at user provided parameters
#>                 Tw         T         X      VoL     VTAW
#> Aggregate 211574.5 174125.93 570818.34 23882.51 -1886.04
#> Average       42.0     34.57    113.33     4.74    -0.37
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
#>   3.2020476626436*1/Phi
#> 
#> Running Delta method computation for user-defined function using robust standard errors
#> 
#>  Expression   Value   s.e. t-ratio (0)
#>         VoL  4.7414 0.1914       24.77
#>        VTAW -0.3744 0.1937       -1.93
#> INFORMATION: The results of the Delta method calculations are returned invisibly as
#>   an output from this function. Calling the function via
#>   result=apollo_deltaMethod(...) will save this output in an object
#>   called result (or otherwise named object).

z = 1.96
delta[,"CI-95%"] <- delta[,"Value"] - 1.96*delta[,"s.e."]
delta[,"CI+95%"] <- delta[,"Value"] + 1.96*delta[,"s.e."]
delta
#>   Expression   Value   s.e. t-ratio (0)    CI-95%   CI+95%
#> 1        VoL  4.7414 0.1914       24.77  4.366256 5.116544
#> 2       VTAW -0.3744 0.1937       -1.93 -0.754052 0.005252
```
