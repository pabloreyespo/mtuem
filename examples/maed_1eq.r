rm(list = ls())
library(apollo)
library(mtuem)
set.seed(42)

# =============================================================================
# Example 1: Single Equation (Work Only)
# =============================================================================
# The simplest model estimates only the work equation.
# No free times or free goods are specified.
# Theta is fixed to 1 for identification.
# =============================================================================

### Initialise
apollo_initialise()
apollo_control <- list(
  modelName = "maed-1eq",
  modelDescr = "Single equation MTUEM (work only)",
  indivID = "PeID",
  outputDirectory = "output"
)

### Load data
database <- mtuem::maed
database <- database[database$EcI > 0, ]

### Parameters: Theta fixed to 1 for identification
apollo_beta <- c(Theta = 1, Phi = 1, thw = 0)
apollo_fixed <- c("Theta")

apollo_inputs <- apollo_validateInputs()

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P <- list()

  work_elasticities <- list(Theta = Theta, Phi = Phi, thw = thw)

  mtuem_settings <- list(
    work_times = c("Tw"),
    work_elasticities = work_elasticities,
    Tc = Tc, Ec = EcI, w = w, tau = 168
  )

  P[["model"]] <- mtuem_likelihood(mtuem_settings, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

model <- apollo_estimate(apollo_beta, apollo_fixed,
                          apollo_probabilities, apollo_inputs, estimate_settings = est_set)
apollo_modelOutput(model)

### Predictions
pred <- apollo_prediction(model, apollo_probabilities, apollo_inputs)
head(pred)

### Values of Time
Tw <- pred$Tw
w  <- database$w
Tc <- database$Tc
Ec <- database$Ec

cteVoL  <- mean((w * Tw - Ec) / (168 - Tw - Tc))
cteVTAW <- mean((w * Tw - Ec) / Tw)

delta <- apollo_deltaMethod(model, list(
  expression = c(
    VoL  = paste0(cteVoL,  "*Theta/Phi"),
    VTAW = paste0(cteVTAW, "*thw/Phi")
  )
))
delta
