# Variables to lavaan model -----------------------------------------------
require(lavaan)

lavaan.model <- function(y,
                         x,
                         m,
                         # Covariates
                         cv = NULL,
                         # Covariates for Y ~ X (and Y ~ M)
                         cv.yx = FALSE,
                         # Covariates for M ~ X
                         cv.mx = FALSE,
                         data,
                         # Number of boostrap resamples
                         B = 5000,
                         process.model = 4) {

  # PROCESS Model 4 (Single Mediator)
  if (process.model == 4) {

    # Direct effect
    # Y ~ X
    direct.eff <- ifelse(
      # Covariates for Y ~ X and Y ~ M automatically:
      !is.null(cv) & isTRUE(cv.yx),
      # Y ~ c*X + b*M + CV1 + ... CVn
      paste0(y, "~", "c*", x, "+", paste0(cv, collapse = "+")),
      # Y ~ c*X
      paste0(y, "~", "c*", x)
    )
    # Mediator
    # M ~ X
    mediator.mx <- ifelse(
      # Covariates for M ~ X
      !is.null(cv) & isTRUE(cv.mx),
      # M ~ a*X + CV1 + ... CVn
      paste0(m, "~", "a*", x, "+", paste0(cv, collapse = "+")),
      # M ~ a*X
      paste0(m, "~", "a*", x)
    )

    # Y ~ M
    mediator.ym <- paste0(y, "~", "b*", m)

    # Indirect effect
    indirect.eff <- "ab := a*b"
    # Total effect
    total.eff <- "total := c + (a*b)"

    # Specify model
    model <- paste(
      direct.eff,
      mediator.mx,
      mediator.ym,
      indirect.eff,
      total.eff,
      sep = "\n"
    )

    # Run mediation analyses
    fit <- sem(model, data = data, se = "bootstrap", bootstrap = B)

  }

  # Output
  cat(
    "\n",
    "\n", "PROCESS Model Specifed: ", process.model,
    "\n",
    "\n",
    model,
    "\n",
    "\n",
    sep ="")

  # Return results
  fit
}
