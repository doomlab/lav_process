# Variables to lavaan model -----------------------------------------------
require(lavaan)

lavaan.model <- function(y,
                         x,
                         # Mediators
                         m,
                         # Categorical variables
                         categorical = NULL,
                         # Covariates
                         cv = NULL,
                         # Covariates for Y ~ X (and Y ~ M)
                         cv.yx = FALSE,
                         # Covariates for M ~ X
                         cv.mx = FALSE,
                         data,
                         # Test statistic
                         # Robust variants:
                         # "Satorra-Bentler" and "Yuan-Bentler"
                         test = "standard",
                         # Number of bootstrap resamples
                         B = 5000,
                         seed = 1234,
                         ci = .95,
                         process.model = 4,
                         print = TRUE) {
  # Set Seed
  set.seed(seed)

  # PROCESS Model 4 (Single Mediator)
  if (process.model == 4) {

    # Calculate SD of Y for partially standardized effect size
    sd.y <- sd(data[, y], na.rm = TRUE)

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

    # Mediator(s)
    # Check length of m
    if (length(m) == 1) {

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
      total.eff <- "total_effect := c + (a*b)"

      # Partially standardized indirect effect size
      a.es <- paste0("ab_es := ab/", sd.y)

      # Specify model
      model <- paste(
        direct.eff,
        mediator.mx,
        mediator.ym,
        indirect.eff,
        total.eff,
        a.es,
        sep = "\n"
      )
    } else {

      # Declare variables for loop
      mediator.mx <- vector("character", length = length(m))
      indirect.eff <- vector("character", length = length(m))
      mediator.ym <- vector("character", length = length(m))
      total.eff <- vector("character", length = length(m))
      indirect.es <- vector("character", length = length(m))

      # Loop through mediators
      for (i in 1:length(m)) {

        # Covariates have been selected for M ~ X
        if (!is.null(cv) & isTRUE(cv.mx)) {
          mediator.mx[i] <- paste0(
            m[i], "~", "a", i, "*", x, "+", paste0(cv, collapse = "+")
          )
        } else {
          mediator.mx[i] <- paste0(
            m[i], "~", "a", i, "*", x
          )
        }

        # Y ~ M
        mediator.ym[i] <- paste0(y, "~", "b", i, "*", m[i])

        # Indirect effect
        indirect.eff[i] <- paste0("a", i, "b", i, ":=", "a", i, "*", "b", i)

        # Total effect
        total.eff[i] <- paste0("(", "a", i, "*", "b", i, ")")

        # Partially standardized indirect effect size
        indirect.es[i] <- paste0(
          paste0("a", i, "b", i, "_es:=", "(a", i, "*", "b", i, ")/", sd.y)
        )
      }

      # Collapse variables
      mediator.mx <- paste0(mediator.mx, collapse = "\n")
      mediator.ym <- paste0(mediator.ym, collapse = "\n")
      indirect.eff <- paste0(indirect.eff, collapse = "\n")
      indirect.es <- paste0(indirect.es, collapse = "\n")

      # Total indirect effect
      total.ind.eff <- paste0(
        "total_indirect:=", paste0(total.eff, collapse = "+")
      )

      # Partially standardized total indirect effect
      total.ind.eff.es <- paste0(
        "total_indirect_es:=", "total_indirect/", sd.y
      )

      # Total direct effect
      total.eff <- paste0("total_effect:=c+", paste0(total.eff, collapse = "+"))

      # Specify model
      model <- paste(
        direct.eff,
        mediator.mx,
        mediator.ym,
        indirect.eff,
        total.ind.eff,
        total.eff,
        indirect.es,
        total.ind.eff.es,
        sep = "\n"
      )
    }

    # Run mediation analyses
    # Categorical variables
    if (!is.null(categorical)) {
      fit <- sem(
        model,
        ordered = categorical,
        data = data,
        # Requires DWLS for bootstrapping
        estimator = "DWLS",
        se = "robust.sem",
        bootstrap = B
      )
    } else {
      fit <- sem(
        model,
        data = data,
        se = "robust.sem",
        missing = "ML",
        test = test,
        bootstrap = B
      )
    }

    # Run bootstrap on coefficients
    # Required for replication to set (iseed)
    RNGkind("L'Ecuyer-CMRG")

    est_boot <- bootstrapLavaan(
      fit,
      FUN = function(x) parameterestimates(x)$est,
      iseed = seed,
      parallel = "multicore",
      R = B
    )

    boot_out <- data.frame(
      param = as.character(parameterestimates(fit)$lab),
      t(
        rbind(
          est = parameterestimates(fit)$est,
          se.boot = apply(est_boot, 2, sd),
          apply(est_boot, 2, quantile, c((1 - ci) / 2, ci + (1 - ci) / 2)),
          p.boot = 2 * pmin(colMeans(est_boot > 0), colMeans(est_boot < 0))
        )
      )
    )

    # Remove empty parameters
    boot_out <- boot_out[boot_out$param != "", ]

    # Set column names
    colnames(boot_out)[c(4:5)] <- c("ci.lb.boot", "ci.ub.boot", "p.value")

    # Set row names
    rownames(boot_out) <- boot_out$param

    # Drop param variable
    boot_out <- boot_out[, -1]

    # Return fit and bootstrap results

    # Show summary results
    if (print) {
      cat(
        "\n",
        "\n", "Summary Results",
        "\n",
        "\n",
        sep = ""
      )

      print(
        parameterEstimates(
          fit,
          zstat = TRUE,
          standardized = FALSE,
          cov.std = FALSE,
          rsquare = TRUE,
          output = "pretty"
        )
      )

      cat(
        "\n",
        "\n", "Bootstrap Results",
        "\n",
        "\n",
        sep = ""
      )

      print(boot_out)
    }

    return(
      list(
        fit = fit,
        boot = boot_out
      )
    )

  }
}
