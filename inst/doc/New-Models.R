## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(bootGOF)
set.seed(1)
X = runif(n = 200, min = 6, max = 14)
d = data.frame(x = X, y = sin(0.5 * X) + rnorm(200, sd = 0.2))
plot(y~x, data = d)

## -----------------------------------------------------------------------------
wrong_model = lm(y ~ I(x^2), data = d)

## -----------------------------------------------------------------------------
mt <- GOF_model(
  data = d,
  model = wrong_model,
  simulator_type = "parametric",
  nmb_boot_samples = 100,
  y_name = "y",
  Rn1_statistic = Rn1_KS$new())
mt$get_pvalue()

## -----------------------------------------------------------------------------
library(minpack.lm)
fit <- minpack.lm::nlsLM(y ~ sin(a * x),
  data = d,
  start = c(a = 0.5),
  control = nls.control(maxiter = 500))
fit

## -----------------------------------------------------------------------------
fit_and_data <- list(fit = fit, data = d)

## -----------------------------------------------------------------------------
library(R6)
my_nls_info_extractor <- R6::R6Class(
  classname = "my_nls_info_extractor",
  inherit = GOF_model_info_extractor,
  public = list(
    yhat = function(model) {
      predict(object = model$fit)
    },
    y_minus_yhat = function(model) {
      residuals(object = model$fit)
    },
    beta_x_covariates = function(model) {
      a_hat <- coef(object = model$fit)
      x <- model$data$x
      ret <- a_hat * x
      return(ret)
    }
  ))
my_info_extractor <- my_nls_info_extractor$new()

## -----------------------------------------------------------------------------
head(d)

## -----------------------------------------------------------------------------
head(my_info_extractor$yhat(model = fit_and_data))

## -----------------------------------------------------------------------------
head(my_info_extractor$y_minus_yhat(model = fit_and_data))

## -----------------------------------------------------------------------------
head(my_info_extractor$beta_x_covariates(model = fit_and_data))

## -----------------------------------------------------------------------------
my_simulator <- GOF_sim_wild_rademacher$new(
  gof_model_info_extractor = my_info_extractor
)

## -----------------------------------------------------------------------------
head(d)

## -----------------------------------------------------------------------------
head(my_simulator$resample_y(model = fit_and_data))

## -----------------------------------------------------------------------------
my_nls_trainer <- R6::R6Class(
  classname = "GOF_nls_trainer",
  inherit = GOF_model_trainer,
  public = list(
    refit = function(model, data) {
      fit <- update(object = model$fit, data = data)
      ret <- list(fit = fit, data = data)
      return(ret)
    }))
my_trainer <- my_nls_trainer$new()

## -----------------------------------------------------------------------------
new_data <- d
new_data$y <- my_simulator$resample_y(model = fit_and_data)
my_trainer$refit(model = fit_and_data, data = new_data)$fit

## -----------------------------------------------------------------------------
set.seed(1)
mt <- GOF_model_test$new(
  model = fit_and_data,
  data = d,
  nmb_boot_samples = 100,
  y_name = "y",
  Rn1_statistic = Rn1_CvM$new(),
  gof_model_info_extractor = my_info_extractor,
  gof_model_resample = GOF_model_resample$new(
    gof_model_simulator = my_simulator,
    gof_model_trainer = my_trainer
  )
)
mt$get_pvalue()

## -----------------------------------------------------------------------------
library(MASS)
set.seed(1)
X1 <- rnorm(100)
X2 <- rnorm(100)
d <- data.frame(
  y = MASS::rnegbin(n = 100, mu = exp(0.2 + X1 * 0.2 + X2 * 0.6), theta = 2),
  x1 = X1,
  x2 = X2)
fit <- MASS::glm.nb(y~x1+x2, data = d)
fit

## -----------------------------------------------------------------------------
my_negbin_trainer <- R6::R6Class(
  classname = "GOF_glmnb_trainer",
  inherit = GOF_model_trainer,
  public = list(
    refit = function(model, data) {
      MASS::glm.nb(formula = formula(model), data = data)
    }))

## -----------------------------------------------------------------------------
set.seed(1)
mt <- GOF_model_test$new(
  model = fit,
  data = d,
  nmb_boot_samples = 100,
  y_name = "y",
  Rn1_statistic = Rn1_CvM$new(),
  gof_model_info_extractor = GOF_glm_info_extractor$new(),
  gof_model_resample = GOF_model_resample$new(
    gof_model_simulator = GOF_glm_sim_param$new(),
    gof_model_trainer = my_negbin_trainer$new()
  )
)
mt$get_pvalue()

## -----------------------------------------------------------------------------
set.seed(1)
mt2 <- GOF_model(
  model = fit,
  data = d,
  nmb_boot_samples = 100,
  simulator_type = "parametric",
  y_name = "y",
  Rn1_statistic = Rn1_CvM$new()
)
mt2$get_pvalue()

