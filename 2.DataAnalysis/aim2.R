
#########################################################
## c-TMLE
#########################################################

fitCTMLE <- function(study.data, outcome.name, cov.names, ...){
  dat.nomiss <- study.data[complete.cases(study.data),]
  cov.df <- dat.nomiss %>% select(all_of(cov.names))
  out.sym <- sym(outcome.name)
  ct.fit <- ctmleDiscrete(Y = dat.nomiss %>% pull(!!out.sym), 
                          A = dat.nomiss %>% pull(primaryTreatmentInt), 
                          W = cov.df, 
                          family='binomial',
                          detailed = TRUE, ...)
  return(ct.fit)
}

#########################################################
## AIPWE
#########################################################

# Use CV fit from Glmnet to pick the covariates for the propensity score model and the outcome model
fitAIPWE <- function(study.data, outcome.name, cov.names, lambda.to.use = "lambda.min", ...){
  ## Use LASSO to select varaibles for the propensity score model
  propensity.candidate.mat <- study.data %>% select(all_of(cov.names)) %>% data.matrix()
  prop.cv <- cv.glmnet(x = propensity.candidate.mat, 
                       y = study.data$primaryTreatmentInt, 
                       family = "binomial", 
                       nfolds = 5)
  prop.coef <- coef(prop.cv, s =lambda.to.use) %>% as.matrix()
  # Get the nonzero var names, remove the intercept
  nonzero.prop.var.names <- (rownames(prop.coef)[prop.coef != 0])[-1]
  
  ## Use LASSO to select varaibles for the outcome model
  main.candidate.mat <- cbind(propensity.candidate.mat, 
                              "primaryTreatmentInt" = study.data$primaryTreatmentInt)
  main.candidate.mat <- model.matrix(~ . + .*primaryTreatmentInt, data.frame(main.candidate.mat))
  main.cv <- cv.glmnet(x = main.candidate.mat, y = study.data$ampFreeSurv2yr, family="binomial")
  main.coef <- coef(main.cv, s = lambda.to.use) %>% as.matrix()
  nonzero.main.var.names <- (rownames(main.coef)[main.coef != 0])[-1]
  
  
  
  ## Construct the component models
  
  prop.mod.formula <- formula(paste(c("~ 1", nonzero.prop.var.names), collapse = " + "))
  main.mod.formula <- formula(paste(c("~ 1", nonzero.main.var.names), sep = "+", collapse = " + "))
  
  mean.mod <- buildModelObj(model = main.mod.formula, 
                            solver.method = "glm", 
                            solver.args = list(family = "binomial"),
                            predict.method = "predict.glm",
                            predict.args = list(type = "response"))
  cont.mod <- buildModelObj(model = ~ 1 + SMOKING + anyAnemia, 
                            solver.method = "glm", 
                            solver.args = list(family = "binomial"),
                            predict.method = "predict.glm",
                            predict.args = list(type = "response"))
  prop.mod <- buildModelObj(model = prop.mod.formula, 
                            solver.method = "glm", 
                            solver.args = list(family = "binomial"),
                            predict.method = "predict.glm",
                            predict.args = list(type = "response"))
  
  ## Run AIPW
  
  #DynTxRegime doesn't place nicely with tibbles sometimes
  study.df <- as.data.frame(study.data)
  out.sym <- sym(outcome.name)
  
  aipw.fit <- earl(moPropen = prop.mod, 
                           moMain = mean.mod,
                           moCont = mean.mod,
                           data = study.df, 
                           response = study.data %>% pull(!!out.sym),  
                           txName = "primaryTreatmentInt",
                           regime =main.mod.formula,
                           iter = 0L,
                           verbose = FALSE)
  return(aipw.fit)
}


#########################################################
## Diagnostics
#########################################################
