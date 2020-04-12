
#########################################################
## c-TMLE
#########################################################

#' Wrapper function for fitting CTMLE using the 
#' 
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

extractCTMLEFinalCovar <- function(ctmle.fit){
  x <- summary(ctmle.fit)
  # This code taken from print.summary.ctmle with minor modication
  # For some reason CTMLE doesn't seem to have an easy way to let you see the covariates in the selected model
  # NF: I believe that you can use: ctmle.fit$candidates$terms[1:ctmle.fit$best_k]
  npercovar <- table(x$covar[1:x$ncand])
  suffix <- paste(x$covar[1:x$ncand], letters[unlist(apply(npercovar, 1,function(x){1:x}))], sep="")
  suffix[cumsum(npercovar)[npercovar==1]] <- names(cumsum(npercovar)[npercovar==1])
  prev_covar <-suffix[cumsum(npercovar[-length(npercovar)])]
  prev_moves <- c(paste(" + epsilon", prev_covar, " * h", prev_covar, sep=""),"")
  current_move <- paste(" + epsilon", suffix, " * h", suffix,sep="")
  TMLEcand <- paste("\tcand ", 1:x$ncand,": Q",1:x$ncand,"(A,W) =", sep="")
  fluctuations <- paste("Q0(A,W)", c(rep("", npercovar[1]),
                                     mapply(function(x){paste(prev_moves[1:(x-1)],collapse="")},x$covar[-(1:npercovar[1])])),
                        current_move, sep="")
  final_update <- c(rep("\n", npercovar[1]), paste("\t                = Q",
                                                   rep(cumsum(npercovar[-length(npercovar)]), times=npercovar[-1]),"(A,W)",
                                                   current_move[(npercovar[1]+1) : x$ncand],",\n", sep=""))
  tx <- c("Intercept only",
          paste(
                mapply(function(x1){paste(x$terms[1:x1], collapse=", ")}, 1:(x$ncand-1)), sep=""))
  return(tx[x$selected])
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
  
  nonzero.var.names <- (rownames(main.coef)[main.coef != 0])[-1]
  nonzero.main.var.names <- nonzero.var.names[str_detect(nonzero.var.names, pattern = ":", negate = TRUE)]
  nonzero.interaction.var.names <- nonzero.var.names[str_detect(nonzero.var.names, pattern = ":")] %>% 
    strsplit(., split = ":") %>% unlist()
  nonzero.interaction.var.names <- nonzero.interaction.var.names[str_detect(nonzero.interaction.var.names, pattern = "primaryTreatmentInt", negate = TRUE)]
  
  ## Construct the component models
  
  prop.mod.formula <- formula(paste(c("~ 1", nonzero.prop.var.names), collapse = " + "))
  main.mod.string <- paste(c("~ 1", nonzero.main.var.names), sep = "+", collapse = " + ")
  main.mod.formula <- formula(main.mod.string)
  cont.mod.string <- paste(c("~ 1", nonzero.interaction.var.names), sep = "+", collapse = " + ")
  
  #print(main.mod.string)
  #print(cont.mod.string)
  
  if(cont.mod.string == "~ 1"){
    # This is kind of goofy but DynTxRegime doesn't seem to like it if Contmod is intercept only and throws an error
    # Solution: take the first nonzero term(s)
    main.lasso <- glmnet(x = main.candidate.mat, y = study.data$ampFreeSurv2yr, family="binomial")
    # Get the coefficient path for the interaction variables
    int.vars <- main.lasso$beta[str_detect(rownames(main.lasso$beta), ":"),]
    cur.index <- 1
    cont.mod.string <- NULL
    while(is.null(cont.mod.string)){
      if(all(int.vars[,cur.index] == 0)){
        cur.index <- cur.index + 1
        next()
      } else{
        nonzero.interaction.var.names <- rownames(int.vars)[int.vars[,cur.index] != 0] %>% 
          strsplit(., split = ":") %>% unlist()
        nonzero.interaction.var.names <- nonzero.interaction.var.names[str_detect(nonzero.interaction.var.names, pattern = "primaryTreatmentInt", negate = TRUE)]
        cont.mod.string <- paste(c("~ 1", nonzero.interaction.var.names), sep = "+", collapse = " + ")
      }
    }
  }
  cont.mod.formula <- formula(cont.mod.string)
  
  mean.mod <- buildModelObj(model = main.mod.formula, 
                            solver.method = "glm", 
                            solver.args = list(family = "binomial"),
                            predict.method = "predict.glm",
                            predict.args = list(type = "response"))
  cont.mod <- buildModelObj(model = cont.mod.formula, 
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
                           moCont = cont.mod,
                           data = study.df, 
                           response = study.data %>% pull(!!out.sym),  
                           txName = "primaryTreatmentInt",
                           regime =cont.mod.formula,
                           iter = 0L,
                           verbose = FALSE)
  return(aipw.fit)
}

calcAIPWate <- function(in.data, outcome.name, cov.names, ...){
  aipw.fit <- fitAIPWE(study.data = in.data, outcome.name = outcome.name, cov.names = cov.names, ...)
  aipw.outcome.mod <- outcome(aipw.fit)$Combined
  # Make counterfactual datasets
  counterfac.data.all.wm <- in.data %>% mutate(primaryTreatmentInt = 1) %>% as.data.frame()
  counterfac.data.all.revasc <- in.data %>% mutate(primaryTreatmentInt = 0) %>% as.data.frame()
  
  #Get predicted outcomes and get predicted probability
  pred.prob.wm <- boot::inv.logit(predict(aipw.outcome.mod, newdata = counterfac.data.all.wm))
  pred.prob.revasc <- boot::inv.logit(predict(aipw.outcome.mod, newdata = counterfac.data.all.revasc))
  
  ate <- mean(pred.prob.wm - pred.prob.revasc)
  return(ate)
}

aipwBootWrapper <- function(in.data, n.boot, outcome.name, cov.names, ...){
  ate.vec <- vector(length = n.boot)
  for(i in 1:n.boot){
    boot.indices <- sample(1:nrow(in.data), size = nrow(in.data), replace = TRUE)
    cur.boot.data <- clti[boot.indices,]
    ate.vec[i] <- calcAIPWate(in.data = cur.boot.data, outcome.name = outcome.name, cov.names = cov.names)
  }
  return(ate.vec)
}

#########################################################
## Diagnostics
#########################################################
