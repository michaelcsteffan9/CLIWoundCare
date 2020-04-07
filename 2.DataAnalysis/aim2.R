
#########################################################
## c-TMLE
#########################################################

fitCTMLE <- function(study.data, outcome.name, cov.names){
  dat.nomiss <- study.data[complete.cases(study.data),]
  cov.df <- dat.nomiss %>% select(all_of(cov.names))
  out.sym <- sym(outcome.name)
  ct.fit <- ctmleDiscrete(Y = dat.nomiss %>% pull(!!out.sym), 
                          A = dat.nomiss$primaryTreatmentInt, 
                          W = cov.df, detailed = TRUE)
  return(ct.fit)
}

#########################################################
## AIPWE
#########################################################

#########################################################
## Diagnostics
#########################################################
