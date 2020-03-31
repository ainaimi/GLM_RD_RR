#' Load Libraries
packages <- c("broom","here","tidyverse","skimr","rlang","sandwich","boot","xtable")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,repos='http://lib.stat.cmu.edu/R/CRAN') 
  }
}

for (package in packages) {
  library(package, character.only=T)
}

#' Define Plot Format
thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.title=element_blank(),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

#' Define the expit (inverse logit) function 
expit <- function(x){1/(1+exp(-x))}

#' Import, look at, and format data
#' This code downloads the data from the website for the book by Hernan and 
#' Robins "Causal Inference." It requires internet access, and will not work if 
#' the link changes.
file_loc <- url("https://cdn1.sph.harvard.edu/wp-content/uploads/sites/1268/1268/20/nhefs.csv")

#' This begins the process of cleaning and formatting the data
nhefs <- read_csv(file_loc) %>% 
  select(qsmk,wt82_71,exercise,sex,age,race,income,marital,school,asthma,bronch) %>% 
  mutate(income=as.numeric(income>15),
         marital=as.numeric(marital>2)) %>% 
  na.omit()

factor_names <- c("exercise","income","marital","sex","race","asthma","bronch")
nhefs[,factor_names] <- lapply(nhefs[,factor_names] , factor)

#' Define outcome
nhefs <- nhefs %>% mutate(wt_delta = as.numeric(wt82_71>median(wt82_71)))
nhefs

#' Quick summary of data
skim(nhefs)

#' Here, we start fitting relevant regression models to the data.
#' modelForm is a regression argument that one can use to regress the 
#' outcome (wt_delta) against the exposure (qsmk) and selected confounders.

formulaVars <- paste(names(nhefs)[-c(2,12)],collapse = "+")
modelForm <- as.formula(paste0("wt_delta ~", formulaVars))
modelForm

#' This model can be used to quantify a conditionally adjusted odds ratio with correct standard error
modelOR <- glm(modelForm,data=nhefs,family = binomial("logit"),control = list(trace = TRUE))
tidy(modelOR)[2,]

#' This model can be used to quantify a conditionally adjusted risk ratio with with correct standard error
#' However, error it returns an error and thus does not provide any results.
modelRR_binom <- glm(modelForm,data=nhefs,family = binomial("log"))

#' This model can be used to quantify a conditionally adjusted risk difference with correct standard error
modelRD_binom <- glm(modelForm,data=nhefs,family = binomial("identity"))

#' This model can be used to quantify a conditionally risk ratio using the Poisson distribuiton and log link function. 
#' However, because the Poisson distribution is used, the model provides incorrect standard error estimates.
modelRR <- glm(modelForm,data=nhefs,family = poisson("log"))
tidy(modelRR)[2,]
#' To obtain the correct variance, we use the "sandwich" function to obtain correct sandwich (robust) standard error estimates.
sqrt(sandwich(modelRR)[2,2])

#' This model can be used to obtain a risk difference with the gaussian distribiton or using ordinary least 
#' squares (OLS, via the lm function). Again, the model based standard error estimates are incorrect. 
modelRD <- glm(modelForm,data=nhefs,family = gaussian("identity"))
modelRD <- lm(modelForm,data=nhefs)
tidy(modelRD)[2,]
#' To obtain the correct variance, we use the "sandwich" function to obtain correct sandwich (robust) standard error estimates.
sqrt(sandwich(modelRD)[2,2])

#' This model shows some potential evidence for interaction between smoking and exercise on the risk difference scale.
formulaVars <- paste(names(nhefs)[-c(2,12)],collapse = "+")
modelForm <- as.formula(paste0("wt_delta ~", formulaVars,"+qsmk*exercise"))
modelForm
summary(lm(modelForm,data=nhefs))

#' Marginal Standardization
##' To avoid assuming no interaction between smoking and any of the other variables
##' in the model, we subset modeling among exposed/unexposed. This code removes smoking from the model,
##' which will allow us to regress the outcome against the confounders among the exposed and 
##' the unexposed searately. Doing so will allow us to account for any potential exposure-covariate interactions
##' that may be present. 
formulaVars <- paste(names(nhefs)[-c(1,2,12)],collapse = "+")
modelForm <- as.formula(paste0("wt_delta ~", formulaVars))
modelForm

#' Regress the outcome against the confounders among the unexposed (model0) and then among the exposed (model1)
model0 <- glm(modelForm,data=subset(nhefs,qsmk==0),family=binomial("logit"))
model1 <- glm(modelForm,data=subset(nhefs,qsmk==1),family=binomial("logit"))
##' Generate predictions for everyone in the sample using the model fit to only the 
##' unexposed (mu0 predictions) and only the exposed (mu1 predictions).
mu1 <- predict(model1,newdata=nhefs,type="response")
mu0 <- predict(model0,newdata=nhefs,type="response")

#' Marginally adjusted risk ratio
marg_stand_RR <- mean(mu1)/mean(mu0)
#' Marginally adjusted risk difference
marg_stand_RD <- mean(mu1)-mean(mu0)

#' Using the bootstrap to obtain confidence intervals for the marginally adjusted 
#' risk ratio and risk difference.
bootfunc <- function(data,index){
  boot_dat <- data[index,]
  model0 <- glm(modelForm,data=subset(boot_dat,qsmk==0),family=binomial("logit"))
  model1 <- glm(modelForm,data=subset(boot_dat,qsmk==1),family=binomial("logit"))
  mu1 <- predict(model1,newdata=boot_dat,type="response")
  mu0 <- predict(model0,newdata=boot_dat,type="response")
  
  marg_stand_RR <- mean(mu1)/mean(mu0)
  marg_stand_RD <- mean(mu1)-mean(mu0)
  res <- c(marg_stand_RD,marg_stand_RR)
  return(res)
}

#' Run the boot function. Set a seed to obtain reproducibility
set.seed(123)
boot_res <- boot(nhefs,bootfunc,R=2000)

boot_RD <- boot.ci(boot_res,index=1)
boot_RR <- boot.ci(boot_res,index=2)

#' Creating Table 1 presented in the manuscript
t11<-as.matrix(tidy(modelRD)[2,2])
t12<-as.matrix(tidy(modelRD)[2,2] - 1.96*sqrt(sandwich(modelRD)[2,2]))
t13<-as.matrix(tidy(modelRD)[2,2] + 1.96*sqrt(sandwich(modelRD)[2,2]))

t21<-as.matrix(exp(tidy(modelRR)[2,2]))
t22<-as.matrix(exp(tidy(modelRR)[2,2] - 1.96*sqrt(sandwich(modelRR)[2,2])))
t23<-as.matrix(exp(tidy(modelRR)[2,2] + 1.96*sqrt(sandwich(modelRR)[2,2])))

t31<-as.matrix(marg_stand_RD)
t32<-as.matrix(boot_RD$bca[4])
t33<-as.matrix(boot_RD$bca[5])

t41<-as.matrix(marg_stand_RR)
t42<-as.matrix(boot_RR$bca[4])
t43<-as.matrix(boot_RR$bca[5])

tab <- data.frame(method=c("GLM","Marginal Standardization"),rd=c(t11,t31),rd.ll=c(t12,t32),rd.ul=c(t13,t33),
                  rr=c(t21,t41),rr.ll=c(t22,t42),rr.ul=c(t23,t43))
xtable(tab)