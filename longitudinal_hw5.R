library(dplyr)
library(haven)
library(MASS)
library(lme4)
library(gee)
library(geepack)
library(ordinal)
library(pbkrtest)


##input data
toenail <- read.csv("/Users/Kin/Downloads/Longitudinal/toenail_month.csv")
glimpse(toenail)

toenail <- toenail %>%
       group_by(ID) %>%
       mutate(visit=row_number())
glimpse(toenail)
anyNA(toenail)

toenail$TRT <- as.factor(toenail$TRT)
levels(toenail$TRT)

##gee model for binary outcome

toe_model1 <- geese(Y~TRT+Month+TRT:Month, id=ID, data=toenail, corstr = "exch", family = binomial)
summary(toe_model1)
toe_model1$alpha

##gee model with random intercept for binary outcome

toe_model2 <- glmer(Y~TRT+Month+TRT:Month+(1|ID), data=toenail, family = binomial)
tt <- summary(toe_model2)
dd <- tt$coefficients
kk <- 1.96*dd[,2]
kk
lower_cf <- dd[,1]-kk
upper_cf <- dd[,1]+kk
exp(dd[,1])
exp(lower_cf)
exp(upper_cf)
tt$sigma

#test random intercept
toe_model3 <- glm(Y~TRT+Month+TRT:Month, data=toenail, family = binomial)
G2 <- -2*logLik(toe_model3)+2*logLik(toe_model2)
pchisq(as.numeric(G2), df=1, lower.tail = F)

