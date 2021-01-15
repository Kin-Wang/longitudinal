library(dplyr)
library(ggplot2)
library(nlme)
library(reshape2)
library(lmeInfo)
library(contrast)
library(AICcmodavg)


dental_ori <- read.csv("/Users/Kin/Downloads/Longitudinal/dentalwide.csv")
glimpse(dental_ori)

dental_long <- melt(dental_ori,id.vars=c("ID","Gender"),
                    measure.vars=c("Y1","Y2","Y3","Y4"),
                    variable.name="Age",
                    value.name="Distance")

levels(dental_long$Age)[levels(dental_long$Age)=="Y1"]='8'
levels(dental_long$Age)[levels(dental_long$Age)=="Y2"]='10'
levels(dental_long$Age)[levels(dental_long$Age)=="Y3"]='12'
levels(dental_long$Age)[levels(dental_long$Age)=="Y4"]='14'
dental_long$Gender <- as.factor(dental_long$Gender)
dental_long$Gender=relevel(dental_long$Gender,ref="M")

dental_model <- dental_long %>% 
  group_by(ID) %>% 
  mutate(Visit=row_number())

glimpse(dental_model)


#unstructured variance model

dental_unst <- gls(Distance~Gender+Age+Gender*Age,data=dental_model,corr=corSymm(form=~Visit|ID),
                   weights=varIdent(form=~1|Visit),
                   na.action=na.omit,method='REML')
summary(dental_unst)
unst_LL <- logLik(dental_unst)

#cs model

dental_cs <- gls(Distance~Gender+Age+Gender*Age,data=dental_model,corr=corCompSymm(form=~Visit|ID),
                 na.action=na.omit,method='REML')
summary(dental_cs)
cs_LL <- logLik(dental_cs)

#hcs model

dental_hcs <- gls(Distance~Gender+Age+Gender*Age,data=dental_model,corr=corCompSymm(form=~Visit|ID),
                  weights=varIdent(form=~1|Visit),
                  na.action=na.omit,method='REML')
summary(dental_hcs)
hcs_LL <- logLik(dental_hcs)

#AR1 model

dental_ar1 <- gls(Distance~Gender+Age+Gender*Age,data=dental_model,corr=corAR1(form=~Visit|ID),
                  na.action=na.omit,method='REML')
summary(dental_ar1)
ar1_LL <- logLik(dental_ar1)

#HAR1 model

dental_har1 <- gls(Distance~Gender+Age+Gender*Age,data=dental_model,corr=corAR1(form=~Visit|ID),
                   weights=varIdent(form=~1|Visit),
                   na.action=na.omit,method='REML')
summary(dental_har1)
har1_LL <- logLik(dental_har1)

#cs+ar1 model

dental_hyb <- lme(Distance~Gender+Age+Gender*Age,data=dental_model,random=~1|ID, 
                  corr=corAR1(form=~Visit|ID),
                  na.action=na.omit,method='REML')
summary(dental_hyb)
hyb_LL <- logLik(dental_hyb)


#-2*logLik
logLik_table <- tibble(Models="-2 Re-LL",unstructured=-2*unst_LL, CS=-2*cs_LL, HCS=-2*hcs_LL, 
                       AR1=-2*ar1_LL, HAR1=-2*har1_LL, CS_AR1=-2*hyb_LL)
logLik_table


#compare AR1 with unst
anova(dental_ar1, dental_unst)


#AIC BIC 
anv_table <- anova(dental_unst, dental_ar1, dental_cs, dental_har1, dental_hcs, dental_hyb)
anv_table
#dental_cs has the smallest AIC and it's a simple model



#a cs model with continuous age

dental_cs <- dental_model
dental_cs$Age <- as.numeric(dental_cs$Age)
dental_cs$Age[dental_cs$Age==1] <- 8
dental_cs$Age[dental_cs$Age==2] <- 10
dental_cs$Age[dental_cs$Age==3] <- 12
dental_cs$Age[dental_cs$Age==4] <- 14
dental_cs$Age
class(dental_cs$Age)

#cs model

cs_model <- gls(Distance~Gender+Age+Gender*Age,data=dental_cs,corr=corCompSymm(form=~Visit|ID),
                na.action=na.omit,method='REML')
summary(cs_model)
#we assume homogeneous variances between subjects
cov_matrix <- getVarCov(cs_model)
cov_matrix
corr_matrix <- cov2cor(cov_matrix)
corr_matrix

#F test, type III test for coefficients
anova(cs_model, type="marginal")

fitted_cs <- fitted(cs_model, glsFit)
fitted_cs
cs_model$sigma
cs_model$residuals

#use AICcmodavg with predictSE to produce fitted values and fittedSE
predict_cs <- predictSE(cs_model, newdata = dental_cs)
class(predict_cs)
predict_cs$fit

#create a new table with fitted values, fittedSE, CI
#df=104, t statistics=1.984
predict_table <- data.frame(ID=dental_cs$ID, Visit=dental_cs$Visit, Gender=dental_cs$Gender, 
                        Age=dental_cs$Age, Distance=dental_cs$Distance,
                        fitted=predict_cs$fit, fitted_se=predict_cs$se.fit,
                        lower=predict_cs$fit-1.984*predict_cs$se.fit,
                        higher=predict_cs$fit+1.984*predict_cs$se.fit)
head(predict_table, n=20L)

#error bar plot with age and fitted mean distance
ggplot(predict_table, aes(x=Age, y=fitted, group=Gender, color=Gender))+
  geom_errorbar(aes(ymin=lower,ymax=higher),width=0.1)+
  geom_point()+geom_line()+
  labs(title="Fitted mean of distances with S.E. over age")+
  xlab("Age / yrs-old")+ylab("Distances / mm")+
  theme(
    plot.title=element_text(hjust=0.5,color="blue",size=12,face="bold")
  )


#contrast different gender at age 14 using CS model
contrast(cs_model, list(Age=14, Gender="F"), list(Age=14, Gender="M"))

#compare two mean models
part_model <- gls(Distance~Age,data=dental_cs,corr=corCompSymm(form=~Visit|ID),
                  na.action=na.omit,method='ML')
full_model <- gls(Distance~Gender+Age+Gender*Age,data=dental_cs,corr=corCompSymm(form=~Visit|ID),
                  na.action=na.omit,method='ML')
anova(part_model, full_model)
#full model is better


#AR1 model with continuous age
AR1_model <- gls(Distance~Gender+Age+Gender*Age,data=dental_cs,corr=corAR1(form=~Visit|ID),
                na.action=na.omit,method='REML')
summary(AR1_model)
AR1_coef <- coef(AR1_model)
class(AR1_coef)
AR1_coef
#F test for coefficients
anova(AR1_model, type = "marginal")

# sandwich estimator
library(clubSandwich)
#robust covariance of betas
robcov.beta <- vcovCR(AR1_model,type='CR0',form='sandwich')
#robust standard error of betas
robse.beta <- sqrt(diag(robcov.beta))
robcov.beta
robse.beta

#T distribution
p_v <- 2*pt(-abs(AR1_coef/robse.beta), df=104)
coef_table <- data.frame(Value=AR1_coef,
                         SE=robse.beta, t=AR1_coef/robse.beta, p_value=p_v)
coef_table
#focus on GenderF







