library(dplyr)
library(haven)
library(MASS)
library(reshape2)
library(gee)
library(geepack)
library(contrast)


resp <- read.csv("/Users/Kin/Downloads/Longitudinal/respir.csv", header = TRUE)
glimpse(resp)

resp_long <- melt(resp, id.vars = c("ID", "Clinic", "Treatment"), 
                  measure.vars = c("Y0", "Y1", "Y2", "Y3", "Y4"),
                  variable.name = "Time", 
                  value.name = "Respiratory")
glimpse(resp_long)
resp_long$Clinic <- as.factor(resp_long$Clinic)
resp_long$Treatment <- as.factor(resp_long$Treatment)
resp_long$Treatment <- relevel(resp_long$Treatment, ref="P")
levels(resp_long$Treatment)
levels(resp_long$Time)

anyNA(resp_long)    #No missing; balanced design;
table(resp_long$Treatment, useNA = "ifany")
table(resp_long$Clinic, useNA = "ifany")
table(resp_long$Respiratory, useNA = "ifany")
table(resp_long$Time, useNA = "ifany")




#binary outcome, longistic, GEE
resp_1 <- gee(Respiratory~Treatment+Time+Treatment:Time, id=ID, data = resp_long,
              corstr = "unstructured", family = binomial)
summary(resp_1)
exp(coef(resp_1))

resp_2 <- gee(Respiratory~Treatment+Time+Treatment:Time, id=ID, data = resp_long,
              corstr = "exchangeable", family = binomial)
summary(resp_2)
exp(coef(resp_2))

#two cov structures are really similar
#we use CS structure

resp_3 <- geeglm(Respiratory~Treatment+Time+Treatment:Time, id=ID, data = resp_long,
                corstr = "exch", family = binomial)
summary(resp_3)

resp_4 <- geese(Respiratory~Treatment+Time+Treatment:Time, id=ID, data = resp_long,
                 corstr = "exch", family = binomial)
summary(resp_4)

coef(resp_3)

#OR for Y3
logor_y3 <- contrast(resp_4, a=list(Treatment="A", Time="Y3"), b=list(Treatment="P", Time="Y3"))
exp(logor_y3$Contrast)    #OR
exp(logor_y3$Lower)
exp(logor_y3$Upper)

anova(resp_3)



#add Clinic in the model

resp_5 <- geeglm(Respiratory~Clinic+Treatment+Time+Clinic*Treatment*Time, id=ID, data = resp_long,
                corstr = "exch", family = binomial)
summary(resp_5)
coef(resp_5)
anova(resp_5)

#drop three terms
resp_6 <- geeglm(Respiratory~Clinic+Treatment+Time+Treatment*Time, id=ID, data = resp_long,
                 corstr = "exch", family = binomial)
summary(resp_6)
anova(resp_6)

#drop interaction
resp_7 <- geeglm(Respiratory~Clinic+Treatment+Time, id=ID, data = resp_long,
                 corstr = "exch", family = binomial)
summary(resp_7)
anova(resp_7)
exp(coef(resp_7))

#extract fitted values
##fitted values=probability
fitted <- resp_7$fitted.values
class(fitted)
fitted_v <- as.data.frame(fitted)
head(fitted_v)

##merge 
fitted_frame <- data.frame(ID=resp_long$ID, Clinic=resp_long$Clinic, 
                           Treatment=resp_long$Treatment, Time=resp_long$Time,
                           Respiratory=resp_long$Respiratory, Fitted=fitted_v$V1)
glimpse(fitted_frame)
summary(fitted_frame$Fitted)


#create a panel variable for plotting
fitted_frame$panel[fitted_frame$Clinic=="1" & fitted_frame$Treatment=="A"] <- 1
fitted_frame$panel[fitted_frame$Clinic=="1" & fitted_frame$Treatment=="P"] <- 2
fitted_frame$panel[fitted_frame$Clinic=="2" & fitted_frame$Treatment=="A"] <- 3
fitted_frame$panel[fitted_frame$Clinic=="2" & fitted_frame$Treatment=="P"] <- 4

library(ggplot2)
#panel plot
ggplot(data = fitted_frame, aes(x=Time, y=Fitted))+
  geom_jitter(width = 0.1, height = 0.02)+ 
  facet_wrap(~panel, labeller = labeller(panel=c("1"="Clinic 1 & Treatment A", "2"="Clinic 1 & Treatment P", "3"="Clinic 2 & Treatment A", "4"="Clinic 2 Treatment P")))


##insomnia data
inso <- read.csv("/Users/Kin/Downloads/Longitudinal/insomnia.csv", header = TRUE)
glimpse(inso)


library(repolr)

inso_model1 <- repolr(Y~Trt+Time+Trt:Time, data=inso, 
                      categories=4, subjects="ID", times=c(1,2),
                      corr.mod="independence",alpha=0.5)
summary(inso_model1)
exp(coef(inso_model1))

#test estimation in repolr whether different from polr
test_model <- polr(factor(Y) ~ factor(Trt)+factor(Time)+factor(Trt):factor(Time),data=inso)
summary(test_model)  #betas are negative



