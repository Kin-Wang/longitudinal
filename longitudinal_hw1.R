library(tidyverse)
library(haven)
library(dplyr)
library(reshape2)
library(ggplot2)


hw1 <- read_sas("topekahw1.sas7bdat")
glimpse(hw1)

#Question 1 a, scatter plot of log transformed FEV1 and age with loess curve
ggplot(hw1, aes(x=AGE, y=LOGFEV1, group=ID))+
  geom_point()+
  stat_smooth(aes(group=1), method="loess")+
  xlab("Age")+ylab("Log FEV1")

#Question 1 b, scatter plot of residuals of height and residuals of FEV1
one_b <- loess(loght~AGE,data=hw1,span=0.5)
resid_1b <- resid(one_b)
newdata <- data.frame(hw1,resid_1b)
one_bb <- loess(LOGFEV1~AGE, data=hw1, span=0.5)
resid_1bb <- resid(one_bb)
new_data <- data.frame(newdata, resid_1bb)

ggplot(data=new_data,aes(x=resid_1b,y=resid_1bb))+geom_point()+
  xlab("Residuals of height")+ylab("Residuals of FEV1")


#Question 1 c, baseline effect VS longitudinal effect
change <- new_data %>%
  group_by(ID) %>%
  mutate(ht.reschange=resid_1b-first(resid_1b),
         FEV1.reschange=resid_1bb-first(resid_1bb))

# remove baseline rows
nobase <- change %>%
  group_by(ID) %>%
  filter(row_number() !=1)

# baseline only groups
baseonly <- change %>%
  group_by(ID) %>%
  filter(row_number() ==1)

# plot for cross-sectional study
ggplot(baseonly,aes(x=resid_1b,y=resid_1bb,group=ID))+geom_point()+
  xlab("Baseline Log Height residuals")+ylab("Baseline Log FEV1 residuals")+
  stat_smooth(aes(group=1),method="loess")

# plot for longitudinal effect
ggplot(nobase,aes(x=ht.reschange,y=FEV1.reschange,group=ID))+geom_point()+
  xlab("Change in Log Height residuals")+ylab("Change in Log FEV1 residuals")+
  stat_smooth(aes(group=1),method="loess")




library(nlme)
library(contrast)

dental <- read.csv("dentalwide.csv")
glimpse(dental)
anyNA(dental)

#melt the original data.frame as a long data.frame
dental_long <- melt(dental,id.vars=c("ID","Gender"),
              measure.vars=c("Y1","Y2","Y3","Y4"),
              variable.name="Age",
              value.name="Distance")

#change level
levels(dental_long$Age)[levels(dental_long$Age)=="Y1"]='8'
levels(dental_long$Age)[levels(dental_long$Age)=="Y2"]='10'
levels(dental_long$Age)[levels(dental_long$Age)=="Y3"]='12'
levels(dental_long$Age)[levels(dental_long$Age)=="Y4"]='14'

#mean, sd, variance for each gender at each age occasion
dental_summary <- dental_long %>%
  group_by(Gender, Age) %>%
  summarize(mean=mean(Distance), sd=sd(Distance), variance=var(Distance))

print.data.frame(dental_summary, digits = 3, row.names = FALSE)

#create two data.frame with different gender
dental_boy <- dental %>%
  filter(Gender=="M") %>%
  select(Y1, Y2, Y3, Y4)
dental_girl <- dental %>%
  filter(Gender=="F") %>%
  select(Y1, Y2, Y3, Y4)

#cov and cor for each gender
cov(dental_boy)
cov(dental_girl)
cor(dental_boy)
cor(dental_girl)


#error bar plot with gender 
#use dental_long
glimpse(dental_long)
anyNA(dental_long)
dental_mean <- dental_long %>%
  group_by(Gender, Age) %>%
  summarize(mean=mean(Distance), sd=sd(Distance), n=n())
dental_mean$se=dental_mean$sd/sqrt(dental_mean$n)
dental_mean

ggplot(dental_mean,aes(x=Age,y=mean,group=Gender,colour=Gender))+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=0.1)+
  geom_point()+geom_line()+
  labs(title="Mean plot of distances with S.E. over time")+
  xlab("Age / yrs-old")+ylab("Distances / mm")+
  theme(
    plot.title=element_text(hjust=0.5,color="blue",size=12,face="bold")
  )



#create a new variable visit
dental_model <- dental_long %>% 
           group_by(ID) %>% 
           mutate(Visit=row_number())
dental_model$Age <- as.numeric(dental_model$Age)
dental_model$Age[dental_model$Age==1] <- 8
dental_model$Age[dental_model$Age==2] <- 10
dental_model$Age[dental_model$Age==3] <- 12
dental_model$Age[dental_model$Age==4] <- 14
dental_model$Gender <- as.factor(dental_model$Gender)
dental_model$Gender=relevel(dental_model$Gender,ref="M")
levels(dental_model$Gender)
dental_model$Age
glimpse(dental_model)

#build GLS models with age and gender respectively
#unstructured model
model_age <- gls(Distance~Age,data=dental_model,corr=corSymm(form=~Visit|ID),
           weights=varIdent(form=~1|Visit),
           na.action=na.omit,method='ML')
summary(model_age)
intervals(model_age, which = "coef")  #CI for beta
vcov(model_age)   
getVarCov(model_age)    #cov matrix
corMatrix(model_age$modelStruct$corStruct)   #correlation matrix; same for each subject


model_gender <- gls(Distance~Gender,data=dental_model,corr=corSymm(form=~Visit|ID),
                 weights=varIdent(form=~1|Visit),
                 na.action=na.omit,method='ML')
summary(model_gender)
intervals(model_gender, which = "coef") 
getVarCov(model_gender)    #cov matrix
corMatrix(model_gender$modelStruct$corStruct)[["1"]]   #correlation matrix; same for each subject

#two predictors: age and gender
model_all <- gls(Distance~Age+Gender,data=dental_model,corr=corSymm(form=~Visit|ID),
                 weights=varIdent(form=~1|Visit),
                 na.action=na.omit,method='ML')
summary(model_all)
intervals(model_all, which="coef")
getVarCov(model_all)
corMatrix(model_all$modelStruct$corStruct)[["1"]]


#hypothesis test with contrast
con <- contrast(model_all,list(Age=0,Gender="F"),list(Age=0, Gender="M"))
summary(con)
con$X
con

#the model with interaction of age and gender

model_inter <- gls(Distance~Age+Gender+Age*Gender, data=dental_model, corr=corSymm(form=~Visit|ID),
                    weights=varIdent(form=~1|Visit),
                    na.action=na.omit, method="ML")
summary(model_inter)
contrast(model_inter, list(Age=10, Gender="F"), list(Age=10, Gender="M"))



