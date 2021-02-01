#Homework 3 for longitudinal 

library(dplyr)
library(ggplot2)
library(nlme)
library(lme4)
library(reshape2)
library(lmeInfo)
library(lmerTest)

#input data
exercise <- read.csv("/Users/Kin/Downloads/Longitudinal/exercise.csv")
glimpse(exercise)
#missing values exist with "."

#melt oringinal dataset
exercise_long <- melt(exercise, id.vars = c("ID", "TRT"),
                      measure.vars = c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7"),
                      variable.name = "Days",
                      value.name = "Strength")
glimpse(exercise_long)
table(exercise_long$Days)
class(exercise_long$Days)


#relevel Days
levels(exercise_long$Days)[levels(exercise_long$Days)=="Y1"] <- "0"
levels(exercise_long$Days)[levels(exercise_long$Days)=="Y2"] <- "2"
levels(exercise_long$Days)[levels(exercise_long$Days)=="Y3"] <- "4"
levels(exercise_long$Days)[levels(exercise_long$Days)=="Y4"] <- "6"
levels(exercise_long$Days)[levels(exercise_long$Days)=="Y5"] <- "8"
levels(exercise_long$Days)[levels(exercise_long$Days)=="Y6"] <- "10"
levels(exercise_long$Days)[levels(exercise_long$Days)=="Y7"] <- "12"
class(exercise_long$Strength)
class(exercise_long$Days)


#Set NA in strength
exercise_long$Strength[exercise_long$Strength=="."] <- NA
table(exercise_long$Strength, useNA = "ifany")
#set Strength as numeric
exercise_long$Strength <- as.numeric(exercise_long$Strength)
#set TRT as factor
exercise_long$TRT <- as.factor(exercise_long$TRT)

#within subject summary stats
exercise_summary_ID <- exercise_long %>%
                  group_by(ID) %>%
                  summarize(mean=mean(Strength, na.rm=TRUE), SD=sd(Strength, na.rm = TRUE), Variance=var(Strength, na.rm = TRUE)) 
                  
exercise_summary_ID
as.data.frame(exercise_summary_ID)
mean(exercise_summary_ID$mean)
sd(exercise_summary_ID$mean)
var(exercise_summary_ID$mean)
#between subjects summary stats for each treatment 
exercise_summary_TRT <- exercise_long %>%
                  group_by(TRT) %>%
                  summarize(mean=mean(Strength, na.rm=TRUE), SD=sd(Strength, na.rm = TRUE), Variance=var(Strength, na.rm = TRUE))
                  
as.data.frame(exercise_summary_TRT)
mean(exercise_long$Strength, na.rm=TRUE)
sd(exercise_long$Strength, na.rm = TRUE)
var(exercise_long$Strength, na.rm = TRUE)

#Spaghetti plot
ggplot(data=exercise_long,aes(x=Days,y=Strength,group=ID, color=ID))+geom_line(size=1)+facet_grid(.~TRT)+
  scale_color_viridis_c(option = "C")+
  ylab("Strength")+xlab("Days")+
  labs(title="Strength vs. Days for Two Treatments")+
  theme(
    plot.title=element_text(hjust=0.5,color="blue",size=12,face="bold")
  )

#error plot for two treatments
exercise_summary <- exercise_long %>%
                group_by(TRT, Days) %>%
                summarize(mean=mean(Strength, na.rm=TRUE), SD=sd(Strength, na.rm = TRUE), Variance=var(Strength, na.rm = TRUE), n=n())
exercise_summary$se <- exercise_summary$SD/sqrt(exercise_summary$n)
exercise_summary

ggplot(exercise_summary,aes(x=Days,y=mean,group=TRT,colour=TRT))+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=0.1)+
  geom_point()+geom_line()+
  labs(title="Error bar plot of strength over time")+
  xlab("Days")+ylab("Strength")+
  theme(
    plot.title=element_text(hjust=0.5,color="blue",size=12,face="bold")
  )


#mixed model with random intercept and slope
#no cov-corr structure in mean model, so no "correlation="
#unstructured in random effects, so set "|" in random
#R same for all subjects, so no "weights"
#interation between strength with days and treatment 
#as.characrter of Days
exercise_long_model <- exercise_long
levels(exercise_long_model$Days)
exercise_long_model$Days <- as.character(exercise_long_model$Days)
exercise_long_model$Days <- as.numeric(exercise_long_model$Days)
table(exercise_long_model$Days)
class(exercise_long_model$Days)
glimpse(exercise_long_model)
head(exercise_long)
head(exercise_long_model)
tail(exercise_long)
tail(exercise_long_model)

#accurate

exercise_model <- lme(Strength~Days+TRT+Days*TRT, random = ~1+Days|ID,
                      data = exercise_long_model, na.action=na.omit, method = "REML")

summary(exercise_model)
G_Matrix <- getVarCov(exercise_model,type="random.effects") #G Matrix
G_Matrix
cov2cor(G_Matrix)
coef(exercise_model)

#There no way to easily test G_Matrix

exercise_model_1 <- lmer(Strength~Days+TRT+Days*TRT+(1+Days|ID), data=exercise_long_model,
                         REML = TRUE)
summary(exercise_model_1)
vc <- VarCorr(exercise_model_1)
print(vc,comp=c("Variance","Std.Dev."),digits=3)
as.data.frame(vc,order="lower.tri") 
coef(exercise_model_1)

#To test random effects
rand(exercise_model_1)
rand_test1 <- lmer(Strength~Days+TRT+Days*TRT+(1+Days|ID), data=exercise_long_model,
     REML = TRUE)
rand_test2 <- lmer(Strength~Days+TRT+Days*TRT+(1|ID)+(0+Days|ID), data=exercise_long_model,
                   REML = TRUE)
rand_test3 <- lmer(Strength~Days+TRT+Days*TRT+(1|ID), data=exercise_long_model,
                   REML = TRUE)
rrt2<- VarCorr(rand_test2)
as.data.frame(rrt2,order="lower.tri") 
VarCorr(rand_test3)
anova(rand_test1, rand_test2)
anova(rand_test2, rand_test3)
anova(rand_test3, rand_test1)
#Yes, you can use LRT to compare covariance models as explained in class.
#There is a R package ‘merDeriv’. Using that, you will get s.e. for each estimator of parameters in G and R matrix. 
#You can use the info to create z-test like given in SAS (covtest)

#cov of Y
getVarCov(exercise_model, type = "marginal")
cov_y1 <- getVarCov(exercise_model,type="marginal", individual=1)
cov_y2 <- getVarCov(exercise_model,type="marginal", individual=2) #cov of Y
getVarCov(exercise_model,type="marginal", individual=9) 
#some subjects have missing
exercise_long_model[is.na(exercise_long_model$Strength),]
#there is no way to easily get corr matrix because SDs are different
S <- matrix(c(1/10.6180, 1/(3.2585*3.2692), 1/(3.2585*3.3215),
              1/(3.2585*3.4134), 1/(3.2585*3.5419), 1/(3.2585*3.7032), 
              1/(3.2585*3.8932),
              1/(3.2585*3.2692), 1/10.688, 1/(3.2692*3.3215),
              1/(3.2692*3.4134), 1/(3.2692*3.5419), 1/(3.2692*3.7032), 
              1/(3.2692*3.8932),
              1/(3.2585*3.3215), 1/(3.2692*3.3215), 1/11.0320,
              1/(3.3215*3.4134), 1/(3.3215*3.5419), 1/(3.3215*3.7032), 
              1/(3.3215*3.8932),
              1/(3.2585*3.4134), 1/(3.2692*3.4134), 1/(3.4134*3.3215),
              1/11.651, 1/(3.4134*3.5419), 1/(3.4134*3.7032), 
              1/(3.4134*3.8932),
              1/(3.2585*3.5419), 1/(3.2692*3.5419), 1/(3.5419*3.3215),
              1/(3.5419*3.4134), 1/12.545, 1/(3.5419*3.7032), 
              1/(3.5419*3.8932),
              1/(3.2585*3.7032), 1/(3.2692*3.7032), 1/(3.7032*3.3215),
              1/(3.7032*3.4134), 1/(3.7032*3.5419), 1/13.7140, 
              1/(3.7032*3.8932),
              1/(3.2585*3.8932), 1/(3.2692*3.8932), 1/(3.8932*3.3215),
              1/(3.8932*3.4134), 1/(3.8932*3.5419), 1/(3.8932*3.7032), 
              1/15.1570), nrow = 7, ncol = 7, byrow = TRUE)
S
tt <- matrix(c(10.6180, 9.9194 , 9.8857 , 9.852 , 9.8183 , 9.7846 , 9.7509,
                9.9194 ,10.6880 ,10.1270 ,10.230, 10.3340 ,10.4380, 10.5410,
                9.8857 ,10.1270, 11.0320 ,10.609, 10.8490 ,11.0900 ,11.3310,
                9.8520 ,10.2300 ,10.6090, 11.651, 11.3650, 11.7430, 12.1210,
                9.8183 ,10.3340 ,10.8490 ,11.365 ,12.5450, 12.3960 ,12.9120,
                9.7846 ,10.4380, 11.0900, 11.743 ,12.3960 ,13.7140 ,13.7020,
                9.7509 ,10.5410 ,11.3310 ,12.121, 12.9120 ,13.7020 ,15.1570),
             nrow = 7, ncol = 7, byrow = FALSE)
tt
corr_mat <- tt*S
corr_mat

#Conditional var(Yi1) and marginal var(Yi1)
#conditional variance is same for each subject at each occassion
#some IDs have missing, so we only use ID 2 with 7 time point
getVarCov(exercise_model, type = "conditional", individuals = 2)  #conditional Yi1 
getVarCov(exercise_model, type = "marginal", individuals = 2)   
#marginal variance at baseline is 10.6180



#model group by TRT
model_TRT <- lme(Strength~Days+TRT+Days*TRT, random = ~1+Days|ID,
                 weights = varIdent(form=~1|TRT),
                 data = exercise_long_model, na.action=na.omit, method = "REML")
summary(model_TRT)
anova(exercise_model, model_TRT)

#conditional
getVarCov(model_TRT, type = "conditional", individuals = 1)   #TRT 1
getVarCov(model_TRT, type = "conditional", individuals = 2)   #TRT 1
getVarCov(model_TRT, type = "marginal", individuals = 1)   #TRT 1
getVarCov(model_TRT, type = "marginal", individuals = 2)   #TRT 1
#ID 9 even has two missing

anova(exercise_model, model_TRT)

#predicted Y for each subject
fitted_y <- fitted(exercise_model,level=0:1) 
fitted_y <- as.data.frame(fitted_y)
fitted_y
nrow(fitted_y)
#missing caused row number is smaller

pred_ex <- predict(exercise_model_1, na.action=na.exclude)
pred_ex
#set missing as NA
pred_mar <- predict(exercise_model_1,re.form=NA, na.action=na.exclude) 
pred_mar
#pred_ex and pred_mar have same length of column with the oringinal data.frame

names(exercise_long_model)

#a data.frame with predictted values
fitted_table <- data_frame(ID=exercise_long_model$ID, TRT=exercise_long_model$TRT,
                           Days=exercise_long_model$Days, Strength=exercise_long_model$Strength,
                           fitted=pred_mar,
                           pred_ID=pred_ex)
fitted_table <- as.data.frame(fitted_table)
glimpse(fitted_table)

#plot predicted values
ggplot(data=fitted_table,aes(x=factor(Days),y=pred_ID,group=ID, color=ID))+
  geom_path(size=1, linetype="dashed")+facet_grid(.~TRT)+
  scale_color_viridis_c(option = "D")+
  ylab("Predicted Strength")+xlab("Days")+
  labs(title="Predicted Strength vs. Days for Two Treatments")+
  theme(
    plot.title=element_text(hjust=0.5,color="blue",size=12,face="bold")
  )

#fixed model without random effects
fixed_model <- exercise_long_model %>% 
  group_by(ID) %>% 
  mutate(Visit=row_number())
fixed_model

#CS model 
exercise_fixed_model_1 <- gls(Strength~Days+TRT+Days*TRT, data = fixed_model, 
                              corr=corCompSymm(form=~Visit|ID), na.action=na.omit,method='REML')
summary(exercise_fixed_model_1)
#no correlation
exercise_fixed_model_2 <- gls(Strength~Days+TRT+Days*TRT, data = fixed_model, 
                              na.action=na.omit,method='REML')
summary(exercise_fixed_model_2)

anova(exercise_model, exercise_fixed_model_1)
anova(exercise_model, exercise_fixed_model_2)

#unstructured 
unst_model <- gls(Strength~Days+TRT+Days*TRT, data = fixed_model, 
    corr=corSymm(form=~Visit|ID), na.action=na.omit,method='REML')
summary(unst_model)

anova(exercise_model, unst_model)
#but the dataset is unbalanced because there are missing values

#no interaction
no_interact1 <- lme(Strength~Days+TRT, random = ~1+Days|ID,
                   data = exercise_long_model, na.action=na.omit, method = "REML")
summary(no_interact1)

no_interact2 <- lmer(Strength~Days+TRT+(1|ID)+(0+Days|ID), data=exercise_long_model,
                     REML = TRUE)
summary(no_interact2)

#mean model comparison
mean_1 <- gls(Strength~Days+TRT, data = fixed_model, 
              corr=corCompSymm(form=~Visit|ID), na.action=na.omit,method='ML')
mean_2 <- gls(Strength~Days+TRT+Days*TRT, data = fixed_model, 
              corr=corCompSymm(form=~Visit|ID), na.action=na.omit,method='ML')

anova(mean_1, mean_2)


#cross-sectional effect and longitudinal effect
cross_data <- exercise_long_model %>%
              group_by(ID) %>%
              mutate(mean_day=mean(Days, na.rm=TRUE))
as.data.frame(cross_data)

cross_model_1 <- lmer(Strength~Days+TRT+mean_day+(1|ID)+(0+Days|ID), data=cross_data,
                    REML = TRUE)
summary(cross_model_1)

cross_model_2 <- lmer(Strength~Days+TRT+mean_day+(1+Days|ID), data=cross_data,
                      REML = TRUE)
summary(cross_model_2)

anova(cross_model_1, cross_model_2)


























