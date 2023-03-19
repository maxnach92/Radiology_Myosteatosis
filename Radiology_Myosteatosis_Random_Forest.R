## Radiology Myosteatosis paper revision ##

library(plyr)
library(dplyr)
library(survex)
library(iBreakDown)
library(survival)
library("DALEX")
library("randomForestSRC")
library("furrr")
library(ranger)
library(reshape)
library(EnvStats)
library("viridis") 

#################################################
# RANDOM FOREST SURVIVAL MODELLING
#################################################


## Generate models accounting for biological sex as covariate (all) or fully-stratified for sex (male/female)

rf_survival_all <- rfsrc(Surv(Time_Death_or_followup, Death_event_incident) ~ ., data = RF_df, 
                         importance = "permute", 
                         nodesize = 4,
                         mtry = 3,
                         ntree = 1000, 
                         nsplit = 50)


rf_survival_male <- rfsrc(Surv(Time_Death_or_followup, Death_event_incident) ~ ., 
                          data = RF_df %>% filter(Sex == "Male") %>% select(-Sex), 
                          importance = "permute")


rf_survival_female <- rfsrc(Surv(Time_Death_or_followup, Death_event_incident) ~ ., 
                            data = RF_df %>% filter(Sex == "Female") %>% select(-Sex), 
                            importance = "permute")

## Extract explainability data from the model and generate partial-dependency plots

rf_explainer_all <- explain(rf_survival_all)
rf_explainer_male <- explain(rf_survival_male)
rf_explainer_female <- explain(rf_survival_female)

rf_profile_all <- model_profile(rf_explainer_all)
rf_profile_male <- model_profile(rf_explainer_male)
rf_profile_female <- model_profile(rf_explainer_female)


## Partial-dependency plots

col <- viridis::plasma(3)

plot(rf_profile_all, facet_ncol = 6, title = "", subtitle = "", 
     variables = c("Age_at_CT","Muscle_HU","BMI",'VAT_Area_cm2', "L3_SMI_cm2m2", "Liver_HU_Median"),
     colors = rev(col))

plot(rf_profile_male, facet_ncol = 6, title = "", subtitle = "", 
     variables = c("Age_at_CT","Muscle_HU","BMI",'VAT_Area_cm2', "L3_SMI_cm2m2", "Liver_HU_Median"),
     colors = rev(col))

plot(rf_profile_female, facet_ncol = 6, title = "", subtitle = "", 
     variables = c("Age_at_CT","Muscle_HU","BMI",'VAT_Area_cm2', "L3_SMI_cm2m2", "Liver_HU_Median"),
     colors = rev(col))

## Compute variable importance in each model

rf_vimp_all <- as.data.frame(rf_survival_all$importance) 
colnames(rf_vimp_all) <- "value"
rf_vimp_all$Sex <- "Both"
rf_vimp_all$variable <- rownames(rf_vimp_all)

rf_vimp_male <- as.data.frame(rf_survival_male$importance)
colnames(rf_vimp_male) <- "value"
rf_vimp_male$variable <- rownames(rf_vimp_male)
rf_vimp_male$Sex <- "Male"

rf_vimp_female <- as.data.frame(rf_survival_female$importance)
colnames(rf_vimp_female) <- "value"
rf_vimp_female$variable <- rownames(rf_vimp_female)
rf_vimp_female$Sex <- "Female"

rf_vimp_df <- rbind(rf_vimp_male,rf_vimp_female,rf_vimp_all)
rownames(rf_vimp_df) <- NULL

