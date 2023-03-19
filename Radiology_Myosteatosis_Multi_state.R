## Radiology Myosteatosis paper revision ##

library("class")
library("dplyr")
library("plyr")
library("RANN")
library(randomForest)
library(gbm)
library(jsonlite)
library(purrr)
library(data.table)
library(ggrepel)
library(rlist)
library(stringr)
library("MSPrep")
library(gtsummary)
library(survival)
library(gt)
library(survminer)
library(mstate)

#################################################
# MULTI STATE EVENT MODELLING
#################################################

### DATA PREPARATION ###

## Create a composite column for the first major cardiovascular event (MACE). 

df_multistate <- df

df_multistate$number_MACE_incident <- df_multistate$CVD_event_incident +
  df_multistate$Heart_failure_event_incident +
  df_multistate$MI_event_incident +
  df_multistate$Ao_Aneur_event_incident 

df_multistate$first_MACE_incident <- ifelse((df$CVD_event_incident + df$Heart_failure_event_incident + 
                                               df$MI_event_incident + df$Ao_Aneur_event_incident) > 0, 1,0)

df_multistate$Time_Death_or_followup <- as.numeric(df_multistate$Time_Death_or_followup)
df_multistate$Death_event <- as.numeric(df_multistate$Death_event)-1

## Compute the time of the first, second, third and fourth MACE time

for(i in 1:nrow(df_multistate)) {
  tmp <- as.numeric(df_multistate[i,c("CVD_event_incident_time", "Heart_failure_event_incident_time",
                                      "MI_event_incident_time", "Ao_Aneur_event_incident_time")])
  
  df_multistate[i,"first_MACE_incident_time"] <- tmp[order(tmp)][1] 
  df_multistate[i,"second_MACE_incident_time"] <- tmp[order(tmp)][2]
  df_multistate[i,"third_MACE_incident_time"] <- tmp[order(tmp)][3]
  df_multistate[i,"fourth_MACE_incident_time"] <- tmp[order(tmp)][4]
}

## For those with NA (meaning the event did not occur), add censoring time (which is death or lost of follow-up)

df_multistate[which(is.na(df_multistate[,"first_MACE_incident_time"])),"first_MACE_incident_time"] <- df_multistate[which(is.na(df_multistate[,"first_MACE_incident_time"])),"Time_Death_or_followup"]
df_multistate[which(is.na(df_multistate[,"second_MACE_incident_time"])),"second_MACE_incident_time"] <- df_multistate[which(is.na(df_multistate[,"second_MACE_incident_time"])),"Time_Death_or_followup"]
df_multistate[which(is.na(df_multistate[,"third_MACE_incident_time"])),"third_MACE_incident_time"] <- df_multistate[which(is.na(df_multistate[,"third_MACE_incident_time"])),"Time_Death_or_followup"]
df_multistate[which(is.na(df_multistate[,"fourth_MACE_incident_time"])),"fourth_MACE_incident_time"] <- df_multistate[which(is.na(df_multistate[,"fourth_MACE_incident_time"])),"Time_Death_or_followup"]

## Add the event column, 0 (did not occur) or 1 if there is an event

df_multistate$first_MACE_event_incident <- 0
df_multistate$second_MACE_event_incident <- 0
df_multistate$third_MACE_event_incident <- 0
df_multistate$fourth_MACE_event_incident <- 0

df_multistate[which(df_multistate[,"number_MACE_incident"] > 0), "first_MACE_event_incident"] <- 1
df_multistate[which(df_multistate[,"number_MACE_incident"] > 1), "second_MACE_event_incident"] <- 1
df_multistate[which(df_multistate[,"number_MACE_incident"] > 2), "third_MACE_event_incident"] <- 1
df_multistate[which(df_multistate[,"number_MACE_incident"] > 3), "fourth_MACE_event_incident"] <- 1

## Quality Check that everything is ok in the number of events and incident column
tmp <- df_multistate[,c("number_MACE_incident","first_MACE_event_incident", "second_MACE_event_incident",
                        "third_MACE_event_incident","fourth_MACE_event_incident")]


tmat <- mstate::transMat(x = list(c(2,6), 
                                  c(3,6), 
                                  c(4,6), 
                                  c(5,6),
                                  c(6),
                                  c()),
                         names = c("CT_scan", "first_MACE", "second_MACE", "third_MACE","fourth_MACE", "Death"))


print(tmat)


## If seems that last contact date is anterior to incident events in a few patients
## Detecting them :

patients_check <- msebmt %>% filter (Tstop ==Inf) %>% select(PT_ID) %>% unique() %>% pull()

tmp_check <- df_multistate[df_multistate$PT_ID %in% patients_check,]


# Problem identified: since none of those died and all had only 1 MACE (and then lost-follow up)
# we will replace censoring time by first MACE time only for those patients


df_multistate_final <- df_multistate
for(id in patients_check) {
  tmp <- df_multistate_final[which(df_multistate_final$PT_ID ==id), 
                             c("number_MACE_incident","first_MACE_incident_time", "second_MACE_incident_time",
                               "third_MACE_incident_time","fourth_MACE_incident_time", 
                               "Time_Death_or_followup")]
  
  position_last_contact <- tmp$number_MACE_incident + 1
  new_last_contact <- tmp[[position_last_contact]]
  
  df_multistate_final[which(df_multistate_final$PT_ID ==id), "Time_Death_or_followup"] <- new_last_contact
  
  if(tmp$number_MACE_incident ==1) {
    
    df_multistate_final[which(df_multistate_final$PT_ID ==id), "second_MACE_incident_time"] <- new_last_contact
    df_multistate_final[which(df_multistate_final$PT_ID ==id), "third_MACE_incident_time"] <- new_last_contact
    df_multistate_final[which(df_multistate_final$PT_ID ==id), "fourth_MACE_incident_time"] <- new_last_contact
  }
  else if(tmp$number_MACE_incident ==2) {
    df_multistate_final[which(df_multistate_final$PT_ID ==id), "third_MACE_incident_time"] <- new_last_contact
    df_multistate_final[which(df_multistate_final$PT_ID ==id), "fourth_MACE_incident_time"] <- new_last_contact
  }
  else if(tmp$number_MACE_incident ==3) {
    df_multistate_final[which(df_multistate_final$PT_ID ==id), "fourth_MACE_incident_time"] <- new_last_contact
  }
  else{print(paste0("There was a problem with patient ID ", id))}
}

# Check everything has been properly adapted 

df_multistate_final_check <- df_multistate_final[df_multistate_final$PT_ID %in% patients_check,] %>% 
  select(c("number_MACE_incident","first_MACE_incident_time", "second_MACE_incident_time",
           "third_MACE_incident_time","fourth_MACE_incident_time", 
           "Time_Death_or_followup"))

msprep_clean <- msprep(data = df_multistate_final, trans = tmat, 
                       time = c(NA, "first_MACE_incident_time", "second_MACE_incident_time","third_MACE_incident_time",
                                "fourth_MACE_incident_time", "Time_Death_or_followup"), 
                       status = c(NA, "first_MACE_event_incident", "second_MACE_event_incident", "third_MACE_event_incident", 
                                  "fourth_MACE_event_incident", "Death_event"), 
                       keep = predictors, id="PT_ID")

mstate::events(msprep_clean) # Check all transitions make sense and compare them with the numbers below

table(df_multistate_final$first_MACE_event_incident)
table(df_multistate_final$second_MACE_event_incident)
table(df_multistate_final$third_MACE_event_incident)
table(df_multistate_final$fourth_MACE_event_incident)


## Last QC before modelling : We add +1 day when 2 event were diagnosed the same day

msprep_clean[which(msprep_clean$Tstart == msprep_clean$Tstop), "Tstop"] <- msprep_clean[which(msprep_clean$Tstart == msprep_clean$Tstop), "Tstop"] + 1


### Multi-state modelling of Death ###

tmat
death_transitions <- c(2,4,6,8,9) # Isolate death transitions

for(transition in death_transitions) {
  
  assign(paste0("cox_model_death_transition_",transition), coxph(Surv(Tstart,Tstop,status) ~ Age_at_CT + Sex + 
                                                                   BMI +
                                                                   Tobacco_yes_no + T2D_at_CT + Sarcopenia +
                                                                   NAFLD_57HU + VAT_Area_cm2 +
                                                                   Muscle_HU_25perc_vs_75_perc_sex_spec +
                                                                   history_of_CV_yes_no ,
                                                                   data = msprep_clean %>% filter(trans == transition)))
}


list_cox_model_death_transition <- list(cox_model_all_death,
                                        cox_model_death_transition_2,cox_model_death_transition_4,
                                        cox_model_death_transition_6, cox_model_death_transition_8,
                                        cox_model_death_transition_9)

### Multi-state modelling of Cardiovascular events ###

tmat
cardiovascular_transitions <- c(1,3,5,7) # Isolate cardiovascular transitions

for(transition in cardiovascular_transitions) {
  
  assign(paste0("cox_model_cardiovascular_transition_",transition), coxph(Surv(Tstart,Tstop,status) ~ Age_at_CT + Sex + BMI +
                                                                            Tobacco_yes_no + T2D_at_CT + Sarcopenia +
                                                                            NAFLD_57HU + VAT_Area_cm2 +
                                                                            Muscle_HU_25perc_vs_75_perc_sex_spec +
                                                                            history_of_CV_yes_no ,
                                                                            data = msprep_clean %>% filter(trans == transition)))
}



list_cox_model_cardiovascular_transition <- list(cox_model_cardiovascular_transition_1,cox_model_cardiovascular_transition_3,
                                                 cox_model_cardiovascular_transition_5, cox_model_cardiovascular_transition_7)
