####----INTRODUCTION----#####
# Manuscript title: ". “Larger territories reduce mortality risk for chimpanzees, wolves, and agents: multiple lines of evidence in a model validation framework"
# Manuscript authors: Kristin N. Crouse, Nisarg P. Desai, Kira A. Cassidy, Erin E. Stahler, Clarence L. Lehman, and Michael L. Wilson
# Code author: Nisarg Desai, desai054[at]umn[dot]edu
# Last update: July 24, 2021

rm(list = ls())
setwd("Directory containing S2_Lethal_Geometry_data.csv and S3_chimp_intergroup_victims.txt")
# Data available at reasonable request to Kristin N. Crouse crou0048[at]umn[dot]edu 

### ANALYSES OF THE DATA FROM THE AGENT BASED MODEL: LETHALGEOMETRY

# Data loading and cleaning steps

geom <- read.csv("Geometry model raw data_stopat1000_ngroups10.csv")
geom$simulation_no <- as.factor(geom$simulation_no)
geom$p_number_of_groups <- as.factor(geom$p_number_of_groups)
geom$p_patch_growth_rate <- as.factor(geom$p_patch_growth_rate)
geom$p_movement_cost <- as.factor(geom$p_movement_cost)
geom$p_aggression_cost <- as.factor(geom$p_aggression_cost)
geom$p_birth_cost <- as.factor(geom$p_birth_cost)
geom$p_stop_at <- as.factor(geom$p_stop_at)

# Define the independent variable: the periphery-territory ratio (to represent percent of territory area that is periphery)

geom$percentperiphery <-geom$group_periphery_count / geom$group_territory_size

# Mixed models to assess the directions of relationships of (i) the mortality rate and (ii) the fertility rate with percent periphery

library(lme4)

# Model of mortality rate with percent periphery with random effects
model.mortality <- glmer(group_death_count ~ percentperiphery + (1|simulation_no) + (1|p_patch_growth_rate) + (1|p_movement_cost) + (1|p_aggression_cost) + (1|p_birth_cost), offset = log(group_population_size), family = poisson(link = "log"), data = geom)
summary(model.mortality) 
CI_mm <- confint(model.mortality)
CI_mm

# Model of fertility rate with percent periphery with random effects
model.fertility <- glmer(group_birth_no ~ percentperiphery + (1|simulation_no) + (1|p_patch_growth_rate) + (1|p_movement_cost) + (1|p_birth_cost), offset = log(group_population_size), family = poisson(link = "log"), data = geom)  # p_aggression_cost removed from random effects as that had 0 variance explained and made the model too complex for the data
summary(model.fertility)
CI_mf <- confint(model.fertility)
CI_mf


### ANALYSES OF THE LETHALGEOMETRY DATA WITH 0 GROUP DEATH COUNTS REMOVED
# the following chunks of code till line 68 assume that the data loading and cleaning steps till the line 22 have been run

# Remove rows with group_death_count = 0

geom_nozero <- geom[which(geom$group_death_count != 0),]


# Model of mortality rate with percent periphery with random effects
model.mortality.nozero <- glmer(group_death_count ~ percentperiphery + (1|simulation_no) + (1|p_patch_growth_rate) + (1|p_movement_cost) + (1|p_aggression_cost) + (1|p_birth_cost), offset = log(group_population_size), family = poisson(link = "log"), data = geom_nozero)
summary(model.mortality.nozero) 
CI_mmno <- confint(model.mortality.nozero)
CI_mmno

# Model of fertility rate with percent periphery with random effects
model.fertility.nozero <- glmer(group_birth_no ~ percentperiphery + (1|simulation_no) + (1|p_patch_growth_rate) + (1|p_movement_cost) + (1|p_birth_cost), offset = log(group_population_size), family = poisson(link = "log"), data = geom_nozero) # p_aggression_cost removed from random effects as that had 0 variance explained and made the model too complex for the data
summary(model.fertility.nozero)
CI_mfno <- confint(model.fertility.nozero)
CI_mfno

### ANALYSES OF THE EMPRIRICAL DATA: CHIMPANZEES

library(MuMIn)

## Models with strongest evidence (observed and inferred cases)
# Fit the 8 possible models
empiricaldata<-read.table("S3_chimp_data.txt")
fm001<-glm(intergroup_victims_obs_inf ~ 1, offset = log(chimp_years), family=poisson, data = empiricaldata)
fm002<-glm(intergroup_victims_obs_inf ~ inverse_r, offset = log(chimp_years), family=poisson, data = empiricaldata)
fm003<-glm(intergroup_victims_obs_inf ~ adult_males, offset = log(chimp_years), family=poisson, data = empiricaldata)
fm004<-glm(intergroup_victims_obs_inf ~ density, offset = log(chimp_years), family=poisson, data = empiricaldata)
fm005<-glm(intergroup_victims_obs_inf ~ inverse_r + adult_males, offset = log(chimp_years), family=poisson, data = empiricaldata)
fm006<-glm(intergroup_victims_obs_inf ~ inverse_r + density, offset = log(chimp_years), family=poisson, data = empiricaldata)
fm007<-glm(intergroup_victims_obs_inf ~ adult_males + density, offset = log(chimp_years), family=poisson, data = empiricaldata)
fm008<-glm(intergroup_victims_obs_inf ~ inverse_r + adult_males + density, offset = log(chimp_years), family=poisson, data = empiricaldata, na.action = "na.fail")
summary.lm(fm008)

# Calculate the dispersion parameter
chat <- deviance(fm008) / df.residual(fm008)

# Create a model selection table by ranking models based on QAICc
ms_tab <- dredge(fm008, rank = "QAICc", chat = chat) 
ms_tab
#write.csv(ms_tab, "Strong evidence with Kahama.csv")

# Perform model-averaging by ranking models based on QAICc
p <- model.avg(fm008, fm007, fm006, fm005, fm004, fm003, fm002, fm001, rank = QAICc, rank.args = list(chat = chat)) # perform model averaging
summary(p)
confint(p)

# Models with strongest evidence WITHOUT Kahama

empiricaldata_nokahama <- empiricaldata[-2, ]
fm001_nk <-glm(intergroup_victims_obs_inf ~ 1, offset = log(chimp_years), family=poisson, data = empiricaldata_nokahama)
fm002_nk <-glm(intergroup_victims_obs_inf ~ inverse_r, offset = log(chimp_years), family=poisson, data = empiricaldata_nokahama)
fm003_nk <-glm(intergroup_victims_obs_inf ~ adult_males, offset = log(chimp_years), family=poisson, data = empiricaldata_nokahama)
fm004_nk <-glm(intergroup_victims_obs_inf ~ density, offset = log(chimp_years), family=poisson, data = empiricaldata_nokahama)
fm005_nk <-glm(intergroup_victims_obs_inf ~ inverse_r + adult_males, offset = log(chimp_years), family=poisson, data = empiricaldata_nokahama)
fm006_nk <-glm(intergroup_victims_obs_inf ~ inverse_r + density, offset = log(chimp_years), family=poisson, data = empiricaldata_nokahama)
fm007_nk <-glm(intergroup_victims_obs_inf ~ adult_males + density, offset = log(chimp_years), family=poisson, data = empiricaldata_nokahama)
fm008_nk <-glm(intergroup_victims_obs_inf ~ inverse_r + adult_males + density, offset = log(chimp_years), family=poisson, data = empiricaldata_nokahama, na.action = "na.fail")

chat_nk <- deviance(fm008_nk) / df.residual(fm008_nk)

ms_tab_nk <- dredge(fm008_nk, rank = "QAICc", chat = chat_nk)
ms_tab_nk
#write.csv(ms_tab_nk, "Strong evidence without Kahama.csv")

p_nk <- model.avg(fm008_nk, fm007_nk, fm006_nk, fm005_nk, fm004_nk, fm003_nk, fm002_nk, fm001_nk, rank = QAICc, rank.args = list(chat = chat_nk))
summary(p_nk)
confint(p_nk)

# Empirical with suspected (all cases)

fm001_ws <-glm(with_suspected_victims ~ 1, offset = log(chimp_years), family=poisson, data = empiricaldata)
fm002_ws <-glm(with_suspected_victims ~ inverse_r, offset = log(chimp_years), family=poisson, data = empiricaldata)
fm003_ws <-glm(with_suspected_victims ~ adult_males, offset = log(chimp_years), family=poisson, data = empiricaldata)
fm004_ws <-glm(with_suspected_victims ~ density, offset = log(chimp_years), family=poisson, data = empiricaldata)
fm005_ws <-glm(with_suspected_victims ~ inverse_r + adult_males, offset = log(chimp_years), family=poisson, data = empiricaldata)
fm006_ws <-glm(with_suspected_victims ~ inverse_r + density, offset = log(chimp_years), family=poisson, data = empiricaldata)
fm007_ws <-glm(with_suspected_victims ~ adult_males + density, offset = log(chimp_years), family=poisson, data = empiricaldata)
fm008_ws <-glm(with_suspected_victims ~ inverse_r + adult_males + density, offset = log(chimp_years), family=poisson, data = empiricaldata, na.action = "na.fail")

chat_ws <- deviance(fm008_ws) / df.residual(fm008_ws)

ms_tab_ws <- dredge(fm008_ws, rank = "QAICc", chat = chat_ws)
ms_tab_ws
#write.csv(ms_tab_ws, "All cases with Kahama.csv")

p_ws <- model.avg(fm008_ws, fm007_ws, fm006_ws, fm005_ws, fm004_ws, fm003_ws, fm002_ws, fm001_ws, rank = QAICc, rank.args = list(chat = chat_ws))
summary(p_ws)
confint(p_ws)

# With suspected without Kahama

fm001_ws_nk <-glm(with_suspected_victims ~ 1, offset = log(chimp_years), family=poisson, data = empiricaldata_nokahama)
fm002_ws_nk <-glm(with_suspected_victims ~ inverse_r, offset = log(chimp_years), family=poisson, data = empiricaldata_nokahama)
fm003_ws_nk <-glm(with_suspected_victims ~ adult_males, offset = log(chimp_years), family=poisson, data = empiricaldata_nokahama)
fm004_ws_nk <-glm(with_suspected_victims ~ density, offset = log(chimp_years), family=poisson, data = empiricaldata_nokahama)
fm005_ws_nk <-glm(with_suspected_victims ~ inverse_r + adult_males, offset = log(chimp_years), family=poisson, data = empiricaldata_nokahama)
fm006_ws_nk <-glm(with_suspected_victims ~ inverse_r + density, offset = log(chimp_years), family=poisson, data = empiricaldata_nokahama)
fm007_ws_nk <-glm(with_suspected_victims ~ adult_males + density, offset = log(chimp_years), family=poisson, data = empiricaldata_nokahama)
fm008_ws_nk <-glm(with_suspected_victims ~ inverse_r + adult_males + density, offset = log(chimp_years), family=poisson, data = empiricaldata_nokahama, na.action = "na.fail")

chat_ws_nk <- deviance(fm008_ws_nk) / df.residual(fm008_ws_nk)

ms_tab_ws_nk <- dredge(fm008_ws_nk, rank = "QAICc", chat = chat_ws_nk)
ms_tab_ws_nk
#write.csv(ms_tab_ws_nk, "All cases without Kahama.csv")

p_ws_nk <- model.avg(fm008_ws_nk, fm007_ws_nk, fm006_ws_nk, fm005_ws_nk, fm004_ws_nk, fm003_ws_nk, fm002_ws_nk, fm001_ws_nk, rank = QAICc, rank.args = list(chat = chat_ws_nk))
summary(p_ws_nk)
confint(p_ws_nk)




### ANALYSES OF THE EMPRIRICAL DATA: WOLVES

# Data loading and cleaning steps

geom_wolf <- read.csv("S4_yellowstone_wolf.csv")
geom_wolf$biological_year <- as.factor(geom_wolf$biological_year)
geom_wolf$pack_name <- as.factor(geom_wolf$pack_name)

# Mixed models to assess the directions of relationships of (i) the mortality rate and (ii) the fertility rate with inverse radius

library(lme4)

# Model of mortality rate with percent periphery with random effects

#add all deaths variable
geom_wolf$all_deaths <- round(geom_wolf$intraspecific_deaths + geom_wolf$intraspecific_deaths_of_litter)

#add a variable to indicate cases with non-zero litter deaths=1. And then a new variable adding 1 litter death
geom_wolf$litter_deaths_indicator <- ifelse(geom_wolf$intraspecific_deaths_of_litter > 0, 1, 0)
geom_wolf$deaths_with_litter <- geom_wolf$intraspecific_deaths + geom_wolf$litter_deaths_indicator

# Model with all deaths
model.mortality.1 <- glmer(all_deaths ~ inverse_radius + no..of.radio.collars.present + (1|biological_year) + (1|pack_name), offset = log(pack_size_avg), family = poisson(link = "log"), weights = territory_size_reliability, data = geom_wolf)
summary(model.mortality.1) 
CI_wm_1 = confint(model.mortality.1)
CI_wm_1
 
# Model with intraspecific + litter deaths as 1
model.mortality.2 <- glmer(deaths_with_litter ~ inverse_radius + no..of.radio.collars.present + (1|biological_year) + (1|pack_name), offset = log(pack_size_avg), family = poisson(link = "log"), weights = territory_size_reliability, data = geom_wolf)
summary(model.mortality.2) 
CI_wm_2 = confint(model.mortality.2)
CI_wm_2

# Model only with intraspecific deaths
model.mortality.3 <- glmer(intraspecific_deaths ~ inverse_radius + no..of.radio.collars.present + (1|biological_year) + (1|pack_name), offset = log(pack_size_avg), family = poisson(link = "log"), weights = territory_size_reliability, data = geom_wolf)
summary(model.mortality.3) 
CI_wm_3 = confint(model.mortality.3)
CI_wm_3

# quick and dirty plot to spot effects
plot(geom_wolf$intraspecific_deaths ~ geom_wolf$inverse_radius)
abline(lm(geom_wolf$intraspecific_deaths ~ geom_wolf$inverse_radius))

### MODELS EXCLUDING CASES WITH 0 DEATHS

# remove rows with 0 deaths
geom_wolf_nozero <- geom_wolf[which(geom_wolf$intraspecific_deaths != 0), ] 
geom_wolf_nozero <- droplevels(geom_wolf_nozero)

# # GLM without random effects
# model.mortality.4 <- glm(intraspecific_deaths ~ inverse_radius, offset = log(pack_size_avg), family = poisson(link = "log"), weights = territory_size_reliability, data = geom_wolf_nozero)
# summary(model.mortality.4)
# plot(model.mortality.4)

# GLMM with all deaths and random effects
model.mortality.4 <- glmer(all_deaths ~ inverse_radius + no..of.radio.collars.present + (1|biological_year) + (1|pack_name), offset = log(pack_size_avg), family = poisson(link = "log"), weights = territory_size_reliability, data = geom_wolf_nozero)
summary(model.mortality.4) 
CI_wmno_1 = confint(model.mortality.4)
CI_wmno_1

# GLMM with intraspecific + litter deaths as 1 and random effects
model.mortality.5<- glmer(deaths_with_litter ~ inverse_radius + no..of.radio.collars.present + (1|biological_year) + (1|pack_name), offset = log(pack_size_avg), family = poisson(link = "log"), weights = territory_size_reliability, data = geom_wolf_nozero)
summary(model.mortality.5) 
CI_wmno_2 = confint(model.mortality.5)
CI_wmno_2

# GLMM with only intraspecific deaths and random effects
model.mortality.6 <- glmer(intraspecific_deaths ~ inverse_radius + no..of.radio.collars.present + (1|biological_year) + (1|pack_name), offset = log(pack_size_avg), family = poisson(link = "log"), weights = territory_size_reliability, data = geom_wolf_nozero)
summary(model.mortality.6) 
CI_wmno_3 = confint(model.mortality.6)
CI_wmno_3

# quick and dirty plot to spot effects
plot(geom_wolf_nozero$intraspecific_deaths ~ geom_wolf_nozero$inverse_radius)
abline(lm(geom_wolf_nozero$intraspecific_deaths ~ geom_wolf_nozero$inverse_radius))