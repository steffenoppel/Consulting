## Little owl survival model
rm(list=ls())


# 1. Importing data ####
ring <- readxl::read_xlsx("data/LO_ringing.xlsx")
ring$year <- as.numeric(format(ring$date, format = "%Y")) # make a year variable
tagg <- readxl::read_xlsx("data/LO_tagging.xlsx")
tagg$ring[tagg$ring == "NA"] <- "no_ring"
tagg$year <- as.numeric(format(tagg$date, format = "%Y")) # make a year variable

# explore
time_tags_active <- aggregate(date ~ ring, data = tagg, FUN = function(x) difftime(max(x), min(x), unit = "days"))
hist(as.numeric(time_tags_active$date)); abline(v = c(30,365), col = c(2,4), lwd = 2) # lines put at 1 month and 1 year - most individuals do not have data over one year of age




# 2. Data wrangling ####

# fix age in ring
ring$DOB # lots of different values
ring$DOB <- gsub(" juv", "-05-20", ring$DOB) # if juv, add May 30th as mean hatching data

## that numeric date is an artefact of Excel - I would fix that by setting an origin, rather than subtracting a number of days as.Date(DOB, origin="1899-12-30")
## see: https://stackoverflow.com/questions/71054152/r-read-excel-reads-dates-as-numbers
## calculate any ages with 'difftime' so you don't need to worry about how many days what year has
ring$DOB <- as.Date(ifelse(grepl("-",ring$DOB), as.character(as.Date(ring$DOB)), # if date format, make date
                            ifelse(ring$DOB == "adult", as.character(as.Date(ring$date)-365), # if adult, take the ringing date and subtract a year
                                   as.character(as.Date(as.numeric(ring$DOB))- 365*70 -17)))) # is number date, make date but fix the issue that all dates seem to be 70 years in the future... subtracting an additional 17 days for the 17 leap years in a 70 period
summary(ring$DOB)
summary(ring$date)
table(ring$DOB < as.Date(ring$date))

ring$age_days <- NULL
for(i in 1:nrow(ring)) {
  ring$age_days[i] <- as.numeric(difftime(as.Date(ring$date[i]), ring$DOB[ring$ring == ring$ring[i]], units = "days"))
}; rm(i)

ind_age <- aggregate(age_days ~ ring, data = ring, FUN = max)
hist(ind_age$age_days/365)
rm(ind_age)

ring$age <- ifelse(ring$age_days < 365, 1, 2) # make categorical variable: 1 = juv, 2 = adult

ring$tag_type <- "ring"
ring <- ring[,c(1:2, 13, 8,10,12)]


# fix age in tagg
tagg$age_days <- NULL
for(i in 1:nrow(tagg)) {
  tagg$age_days[i] <- as.numeric(difftime(tagg$date[i], tagg$DOB[tagg$ring == tagg$ring[i]], units = "days"))
}; rm(i)

ind_age <- aggregate(age_days ~ ring, data = tagg, FUN = max)
hist(ind_age$age_days/365)
rm(ind_age)

tagg$age <- ifelse(tagg$age_days < 365, 1, 2)


# Easy version: transform tagging data into yearly cap/recap data
tagg[tagg$ring == "K84919",] # although less than 2 years old, this individual was seen in 3 years, but the record of the middle year is not included in the dataset. Here, we fill this missing year in.
tagg <- tagg[,c(1:2,4,9,12,14)]
for(i in unique(tagg$ring)) { # some birds lived more than one year, but that data is omitted, so we add missing years between tagging and last observation
  x <- unique(range(tagg$year[tagg$ring == i]))
  x <- min(x):max(x)
  x <- x[!(x %in% tagg$year[tagg$ring == i])]
  for(j in x) {
    tagg <- rbind(tagg, tagg[tagg$ring == i,][1,])
    tagg[nrow(tagg), "year"] <- j
    tagg[nrow(tagg), "age"] <- 2 # if we have to ad a year, then the individual is definitely an adult in this year
    tagg[nrow(tagg), "tag_type"] <- tagg[tagg$ring == i, "tag_type"][nrow(tagg[tagg$ring == i,])-1,] # the tag type is that of the last observation before the added one
  }
}; rm(i,j,x)



# 3. Make CH ####
library(tidyr)
library(dplyr)

# merge ring and tagging data
df <- rbind(ring, tagg)  

# dealing with recoveries in the same year as alive observations
for(i in unique(df$ring)) {
  if(!(0 %in% df$alive[df$ring == i])) next # if individual is not recovered, skip to next
  if(!(df$year[df$ring == i & df$alive == 0] %in% df$year[df$ring == i & df$alive == 1])) next # if there is no conflict between recovery and observation within the same year, skip to next
  df <- rbind(df, df[df$ring == i & df$alive == 0,]) # add a new row of data with the recovery
  df$year[nrow(df)] <- df$year[nrow(df)] + 1 # in that new row, change the year of recovery to the next one
  df[df$ring == i & df$alive == 0,][1,]$alive <- 1 # in the original recovery observation, change this recovery to alive
}; rm(i) # note: this overestimates survival for those individuals that were ringed once and recovered in the next year, as now they are recovered in the 3rd year - we could see how much this impacts the estimates alter

# wrangle data
df <- df %>%
  add_row(year = 2011) %>% # no data in 2011 - we add this manually so that the Cap His matrix has a column for every year
  arrange(ring, year)  %>%
  group_by(ring, year) %>%
  summarise(sex = first(sex), tag_type = paste0(unique(tag_type), collapse = "_"), alive = max(alive), age = first(age)) %>%
  group_by(ring)

## create capture history matrix
CH <- df %>% 
  arrange(year) %>% 
  pivot_wider(id_cols = ring, names_from = year, values_from = alive) %>% 
  dplyr::filter(!is.na(ring)) %>% # that extra year we added earlier also created an NA in ring, which we get rid of here
  tibble::column_to_rownames("ring") # the CH can only contain the states for each individual for each year
View(CH) # 1 = seen alive, 0 = recovered dead, NA = not seen. However, we must code these observations differently for JAGS to understand it

table(as.matrix(CH))
CH[CH == 0] <- 2 # in JAGS, observation 2 means 'recovered dead'  - but how can this then be lumped with 'not seen' in L. 183?
CH[is.na(CH)] <- 3 # in JAGS, observation 3 means 'not seen'

## create vector of first encounter
first <- apply(CH, 1, function(x) min(which(x == 1)))

## create age matrix: each individual has an age attributed to it (juv or adult) - this matrix will have the same dimensions as the CH
AM <- df %>% 
  arrange(year) %>% 
  pivot_wider(id_cols = ring, names_from = year, values_from = age) %>% 
  dplyr::filter(!is.na(ring)) %>% # that extra year we added earlier also created an NA in ring, which we get rid of here
  tibble::column_to_rownames("ring") # the CH can only contain the states for each individual for each year
AM # few individuals were caught as adults: their first record is 2; some individuals are recorded as juvenile 2 years in a row, this needs to be fixed

# make sure the that AM is filled after the first occasion
for(i in 1:nrow(AM)) {
  if(first[i] == ncol(AM)) next
  AM[i,(first[i]+1):ncol(AM)] <- 2
}; rm(i)
AM

## create the tag type matrix
TM <- df %>% 
  mutate(tag_type = ifelse(tag_type == "NA", NA, tag_type),
         tag_type = factor(gsub("ring_|_VHF","",tag_type)), # important: we assume that GPS tags have the highest priority, then VHF, then ring: if an individual has a GPS and VHF and ring, its TM data will show GPS only
         tag_type = as.numeric(tag_type)) %>% # levels: 1 = ring, 2 = GPS, 3 = VHF
  arrange(year) %>% 
  pivot_wider(id_cols = ring, names_from = year, values_from = tag_type) %>% 
  dplyr::filter(!is.na(ring)) %>% # that extra year we added earlier also created an NA in ring, which we get rid of here
  tibble::column_to_rownames("ring") # the CH can only contain the states for each individual for each year
View(TM) # this looks good, though again we have to fill this matrix

table(apply(TM, 1, function(x) length(unique(na.exclude(x))))) # Some individuals have two different tags over time

for(i in 1:nrow(TM)) {
  for(t in first[i]:ncol(TM)) {
    if(is.na(TM[i,t])) TM[i,t] <- TM[i,t-1] # at the first occasion all individuals have a tag type. After, if a cell is NA, fill it with the tag type from the previous cell
  }
}; rm(i,t)

TM2 <- TM
TM2[TM == 3] <- 2

TM3 <- NULL
for(i in 1:dim(TM)[1]) {TM3[i] <- TM[i,first[i]]}; rm(i)
table(TM3)


## preparing sex vector
df$sex <- ifelse(df$sex == "f", 1, ifelse(df$sex == "m", 2, NA))
sex <- aggregate(sex ~ ring, data = df, FUN = dplyr::first, na.action = na.pass)$sex
table(sex, useNA = "i") # equal number of females and males, as many combined NAs
sex.ratio <- prop.table(table(sex))


## Save data
save.image("data/jags_prep.RData")





# 4a. Running CJS survival model in JAGS ####
library(jagsUI)
rm(list=ls())

# load data
load("data/jags_prep.RData")

# Making 3-states into 2-states
CH[CH == 3] <- 2 # 1 = seen, 2 = not seen  ## this seems to lump the 'recovered dead' from L. 113 together with 'not seen', which seems strange??
table(as.matrix(CH))

# Bundle data
data_jags <- list(CH = CH, n_ind = dim(CH)[1], n_occ = dim(CH)[2], first = first, # basic data
                  AM = AM, sex = sex, sex.ratio = sex.ratio, TM = TM, n_tag = max(TM, na.rm = T), # covariate data 
                  FR = rep(1,dim(CH)[1]), ones = rep(1,dim(CH)[1])) # for likelihood marginalization

# Parameters
parameters <- c("s_age","s_age_sex", # survival per age and/or sex category
                "p_age", "p_tag", "p_age_tag", # detection per age and/or tag type
                "r_age", "r_tag", "r_age_tag") # recovery rate per age and/or tag type

# MCMC settings
{ni <- 1000
  nt <- 4
  nb <- ni*0.5
  nc <- 4}

## cannot run this because did not get model files
# #  Model selection 
# m1 <- jags(data_jags, inits = NULL, parameters, "models/model_CJS_v0.R", n.chains = nc, n.iter = ni, n.thin = nt, n.burnin = nb, parallel = T)
# m2 <- jags(data_jags, inits = NULL, parameters, "models/model_CJS_v0_Sex.R", n.chains = nc, n.iter = ni, n.thin = nt, n.burnin = nb, parallel = T)
# m3 <- jags(data_jags, inits = NULL, parameters, "models/model_CJS_v1_AM_TM.R", n.chains = nc, n.iter = ni, n.thin = nt, n.burnin = nb, parallel = T)
# m4 <- jags(data_jags, inits = NULL, parameters, "models/model_CJS_v2_TM.R", n.chains = nc, n.iter = ni, n.thin = nt, n.burnin = nb, parallel = T)
# m5 <- jags(data_jags, inits = NULL, parameters, "models/model_CJS_v3_Sex_AM_TM.R", n.chains = nc, n.iter = ni, n.thin = nt, n.burnin = nb, parallel = T)
# 
# DIC <- sapply(list(m1,m2,m3,m4,m5), function(x) x$DIC)
# dfDIC <- data.frame(model = paste0("m",1:length(DIC)),
#                     file = sapply(list(m1,m2,m3,m4,m5), function(x) x$modfile),
#                     DIC = DIC,
#                     dDIC = DIC - min(DIC)) %>% arrange(dDIC)
# dfDIC
# 
# # final model 
# m3
# plot(m3)
# 
# save(m, data_jags, file = "data/jags_CJS_output.RData")
# 
# 
# # quick plot for parameter estimates
# whiskerplot(m2, "s_age_sex")
# whiskerplot(m3, "s_age")
# whiskerplot(m3, "p_age_tag")
# whiskerplot(m4, "p_tag")







# 4b. Running Recapture-Recovery survival model in JAGS ####
library(jagsUI)
rm(list=ls())

# load data
load("data/jags_prep.RData")

# Bundle data
data_jags <- list(CH = CH, n_ind = dim(CH)[1], n_occ = dim(CH)[2], first = first, # basic data
                  AM = AM, sex = sex, sex.ratio = sex.ratio, TM = TM, n_tag = max(TM, na.rm = T), # covariate data 
                  FR = rep(1,dim(CH)[1]), ones = rep(1,dim(CH)[1])) # for likelihood marginalization

# Parameters
parameters <- c("s_age","s_age_sex", # survival per age and/or sex category
                "p_age", "p_tag", "p_age_tag", # detection per age and/or tag type
                "r_age", "r_tag", "r_age_tag") # recovery rate per age and/or tag type

# MCMC settings
{ni <- 5000
  nt <- 4
  nb <- ni*0.5
  nc <- 4}


#  Model selection 
# m1 <- jags(data_jags, inits = NULL, parameters, "models/model_RR_v0.R", n.chains = nc, n.iter = ni, n.thin = nt, n.burnin = nb, parallel = T)
# m2 <- jags(data_jags, inits = NULL, parameters, "models/model_RR_v0_Sex.R", n.chains = nc, n.iter = ni, n.thin = nt, n.burnin = nb, parallel = T)
# m3 <- jags(data_jags, inits = NULL, parameters, "models/model_RR_v1_AM_TM.R", n.chains = nc, n.iter = ni, n.thin = nt, n.burnin = nb, parallel = T)
m4 <- jags(data_jags, inits = NULL, parameters, "models/model_RR_v2_TM.R", n.chains = nc, n.iter = ni, n.thin = nt, n.burnin = nb, parallel = T)
# m5 <- jags(data_jags, inits = NULL, parameters, "models/model_RR_v3_Sex_AM_TM.R", n.chains = nc, n.iter = ni, n.thin = nt, n.burnin = nb, parallel = T)

DIC <- sapply(list(m1,m2,m3,m4,m5), function(x) x$DIC)
dfDIC <- data.frame(model = paste0("m",1:length(DIC)),
                    file = sapply(list(m1,m2,m3,m4,m5), function(x) x$modfile),
                    DIC = DIC,
                    dDIC = DIC - min(DIC))
dfDIC <- dfDIC[order(dfDIC$dDIC),]
dfDIC

# final model 
m3
m4
plot(m4)

save(m4, data_jags, file = "data/jags_RR_output.RData")


# quick plot for parameter estimates
whiskerplot(m2, "s_age_sex")
whiskerplot(m3, "s_age")
whiskerplot(m3, "p_age_tag")
whiskerplot(m3, "r_age_tag")

par(mfrow=c(3,1), cex = 1, bty = "l")
whiskerplot(m4, "s_age")
whiskerplot(m4, "p_tag")
whiskerplot(m4, "r_tag")


# Survival plot
load("data/jags_RR_output.RData")
m <- m4
jpeg("mean_survival_age.jpeg", quality = 100, res = 300, width = 5, height = 6.5, units = "in")
plot(1:2, m$mean$s_age, bty = "l", las = 1, xaxt = "n", xlab = "Age", ylab = "Annual survival", xlim = c(.5,2.5), ylim = c(0,1), pch = 16, cex = 1.4, col = c(2,4), cex.lab = 1.2)
segments(1:2, m$q97.5$s_age, 1:2, m$q2.5$s_age, col = c(2,4), lwd = 2)
axis(1, 1:2, c("Juv", "Adult"))
dev.off()



# 4c. Running Recapture-Recovery survival model in JAGS - one tag type per individual ####
library(jagsUI)
rm(list=ls())

# load data
load("data/jags_prep.RData")

# Bundle data
data_jags <- list(CH = CH, n_ind = dim(CH)[1], n_occ = dim(CH)[2], first = first, # basic data
                  AM = AM, sex = sex, sex.ratio = sex.ratio, TM = TM3, n_tag = max(TM3, na.rm = T), # covariate data 
                  FR = rep(1,dim(CH)[1]), ones = rep(1,dim(CH)[1])) # for likelihood marginalization

# Parameters
parameters <- c("s_age","s_age_sex", # survival per age and/or sex category
                "p_age", "p_tag", "p_age_tag", # detection per age and/or tag type
                "r_age", "r_tag", "r_age_tag") # recovery rate per age and/or tag type

# MCMC settings
{ni <- 1000
  nt <- 4
  nb <- ni*0.5
  nc <- 4}


#  Model selection 
m1 <- jags(data_jags, inits = NULL, parameters, "models/model_RR_v0.R", n.chains = nc, n.iter = ni, n.thin = nt, n.burnin = nb, parallel = T)
m2 <- jags(data_jags, inits = NULL, parameters, "models/model_RR_v0_Sex.R", n.chains = nc, n.iter = ni, n.thin = nt, n.burnin = nb, parallel = T)
m3 <- jags(data_jags, inits = NULL, parameters, "models/model_RR_v1_AM_TM3.R", n.chains = nc, n.iter = ni, n.thin = nt, n.burnin = nb, parallel = T)
m4 <- jags(data_jags, inits = NULL, parameters, "models/model_RR_v2_TM3.R", n.chains = nc, n.iter = ni, n.thin = nt, n.burnin = nb, parallel = T)
m5 <- jags(data_jags, inits = NULL, parameters, "models/model_RR_v3_Sex_AM_TM3.R", n.chains = nc, n.iter = ni, n.thin = nt, n.burnin = nb, parallel = T)

DIC <- sapply(list(m1,m2,m3,m4,m5), function(x) x$DIC)
dfDIC <- data.frame(model = paste0("m",1:length(DIC)),
                    file = sapply(list(m1,m2,m3,m4,m5), function(x) x$modfile),
                    DIC = DIC,
                    dDIC = DIC - min(DIC)) %>% arrange(dDIC)
dfDIC

# final model 
m3
plot(m3)

save(m4, data_jags, file = "data/jags_RR_output_TM3.RData")


# quick plot for parameter estimates
whiskerplot(m2, "s_age_sex")
whiskerplot(m3, "s_age")
whiskerplot(m3, "p_age_tag")
whiskerplot(m3, "r_age_tag")
whiskerplot(m4, "p_tag")
whiskerplot(m4, "r_tag")
