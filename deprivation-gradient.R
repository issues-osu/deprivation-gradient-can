
library(rgdal)
library(pbugs)
library(sociome)
library(sp)
library(sf)
library(spdep)
library(dplyr)
library(tidyverse)
library(tigris)
library(tidycensus)

#Remove scientific notation
options(scipen = 999)

#Points to Winbugs Directory
bugs_dir <- "C:/Users/barboza-salerno.1/OneDrive - The Ohio State University/Documents/winbugs/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14/"
setwd("C:/Users/barboza-salerno.1/OneDrive - The Ohio State University/Documents/winbugs/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14")

#Download ADI
BC_adi_ind <- get_adi("tract", state = "NM",  county = "Bernalillo", keep_indicators = TRUE, year = 2010, geometry = FALSE)

#Get factor loadings for ADI
attr(BC_adi_ind$ADI, "loadings")

#Summarize indicators
BC_adi_ind %>% select(c(7:21)) %>%
  dplyr::summarise(across(
    .cols = is.numeric, 
    .fns = list(Mean = mean, SD = sd), na.rm = TRUE, 
    .names = "{col}_{fn}"
  )) %>% t(.)

#Create quintiles
Q <- quantile(BC_adi_ind$ADI, c(.20, .40, .60, .80, 1), na.rm=TRUE) 

#Recode ADI into quintiles
(BC_adi_ind <- BC_adi_ind %>%
    mutate(ADI_Q = case_when(
      ADI <= Q[1]    ~ "Q20",
      (ADI > Q[1]    & ADI <= Q[2] )  ~ "Q40",
      ADI > Q[2]  & ADI <= Q[3]  ~ "Q60",
      ADI > Q[3]  & ADI <= Q[4]  ~ "Q80",
      ADI > Q[4]  & ADI <= Q[5]     ~ "Q100"
    ))) 

#Summarize ADI
BC_adi_ind %>% 
  summarise(mean = mean(ADI), sd = sd(ADI), Q2 = median(ADI), min = min(ADI), max = max(ADI), IQR = IQR(ADI)) 

#Download census tracts for county
BC_sf <- tracts(state = "NM", county = "Bernalillo", year = 2010)

#Order the tracts spatially, all data must be ordered in same way
BC_sf <- BC_sf[order(BC_sf$GEOID10),]

#Create spatialpolygon object
BC_sp <- as(BC_sf, "Spatial")

#create contiguity neighbors
BC_nb <- poly2nb(BC_sp)
BC_wb <- nb2WB(BC_nb)

#load in observed and expected counts
obs <- read.csv("obs.csv") %>% select(pn) %>% as.matrix()
exp <- read.csv("exp.csv")%>% select(pn)%>% as.matrix() 

#Set rownames of matrix to census tract geography
rownames(obs) <-BC_sp@data$GEOID
rownames(exp) <-BC_sp@data$GEOID 

#First column has counts of physical abuse, second has counts of neglect
obs <- obs[,1]
exp <- exp[,1] 

#Subset ADI
BC_adi_sub <-  BC_adi_ind  %>% 
  select(ADI) %>%  
  rename(Deprivation  = ADI ) %>% 
  as.matrix() 

Deprivation <- BC_adi_sub[,1] 

### Model code provided by Martinez-Benito & Botella-Rocamora (2018) Disease Mapping: From Foundations to Multidimensional Modeling (see Book for code to implement stepwise models)

DICPoisson = function(Simu.sSMRs, O, E) {
  mu = t(apply(Simu.sSMRs/100, 1, function(x) {
    x * E
  }))
  D = apply(mu, 1, function(x) {
    -2 * sum(O * log(x) - x - lfactorial(O))
  })
  Dmean = mean(D)
  mumean = apply(Simu.sSMRs/100, 2, mean) * E
  DinMean = -2 * sum(O * log(mumean) - mumean - lfactorial(O))
  # if(save==TRUE){return(c(Dmedia,Dmedia-DenMedia,2*Dmedia-DenMedia))}
  cat("D=", Dmean, "pD=", Dmean - DinMean, "DIC=", 2 * Dmean - DinMean, "\n")
}

cuts.5 <- quantile(Deprivation, prob = c(0, 0.2, 0.4, 0.6, 0.8, 1))
Deprivation.5 <- as.numeric(cut(Deprivation, cuts.5, include.lowest = TRUE)) 

RegEcoLinear = function(){
  for(i in 1:n){
    O[i] ~ dpois(lambda[i])
    log(lambda[i]) <- log(E[i])+log.theta[i]
    log.theta[i] <- mu+beta*Deprivation[i]+sd.sp*sp[i]+sd.h*het[i]
    sSMR.withoutDep[i] <- 100*exp(mu+sd.sp*sp[i]+sd.h*het[i])
    sSMR.withDep[i] <- 100*exp(log.theta[i])
    het[i] ~ dnorm(0,1)
  }
  sp[1:n] ~ car.normal(adj[],weights[],num[],1)  
  
  for (i in 1:nVec){
    weights[i] <- 1
  }
  
  #Prior distributions
  mu ~ dflat()
  sd.sp ~ dunif(0,5)
  sd.h ~ dunif(0,5)
  beta ~ dflat()
}

data = list(O = obs, E = exp, Deprivation = Deprivation, n = length(obs), 
            nVec = length(BC_wb$adj), 
            num = BC_wb$num, adj = BC_wb$adj)
inits = function() {
  list(mu = rnorm(1), beta = rnorm(1), sd.sp = runif(1), 
       sd.h = runif(1), sp = rnorm(length(obs)), 
       het = rnorm(length(obs)))
}
param = c("sSMR.withDep", "sSMR.withoutDep", "mu", "beta", "sd.sp", "sd.h")
ResulLin = pbugs(data = data, inits = inits, parameters = param, model.file = RegEcoLinear, n.iter = 31000, n.burnin = 30000,  DIC = F, bugs.seed = 10, bugs.dir = bugs_dir) 

# Once model finishes the parameters can be extracted as follows
exp(ResulLin$mean$beta*Q[1]) # etc. for each ADI quintile
exp(round(ResulLin$summary["beta", ], 3))
round(ResulLin$summary["sd.sp", ], 3)
round(ResulLin$summary["sd.h", ], 3) 


