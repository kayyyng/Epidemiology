#library
library(dplyr)
library(reshape2)
library(pracma)
library(NlcOptim)

#source other R scripts
source("deathDeconv.R")
source("totalLogLikelihoodRelTransWeeklyAllCountry.R")
source("negTotalLogLikelihoodAllCountry.R")
source("mcmcParallel.R")

# Relative fitness of D614G

# 1. Data
countryList = c("Australia", "Belgium", "Denmark", "India",
                "Netherlands", "Spain", "Portugal", "United-Kingdom", "USA")
countryListConfirmed = c("Australia", "Belgium", "Denmark", "India", 
                         "Netherlands", "Spain", "Portugal", "United Kingdom", "United States")

# Fixed generation time distribution
meanGenTime = 5.4
sdGenTime = 3.8
scaleGenTime = sdGenTime^2/meanGenTime
shapeGenTime = meanGenTime/scaleGenTime

# Fixed infection-to-death distribution
meanInfToDeath = 28
sdInfToDeath = 8.4

# Date zero of the entire database
dateZero = as.Date("2019/12/21")

# phylo methods
phyloMethodsArr = "strict"

# sensitivity analysis for subclade size
subCladeSizeArr = c("equal_more_2", "equal_more_3", "equal_more_5")

# Inference
for (iiPhylo in 1:length(phyloMethodsArr)) {
  iiPhyloMethods = phyloMethodsArr[iiPhylo]
  for (iiSubSize in 1:length(subCladeSizeArr)) {
    iiSubCladeSize = subCladeSizeArr[iiSubSize]
    genTimeDistrD614 = numeric()
    startWeekInf = numeric()
    allObsIncidence = data.frame(matrix(ncol=0,nrow=0)) 
    allObsTest = data.frame(matrix(ncol=0,nrow=0)) 
    allObsMutant = data.frame(matrix(ncol=0,nrow=0)) 
    
    for (iiCountry in 1:length(countryList)) {
      # Country string
      countryStr = countryList[iiCountry]
      
      # Data
      ## Proxy by confirmed death
      obsDailyIncidence = read.csv(paste("../0. data/confirmed_by_country/owid_death_",
                                         countryListConfirmed[iiCountry], ".csv", sep = ""))
      obsDailyIncidence = deathDeconv(obsDailyIncidence, meanInfToDeath, sdInfToDeath)
      #obsDailyIncidence$week_num = obsDailyIncidence$week_num - (1 + 4)
      
      ## Proxy by confirmed case
      #obsDailyIncidence = read.csv(paste("../0. data/confirmed_by_country/owid_", countryListConfirmed[iiCountry], ".csv", sep = ""))
      
      # D and G data
      if (iiPhyloMethods == "relaxed") {
        obsLabData = read.csv("../0. data/2020_08_17_cladeABC/classification_ABC.csv")
      } else if (iiPhyloMethods == "strict") {
        obsLabData = read.csv("../0. data/2020_08_25_cladeABC/classification_allcountries_strict.csv")
      }
      geoData = read.csv("../0. data/geo2_sequence_data_d614g.csv")
      obsLabData = obsLabData %>% merge(geoData, by = c("sequence", "sample_date", "mutant", "country")) %>%
        filter(country == countryStr)
      
      # select "B" and "C"
      obsLabData = obsLabData %>% filter(!!as.symbol(iiSubCladeSize) %in% c("B", "C"))
      #obsLabData = obsLabData %>% filter(!!as.symbol(iiSubCladeSize) == "C")
      
      # adjust time frame 
      obsLabData$date_num = as.integer(as.Date(obsLabData$sample_date, "%m/%d/%Y") - dateZero)
      obsLabData$week_num = floor(obsLabData$date_num/7)+1
      
      # count incidence of D614 and G614 by week
      obsWeeklyLabData = data.frame(week_num = min(obsLabData$week_num):max(obsLabData$week_num))
      obsCountLabData = dcast(obsLabData, week_num~mutant, fun.aggregate = length)
      obsWeeklyLabData = obsWeeklyLabData %>% merge(obsCountLabData, by = "week_num", all.x = TRUE)
      obsWeeklyLabData[is.na(obsWeeklyLabData)] = 0
      
      # correct reporting delay
      obsWeeklyLabData$week_num = obsWeeklyLabData$week_num-1
      
      startWeekInf[iiCountry] = min(obsWeeklyLabData[obsWeeklyLabData$G614 >0,"week_num"])
      
      # Start week
      startWeek = startWeekInf[iiCountry]
      obsIncidence = obsDailyIncidence[obsDailyIncidence$week_num >= startWeek, "new_incidence"]
      #obsIncidence = obsDailyIncidence[obsDailyIncidence$week_num >= startWeek, "new_cases"]
      obsIncidence[obsIncidence == 0] = 1e-6
      obsMutant = obsWeeklyLabData[obsWeeklyLabData$week_num >= startWeek, c("D614","G614")]
      obsMutant$D614 = obsMutant$D614 + obsMutant$G614
      allObsIncidence[1:length(obsIncidence), iiCountry] = obsIncidence
      allObsIncidence[is.na(allObsIncidence)] = 0
      allObsTest[1:(nrow(obsMutant)-3), iiCountry] = obsMutant[1:(nrow(obsMutant)-3), 1]
      allObsTest[is.na(allObsTest)] = 0
      allObsMutant[1:(nrow(obsMutant)-3), iiCountry] = obsMutant[1:(nrow(obsMutant)-3), 2]
      allObsMutant[is.na(allObsMutant)] = 0
      
      # Fixed generation time distribution
      for (ii in 1:nrow(obsDailyIncidence)) {
        if (ii == 1) {
          genTimeDistrD614[ii] = pgamma(ii+0.5, shape = shapeGenTime, scale = scaleGenTime) - 
            pgamma(ii-1, shape = shapeGenTime, scale = scaleGenTime)
        } else {
          genTimeDistrD614[ii] = pgamma(ii+0.5, shape = shapeGenTime, scale = scaleGenTime) - 
            pgamma(ii-0.5, shape = shapeGenTime, scale = scaleGenTime)
        }
      }
    }
    
    # Log likelihood function
    relTransmiss = 1.1
    relGentime = 0.95
    rho_1 = rep(0.1, length(countryList))
    x0 = c(relTransmiss,relGentime,rho_1)
    x0LowerBound = c(0.01, 1e-10, rep(1e-10,length(countryList)))
    x0UpperBound = c(100, 1, rep(1,length(countryList)))
    
    tic()
    print(totalLogLikelihoodRelTransWeeklyAllCountry(
      min(startWeekInf), startWeekInf,allObsIncidence,allObsTest,allObsMutant,genTimeDistrD614,
      meanGenTime,sdGenTime,relTransmiss,relGentime,rho_1))
    print(negTotalLogLikelihoodAllCountry(
      x0, min(startWeekInf), startWeekInf, allObsIncidence, allObsTest, allObsMutant, genTimeDistrD614, 
      meanGenTime, sdGenTime))
    toc()
    
    #MLE
    #redeffun = function(x) negTotalLogLikelihoodAllCountry(x, min(startWeekInf), startWeekInf, allObsIncidence, allObsTest, allObsMutant, genTimeDistrD614, meanGenTime, sdGenTime)
    #fmin = fmincon(x0, redeffun, lb = x0LowerBound, ub = x0UpperBound, maxfeval = 30000)
    fminParameter = x0
    
    #MCMC
    dir.create(file.path("mcmc_result", iiPhyloMethods, iiSubCladeSize), showWarnings = FALSE)
    mcSteps = 30000
    stepSize = bsxfun(min, 0.05*(x0UpperBound-x0LowerBound), 0.05)
    mcmcParallel('all_country',iiPhyloMethods,iiSubCladeSize,mcSteps, min(startWeekInf),
                 startWeekInf,allObsIncidence,allObsTest,allObsMutant,genTimeDistrD614,
                 meanGenTime,sdGenTime,fminParameter,stepSize,x0LowerBound,x0UpperBound)
  }
}
