# Deconvolution function
# To deconvolute a death curve into an incidence curve
# Refernce: https://doi.org/10.1073/pnas.0902958106
# Input:
  # obsDeath: observed death curve
  # meanDistr, sdDistr: mean and standard deviation of gamma distribution of the time between infection and death
# Output:
  # incidence curve of the number of deaths result from an infection of a given day, NOT the total incidence
deathDeconv = function (obsDeath, meanDistr, sdDistr) {
  #Trim observed data
  obsDeath = obsDeath[complete.cases(obsDeath),]
  firstDeathDate = obsDeath[obsDeath$new_deaths !=0, ][1,"day_num"]
  obsDeath = obsDeath[which(obsDeath$day_num == firstDeathDate):nrow(obsDeath),c("new_deaths","day_num")]
  numOfDays = nrow(obsDeath)
  
  # Time-to-death distribution
  scaleDistr = sdDistr^2/meanDistr
  shapeDistr = meanDistr/scaleDistr
  toDeathDistr = numeric()
  day_prior = round(qgamma(0.95, shape = shapeDistr, scale = scaleDistr)) # estimated number of days from first incidence to available data
  for (ii in 1:(numOfDays + day_prior - 1)) {
    if (ii == 1) {
      toDeathDistr[ii] = pgamma(ii+0.5, shape = shapeDistr, scale = scaleDistr) - 
        pgamma(ii-1, shape = shapeDistr, scale = scaleDistr)
    } else {
      toDeathDistr[ii] = pgamma(ii+0.5, shape = shapeDistr, scale = scaleDistr) - 
        pgamma(ii-0.5, shape = shapeDistr, scale = scaleDistr)
    }
  }
  toDeathDistr[toDeathDistr == 0] = 1e-6
  
  # Initial values: shift death curve by meanDistr
  initialParam = obsDeath
  initialParam$day_num = initialParam$day_num-meanDistr
  deconvIncidence = data.frame(day_num = (obsDeath$day_num[1]-day_prior):(obsDeath$day_num[numOfDays]-1))
  deconvIncidence = merge(deconvIncidence, initialParam, by = "day_num", all.x = TRUE)
  colnames(deconvIncidence)[2] = "new_incidence"
  deconvIncidence[is.na(deconvIncidence)] = 1e-6
  
  # Fraction of deaths observed 
  fracObserved = numeric()
  for (jj in 1:nrow(deconvIncidence)) {
    fracObserved[jj] = 1 - pgamma(day_prior - jj + 0.5, shape = shapeDistr, scale = scaleDistr) -
      (1 - pgamma(nrow(deconvIncidence) +1 - jj + 0.5, shape = shapeDistr, scale = scaleDistr))
  }
  fracObserved[fracObserved == 0] = 1e-6
  
  # EM-algorithm
  #chi_sqr = 10000
  expDeath = numeric()
  nextParam = numeric()
  for (k in 1:10) { 
  #while(chi_sqr >= 1) {
    # Expected deaths
    for (ii in 1:numOfDays) {
      expDeath[ii] = t(rev(toDeathDistr[1:(day_prior-1+ii)])) %*% deconvIncidence[1:(day_prior-1+ii),2]
    }
    
    # Chi-square value
    #chi_sqr = sum((expDeath - obsDeath$new_deaths)**2/expDeath)/numOfDays
    
    # New set of parameters
    for (jj in 1:nrow(deconvIncidence)) {
      nextParam[jj] = (deconvIncidence[jj,2] / fracObserved[jj]) *
        sum(toDeathDistr[max(day_prior+1-jj, 1):(nrow(deconvIncidence)-jj+1)] * 
              obsDeath$new_deaths[max(jj-day_prior+1,1):numOfDays] / 
              expDeath[max(jj-day_prior+1,1):numOfDays])
    }
    deconvIncidence[,2] = nextParam
  }
  
  deconvIncidence$week_num = deconvIncidence$day_num %/% 7 + 1
  return(deconvIncidence)
}
