totalLogLikelihoodRelTransWeeklyAllCountry = function(
  startWeek, startWeekInf, allObsIncidence, allObsTest, allObsMutant, genTimeDistrD614, 
  meanGenTimeD614, sdGenTimeD614, relTransmiss, relGentime, rho_1
) {
  # gentime distr of mutant
  meanGenTimeG614 = relGentime*meanGenTimeD614
  sdGenTimeG614 = sdGenTimeD614/meanGenTimeD614*meanGenTimeG614
  scaleMutant = sdGenTimeG614^2/meanGenTimeG614
  shapeMutant = meanGenTimeG614/scaleMutant
  
  genTimeDistrG614 = numeric()
  for (ii in 1:length(genTimeDistrD614)) {
    if (ii == 1) {
      genTimeDistrG614[ii] = pgamma(ii+0.5, shape = shapeMutant, scale = scaleMutant) - 
        pgamma(ii-1, shape = shapeMutant, scale = scaleMutant)
    } else {
      genTimeDistrG614[ii] = pgamma(ii+0.5, shape = shapeMutant, scale = scaleMutant) - 
        pgamma(ii-0.5, shape = shapeMutant, scale = scaleMutant)
    }
  }
  
  totalLogL = numeric()
  for (iiCountry in 1:length(allObsIncidence[1,])) {
    iiStartWeekInf = startWeekInf[iiCountry]
    startIdx = iiStartWeekInf-startWeek+1
    obsIncidence = allObsIncidence[(7*(startIdx-1)+1):nrow(allObsIncidence), iiCountry]
    obsLabData = cbind(allObsTest[startIdx:nrow(allObsTest), iiCountry], 
                       allObsMutant[startIdx:nrow(allObsTest), iiCountry])
    # Relative transmissibility p(t)
    rho_d = numeric()
    for (ii in 1:length(obsIncidence)) {
      if (ii == 1) {
        rho_d[ii] = rho_1[iiCountry]
      } else {
        XG614 = relTransmiss*(t(rev(genTimeDistrG614[1:(ii-1)])) %*% (obsIncidence[1:(ii-1)]*rho_d[1:(ii-1)]))
        XD614 = t(rev(genTimeDistrD614[1:(ii-1)])) %*% (obsIncidence[1:(ii-1)]*(1-rho_d[1:(ii-1)]))
        rho_d[ii] = XG614/(XG614 + XD614)
      }
    }
    # Weekly likelihood
    sumTest = numeric() ; sumRes = numeric()
    for (jj in 1:(length(obsIncidence)/7)) {
      sumTest[jj] = sum(obsIncidence[((jj-1)*7+1):(jj*7)])
      sumRes[jj] = t(obsIncidence[((jj-1)*7+1):(jj*7)])%*%(rho_d[((jj-1)*7+1):(jj*7)])
    }
    weeklyAvgRho = sumRes/sumTest
    weeklyAvgRho = weeklyAvgRho[1:min(length(weeklyAvgRho),nrow(obsLabData))]
    obsLabData = obsLabData[1:min(length(weeklyAvgRho),nrow(obsLabData)),]
    iiLikelihood = obsLabData[obsLabData[,1]>0,2]*log(weeklyAvgRho[obsLabData[,1]>0]) +
      (obsLabData[obsLabData[,1]>0,1]-obsLabData[obsLabData[,1]>0,2])*log(1-(weeklyAvgRho[obsLabData[,1]>0]))
    totalLogL[iiCountry] = sum(iiLikelihood)
  }
  return(sum(totalLogL))
}
