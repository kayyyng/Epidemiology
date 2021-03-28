calcMCMCTotalLogLikelihood = function(
  x0, startWeek, startWeekInf, allObsIncidence, allObsTest, allObsMutant, genTimeDistrD614, 
  meanGenTimeD614, sdGenTimeD614
) {
  relTransmiss = x0[1]
  relGentime = x0[2]
  rho_1 = x0[3:length(x0)]
  
  totalLogLikelihoodRelTransWeeklyAllCountry(
    startWeek, startWeekInf, allObsIncidence, allObsTest, allObsMutant, genTimeDistrD614, 
    meanGenTimeD614, sdGenTimeD614, relTransmiss, relGentime, rho_1)
}
