mcmcParallel = function(countryText, phyloMethod, subCladeSize, mcSteps, startWeek, startWeekInf,
                        allObsIncidence,allObsTest,allObsMutant,genTimeDistrD614,meanGenTimeD614,
                        sdGenTimeD614, parameters,parameterSteps,lowerLimit, upperLimit) {
  goodMC = FALSE
  pAccept = rep(0.5, length(parameters))
  numWhileLoop = 0
  numPars = length(parameters)
  
  source("calcMCMCTotalLogLikelihood.R")
  source("mcmcProposal.R")
  
  while (goodMC == FALSE) {
    numWhileLoop = numWhileLoop + 1
    iPar = zeros(mcSteps, numPars)
    iiLogL = zeros(mcSteps, 1)
    nAccept = zeros(1, numPars)
    parameterSteps[pAccept > 0.7] = parameterSteps[pAccept > 0.7]*1.2
    parameterSteps[pAccept < 0.3] = parameterSteps[pAccept < 0.3]/1.25
    parameterSteps = bsxfun(min, parameterSteps, upperLimit - lowerLimit)
    
    # Current loglikelihood
    parameters_c = parameters
    logLikelihood = calcMCMCTotalLogLikelihood(parameters_c,startWeek, startWeekInf,
                                             allObsIncidence,allObsTest,allObsMutant,
                                             genTimeDistrD614,meanGenTimeD614,sdGenTimeD614)
    for (tt in 1:mcSteps) {
      tempNew = mcmcProposal(parameters_c, parameterSteps, lowerLimit, upperLimit)
      for (ii in 1:numPars) {
        parameters_mc = parameters_c
        parameters_mc[ii] = tempNew[ii]
        logLikelihood_new = calcMCMCTotalLogLikelihood(parameters_mc,startWeek, startWeekInf,allObsIncidence,
                                                       allObsTest,allObsMutant,genTimeDistrD614,meanGenTimeD614,sdGenTimeD614)
        alpha = min(1,exp(logLikelihood_new - logLikelihood))
        uu = runif(1)
        if (uu <= alpha) {
          parameters_c = parameters_mc
          logLikelihood = logLikelihood_new
          nAccept[ii] = nAccept[ii] + 1
        }
      }
      iiLogL[tt,1] = logLikelihood
      iPar[tt,] = parameters_c
      # Check MCMC
      if (tt %% (mcSteps/100) == 0) {
        pAccept = nAccept/tt
        write.table(as.matrix(t(parameterSteps)), file = paste("mcmc_result/",phyloMethod,"/",subCladeSize,"/",countryText,"_parameter_step.csv", sep = ""),
                    col.names = FALSE, row.names = FALSE, sep = ",")
        write.table(iPar, file = paste("mcmc_result/",phyloMethod,"/",subCladeSize,"/",countryText,"_mcmc_res.csv", sep = ""),
                    col.names = FALSE, row.names = FALSE, sep = ",")
        write.table(iiLogL, file = paste("mcmc_result/",phyloMethod,"/",subCladeSize,"/",countryText,"_log_likelihood.csv", sep = ""),
                    col.names = FALSE, row.names = FALSE, sep = ",")
        print("MCMC")
        print(tt/mcSteps)
        print("Acceptance probability")
        print(pAccept)
        if (any(pAccept > 0.7) || any(pAccept < 0.3)){
          goodMC = FALSE
          if (numWhileLoop < 20) {
            if (tt/mcSteps >= 0.02) break
          } else goodMC = TRUE
        } else goodMC = TRUE
      }
    }
  }
}
