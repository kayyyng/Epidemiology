mcmcProposal = function(parameters_c,parameterSteps,lowerRange,upperRange){
  parameters_new = rep(-999, length(parameters_c))
  while (any(parameters_new < lowerRange) || any(parameters_new > upperRange)){
    parameters_new = parameters_c + parameterSteps*(runif(length(parameters_c))-0.5)*2
    parameters_new[parameters_new<lowerRange] = lowerRange[parameters_new<lowerRange] + (lowerRange[parameters_new<lowerRange] - parameters_new[parameters_new<lowerRange])
    parameters_new[parameters_new>upperRange] = upperRange[parameters_new>upperRange] - (parameters_new[parameters_new>upperRange] - upperRange[parameters_new>upperRange])
  }
  return(parameters_new)
}
