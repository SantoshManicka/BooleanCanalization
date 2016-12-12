# This program computes the Keff of a Boolean function when its schemata are supplied rather than the
# function itself (specified by its output vector); compare with 'ComputeKeff.R' which requires the latter.

# Basic idea:
# Decompose or split the schemata going to each output into a set of non-overlapping schemata while remembering the
# original number of their wildcard symbols (WCs). It is then straightforward to use this information to compute Keff.

# Input:
# 'SchemataParts' -- a list containing two items; each item contains a matrix of schemata mapping to an output, where 
# the first matrix maps to output 0, and the second matrix maps to output 1.

# Output:
# A scalar Keff value of the Boolean function that the 'SchemataParts' defines.

# USAGE EXAMPLE:
# K = 4
# Func = sample(c(0,1), 2^K, replace = T)
# SchemataParts = list()
# SchemataOutputs = ComputePrimeImplicants(Func, K)
# SchemataMatrix = SchemataOutputs$PrimeImplicants
# Outputs = SchemataOutputs$Outputs
# SchemataParts[[1]] = matrix(SchemataMatrix[Outputs == 0, ], ncol = K)
# SchemataParts[[2]] = matrix(SchemataMatrix[Outputs == 1, ], ncol = K)
# Keff = ComputeKeffFromSchemata(SchemataParts)
# print(Keff)
# Keff = ComputeKeff(Func, K)  # Compare with Keff calculated from the function itself
# print(Keff)

source('ComputeDecomposeSchemataNoOverlaps.R')

# Note that the Keff computed from a set of schemata that is not minimal is also not minimal.


ComputeKeffFromSchemata <- function(SchemataParts) {
  
  N = ncol(SchemataParts[[1]])
  
  Sch0 = SchemataParts[[1]]  # output = 0
  
  if ( !anyNA(Sch0) ) {
    
    DecomposedSchemataWC = DecomposeSchemataNoOverlaps(Sch0, ReturnOriginalNumWCs = T, NumCores = 1)
    
    DecomposedSchemata0 = DecomposedSchemataWC$DecomposedSchemata
    
    OrigNumWCs0 = DecomposedSchemataWC$OriginalNumWCs
    
  } else {
    
    DecomposedSchemata0 = NULL
    OrigNumWCs0 = NULL
    
  }
  
  Sch1 = SchemataParts[[2]]  # output = 1
  
  if ( !anyNA(Sch1) ) {
    
    DecomposedSchemataWC = DecomposeSchemataNoOverlaps(Sch1, ReturnOriginalNumWCs = T, NumCores = 1)
    
    DecomposedSchemata1 = DecomposedSchemataWC$DecomposedSchemata
    
    OrigNumWCs1 = DecomposedSchemataWC$OriginalNumWCs
    
  } else {
    
    DecomposedSchemata1 = NULL
    OrigNumWCs1 = NULL
    
  }
  
  DecomposedSchemata = rbind(DecomposedSchemata0, DecomposedSchemata1)
  
  OrigNumWCs = c(OrigNumWCs0, OrigNumWCs1)
  
  CurrNumWCs = apply(DecomposedSchemata, 1, function(s) sum(s == 2))
  
  DecomposedSchemataCardinality = 2^CurrNumWCs
  
  OrigSchemataKeffVals = N - OrigNumWCs
  
  TotalSchemataKeff = sum(OrigSchemataKeffVals * DecomposedSchemataCardinality)
  
  # Note that 'TotalCardinality' below must be equal to 2^K (K is the in-degree), if the schemata cover all of 
  # the LUT entries.
  
  TotalSchemataCardinality = sum(DecomposedSchemataCardinality)
  
  AvgKeff = TotalSchemataKeff / TotalSchemataCardinality
  
  return(AvgKeff)
  
}



















