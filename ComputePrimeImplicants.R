# This program computes the prime implicants of a Boolean function, also providing the option of all computing 
# all prime implicants.
# It proceeds by first computing the prime implicants that are essential to compute the Keff just as in 'ComputerKeff.R'.
# This set consists of the essential prime implicants, and *may* contain others (the generation of the latter is not 
# systematically controlled). If all prime implicants are requested, then this set is passed on to a function known as
# 'UnionSchemataSet' defined in 'ComputeSchemaSetOperations.R' that performs logical disjunction operations on it and in
# the process identifies the other non-essential prime implicants. In essence, the whole process may be thought of as the
# Quine-McCluskey procedure done in reverse.

# USAGE EXAMPLES:

# Note: Make sure you also have 'ComputeDetectCubes.R' and 'ComputeSchemaSetOperations.R' in the same folder as this file.

# K = 3
# Func = c(0,1,1,1,1,1,1,0)

# > ComputePrimeImplicants(Func, K, AllPrimeImplicants = F)
# $PrimeImplicants
# [,1] [,2] [,3]
# [1,]    0    0    0
# [2,]    1    1    1
# [3,]    0    2    1
# [4,]    2    0    1
# [5,]    0    1    2
# [6,]    2    1    0
# [7,]    1    0    2
# 
# $Outputs
# [1] 0 0 1 1 1 1 1

# > ComputePrimeImplicants(Func, K, AllPrimeImplicants = T)
# $PrimeImplicants
# [,1] [,2] [,3]
# [1,]    0    0    0
# [2,]    1    1    1
# [3,]    1    2    0
# [4,]    0    2    1
# [5,]    2    0    1
# [6,]    0    1    2
# [7,]    2    1    0
# [8,]    1    0    2
# 
# $Outputs
# [1] 0 0 1 1 1 1 1 1

# K = 5
# Func = c(1,0,1,0,0,0,0,0, 1,0,1,0,0,0,1,1, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,1,1)

# > ComputePrimeImplicants(Func, K, AllPrimeImplicants = F)
# $PrimeImplicants
# [,1] [,2] [,3] [,4] [,5]
# [1,]    2    0    2    2    1
# [2,]    2    2    0    2    1
# [3,]    2    2    2    0    1
# [4,]    2    0    1    2    2
# [5,]    2    2    1    0    2
# [6,]    1    0    2    2    2
# [7,]    1    2    0    2    2
# [8,]    0    2    0    2    0
# [9,]    2    1    1    1    2
# 
# $Outputs
# [1] 0 0 0 0 0 0 0 1 1

# > ComputePrimeImplicants(Func, K, AllPrimeImplicants = T)
# $PrimeImplicants
# [,1] [,2] [,3] [,4] [,5]
# [1,]    1    2    2    0    2
# [2,]    2    0    2    2    1
# [3,]    2    2    0    2    1
# [4,]    2    2    2    0    1
# [5,]    2    0    1    2    2
# [6,]    2    2    1    0    2
# [7,]    1    0    2    2    2
# [8,]    1    2    0    2    2
# [9,]    0    1    2    1    0
# [10,]    0    2    0    2    0
# [11,]    2    1    1    1    2
# 
# $Outputs
# [1] 0 0 0 0 0 0 0 0 1 1 1

# The schema symbols are encoded according to the following mapping: '0' -> 0, '1' -> 1, '#' -> 2

#--------------------------------------------------------------------------------------------------------------------

source('ComputeDetectCubes.R')
source('ComputeSchemaSetOperations.R')

ComputePrimeImplicants <- function(Func, K, AllPrimeImplicants = F) {
  
  # Convert a decimal to a binary of a specified length. Max allowed value for len = 32.
  
  dec2bin <- function(x, len) {
    
    sfx = as.integer(strsplit(intToBin(x),'')[[1]])
    
    c(rep(0, len - length(sfx)), sfx)  
    
  }
  
  # LUT = t(sapply(0: ((2^K) - 1), function(i) rev(dec2bin(i, K))))  # using 'dec2bin' defined in package 'BoolNet'
  LUT = t(sapply(0: ((2^K) - 1), function(i) dec2bin(i, K)))
  
  if (K == 1) LUT = matrix(LUT, ncol = K)
  
  AllCvrDims = c()
  
  SubLUT = matrix(LUT[Func == 1, ], ncol = K)
  
  if (nrow(SubLUT) > 0) {
    
    Cubes = DetectCubes(SubLUT, K)
    
    Sch1 = Cubes$Schemata
    
    Sch1 = matrix(unlist(Sch1), byrow = T, ncol = K)
    
    if (AllPrimeImplicants) {
      
      SchSet1 = list()
      SchSet1[[2]] = Sch1
      
      # Comments above the definition of 'UnionSchemataSet' describe the meaning of the parameters
      
      Sch1 = UnionSchemataSet(SchemataSet = SchSet1, MaxWCLossProp = 1, MaxWCLoss = Inf, InitNewSet = T)
      
    }
    
    NumSch1 = nrow(Sch1)
    
  } else {
    
    Sch1 = NULL
    
    NumSch1 = 0
    
  }
  
  SubLUT = matrix(LUT[Func == 0, ], ncol = K)
  
  if (nrow(SubLUT) > 0) {
    
    Cubes = DetectCubes(SubLUT, K)
    
    Sch0 = Cubes$Schemata
    
    Sch0 = matrix(unlist(Sch0), byrow = T, ncol = K)
    
    if (AllPrimeImplicants) {
      
      SchSet0 = list()
      SchSet0[[2]] = Sch0
      
      # Comments above the definition of 'UnionSchemataSet' describe the meaning of the parameters
      
      Sch0 = UnionSchemataSet(SchemataSet = SchSet0, MaxWCLossProp = 1, MaxWCLoss = Inf, InitNewSet = T)
      
    }
    
    NumSch0 = nrow(Sch0)
    
  } else {
    
    Sch0 = NULL
    
    NumSch0 = 0
    
  }
  
  PrimeImplicants = rbind(Sch0, Sch1)
  
  Outputs = c(rep(0, NumSch0), rep(1, NumSch1))
  
  return(list(PrimeImplicants = PrimeImplicants, Outputs = Outputs)) 
  
}





