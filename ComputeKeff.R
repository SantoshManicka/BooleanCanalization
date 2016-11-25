# This program computes the effective connectivity of a Boolean function.

# USAGE EXAMPLE:

# Note: Make sure you also have 'ComputeDetectCubes.R' in the same folder as this file.

# k = 3
# func = c(1,1,1,1,0,1,1,1)
# ComputeKeff(func, k, ReturnSchemata = T)

# The following output is generated:
# $EssentialSchemata
# [,1] [,2] [,3]
# [1,]    1    0    0
# [2,]    0    2    2
# [3,]    2    2    1
# [4,]    2    1    2
# 
# $SchemaOutputs
# [1] 0 1 1 1
# 
# $Keff
# [1] 1.25

# The schema symbols are encoded according to the following mapping: '0' -> 0, '1' -> 1, '#' -> 2

#--------------------------------------------------------------------------------------------------------------------

source('ComputeDetectCubes.R')

ComputeKeff <- function(Func, K, ReturnSchemata = F) {
  
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
	  
	  AllCvrDims = c(AllCvrDims, Cubes$CvrDims)
	  
	} else Sch1 = NULL
	
	SubLUT = matrix(LUT[Func == 0, ], ncol = K)
	
	if (nrow(SubLUT) > 0) {
	
	  Cubes = DetectCubes(SubLUT, K)
	  
	  Sch0 = Cubes$Schemata
	  
	  AllCvrDims = c(AllCvrDims, Cubes$CvrDims)
	  
	} else Sch0 = NULL
		
	Keff = mean(K - AllCvrDims)
	
	if (ReturnSchemata) {
		
		AllSchemata = list(Sch0, Sch1)
		
		AllSchemata = matrix(unlist(AllSchemata), byrow = T, ncol = K)
		
		SchemaOutputs = c(rep(0, length(Sch0)), rep(1, length(Sch1)))
		
		return(list(EssentialSchemata = AllSchemata, SchemaOutputs = SchemaOutputs, Keff = Keff)) 
	
	} else return(Keff)

}





