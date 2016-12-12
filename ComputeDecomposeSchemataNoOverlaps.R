# This program decomposes a set of overlapping schemata into one of non-overlapping schemata, while preserving 
# the overall cover.

# Basic idea: 
# Iteratively subtract larger schemata from the smaller ones; any overlap is kept with the larger schema.

# Outline of the method:
# First arrange the schemata in decreasing order of the number of their wildcards (NumWCs). Then, go through the list
# and subtract each schema from the schemata in the rest of the list. Each subtraction step creates a new set of schemata
# that the list is updated with. Repeat the process until the (dynamically updated) list is exhausted, leaving a list with
# a list of schemata where no schema overlaps with any other schema. It is important to note that the LUT entries covered
# by the original set of schemata and the decomposed set are one and the same.

# Optional parameters:
# (1) Preserving the original NumWCs: The parameter 'ReturnOriginalNumWCs' is used for this purpose. The subtraction operations 
# decompose schemata into smaller parts; thus they lose their WCs. However, there may be a need to preserve their original 
# NumWCs, as in the case computing Keff from a set of schemata (see 'ComputeKeffFromSchemata.R'). Setting ReturnOriginalNumWCs
# to True returns the original NumWCs of each schema in the decomposed set.
# (2) Parallelization: The parameter 'NumCores' is used for this purpose. The subtraction operations can be parallelized 
# particularly when a large number of schemata are involved. Set 'NumCores' to the number of cores you want to run this 
# program simultaneously on.

# Input:
# A schemata matrix (rows contain schemata).

# Output:
# A list containing two items:
# (1) DecomposedSchemata: A matrix with the same number of columns as the original Schemata matrix, containing the decomposed
# schemata.
# (2) OrigNumWCs: A list containing the original NumWCs of the schemata in DecomposedSchemata.

# Applications: 
# (1) To compute the cardinality of the set of LUT entries covered by a set of overlapping schemata; we propose that this 
# method may be more efficient than the traditional inclusion-exclusion based methods.

# Remarks:
# Non-essential schemata will be lost by definition, since they are the "bridges" that overlap among the essential schemata.

# USAGE EXAMPLE:
# k = 3
# func = c(1,1,1,1,1,0,0,0)
# SchemataOutputs = ComputePrimeImplicants(func, k)
# Schemata = SchemataOutputs$PrimeImplicants[SchemataOutputs$Outputs == 1, ]
# DecomposedSchemata = DecomposeSchemataNoOverlaps(Schemata, ReturnOriginalNumWCs = F)  # don't return orignal NumWCs
# Output:
# > Schemata
# [,1] [,2] [,3]
# [1,]    0    2    2
# [2,]    2    0    0
# > DecomposedSchemata
# $DecomposedSchemata
# [,1] [,2] [,3]
# [1,]    0    2    2
# [2,]    1    0    0
# > DecomposeSchemataNoOverlaps(Schemata, ReturnOriginalNumWCs = T)  # return orignal NumWCs
# $DecomposedSchemata
# [,1] [,2] [,3]
# [1,]    0    2    2
# [2,]    1    0    0
# 
# $OrigNumWCs
# [1] 2 1

require(parallel)

source('ComputeSchemaSetOperations.R')

DecomposeSchemataNoOverlaps <- function(Schemata, MaxIters = Inf, ReturnOriginalNumWCs = F, NumCores = 1) {
  
  Subtract <- function(j, SchemataList, Schema) {
    
    WC = OrigWCList[j]
    
    FromSchema = SchemataList[[j]]
    
    RemSchList = SubtractSchema(FromSchema, Schema)
    
    if ( !is.null(RemSchList) ) {
      
      if ( !anyNA(RemSchList) ) {
        
        return(RemSchList)
        
      } else return(NULL)
      
    } else {
      
      return(list(FromSchema))
      
    }
    
    # NULL is returned if FromSchema is duplicate/subsumed 
    
  }
  
  # Rearrage schemata and the corresponding list of WCs in decreasing order of their size (tallied by the number of WCs)
  
  SchemataList = unlist(apply(Schemata, 1, function(v) list(v)), recursive = F)
  
  OrigWCList = sapply(SchemataList, function(s) sum(s==2))
  
  SortOrder = order(OrigWCList, decreasing = T)
  
  SchemataList = SchemataList[SortOrder]
  OrigWCList = OrigWCList[SortOrder]
  
  i = 1
  
  SchemataListLen = length(SchemataList)	
  
  while ( (i < SchemataListLen) & (i < MaxIters) ) {  # 'SchemataList' changes dynamically, so use 'while' instead of 'for'
    
    Schema = SchemataList[[i]]
    
    NewSchemataList = SchemataList[1 : i]
    
    if ( ReturnOriginalNumWCs ) NewOrigWCList = OrigWCList[1 : i]
    
    # Subtract schema at 'i' from all the schemata from position 'i'+1 onwards
    
    # If OriginalNumWCs is to be remembered, then the original order of schemata is important and must be retained
    # (see note below). Parallelizing the subtractions could meddle with that ordering, so it's better to avoid it.
    
    if ( ReturnOriginalNumWCs ) NumCores = 1
    
    RemSchLists = mclapply ((i + 1) : SchemataListLen, Subtract, SchemataList, Schema, mc.cores = NumCores)

    RemLen = sapply(RemSchLists, function(x) ifelse(is.null(x), 0, length(x)))

    IncludeIdx = (RemLen > 0)
    
    # Process only the non-NULL remainders from the subtraction above

    if (sum(IncludeIdx) > 0) {
      
      # Since a single subtraction operation could leave a remainder containing multiple schemata, the original NumWCs
      # list needs to be updated accordingly.
      
      if ( ReturnOriginalNumWCs ) {
        
        RemOrigWCList = mclapply(which(IncludeIdx), function(j) rep(OrigWCList[j + i], RemLen[j]), mc.cores = NumCores)
        RemOrigWCList = unlist(RemOrigWCList)
        
        NewOrigWCList = c(NewOrigWCList, RemOrigWCList)
        
      }
      
      # Make a list of the remainder schemata and append it to SchemataList

      RemSchLists = RemSchLists[IncludeIdx]
      
      RemSchListFlat = unlist(RemSchLists, recursive = F)
      
      NewSchemataList = append(NewSchemataList, RemSchListFlat)		
      
    }
    
    SchemataList = NewSchemataList
    
    if ( ReturnOriginalNumWCs ) OrigWCList = NewOrigWCList
    
    SchemataListLen = length(SchemataList)	
    
    if (i == SchemataListLen) break else i = i + 1   # 'i' cannot be greater than SchemataListLen
    
    # Pick a schema of the largest size from among those at the new position 'i' onwards and place it at 'i',
    # so that it becomes the schema to be subtracted from the schemata in the rest of the list.
    # NOTE: The motivation behind this logic is that preserving the larger schemata in the subtractions lowers
    # the computational expense. However, if the OriginalNumWCs are to be remembered, then this logic won't 
    # work, as it could meddle with the same. For instance, consider the schemata {##00, 1#0#, 11#1}. Start with
    # ##00 for the subtractions. Subtract(1#0#, ##00) -> 1#01; Subtract(11#1, ##00) -> NULL. This leaves a new
    # set: {##00, 1#01, 11#1}. OrigNumWCs = {2,2,1}. In the second round, we pick one of two schemata with only 
    # one '#'. If for reason, we happened to pick 11#1 and subtract it from 1#01, then that would leave the set
    # {##00, 11#1, 1001} and OrigNumWCs = {2,1,2}. Whereas, what we really want to the result to be is the following:
    # {##00, 1#01, 1111} and OrigNumWCs = {2,2,1}. In summary, ReturnOriginalNumWCs = T may lower the efficiency of the
    # decomposition, but it is a necessary expense.
    
    if ( !ReturnOriginalNumWCs ) {
      
      SchemataSubList = SchemataList[i : SchemataListLen]
      
      CurrWCSubList = sapply(SchemataSubList, function(s) sum(s==2))
      
      MaxCurrWCIdx = match(max(CurrWCSubList), CurrWCSubList)
      
      MaxCurrWCSchema = SchemataSubList[[MaxCurrWCIdx]]
      temp = SchemataList[[i]]
      SchemataList[[i]] = MaxCurrWCSchema
      SchemataList[[MaxCurrWCIdx + i - 1]] = temp
    
      rm(SchemataSubList)
      rm(CurrWCSubList)
      
    }
      
  }
  
  SchemataList = SchemataList[1 : i]
  
  K = length(SchemataList[[1]])
  
  DecomposedSchemata = matrix(unlist(SchemataList), ncol = K, byrow = T)
  
  OrigNumWCs = OrigWCList[1 : i]
  
  if ( ReturnOriginalNumWCs ) return (list(DecomposedSchemata = DecomposedSchemata, OriginalNumWCs = OrigNumWCs))
  
  return (list(DecomposedSchemata = DecomposedSchemata))
  
}





















