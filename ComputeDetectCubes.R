# BASIC IDEA:
# Rather than consider all possible subcubes and verify if each one qualifies as a prime implicant (PI), we
# first consider sets of LUT entries that could potentially contain one or more subcubes. We identify such 
# sets based on the number of 1s of each input value of that set. 
# For a set of LUT entries to constitute a D-dimensional subcube, the distribution of the numbers of 1s associated 
# with all of the 2^D corners of the subcube must match a particular distribution. Let's say that the bottom left (BL) 
# corner of the subcube contains an LUT entry with number of 1s equal to X. If it is the BL of a subcube (PI), then it must
# be the case that the top right (TR) corner of the same subcube constains X + D number of 1s. Further, there must be 
# exactly C(D, X+Ni) number of corners with (X+Ni) number of 1s. Thus, there is a specific (binomial) distribution of the number of 1s 
# in the LUT entries associated with the corners of a subcube. In other words, it is a necessary condition for a set of LUT
# entries to constitue a subcube.

# DETAILS OF THE METHOD: 
# We first identify a set of LUT entries that potentially contains one or more D-dimensional subcubes. The associated necessary condition 
# is that the set contains a net distribution of the number of 1s must be equal to or greater than the one required (see above).
# Within such a set, we then proceed to identify every D-dimensional cube that exists. The procedure is as follows.
# We first identify all possible BL-TR corners for a given D. For an LUT entry to be a potential the BL of a D-subcube, it must
# contain at least D 0s. Likewise, a potential TR of a D-subcube must contain at least D 1s. For example, let D=3. Then, (01000)
# and (10100) are potential BL, whereas (11100) is not. Likewise, (01111) and (11111) are potential TR, whereas (11000) is not.
# For a particular subcube, there can be only one BL-TR pair. Therefore, the next step is to identify such pairs. The necessary
# and sufficient conditions for such a pair to exist are: 
# (1) the difference in the number of 1s in the TR and that in the BL should be D;
# (2) the number of mismatches between the TR and BL should be equal to D.
# REMARK 1 below explains this in a little more detail.
# For every valid BL-TR pair identified above, we identify the 'fixed' variables and their values (the variables whose
# values match in the BL-TR pair) of the potential subcube. For example, if (10000) and (11011) constitute the BL-TR pair, then
# the 'fixed' variables are i1 and i3 with value 1 and 0 respectively. Then, we extract all other LUT entries whose values corresponding
# to the fixed variables inferred above match the inferred values. For the joint set consisting of BL, TR and the set identified in this
# step to constitute a subcube, it must be the case that its cardinality is equal to 2^D. 
# This completes the series of steps required to identify a single D-dimensional subcube. In a similar manner, we proceed to 
# identify the other subcubes, completing the series of steps required to identify all subcubes with dimension D.
# We repeat the above steps starting from the maximum D possible down to D = 1. If all parts of a subcube are already covered by
# larger subcubes or those of the same dimension, then we ignore it. That is, we only identify the essential prime implicants.

# REMARK 1: 
# Consider two binary vectors, v1 and v2, with the number of 1s in v2 equal number of 1s in v1 + D,
# where D is the dimension of the subcube. We will represent it as S(v2) = S(v1) + D; 'S' standing for sum.
# Assumption: S(v2) = S(v1) + D
# This implies that:
# (1) There must be at least D mismatches between v1 and v2. 
# (2) Those minimal mismatches must correspond to some D positions of v1 containing 0s and the same D positions of v2 containing 1s.
# Proof of (1) follows from the assumption.
# Proof of (2) follows from the assumption and (1): first, maximally match every bit of v1 with a bit
# from v2; what's left is necessarily D mismatching 0s of v1 and D 1s of v2.
# Thus, if the number of mismatches between v1 and v2 is exactly equal to D, then that implies that v1 and
# v2 are potential BL and TR corners of a subcube. Any other number of mismatches means otherwise.
# Examples: D = 2; v1 = 11000; minimal mismatch: v2 = 11011; maximal mismatch: v2 = 10111.

DetectCubes <- function(SubLUT, K) {
  
  # For a given fixed variable NumOnes (numfxdones), identify all possible BL and TR, match them and detect subcubes.
  
  MatchCornersDetectCubes <- function(SubLUT, numfxdones, LUTNumOnes, 
                                      NumCubeVars, K, CvrDims) {
    
    # Identify BL and TR based on their NumOnes, given numfxdones and subcube Dim (NumCubeVars)
    
    LUTBLIdx = which(LUTNumOnes == numfxdones)
    LUTTRIdx = which(LUTNumOnes == (numfxdones + NumCubeVars))
    
    PsblCubeBottLeftCorners = matrix(SubLUT[LUTBLIdx, ], ncol = K)
    PsblCubeTopRightCorners = matrix(SubLUT[LUTTRIdx, ], ncol = K)
    
    NumBL = nrow(PsblCubeBottLeftCorners)
    NumTR = nrow(PsblCubeTopRightCorners)
    
    bl = 1
    
    AllCubes = list()
    AllCvrdStsIndices = list()
    
    ReqdNumOtherCorners = (2 ^ NumCubeVars) - 2
    
    CubeListIdx = 1
    
    # A single BL could be associated with multiple subcubes; so is a single TR.
    
    while (bl <= NumBL & NumTR > 0) {
      
      BLInput = PsblCubeBottLeftCorners[bl, ]
      
      CubeExists = T
      
      tr = 1
      
      while (tr > 0 & tr <= NumTR) {
        
        TRInput = PsblCubeTopRightCorners[tr, ]
        
        # Compute mismatch between BL and TR
        
        MtchIdx = !xor(BLInput, TRInput)  # xor(0, 1) -> T; xor(0, 0) -> F
        
        cube = c()
        
        if (sum(MtchIdx) == (K - NumCubeVars)) {  # MATCH: the only elements that don't match are 0s of BL and 1s of TR (REMARK 1)
          
          cube = c(cube, LUTBLIdx[bl], LUTTRIdx[tr])
          
          if (NumCubeVars < K) {
            
            FxdVars = which(MtchIdx)
            FxdVals = BLInput[MtchIdx]
            
          }
          
          if (NumCubeVars > 1) {
            
            OtherNumOnesLB = numfxdones + 1
            OtherNumOnesUB = numfxdones + NumCubeVars - 1
            
            LUTOtherCornerIdx = which(LUTNumOnes >= OtherNumOnesLB & 
                                        LUTNumOnes <= OtherNumOnesUB)
            
            PsblCubeOtherCorners = matrix(SubLUT[LUTOtherCornerIdx, ], ncol = K)
            
            if (NumCubeVars < K) {
              
              FxdColMtch = sapply(1 : length(FxdVars), function(i) 
                !xor(PsblCubeOtherCorners[, FxdVars[i]], FxdVals[i]))
              
              
              NetFxdColMtch = rep(T, nrow(FxdColMtch))
              
              for (fxcol in 1 : ncol(FxdColMtch)) NetFxdColMtch = NetFxdColMtch & FxdColMtch[, fxcol]
              
            } else NetFxdColMtch = rep(T, nrow(PsblCubeOtherCorners))
            
            OtherLUTNumOnesIdx = LUTOtherCornerIdx[NetFxdColMtch]
            
            if (length(OtherLUTNumOnesIdx) == ReqdNumOtherCorners) 
              cube = c(cube, LUTOtherCornerIdx[NetFxdColMtch])
            else CubeExists = F
            
          } # else a 1D cube exists (that contains only the BL and TR)
          
          if (CubeExists) {
            
            NewCvrs = (CvrDims[cube] < NumCubeVars)
            
            if (sum(NewCvrs) > 0) {
              
              CvrDims[cube][NewCvrs] = NumCubeVars
              
              Schema = rep(2, K)
              
              if (NumCubeVars < K) Schema[FxdVars] = FxdVals
              
              AllCubes[[CubeListIdx]] = Schema
              AllCvrdStsIndices[[CubeListIdx]] = cube
              
              CubeListIdx = CubeListIdx + 1
              
            }
            
            # Prune search (remove BL and TR that cannot be candidates anymore)
            
            if (sum(TRInput) == NumCubeVars) {  # this TR corner won't match with any other BL corner
              
              LUTTRIdx = LUTTRIdx[-tr]
              
              if (NumTR <= 2) {  # Removing one row from a 2-row matrix will return a vector, not a matrix
                
                PsblCubeTopRightCorners = matrix(PsblCubeTopRightCorners[-tr, ], ncol = K)
                
              } else {
                
                PsblCubeTopRightCorners = PsblCubeTopRightCorners[-tr, ]
                
              }
              
              NumTR = NumTR - 1
              
            } else if (sum(BLInput) == (K - NumCubeVars)) {  # this BL corner won't match with any other TR corner
              
              break  # no need to match other TR corners
              
            } else tr = tr + 1
            
          } else {  # Cube not found
            
            CubeExists = T
            
            tr = tr + 1
            
          }
          
        } else tr = tr + 1  # NO MATCH: move to the next tr if the current topRight does not match the current bottomLeft
        
      }  # end of search for suitable topRight corner for the current bottomLeft corner
      
      bl = bl + 1
      
    }  # end of scanning through all bottomLeft corners
    
    return(list(CvrDims = CvrDims, AllCubes = AllCubes, AllCvrdStsIndices = AllCvrdStsIndices))
    
  }  # end of function definition
  
  LUTNumOnes = apply(SubLUT, 1, sum)
  
  # Compute distribution of the number of ones in the LUT
  
  LUTNumOnesCnts = hist(LUTNumOnes, breaks = 0 : (K+1), right = F, plot = F)$counts
  
  # Compute maximum possible subcube dimension
  
  MaxCubeSize = nrow(SubLUT)
  
  if (MaxCubeSize == 0) return(list(Schemata = NULL, CvrdStsIndices = NULL, CvrDims = NULL))
  
  MaxCubeDim = floor(log(MaxCubeSize, 2))
  
  CvrDims = rep(0, nrow(SubLUT))
  
  AllSchemata = list()
  AllCvrdStsIdx = list()
  
  # Consider every possible D from max down to 1.
  
  if (MaxCubeDim >= 1) {
    
    PsblCubeDims = MaxCubeDim : 1
    
    for (cbdim in PsblCubeDims) {
      
      NumCubeVars = cbdim
      
      # For every subcube dimension D, there are K-D 'fixed' variables with NumOnes ranging in [0, K-D].
      
      NumFixedVars = K - NumCubeVars
      
      PsblNumOnesCnts = sapply(0 : NumCubeVars, function(x) choose(NumCubeVars, x))
      
      # For every possible fixed variable NumOnes, there is a potential subcube.
      
      for (numfxdones in 0 : NumFixedVars) {
        
        ShftNumOnesCnts = LUTNumOnesCnts[(0 : NumCubeVars) + 1 + numfxdones]
        
        # For a given fixed variable NumOnes, there is a binomial distribution that every potential subcube must match.
        
        if (all(ShftNumOnesCnts >= PsblNumOnesCnts)) {
          
          PsblCubeIdx = (LUTNumOnes >= numfxdones & 
                           LUTNumOnes <= (numfxdones + NumCubeVars))
          
          if (min(CvrDims[PsblCubeIdx]) < NumCubeVars) {
            
            Results = MatchCornersDetectCubes(SubLUT, numfxdones, LUTNumOnes, 
                                              NumCubeVars, K, CvrDims)
            
            CvrDims = Results$CvrDims
            
            AllSchemata = append(AllSchemata, Results$AllCubes)
            AllCvrdStsIdx = append(AllCvrdStsIdx, Results$AllCvrdStsIndices)
            
          }
          
        }
        
      }
      
    }
    
  }
  
  StandAlones = (CvrDims == 0)
  
  if (sum(StandAlones) > 0) {
    
    sch = unlist(apply(matrix(SubLUT[StandAlones, ], ncol = K), 1, list), recursive = F)
    sts = unlist(lapply(which(StandAlones), list), recursive = F)
    
    AllSchemata = append(AllSchemata, sch)
    AllCvrdStsIdx = append(AllCvrdStsIdx, sts)
    
  }
  
  return(list(Schemata = AllSchemata, CvrdStsIndices = AllCvrdStsIdx, CvrDims = CvrDims))
  
}

















