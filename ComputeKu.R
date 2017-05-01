# This program computes the unified canalization (k_u) of a Boolean function.

# Appendix A of the documentation describes the main concepts in detail with the help of an example.

# USAGE EXAMPLE:

# Note: Make sure you also have 'ComputeDetectCubes.R' in the same folder as this file.

# k = 3
# func = c(1,1,1,1,0,1,1,1)
# ComputeKu(func, k, ComputeKeff = T, ReturnSymmSchs = T)

# The following output is generated:
# $EssentialSchemata
# [,1] [,2] [,3]
# [1,]    0    5    5
# [2,]    2    4    5
# [3,]    1    0    0
# 
# $SchemaOutputs
# [1] 1 1 0
# 
# $Keff
# [1] 1.25
# 
# $Ku
# [1] 0.8112781

# Optional program parameters: 
# 'ReturnSymmSchs = T' will return the essential two-symbol schemata including ordinary schemata.
# 'ComputeKeff = T' will return the Keff value as well.

# Description of output:
# If 'ReturnSymmSchs = T' then the essential two-symbol schemata and the corresponding output values are returned.
# The schema symbols are encoded according to the following mapping: '0' -> 0, '1' -> 1, '#' -> 2, '0'' -> 3, '1'' -> 4, '#'' -> 5
# Non-permuting '#' are also represented with 5 when they are the only kind of '#' in the schema.
# So, (02345) mean (0#0'1'#'); (015) means (01#).

#-------------------------------------------------------------------------------------------------------------------------------------

# BASIC IDEA: 
# Identify "two-symbol schemata" on the basis of the number of 1s in the input values. We may think of this as a
# generalization of identifying input schemata where the number of 1s must cover the entire range of possible values
# in the interval [0, D] where D is the dimension of the cube. Two-symbol schemata, on the other hand, may comprise input values
# with number of 1s occupying partial sub-intervals of [0, D].

# TWO-SYMBOL SCHEMATA: 
# Following is a method for identifying two-symbol schemata. For each cube and for a given output value (0 or 1), extract the 
# set of input vectors corresponding to that cube and that go to the given output. Count the number of 1s in each input value
# of that set. For a set of input values with a given number of 1s to form a two-symbol schema, its cardinality must be equal
# to C(D, NumOnes) where D is the dimension of the cube. For example, |(0'1'1')| = C(3, 2) = 3. We only identify two-symbol
# schemata that consist of at least one '#'. For example, we consider (0'1'#'), which is a combination of (0'1'0') and (0'1'1'), 
# but we don't consider (0'1'1') or (0'1'0') on their own, if they are stand-alone schemata. Further, a cube may contain 
# more than one two-symbol schema. For example, (0'0'0'1'1'#') and (1'1'1'1'1'#') are two non-overlapping schemata in a 
# 6-dimensional cube.

# COMPOSITE SCHEMATA:
# A composite schema is a two-symbol schema with at least one fixed (non-permuting) '#' and at least one permuting permuting 
# '#'; it is named so because it contains the characteristic features of both an ordinary schema and a two-symbol schema.
# We identify 'composite' schemata by combining two or more two-symbol schemata, as identified above, on the basis that
# one or more non-permuting variables in the combined schema may be redundant. For example, we can combine (011'#') and
# (001'#') to form (0#1'#'). This is the only kind of composite schemata that we identify. There are others, for example, (1'#')(0'#') 
# and even more complex ones like (1'#')'(0'#')' that may exist but are not identified here. 
# Instead of combing through the schemata set and combining them, we enumerate the 'signatures' and then combine schemata that match 
# a given signature. For example, given NumPermutVars = 2, the possible signatures are (#'#'), (0'#'), and (1'#').
# A given schema does not have to match a signature exactly -- partial matches also work. For example, the schema (1'#'#') 
# partially matches the signature (0'1'#') since the former fully covers the latter. In fact, partial matches are necessary 
# to identify composite schemata that bridge portions of parallel schemata. 
# We identify composite from sets of 'opposite-facing' or 'parallel' cubes. For example, (00'#') and (10'#') are schemata 
# from a pair of parallel cubes, identified by i1=0 and i1=1 respectively, that combine to give the composite schema (#0'#').
# Schemata that are not parallel to each other can also combine to give composite schemata. A simple example is: 
# (#1#) and (##1) combine to give (#1'#'). However, it is straightforward to show that any composite schema can be obtained
# by combining two or more parallel two-symbol schemata. We also identify composite schemata that bridge *portions* of two 
# parallel schemata (by partially matching signatures). For example, (00'1'#') and (1###) are parallel schemata but they combine 
# to give (#0'1'#'). Likewise, (00'#'#') and (11'#'#') can combine to give (#0'1'#').

# OUTLINE OF THE ALGORITHM:
# Repeat steps 1 to 4 below for every possible dimension $D$ in decreasing order from $k$ to 1, and for every possible set 
# of parallel cubes of a given dimension. 
# Step 1: Consider a single set of parallel subcubes of a given dimension.
# Step 2: Identify all ordinary two-symbol schemata in each subcube in the set obtained in the previous step and for each 
# output value.
# Step 3: Identify all composite schemata from the set of ordinary two-symbol schemata obtained in the previous step and for 
# each output value.
# Step 4: Record the dimension of the two-symbol schema against every input vector it covers if and only if the current largest 
# covering two-symbol schema's dimension is smaller. 
# Step 5: Compute $k_r^*$ by averaging over the covering dimensions of all $2^k$ input vectors. Compute $k_e^* = k - k_r^*$.

#-------------------------------------------------------------------------------------------------------------------------------------

# Identifying SymmSch without any constraints on the number of permuting wildcards:
# As mentioned above, this program identifies only those two-symbol schemata that contain at least one permuting `#'.
# The places in the code where that constraint is enforced are marked as 'MinNumWC(A)', 'MinNumWC(B)' and 'MinNumWC(C)'.
# Modifying the conditions specified in all of those places will enable removing the constraint.

# Efficiency considerations: Unlike the Keff calculation method where if a higher-dimensional schema covering a certain region 
# of the hypercube exists that region wouldn't have to be combed further, the Ku calculation method here would have to
# continue combing. For example, at D = 4, we might have found (1'1'#'#') whose cardinality is 11. However, it's possible that 
# a schema like (##1'#') exists whose cardinality is higher: 12. Here, even though the last two variables is fully covered by a higher-
# dimensional schema (former), those two variables would still have to be combed to identify composite schemata like the (##1'#'). We
# perhaps need an analytical way to prune some of these searches. We already have such a 'Pruning step' (see note in the code below).
# We need more like these to improve efficiency.

# Note of reassurance: Before coding this version, we had two versions that identified two-symbol schemata in an inefficient 
# way -- by enumerating all possible two-symbol schemata for a given dimension (number of permuting variables), and then matching 
# them with the given cubes. They are called 'ComputeKeffGeneral(CvrSzOrder).R' and 'ComputeKeffGeneral(TreeLvlOrder).R'.
# These turned out to be useful because the output of this version was validated against the output of the above alternatives 
# which are sure-shot methods.

# Misc note: In an earlier version (ComputeKeffGeneral(Rotate).R), we rotated every cube once, and considered only those input
# values whose outputs remained invariant after the rotation. Turns out that that step is not necessary. Accordingly, this version
# is rid of the rotation step.

#-------------------------------------------------------------------------------------------------------------------------------------

source('ComputeDetectCubes.R')

require(parallel)
require(R.utils)  # for 'intToBin'

ComputeKu <- function(Func, K, ComputeKeff = F, ReturnSymmSchs = F, NumCores = 1) {
  
  ## Miscellaneous functions 
  
  # Convert a decimal to a binary of a specified length. Max allowed value for len = 32.
  
  dec2bin <- function(x, len) {
    
    sfx = as.integer(strsplit(intToBin(x),'')[[1]])
    
    c(rep(0, len - length(sfx)), sfx)  
    
  }
  
  # Compute the number of input states that a given symemtric group covers.
  # K_p = Number of permuting dimensions
  # W = Number of permuting wildcards
  # O = Number of permuting ones
  # K_f = Number of frozen dimensions with wildcards
  
  ComputeNumCvrdInputStates <- function(K_p, W, O, K_f) {
    
    O_LB = O
    O_UB = O + W
    
    (2 ^ K_f) * sum(sapply(O_LB: O_UB, function(ones) choose(K_p, ones)))
    
  }
  
  ## End of miscellaneous functions 
  
	# Identifies the two-symbol schemata present in a cube of a given Dim. For Dim = 2, for example, 
	# the possible two-symbol schemata are (#'#'), (0'#'), (1'#'). The following are examples that are NOT
	# two-symbol schemata and therefore won't be identified by this method: (0'0'), (1'0'), (1#). 
	# Note that even though a schema like (1#) won't be identified for Dim = 2, it will be identified
	# when Dim = 1, since it's equivalent to (1#').
	
	ScanCube <- function(LUT, Func, fxdvars, fxdvals, Dim, K) {

		if (!is.null(fxdvars)) {
			
			FixedLUT = as.matrix(LUT[, fxdvars], ncol = length(fxdvars))
			
			FixedCubeIdx = apply(FixedLUT, 1, function(r) {x = (r == fxdvals); all(x == T)})
		
		} else {
	
			FixedCubeIdx = rep(T, nrow(LUT))
			
		}
		
		permutvars = setdiff(1 : K, fxdvars)
		
		PermutLUT = LUT[, permutvars]
		if (length(permutvars) == 1) PermutLUT = matrix(PermutLUT, ncol = 1)
		
		PermutLUTNumOnes = apply(PermutLUT, 1, sum)
		
		DesiredNumOnes = sapply(0 : Dim, function(x) choose(Dim, x))  # the binomial distribution
		
		Schemata = list()
		CvrdStates = list()
		OutputVals = c()
		
		ScanResIdx = 1
		
		for (OutputVal in 0 : 1) {
			
			OutIdx = ((Func == OutputVal) & FixedCubeIdx)
		
			NumStates = sum(OutIdx)
			
			# Either the entire cube of Dim is a SymmSch with only permuting '#' or there may be one or more SymmSch
			
			if (NumStates == 2^Dim) {  # the entire cube is a single two-symbol schema with only permuting '#'
				
				Schema = rep(-1, K)
				
				if (!is.null(fxdvars)) Schema[fxdvars] = fxdvals
				
				Schema[permutvars] = rep(5, length(permutvars))
				
				CvrdStates = list(which(OutIdx))
				
				Schemata[[ScanResIdx]] = Schema
						
				CvrdStates[[ScanResIdx]] = which(OutIdx)
				
				OutputVals  = c(OutputVals, OutputVal)
								
			} else {  # may contain (multiple) two-symbol schemata
				
				if (NumStates < (Dim + 1)) next  # MinNumWC(A): if true, a SymmSch with at least one permuting `#' can't exist
				
				NumOnes = PermutLUTNumOnes[OutIdx]
				
				# The following extracts from the cube all possible SymmSch with at least one permuting '#'.
				# Basic idea: As mentioned earlier, every SymmSch in a cube of dimension D covers a set of LUT entries whose
				# NumOnes range over some sub-interval of [0,D] and whose distribution in that range is binomial ('DesiredNumOnes').
				# SymmSch with at least one permuting '#' has an associated NumOnes range whose width is >= 1. For example,
				# the NumOnes range of (1'#') is [1,2] whose width is 1, whereas that of (0'1') is 0 since its range is [1,1].
				# A cube of dimension D may contain multiple SymmSch with non-overlapping NumOnes ranges that jointly don't cover
				# the full range [0,D]. For example, in a cube with D=4, the following two SymmSch could both exist: (0'0'0'#') and
				# (1'1'1'#') with NumOnes ranges [0,1] and [3,4] respectively. In a cube with D=3, on the other hand, the SymmSch 
				# (0'0'#') and (1'1'#') can't both exist since together they cover the entire cube; their NumOnes ranges [0,1] and
				# [1,2] together cover the range [0,2].
				
				# Extract LUT entries whose NumOnes distributions are binomial
				
				IsSymmSch = sapply(0 : Dim, function(ones) sum(NumOnes == ones) == DesiredNumOnes[ones + 1])
				IsSymmSch = c(IsSymmSch, F)  # add a F at the end to pick up any trailing SymmSchs
				
				GapInds = which(IsSymmSch == F)  # specify gaps between multiple SymmSch in the same cube
				
				# Extract all SymmSch in the cube
				
				st = -1
				for (i in 1 : length(GapInds)) {
					
					nd = GapInds[i] - 1
					
					NumPermutW = nd - st - 2
					
					if (NumPermutW > 0) {  # MinNumWC(B): guarantees SymmSch with at least one permuting `#'
						
						NumPermutO = st + 1
						NumPermutZ = Dim - (NumPermutO + NumPermutW)
						
						SymmSch = c(rep(3, NumPermutZ), rep(4, NumPermutO), rep(5, NumPermutW))
						
						Schema = rep(-1, K)
				
						if (!is.null(fxdvars)) Schema[fxdvars] = fxdvals
						
						Schema[permutvars] = SymmSch
						
						Schemata[[ScanResIdx]] = Schema
						
						CvrdStates[[ScanResIdx]] = which(OutIdx & 
																			(PermutLUTNumOnes >= st + 1 & PermutLUTNumOnes <= nd - 1))
						
						OutputVals  = c(OutputVals, OutputVal)
						
						ScanResIdx = ScanResIdx + 1
						
					}
					
					st = nd
					
				}  # end of 'GapInds' loop
				
			}  # end of 'SymmSch' if-else
		
		}  # end of 'OutputVal' loop
		
		ScanResults = list(Schemata = Schemata, CvrdStates = CvrdStates, OutputVals = OutputVals)
		
		return(ScanResults)
			 		
	} # end of function ScanCube
	
	# Scan opposite-facing or 'parallel' cubes. Ex: {0##, 1##} are such a pair of parallel cubes; 
	# {00#, 01#, 10#, 11#} and {0#0#, 0#1#, 1#0#, 1#1#} are sets of parallel cubes (visualize them).
	# The only motivation for scanning parallel cubes: to identify composite schemata that can only form
	# through combinations of subsets of schemata in parallel cubes.
	
	ScanParallelCubes <- function(idx, FixedVars, FixedVals, LUT, NodeKred, NodeKredGeneral, Func, Dim, K, 
																ComputeKeff, RetSS) {
		
	  # Step 1: Consider a single set of parallel subcubes of a given dimension
	  
		# First scan each cube in the set of parallel cubes
		
		if (! is.null(FixedVars)) {
		
			fxdvars = FixedVars[ , idx]
			NumFixedVars = length(fxdvars)
			NumFixedVals = ncol(FixedVals)
		
		} else {
			
			fxdvars = NULL
			NumFixedVars = 0
			NumFixedVals = 0
			
		}
							
		permutvars = setdiff(1 : K, fxdvars)
		
		NumPermutVars = length(permutvars)
		
		# Step 2: Identify all ordinary two-symbol schemata in each subcube in the set obtained in the previous step and for each 
		# output value.
		
		if (NumFixedVars > 0) {  # Dim < K

			AllSchemata = list()
			AllCvrdStates = list()
			AllOutputVals = c()
					
			for (nfi in 1 : NumFixedVals) {
				
				fxdvals = FixedVals[, nfi]
				
				ScanResults = ScanCube(LUT, Func, fxdvars, fxdvals, Dim, K)
				
				if (length(ScanResults$Schemata) > 0) {
					
					AllSchemata = append(AllSchemata, ScanResults$Schemata)
					AllCvrdStates = append(AllCvrdStates, ScanResults$CvrdStates)
					AllOutputVals = c(AllOutputVals, ScanResults$OutputVals)
				
				}
								
			}
		
		} else {  # Dim = K; the set of parallel cubes consists of only one cube
		  
			ScanResults = ScanCube(LUT, Func, NULL, NULL, Dim, K)
			
			if (length(ScanResults$Schemata) > 0) {
				
				AllSchemata = ScanResults$Schemata
				AllCvrdStates = ScanResults$CvrdStates
				AllOutputVals = ScanResults$OutputVals
			
			} else {

				if (RetSS) return(list(EssentialSchemata = NULL, NodeKred = NodeKred, NodeKredGeneral = NodeKredGeneral))
				else return(list(NodeKred = NodeKred, NodeKredGeneral = NodeKredGeneral))
				
			}
			
		}
		
		if (length(AllSchemata) == 0) {
			
			if (RetSS) return(list(EssentialSchemata = NULL, NodeKred = NodeKred, NodeKredGeneral = NodeKredGeneral))
			else return(list(NodeKred = NodeKred, NodeKredGeneral = NodeKredGeneral))
			
		}
		
		# Then, process the above-identified schemata to form composite schemata.
		# Also, record dimensions of covering schemata.
		
		AllSchemata = matrix(unlist(AllSchemata), ncol = K, byrow = T)
		
		if (RetSS) {
		  EssentialSchemata = list()
		  EssSchCvrdStates = list()
		  EssSchOutputs = list()
		  EssSchIdx = 1
		}
		
		# Step 3: Identify all composite schemata from the set of ordinary two-symbol schemata obtained in the previous step 
		# and for each output value.
		
		# For each OutputVal, process schemata set obtained above from the set of parallel cubes 
		
		for (OutputVal in 0 : 1) {  
			
			schidx = (AllOutputVals == OutputVal)
			
			if (sum(schidx) == 0) next
			
			Schemata = matrix(AllSchemata[schidx, ], ncol = K)
			CvrdStates = AllCvrdStates[schidx]
			
			NumSchemata = nrow(Schemata)
			
			SymmSchs = Schemata[ , permutvars]
			if (NumSchemata == 1 | NumPermutVars == 1) SymmSchs = matrix(SymmSchs, ncol = NumPermutVars)
			
			# Create a matrix contain only the FxdVals of the SymmSch; NumFixedVars = K - NumPermutVars (or D)
			
			if (NumFixedVars > 0) {
				
				SymmSchsFxdVals = Schemata[ , fxdvars]
				if (NumSchemata == 1| NumFixedVars == 1) SymmSchsFxdVals = 
																							matrix(SymmSchsFxdVals, ncol = NumFixedVars)
			
			}
			
			# Create an 'indicator' matrix that indicates the NumOnes associated with each SymmSch in the range [0,D]
			
			SchRowIndices = 1 : NumSchemata
			
			SymmSchIndMat = matrix(0, NumSchemata, NumPermutVars + 1)
			
			for (r in 1 : NumSchemata) {
				
				SS = SymmSchs[r, ]
				
				MinNumOnes = sum(SS == 4)
				MaxNumOnes = MinNumOnes + sum(SS == 5)
				
				SymmSchIndMat[r, MinNumOnes : MaxNumOnes + 1] = 1
				
			}
			
			# Instead of combing through the schemata set and combining them, we enumerate the 'signatures' and then
			# combine schemata that match a given signature. For example, given NumPermutVars = 2, the possible signatures
			# are (#'#'), (0'#'), and (1'#'). A given schema does not have to match a signature exactly -- partial matches also work.
			# For example, the schema (1'#'#') partially matches the signature (0'1'#') since the former covers by the latter.
			# In fact, partial matches are necessary to identify composite schemata that bridge portions of parallel schemata; 
			# example mentioned below.
			# Furthermore, whether a composite SymmSch is detected or not, below is where the dimension of the largest SymmSch
			# covering each input value (state) is recorded.
			
			ProcessedSymmSch = rep(F, NumSchemata)
			
			# "Processing" the schemata entails the following:
			# Identify composite schemata, which is possible only when NumPermutVars > 1 & NumPermutVars < K.
			# Record the dimensions of the largest covering SymmSch against each covered LUT entry.
			# Record the covering SymmSch if requested.
			
			if (NumPermutVars > 1 & NumPermutVars < K) {  # Composite schemata possible
				
				# Note (01/28/2016): We tried to optimize the below for-for loop by narrowing down the iteration ranges by making 
				# them 'data driven'. However, it seemed difficult to come up with those. We tested a few, but they failed.
			  
			  # Difficulties with pruning in general (11/10/2016): 
			  # The following example illustrates one of the obstacles to pruning. Suppose that K=3. In this case, a pruning
			  # procedure, upon discovering the SymmSch (0'#'), say, would stop there. This is because no other SymmSch (e.g.,1'#')
			  # can exist, lest the whole cube be covered in which case the valid SymmSch would be (#'#'). The basic idea there
			  # is that the first SymmSch covers a range [0,1] of NumOnes, whereas the latter covers [1,2]; together they cover [0,2].
			  # Now suppose K=4 and we have the following set of SymmSch: {(1#0'#', #11'#', 11#'#')}. The above pruning won't work for
			  # the last two permuting variables, since (1#0'#') and (#11'#') can co-exist without coalescing into one.
			  # Thus, a more complex pruning procedure would be needed when not all K variables are permuting.
			  # Even when all K variables permute, it might seem that a SymmSch with fewer permuting '#' would always be covered
			  # by one with more permuting '#'. An example that violates this seeming constraint is: (0'0'#'#') and (0'1'1'#').
			  # In this case, K=4 and the corresponding NumOnes ranges [0,2] and [2,3] (that is, the case with all 4 ones is not
			  # covered by either SymmSch). This means that for K=6, say, the following SymmSch are possible: (##0'0'#'#') and 
			  # (1#0'1'1'#') even though the '##' of the first SymmSch fully covers the `1#' of the other. 
			  
			  # The other important point is that even though subsequent SymmSch in the loop (sequence) below need not be fully
			  # covered by those found earlier in the sequence, it's guaranteed that later SymmSch cannnot fully cover earlier ones.
			  # E.g., (0'#'#') covers (0'1'#'), and (0'0'#'#') does not cover (0'1'1'#'), but in neither case the latter covers the former.
			  # The basic idea behind the above fact is that subsequent SymmSch necessarily cover smaller ranges of NumOnes. Therefore, 
			  # it's possible that a smaller NumOnes range is either fully or partially covered by a larger range, the former can 
			  # never fully cover the latter. Thus, (###1'1'#') does not cover (11#1'#'#') even though the cardinality of the former
			  # (32) is quite larger than that of the latter (14). 

				for (NumWC in NumPermutVars : 1) {  # MinNumWC(C): guarantees SymmSch with at least one permuting `#'
					
					for (MinNumOnes in 0 : (NumPermutVars - NumWC)) {

						MaxNumOnes = MinNumOnes + NumWC
						
						FiltSSIndMat = matrix(SymmSchIndMat[, MinNumOnes : MaxNumOnes + 1], ncol = NumWC + 1)
						
						# Here is where we effectively combine schemata like (00'1'#') and (1#'#'#') to get (#0'1'#'), or
						# (00'#'#') and (11'#'#') to get (#0'1'#').
						
						# Any SymmSch in a cube of dimension D is fully identified by the associated min NumOnes and the max NumOnes,
						# indicated here by MinNumOnes and MaxNumOnes respectively (signature). We run these signatures through the 
						# indicator matrix and pick up any matches.
						
						FitsSSSign = T
						for (SScol in 1 : ncol(FiltSSIndMat)) FitsSSSign = FitsSSSign  & FiltSSIndMat[, SScol]
						
						# The number of matches could range from 0 to NumSchemata.
						
						NumMtchngSS = sum(FitsSSSign)
						
						if (NumMtchngSS == 0) next  # no matching SymmSchs
												
						if (NumMtchngSS > 1 & NumWC < NumPermutVars) {  # form Composite SymmSchs
							
							# Pruning steps: no need to identify composite SymmSchs that don't contribute to K_e_g
							
							# CAUTION: The following line of code may seem like a valid pruning step, but it's not.
							# An example of where it will fail: Let k=4 and SymmSch set going to 1 = {00'#'#',11'#'#',#0'1'#'}.
							# In this case, after detecting 00'#'#' and 11'#'#', #0'1'#' would be skipped if the below condition 
							# is included. This is because 00'#'#' and 11'#'#' are both *separately* processed and #0'1'#' is a
							# bridge that runs between the two.
							# if (all(ProcessedSymmSch[FitsSSSign])) next  # DO NOT INCLUDE
							
							# Step 1: no need to identify SymmSch if the max possible schema dim is less than that already exists.
							# Note that if Pruning step 1 checks out, Pruning step 2 checks out as well; by separating them we reduce
							# computational expense.
							
							MaxPsblSchDim = floor(log(NumMtchngSS, 2)) +   # contribution from fixed '#'
							  log(ComputeNumCvrdInputStates(NumPermutVars, NumWC, MinNumOnes, 0), 2)   # contribution from permuting symbols
							
							MinKredGeneralSoFar = min(NodeKredGeneral[unlist(CvrdStates[FitsSSSign])])

							if ((MinKredGeneralSoFar >= MaxPsblSchDim)) next
							
							# If pruning conditions are not satisfied, then it is possible that a valid composite schema exists.
							# Extract the LUT columns corresponding to the FixedVars columns and rows that match the current SymmSch signature.
							# Then pass it to your favorite prime implicant detection function, which here is called 'DetectCubes'.
							# 'DetectCubes' is separately defined in 'ComputeDetectCubes.R'. Note that the prime implicants here constitute
							# the non-permuting parts of the composite schema. For example, the first two symbols in 1#0'#' represent a prime
							# implicant and would be identified by combining 100'#' and 110'#' (the first two 10 and 11 combine to form 1#).

							FixedVarSubLUT = matrix(SymmSchsFxdVals[FitsSSSign, ], ncol = NumFixedVars)
							
							Cubes = DetectCubes(FixedVarSubLUT, NumFixedVars)
							
							# The number of 'Cubes' or prime implicants is at least one and at most equal to the number of SymmSch that
							# match the current signature or nrow(FixedVarSubLUT).
							
							for (cbidx in 1 : length(Cubes$Schemata)) {  # this is assuredly non-empty
							
								SymmSchSubIdx = SchRowIndices[FitsSSSign][Cubes$CvrdStsIndices[[cbidx]]]  # e.g., (1:5)[c(T,T,F,T,F)][c(2,3)] gives c(2,4)
								
								if (length(SymmSchSubIdx) > 1) {  # Composite SymmSch
								  
								  # Since composite schemata are a bridge between one or more SymmSch, the set of covered LUT entries corresponding
								  # to the newly formed composite schema is not simply the union of the covered LUT entries of the SymmSch from which
								  # it is obtained. Instead, it is the subset whose NumOnes match the current SymmSch signature. For example,
								  # 0## and 11'#' combine to form #1'#' whose covered set of LUT entries is a subset of the union of the two.
								
									FilterLUTIdx = unlist(CvrdStates[SymmSchSubIdx])  # Covered states by parallel cubes are guaranteed to be non-overlapping; no need to 'unique'
									
									FilterLUT = LUT[FilterLUTIdx, permutvars]
									
									NumOnes = apply(FilterLUT, 1, sum)
									
									# Filter LUT entries that match current signature
									
									cvrdsts = FilterLUTIdx[NumOnes >= MinNumOnes & NumOnes <= MaxNumOnes]  
																		
									ProcessedSymmSch[SymmSchSubIdx] = T
										
								} else {  # Regular SymmSch; e.g., {110'#',000'#',010'#'} -> {110'#',000'#',010'#'}; here 110'#' is a lone regular SymmSch 
									
									if (ProcessedSymmSch[SymmSchSubIdx]) next   # subsequent matches can only cover LUT entries that are a subset of those from the first match
									
									else {
										
										cvrdsts = CvrdStates[[SymmSchSubIdx]]
										
										ProcessedSymmSch[SymmSchSubIdx] = T
									
									}
									
								}
								
								# Step 4: Record the dimension of the two-symbol schema against every input vector it covers if and only if 
								# the current largest covering two-symbol schema's dimension is smaller. 
								
								# Record dim of this SymmSch against covered LUT entries if their covering dimensions are smaller
															
								SchSize = length(cvrdsts)
								
								SchDim = log(SchSize, 2)
								
								NewCoverIdx = which(NodeKredGeneral[cvrdsts] < SchDim)
							
								if (length(NewCoverIdx) > 0) {
									
									if (RetSS) {
										
										sch = rep(-1, K)
										
										MinNumZeros = NumPermutVars - MaxNumOnes
										sch[permutvars] = c(rep(3, MinNumZeros), rep(4, MinNumOnes), rep(5, NumWC))
										
										FxdVarSchema = Cubes$Schemata[[cbidx]]
										sch[fxdvars] = FxdVarSchema
										
										EssentialSchemata[[EssSchIdx]] = sch
										EssSchCvrdStates[[EssSchIdx]] = cvrdsts
										EssSchOutputs[[EssSchIdx]] = OutputVal
										
										EssSchIdx = EssSchIdx + 1
										
									}
									
									NodeKredGeneral[cvrdsts][NewCoverIdx] = SchDim
								
								}
								
								if (ComputeKeff) {
									
									NumFixedWC = sum(Cubes$Schemata[[cbidx]] == 2)
									
									SchDimKeff = NumFixedWC + NumWC
									
									NewCoverIdxKeff = which(NodeKred[cvrdsts] < SchDimKeff)
									
									if (length(NewCoverIdxKeff) > 0) NodeKred[cvrdsts][NewCoverIdxKeff] = SchDimKeff
									
								}
								
							}  # end of (for-loop) processing cubes returned by DetectCubes
						
						} else {  
						  
						  # One or more of the following must be true:
						  # a lone SymmSch (NumMtchngSS = 1, NumWC < NumPermutVars) 
						  # a set of parallel ordinary schemata (NumMtchngSS > 1, NumWC = NumPermutVars) 
						  # a lone ordinary schema (NumMtchngSS = 1, NumWC = NumPermutVars)
							
							for (r in which(FitsSSSign)) {
							  
							  # Record dim of this SymmSch against covered LUT entries if their covering dimensions are smaller
								
								cvrdsts = CvrdStates[[r]]
								
								SchSize = length(cvrdsts)
								
								SchDim = log(SchSize, 2)								
						
								NewCoverIdx = which(NodeKredGeneral[cvrdsts] < SchDim)
								
								if (length(NewCoverIdx) > 0) {
									
									if (RetSS) {
										
										sch = rep(-1, K)
										
										MinNumZeros = NumPermutVars - MaxNumOnes
										sch[permutvars] = c(rep(3, MinNumZeros), rep(4, MinNumOnes), rep(5, NumWC))
										
										FxdVarSchema = SymmSchsFxdVals[r, ]
										sch[fxdvars] = FxdVarSchema
										
										EssentialSchemata[[EssSchIdx]] = sch
										EssSchCvrdStates[[EssSchIdx]] = cvrdsts
										EssSchOutputs[[EssSchIdx]] = OutputVal
										
									  EssSchIdx = EssSchIdx + 1
										
									}
								
								}
								
								# Step 4: Record the dimension of the two-symbol schema against every input vector it covers if and only if 
								# the current largest covering two-symbol schema's dimension is smaller.
									
								if (ComputeKeff) {
								
									SchDimKeff = NumWC
									
									NewCoverIdxKeff = which(NodeKred[cvrdsts] < SchDimKeff)
									
									if (length(NewCoverIdxKeff) > 0) NodeKred[cvrdsts][NewCoverIdxKeff] = SchDimKeff
									
								}
								
								NodeKredGeneral[cvrdsts][NewCoverIdx] = SchDim
																
								ProcessedSymmSch[r] = T
							
							}  # end of lone SymmSchs or ordinary Schs loop
							
						}  # end of if-else condn of forming Composite SymmSchs or processing lone SymmSchs/ordinary Schs
						
					}  # end of MinNumOnes loop
					
				}  # end of NumWC loop
			
			} else {  # NumPermutVars = 1 (only ordinary Schs possible here) or NumPermutVars = K; in either case, composite schemata not possible
				
				for (r in 1 : NumSchemata) {
				  
				  # Record dim of this SymmSch against covered LUT entries if their covering dimensions are smaller
								
					cvrdsts = CvrdStates[[r]]
					
					SchSize = length(cvrdsts)
					
					SchDim = log(SchSize, 2)								
			
					NewCoverIdx = which(NodeKredGeneral[cvrdsts] < SchDim)
					
					# Step 4: Record the dimension of the two-symbol schema against every input vector it covers if and only if 
					# the current largest covering two-symbol schema's dimension is smaller.
										
					if (length(NewCoverIdx) > 0) {
						
						if (RetSS) {
							
							sch = Schemata[r, ]
							
							EssentialSchemata[[EssSchIdx]] = sch
							EssSchCvrdStates[[EssSchIdx]] = cvrdsts
							EssSchOutputs[[EssSchIdx]] = OutputVal
							
							EssSchIdx = EssSchIdx + 1
							
						}
						
						NodeKredGeneral[cvrdsts][NewCoverIdx] = SchDim
					
					}
					
					if (ComputeKeff) {
						
						sch = Schemata[r, ]
									
						SchDimKeff = sum(sch == 5)
						
						NewCoverIdxKeff = which(NodeKred[cvrdsts] < SchDimKeff)
						
						if (length(NewCoverIdxKeff) > 0) NodeKred[cvrdsts][NewCoverIdxKeff] = SchDimKeff
						
					}
					
				}  # end of Processing Dim = K or Dim = 1 SymmSchs
				
			}
			
		}  # end of 'OutputVal' loop
		
		# remove variables
		rm(AllSchemata)
		rm(AllCvrdStates)
		rm(AllOutputVals)
		
		if (RetSS){
			
			if (length(EssentialSchemata) > 0) {
			  
			  EssSchInfo = list(EssentialSchemata = EssentialSchemata, EssSchCvrdStates = EssSchCvrdStates, EssSchOutputs = EssSchOutputs)
			
			} else EssSchInfo = list(EssentialSchemata = NULL, EssSchCvrdStates = NULL, EssSchOutputs = NULL)
						
			return(list(EssSchInfo = EssSchInfo, NodeKred = NodeKred, NodeKredGeneral = NodeKredGeneral))
		
		} else return(list(NodeKred = NodeKred, NodeKredGeneral = NodeKredGeneral))
		
	}  # end of function ScanCubeSet
	
	
	### MAIN PROCESS ###
	
	# LUT = t(sapply(0: ((2^K) - 1), function(i) rev(dec2bin(i, K))))  # using 'dec2bin' defined in package 'BoolNet'
	LUT = t(sapply(0: ((2^K) - 1), function(i) dec2bin(i, K)))
	
	if (K == 1) LUT = matrix(LUT, ncol = K)
	
	NodeKred = rep(0, 2^K)
	
	NodeKredGeneral = rep(0, 2^K)
	
	RetSS = ReturnSymmSchs
	
	if (RetSS) {
	  
	  EssentialSchemata = list()
	  EssSchCvrdStates = list()
	  EssSchOutputs = list()
	  
	}
	
	for (Dim in K : 1) {  # Dim = 1 must be included because ordinary schemata with no permuting symbols are SymmSch as well
		
		NumFixedVars = K - Dim
		
		if (NumFixedVars > 0) {	
	
			FixedVars = combn(K, NumFixedVars)
			if (NumFixedVars == 1) FixedVars = matrix(FixedVars, nrow = 1)
			
			# FixedVals = sapply(0 : (2^NumFixedVars - 1), function(x) rev(dec2bin(x, NumFixedVars)))  # using 'dec2bin' defined in package 'BoolNet'
			FixedVals = sapply(0 : (2^NumFixedVars - 1), function(x) dec2bin(x, NumFixedVars))
			
			if (NumFixedVars == 1) FixedVals = matrix(FixedVals, nrow = 1)
			
			NumFixedVarCombos = ncol(FixedVars)
			
			ScanResults = mclapply(1 : NumFixedVarCombos, ScanParallelCubes, FixedVars, FixedVals, LUT, 
																NodeKred, NodeKredGeneral, Func, Dim, K, ComputeKeff, RetSS, 
																mc.cores = min(NumFixedVarCombos, NumCores))
			
			if (ComputeKeff) {
				
				NodeKred = c()
				
				for (scn in 1 : length(ScanResults)) NodeKred = c(NodeKred, ScanResults[[scn]]$NodeKred)
				
				NodeKred = matrix(NodeKred, ncol = 2^K, byrow = T)
				
				NodeKred = apply(NodeKred, 2, max)
			
			}
			
			NodeKredGeneral = c()
			
			for (scn in 1 : length(ScanResults)) NodeKredGeneral = c(NodeKredGeneral, ScanResults[[scn]]$NodeKredGeneral)
			
			NodeKredGeneral = matrix(NodeKredGeneral, ncol = 2^K, byrow = T)
			
			NodeKredGeneral = apply(NodeKredGeneral, 2, max)
			
			if (RetSS) {
			  
				for (scn in 1 : length(ScanResults)) {
				  
				  ResEssSch = ScanResults[[scn]]$EssSchInfo$EssentialSchemata
				  
				  if (!is.null(ResEssSch)) {
				    
				    EssentialSchemata = c(EssentialSchemata, ResEssSch)
				    EssSchCvrdStates = c(EssSchCvrdStates, ScanResults[[scn]]$EssSchInfo$EssSchCvrdStates)
				    EssSchOutputs = c(EssSchOutputs, ScanResults[[scn]]$EssSchInfo$EssSchOutputs)
				    
				  }
				  
				}
			  
			}
			
		} else {  # NumFixedVars = 0, that is, Dim = K (there is only one cube at this Dim)
			
			ScanResults = ScanParallelCubes(-1, NULL, NULL, LUT, NodeKred, NodeKredGeneral, Func, Dim, K, ComputeKeff, RetSS)

			if (ComputeKeff) NodeKred = ScanResults$NodeKred
			
			NodeKredGeneral = ScanResults$NodeKredGeneral
			
			if (RetSS) {
				
				EssentialSchemata = c(EssentialSchemata, ScanResults$EssSchInfo$EssentialSchemata)
				EssSchCvrdStates = c(EssSchCvrdStates, ScanResults$EssSchInfo$EssSchCvrdStates)
				EssSchOutputs = c(EssSchOutputs, ScanResults$EssSchInfo$EssSchOutputs)
				
			}
			
		}
		
	} # end of 'Dim' loop
	
	# Eliminate redundant schemata
	
	RmIdx = c()
	
	for (scn in 1 : length(EssSchCvrdStates)) {
	  
	  CvrdStates = EssSchCvrdStates[[scn]]
	  
	  SchDim = log(length(CvrdStates), 2)
	  
	  CvrdStatesDim = NodeKredGeneral[CvrdStates]
	  
	  if(all(CvrdStatesDim > SchDim)) RmIdx = c(RmIdx, scn)
	  
	}
	
	EssentialSchemata = EssentialSchemata[ -RmIdx ]
	EssSchCvrdStates = EssSchCvrdStates[ -RmIdx ]
	EssSchOutputs = EssSchOutputs[ -RmIdx ]
	
	# Step 5: Compute $k_r^*$ by averaging over the covering dimensions of all $2^k$ input vectors. 
	# Compute $k_e^* = k - k_r^*$.
	
	if (ComputeKeff) Keff = K - mean(NodeKred)
		
	Ku = K - mean(NodeKredGeneral)
	
	if (ComputeKeff) {

		if (RetSS) {
		  
		  EssentialSchemata = unlist(EssentialSchemata)
		  
		  if (!is.null(EssentialSchemata)) EssentialSchemata = matrix(EssentialSchemata, byrow = T, ncol = K) 
		  
		  EssentialSchemata = rbind(EssentialSchemata, LUT[which(NodeKredGeneral == 0), ])  # include uncovered LUT entries (Dim = 0)
		  
		  if (K == 1) EssentialSchemata = matrix(EssentialSchemata, ncol = K)
		  
		  SchemaOutputs = unlist(EssSchOutputs)
		  
		  SchemaOutputs = c(SchemaOutputs, Func[which(NodeKredGeneral == 0)])  # include outputs with uncovered LUT entries (Dim = 0)
		  
		  return(list(EssentialSchemata = EssentialSchemata, SchemaOutputs = SchemaOutputs, Keff = Keff, Ku = Ku))
		  
		} else return(list(Keff = Keff, Ku = Ku))
	
	} else {
		
		if (RetSS) {
		  
		  EssentialSchemata = unlist(EssentialSchemata)
		  
		  if (!is.null(EssentialSchemata)) EssentialSchemata = matrix(EssentialSchemata, byrow = T, ncol = K) 
		  
		  EssentialSchemata = rbind(EssentialSchemata, LUT[which(NodeKredGeneral == 0), ])  # include uncovered LUT entries (Dim = 0)
		  
		  if (K == 1) EssentialSchemata = matrix(EssentialSchemata, ncol = K)
		  
		  SchemaOutputs = unlist(EssSchOutputs)
		  
		  SchemaOutputs = c(SchemaOutputs, Func[which(NodeKredGeneral == 0)])  # include outputs with uncovered LUT entries (Dim = 0)
		  
		  return(list(EssentialSchemata = EssentialSchemata, SchemaOutputs = SchemaOutputs, Ku = Ku))
		  
		} else return(Ku)
		
	}

}
	

















