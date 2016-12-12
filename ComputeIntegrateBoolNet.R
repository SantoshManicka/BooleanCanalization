# This program integrates a Boolean network for a specified number of time steps.

# What is meant by 'integration'?
# 'Integration' here means the same as integration of a dynamical system. Integration is the method used to solve a 
# dynamical system specified by ODE, say, that can be solved. The outcome is an analytic expression that computes the
# state of the dynamical system after a given time t from an initial condition. If an analytic solution is not possible,
# which is typically the case, numerical integration methods (e.g., Euler method) can still be used to compute phase portraits,
# which are graphical solutions, from which trajectories following initial conditions can be determined.
# For discrete dynamical systems, like Boolean networks, the state transition graph (STG) is the phase portrait obtained by
# iteratively updating the state of the system starting from all possible initial states. The STG thus helps determine the
# state of the network after any given number of discrete time steps starting from an initial state. Thus, we can think of the 
# iterative process that constructs the STG as a kind of integration (akin to Euler integration). 
# We define a more general method of Boolean network integration that "solves" for each node of the network. It constructs a 
# 'schema causation chain' (SCC) -- a chain of a set of Boolean functions that determines the states of the corresponding nodes 
# after a discrete number of time steps t from any initial condition (t is the position of a given set in the chain). In other 
# words, a given Boolean function in the set at a given position in the SCC computes the state of the associated node after t 
# time steps from any initial state of the network. The SCC is more general than the STG in the sense that the latter can be 
# deduced from the former since the state of the entire network after t steps can be computed simply by computing the states of
# all the nodes from a single initial condition off of the SCC.

# Input:
# A Boolean network data structure whose format is specified in package 'BoolNet'. 

# Output:
# A time series called 'SchemaCausationChain' in the form of a list of list of list of schemata matrices; 
# a given matrix is associated with a particular node, time step and output, indexed in that order.

# Optional parameters:
# (1) Parallelization: The parameter 'NumCores' is used for that purpose. The intersection and union steps for a single node 
# may have to deal with a large number of schemata. You can avail the parallel processing option for that purpose. If you are 
# running this program on a multicore computer, you can set 'NumCores' option to set the number of cores you want to use. 
# Parallel processing is particularly useful when you are dealing with a large number of schemata for any node; it could be an 
# inefficient method to process just a few schemata (due to communication overhead).
# (2) Randomization: The parameters 'SetMaxPI' and 'MaxNumPIs' are used for that purpose, and are utilized by the 'UnionSchemataSet' 
# function (see related documentation for more details). 
# (3) Compression nonlinearity: The parameters 'MaxWCLossProp' and 'MaxWCLoss' are used for that purpose, and are utilized by 
# the 'UnionSchemataSet' function (see related documentation for more details). 

# USAGE EXAMPLE:
# require(BoolNet)
# Net <- loadNetwork("./Data/ExampleBooleanNetwork.bn", symbolic = F)  # the .bn file must be in the format specified in the BoolNet package
# SchemaCausationChain = IntegrateBoolNet(Net, NumIters = 2)

# An example of 'Net' (used in the dissertation and paper):
# > Net
# Boolean network with 3 genes
# 
# Involved genes:
# x1 x2 x3
# 
# Transition functions:
# x1 = x1
# x2 = (x2 & x3) | (x1 & !x2 & !x3)
# x3 = (x1 & x2)

# Output for the above 'Net' (the list is indexed as [[node]][[timeStep]][[output+1]]):
# > SchemaCausationChain
# [[1]]
# [[1]][[1]]
# [[1]][[1]][[1]]   # node = 1, time = 1, output = 0
# [,1] [,2] [,3]
# [1,]    0    2    2
# 
# [[1]][[1]][[2]]   # node = 1, time = 1, output = 1
# [,1] [,2] [,3]
# [1,]    1    2    2
# 
# 
# [[1]][[2]]
# [[1]][[2]][[1]]   # node = 1, time = 2, output = 0
# [,1] [,2] [,3]
# [1,]    0    2    2
# 
# [[1]][[2]][[2]]   # node = 1, time = 2, output = 1
# [,1] [,2] [,3]
# [1,]    1    2    2
# 
# 
# 
# [[2]]
# [[2]][[1]]
# [[2]][[1]][[1]]   # node = 2, time = 1, output = 0
# [,1] [,2] [,3]
# [1,]    0    0    2
# [2,]    0    2    0
# [3,]    2    0    1
# [4,]    2    1    0
# 
# [[2]][[1]][[2]]   # node = 2, time = 1, output = 1
# [,1] [,2] [,3]
# [1,]    2    1    1
# [2,]    1    0    0
# 
# 
# [[2]][[2]]
# [[2]][[2]][[1]]   # node = 2, time = 2, output = 0
# [,1] [,2] [,3]
# [1,]    2    2    0
# [2,]    0    2    2
# 
# [[2]][[2]][[2]]   # node = 2, time = 2, output = 1
# [,1] [,2] [,3]
# [1,]    1    2    1
# 
# 
# 
# [[3]]
# [[3]][[1]]
# [[3]][[1]][[1]]   # node = 3, time = 1, output = 0
# [,1] [,2] [,3]
# [1,]    0    2    2
# [2,]    2    0    2
# 
# [[3]][[1]][[2]]   # node = 3, time = 1, output = 1
# [,1] [,2] [,3]
# [1,]    1    1    2
# 
# 
# [[3]][[2]]
# [[3]][[2]][[1]]   # node = 3, time = 2, output = 0
# [,1] [,2] [,3]
# [1,]    0    2    2
# [2,]    2    0    1
# [3,]    2    1    0
# 
# [[3]][[2]][[2]]   # node = 3, time = 2, output = 1
# [,1] [,2] [,3]
# [1,]    1    0    0
# [2,]    1    1    1


require(BoolNet)
require(parallel)

source('ComputePrimeImplicants.R')
source('ComputeSchemaSetOperations.R')

IntegrateBoolNet <- function(Net, NodeNums = 1 : length(Net$genes), NumIters=10, SetMaxPI=T, 
											  	MaxNumPIs=500, MaxWCLossProp=0.0, MaxWCLoss=1, NumCores = 1) {
	
  # For every automaton in the Boolean network, compute the prime implicants or schemata.
  
	LoadFuncPISet <- function() {
		
		FuncSet = list()
		
		for (n in 1 : N) FuncSet = append(FuncSet, list(Net$interactions[[n]]$func))
				
		PISet = lapply(FuncSet, function(Func) {
			
			K = log(length(Func), 2)
			
			if (all(Func == 1) | all(Func == 0)) {
				
				PIs = list()
				
				if (K > 0) PIs[[1]] = rep(2, K) else PIs[[1]] = 2
				
				PIOutputs = Func[1]
				
				list(PIs = matrix(unlist(PIs), nrow = 1, ncol = 1, byrow = T), Outputs = PIOutputs, NumPIs = length(PIs))
				
			} else {
			
				PIsData = ComputePrimeImplicants(Func, K, AllPrimeImplicants = T)
				
				PIs = PIsData$PrimeImplicants
				
				PIOutputs = PIsData$Outputs
				
				NumPIs = nrow(PIs)
				
				list(PIs = PIs, Outputs = PIOutputs, NumPIs = NumPIs)
						
			}
		
		})
		
		return(list(FuncSet = FuncSet, PISet = PISet))
	
	}
	
	# 'SchemataOutputSet' is a list containing two items namely 'SchemataSet' and 'OutputsSet'. 
	# 'SchemataSet' is a list of N items; each item is a 'schemata matrix' associated with a single node. 
	# A schemata matrix is a N-column matrix whose rows contain the schemata.
	# 'OutputsSet' is a list of N items containing the set of output vectors corresponding to each node, where
	# each output value corresponds to a schema in the associates Schemata matrix.
	# 'RefSchemataOutputSet' is a referece 'SchemataOutputSet' that contains the original schemata and the corresponding
	# outputs that define the Boolean network.
	
	ComputeRefSchemataOutputSet <- function(PISet, InputsPerNode, NumPIsPerNode) {
		
		SchemataSet = list()
		
		OutputsSet = list()
		
		for (node in 1 : N) {
			
			SchemataMatrix = matrix(2, nrow = NumPIsPerNode[node], ncol = N)
			
			inputs = InputsPerNode[[node]]
			
			SchemataMatrix[ , inputs] = PISet[[node]]$PIs
			
			SchemataSet[[node]] = SchemataMatrix
			
			OutputsSet[[node]] = PISet[[node]]$Outputs
			
		}
		
		SchemataOutputSet = list(SchemataSet = SchemataSet, OutputsSet = OutputsSet)
		
		return(SchemataOutputSet)
	
	}

	# Propagate constant states through the network, and simplify 'SchemataOutputSet' in the process.
	
	# Basic idea: 
	# If a k-input function has m constant-state inputs, then it is effectively a function with (k - m) inputs.
  
	# Method:
	# Remove schemata with constant enput nodes and whose values contradict the corresponding constant enput states.
	
	SimplifySchemataOutputSet <- function(SchemataOutputSet) {
	  
	  OutputsSet = SchemataOutputSet$OutputsSet
	  
	  ConstNodes = which(sapply(OutputsSet, function(v) all(v == 0) | all(v == 1)))
	  
	  ChangeDetected = T
	  
	  # This loop takes into account the possibility that nodes could be iteratively rendered constant by other 
	  # constant nodes.
	  
	  while ((length(ConstNodes) > 0) & ChangeDetected) {  
	    
	    ChangeDetected = F
	    
	    ConstNodeOutputs = sapply(ConstNodes, function(c) OutputsSet[[c]][1])
	    
	    NewSchemataOutputSet = SchemataOutputSet
	    
	    for (node in 1 : N) {
	      
	      SchSet = NewSchemataOutputSet$SchemataSet[[node]]
	      OutSet = NewSchemataOutputSet$OutputsSet[[node]]
	      
	      ImpslSchIdx = apply(SchSet, 1, function(v) any((v[ConstNodes] != ConstNodeOutputs) & (v[ConstNodes] != 2)))
	      
	      # Remove "impossible" schemata with enput values contradicting the constant state of that input
	      
	      SchSet = matrix(SchSet[ !ImpslSchIdx, ], ncol = N)
	      OutSet = OutSet[ !ImpslSchIdx ]
	      
	      # The above leads to one of two possible outcomes:
	      # (1) The current node is forced to be a constant due to upstream constant nodes, in which case all of its
	      # inputs could be severed, thus  maximally simplifying or compressing the corresponding schemata.
	      # (2) The current node is not forced to be a constant, in which case any depending on the inputs with constant
	      # states could be severed, thus simplifying the schemata to some extent.

	      if (all(OutSet == 0) | all(OutSet == 1)) {  # a forced-constant node
	        
	        ConstOut = OutSet[1]
	        
	        NewSchemataOutputSet$SchemataSet[[node]] = matrix(rep(2, N), ncol = N)
	        NewSchemataOutputSet$OutputsSet[[node]] = ConstOut
	        
	      } else {  # both output values have at least one schema mapping to it => 
	                # at least two rows in SchSet => at least one non-constant input with two different values
	        
	        NonConstInputVars = setdiff(1 : N, ConstNodes)
	        
	        NumNonConstInputs = length(NonConstInputVars)
	        
	        # Remove any dependency on the constant inputs
	        
	        SchSet = matrix(SchSet[ , NonConstInputVars], ncol = NumNonConstInputs)
	        
	        # Simplify or compress the schemata resulting from the last step. 
	        # Here is a simple example to understand this step. Consider the following set of schemata: {#01, 10#}.
	        # Suppose that node 1 is a constant node. So, we remove the dependency on that node resulting in a new
	        # schemata set: {#01, #0#}. Clearly this set can be further compressed to obtain: {#0#}.
	        
	        # Note that since there is at least one schema mapping to each output, SchemataSet0 and SchemataSet1 
	        # below can never be empty.
	        
	        SchemataSet0 = list()
	        SchemataSet0[[2]] = matrix(SchSet[OutSet == 0, ], ncol = NumNonConstInputs)
	        
	        Sch0 = UnionSchemataSet(SchemataSet = SchemataSet0, MaxWCLossProp = 1, MaxWCLoss = Inf, InitNewSet = T)
	        
	        SchemataSet1 = list()
	        SchemataSet1[[2]] = matrix(SchSet[OutSet == 1, ], ncol = NumNonConstInputs)
	        
	        Sch1 = UnionSchemataSet(SchemataSet = SchemataSet1, MaxWCLossProp = 1, MaxWCLoss = Inf, InitNewSet = T)
	        
	        SchSubSet = rbind(Sch0, Sch1)
	        SchSet = matrix(2, nrow(SchSubSet), N)
	        SchSet[ , NonConstInputVars] = SchSubSet
	        
	        OutSet = c(rep(0, nrow(Sch0)), rep(1, nrow(Sch1)))
	        
	        NewSchemataOutputSet$SchemataSet[[node]] = SchSet
	        NewSchemataOutputSet$OutputsSet[[node]] = OutSet
	        
	      }

	      if (any(ImpslSchIdx)) ChangeDetected = T
	      
	    }

	    SchemataOutputSet = NewSchemataOutputSet
	    
	    OutputsSet = SchemataOutputSet$OutputsSet
	    
	    ConstNodes = which(sapply(OutputsSet, function(v) all(v == 0) | all(v == 1)))
	    
	  }
	  
	  return(SchemataOutputSet)
	  
	}
	
	# Compute the predecessors of a schema:
	# A schema is nothing but a partial state where the states of a subset of nodes is known and the others unknown
	# (represented as a wildcard symbol '#'). The predecessors of a given schema are the schemata that guarantee the
	# the truth of that schema in the following time step. In other words, the predecessors are the schemata representing
	# a set of logical conditions that guarantee the logical condition specified by the successor schema in the subsequent
	# time step.
	
	# Method:
	# In a schema replace every literal (0 or 1) with the schemata of the corresponding node that outputs that value
	# (fetched from RefSchemataOutputSet) and intersect (logical conjunction) the same.
	
	ComputeSchemaPredecessors <- function(Idx, Schemata, RefSchemataOutputSet) {
		
		Schema = Schemata[Idx, ]
		
		EnputIdx = which(Schema == 0 | Schema == 1)
		
		SchemataSet = list()
		
		ImpsblSchema = F
		
		for (i in 1 : length(EnputIdx)) {
			
			EnputNode = EnputIdx[i]
			
			EnputValue = Schema[EnputNode]
			
			SchIdx = which(RefSchemataOutputSet$OutputsSet[[EnputNode]] == EnputValue)
			
			# Note that the following condition should not be satisfied if 'SimplifySchemataOutputSet' had been called.
			# We still have it for clarity.
			
			if (length(SchIdx) == 0) {  # either contradiction (EnputValue = 1) or tautology (EnputValue = 0)
				
				ImpsblSchema = T
				
				break
			
			}
			
			Schemata = RefSchemataOutputSet$SchemataSet[[EnputNode]][SchIdx, ]
			
			Schemata = matrix(Schemata, ncol = N)
			
			SchemataSet[[i]] = Schemata
			
		}
		
		# Note that if SchemataSet contains only one schemata matrix, IntersectSchemataSet will automatically
		# return it without any processing (so it's ok not to handle that condition here).
		
		if ( !ImpsblSchema ) PredecessorSchemata = IntersectSchemataSet(SchemataSet, MaxWCLossProp, MaxWCLoss, 
		                                                                SetMaxPI, MaxNumPIs)
		else PredecessorSchemata = NULL
		
		return(PredecessorSchemata)
		
	}
	
	N = length(Net$genes)
	
	FuncPISet = LoadFuncPISet()
	
	FuncSet  = FuncPISet$FuncSet
	
	PISet = FuncPISet$PISet
	
	# NOTE: The ordering of Net$interactions[[n]]$input absolutely matters
	
	InputsPerNode = lapply(1 : N, function(n) Net$interactions[[n]]$input)
	
	NumPIsPerNode = sapply(PISet, function(PI) PI$NumPIs)
	
	RefSchemataOutputSet = ComputeRefSchemataOutputSet(PISet, InputsPerNode, NumPIsPerNode)
	
	# Propagate constant node states through the network
	
	RefSchemataOutputSet = SimplifySchemataOutputSet(RefSchemataOutputSet)
	
	# Main loop constituting the Boolean network integration method:
	# For every node and for each output, do the following:
	# (1) Compute predecessors (Intersection step): Compute the predecessors of every schema that go to the given output, 
	#     to obtain a set of predecessor schemata. This step is named 'intersection' because a schema could contain multiple 
	#     literals or enput values each specifying a single logical condition all of which must be satisfied by the predecessor
	#     schemata; thus the predecessor schemata associated with each literal in the given schema must be combined
	#     by logical conjunction or intersection.
	# (2) Compress predecessors (Union step): Compress the set of predecessor schemata obtained in step (1) by logical disjunction 
	#     or union.
	
	# Interpretation: 
	# If the set of schemata at the beginning of the iteration maps to output 'x' of node 'i' after 't' time steps,
	# then the set of schemata obtained at the end of iteration maps to the same after 't+1' time steps.
	# The original schemat step contained in 'RefSchemataOutputSet' has a mapping time scale equal to 1.

	SchemaCausationChain = list()
	
	for (NodeNum in NodeNums) {
	  
	  # First fill iter = 1 data, which is nothing but the data contained in RefSchemataOutputSet.
		
	  ConstantNode = F

	  # 'SchemataParts' is a list containing two items; each item contains a matrix of schemata mapping
	  # to an output, where the first matrix maps to output 0, and the second matrix maps to output 1.
	  
		PrevSchemataParts = list()  # associated with a single node NodeNum
		
		for (output in 0 : 1) {
		  
		  SchIdx = which(RefSchemataOutputSet$OutputsSet[[NodeNum]] == output)
		  
		  if (length(SchIdx) > 0) {
		  
		    Schemata = RefSchemataOutputSet$SchemataSet[[NodeNum]][SchIdx, ]
		    
		    Schemata = matrix(Schemata, ncol = N)
		    
		    PrevSchemataParts[[output + 1]] = Schemata
		    
		  } else {
		    
		    PrevSchemataParts[[output + 1]] = NA
		    
		    ConstantNode = T
		    
		  }
		  
		}
		
		SchemaCausationChain[[NodeNum]] = list()
		
		SchemaCausationChain[[NodeNum]][[1]] = PrevSchemataParts
		
		if (ConstantNode) next  # constant state is an "attractor"
		
		# Now begin the integration process from iter = 2 onwards.

		AttrctRchd = F
		
		for (iter in 2 : NumIters) {
			
			SchemataParts = list()
			
			for (output in 0 : 1) {
			
				Schemata = PrevSchemataParts[[output + 1]]
				
				# Step 1: Compute predecessor schemata
				
				SchemataSet = mclapply(1 : nrow(Schemata), ComputeSchemaPredecessors, Schemata, RefSchemataOutputSet,
				                       mc.cores = NumCores)
				
				NullIdx = unlist(lapply(SchemataSet, is.null))
				
				SchemataSet = SchemataSet[ !NullIdx ]
				
				SchemataSetLength = length(SchemataSet)
				
				# Step 2: Compress predecessor schemata
				
				while (SchemataSetLength > 1) {
				  
				  # The NumCores parameter specifies the number of parallel processing cores over which to run
				  # the union of the schemata matrices in SchemataSet. However, since at least two schemata matrices 
				  # are needed to pass to a single instance of the union function, the number of cores needed is
				  # dynamically determined, and specified in the parameter NC, depending on the number of available 
				  # schemata matrices in the set.
					
					NC = NumCores
					
					NumSchPerCore = floor(SchemataSetLength / NC)
					
					if (NumSchPerCore < 2) {  # at least two schemata matrices are needed for union
						
						NumSchPerCore = 2
						
						NC = floor(SchemataSetLength / NumSchPerCore)	
						
						NumSchLastCore = NumSchPerCore + (SchemataSetLength %% NumSchPerCore)
						
					} else {
				
						NumSchLastCore = NumSchPerCore + (SchemataSetLength %% NC)
					
					}
						
					NumSchPerCoreList = c(rep(NumSchPerCore, NC - 1), NumSchLastCore)
					
					SchemataSet = mclapply(1 : NC, UnionSchemataSet, SchemataSet, NumSchPerCoreList,
					                       MaxWCLossProp, MaxWCLoss, F, SetMaxPI, MaxNumPIs,
					                       mc.cores = NC)
	
					SchemataSetLength = length(SchemataSet)
					
				}
				
				# Note that the union of a non-empty set of schemata can never be empty.
				
				if (SchemataSetLength == 1) {
			
					Schemata = SchemataSet[[1]]
													
					SchemataParts[[output + 1]] = Schemata
					
					# MacroKeffEst = ComputeKeffIgnoreOverlaps(Schemata)
					
					# MacroKeffEst = round(MacroKeffEst, 2)
					
					if (AttrctRchd) break  

				} else {  # SchemataSetLength = 0; intersection of Schemata is empty or NULL
					
					SchemataParts[[output + 1]] = NA
					
					AttrctRchd = T
					
					if (output == 1) break
									
				}  # end of 'if (SchemataSetLength == 1)'
				
			}  # end of 'for (output in 0 : 1)'
			
			SchemaCausationChain[[NodeNum]][[iter]] = SchemataParts
								
			PrevSchemataParts = SchemataParts
			
			if (AttrctRchd) break  # proceed to the next node
		
		}  # end of 'for (iter in 2 : NumIters)'
		
	}  # end of 'for (NodeNum in 1 : N)'
	
	return(SchemaCausationChain)

}
																		



















