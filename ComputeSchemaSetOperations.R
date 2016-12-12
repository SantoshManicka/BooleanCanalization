# Disjunct a pair of schemata:

# A schema represents a logical condition, in particular a conjunctive clause, also known as implicant.
# A union of two schemata is thus a union of two logical conditions. The goal of the following function is to identify any
# essential schema that subsumes both of the uniting schemata or an additional 'non-essential' schema that is not fully covered 
# by either of the uniting schemata; such schemata "bridge" portions of the uniting schemata. In other words, given a set of 
# implicants or prime implicants, that the schemata represent, this function identifies any encompassing essential or additional 
# non-essential prime implicant that does not contradict the logical conditions specified by the uniting schemata. For example,
# (#10) is a bridge and a non-essential PI uniting the schemata (0#0) and (11#).

# The function works by sequentially considering each position in the schemata and combining the symbols located there. 
# The rules of the combination are: {0, 0} => 0; {1, 1} => 1; {0, 1} => # (allowed only once); {0, #} => 0; {1, #} => 1.
# These rules answer the question: what pieces of the logical conditions are common to, or are irrelevant in, the schemata? So,
# they involve both logical conjunction and disjunction operations (intuitively relates to the notion of the "bridge"). Perhaps, 
# there is a way to deduce these rules from the axioms of Boolean logic, but for now let them just appeal to the intuition.

# Parameters 'MaxWCLossProp' and 'MaxWCLoss' control the proportion of wildcards (WC) and the absolute number of WCs respectively
# that can be "lost" by the schema with fewer WCs, and as a consequence, the number of WCs in the resulting schema. The key point 
# is that when a non-esential PI is produced in a union, both of the uniting schemata necessarily "lose" one or more WCs. 
# For example, 00#0#0 U ##111# => 001#10; the first schema loses one WC and the second schema loses two WCs to produce a schema, 
# a non-essential PI, containing a single WC gained from a piece of both schemata.

# Setting MaxWCLossProp = 1 and/or MaxWCLoss = Inf will remove all constraints on the number of wildcard symbols contained
# in the resultant schema; that is, the number of WCs in the non-essential PI produced can be any in number. Note that even though 
# these parameter values suggest that all wildcard symbols could be "lost" by the smallest of the combining schemata, a schema 
# with no wilcard symbol will never be produced since such a "schema" would be covered by both of the combining schemata, meaning 
# it can't be a prime implicant; they rather suggest that the resultant schema could contain just one wildcard regardless of the 
# number of those in the combining schemata (see last example in the list below).
# On the other hand, setting MaxWCLossProp = 0 and/or MaxWCLoss = 0 will still produce a non-essential PI with the same number
# of #s as in the schema with the smallest number of those, if any.

# If no schema could be produced either because none exists or because the parameter boundaries are violated, NULL is returned. 
# If one of the uniting schemata is covered by the other it is marked "duplicate"; if covered by the resultant schema, it is 
# marked "subsumed".

# Conceptual examples:
# (00#) U (01#) => (0##); resultant is an (encompassing) essential PI
# (0#0) U (11#) => (#10); resultant is a non-essential PI
# (10#0) U (#11#) => (1#10); resultant is a non-essential PI
# (1#) U (#1) => NULL; neither an essential PI nor an additional non-essential PI exists
# (1##0) U (#11#) => NULL; neither an essential PI nor an additional non-essential PI exists
# (#1) U (#1) => NULL; first schema is a "duplicate" as it's the same as the first schema
# (#1) U (01) => NULL; second schema is a "duplicate" as it's covered by the first schema
# (#1) U (10) => (1#); second schema is a "subsumed" schema as it's covered by the resultant schema
# (01#) U (1##) => (#1#); first schema is a "subsumed" schema as it's covered by the resultant schema
# (10) U (11) => (1#); both uniting schemata are "subsumed" as they are both covered by the resultant schema
# (0#) U (1#) => (##); both uniting schemata are "subsumed" as they are both covered by the resultant schema
# (0#0#0) U (#111#) => (01#10); notice that the resultant schema has fewer WC than either of the uniting schemata

# Examples of output:
# (00#) U (01#) => Mismatches = 1, SubsumedSch = 3, UniSchema = (0##), DuplSch = F.	
# (01#) U (01#) => Mismatches = 0, SubsumedSch = -1, UniSchema = NULL, DuplSch = T.	
# (01#) U (0##) => Mismatches = 0, SubsumedSch = 1, UniSchema = NULL, DuplSch = T.	
# (0#0) U (11#) => Mismatches = 1, SubsumedSch = 0, UniSchema = (#10), DuplSch = F.	
# (0#1##1) U (#0#100) => Mismatches = 1, SubsumedSch = 0, UniSchema = 00100#, DuplSch = F.	
# (0##) U (##1) => Mismatches = 0, SubsumedSch = 0, UniSchema = NULL, DuplSch = F (such cases will be ignored).

UnionSchemata <- function(sch1, sch2, MaxWCLossProp = 0, MaxWCLoss = 1) {
  
  N = length(sch1)
  
  UniSchema = rep(-1, N)
  
  NumWC = c(sum(sch1 == 2), sum(sch2 == 2))
  
  MinWCSch = ifelse(NumWC[1] < NumWC[2], 1, 2)
  
  MaxWCLossP = floor(MaxWCLossProp * NumWC[MinWCSch])
  
  MaxWCLoss = max(MaxWCLossP, MaxWCLoss)
  
  NumWCLost = c(0, 0)
  
  SubsumedSch = -1
  Mismatches = 0
  DuplSch = F
  MoreWCLost = F
  
  # SubsumedSch = -1 if and only if the uniting schemata are exactly the same.
  # SubsumedSch = 0 if neither schemata covers the other, and thus remains 0 once set to 0.
  
  for (i in 1 : N) {
    
    s1 = sch1[i]
    s2 = sch2[i]
    
    if (s1 != s2) {
      
      if (s1 == 2) {
        
        if (SubsumedSch == 1 | SubsumedSch == 0) SubsumedSch = 0 else SubsumedSch = 2
        
        NumWCLost[1] = NumWCLost[1] + 1
        
        UniSchema[i] = s2
        
      } else if (s2 == 2) {
        
        if (SubsumedSch == 2 | SubsumedSch == 0) SubsumedSch = 0 else SubsumedSch = 1
        
        NumWCLost[2] = NumWCLost[2] + 1
        
        UniSchema[i] = s1
        
      } else {
        
        if (SubsumedSch == -1) SubsumedSch = 3
        
        UniSchema[i] = 2
        
        Mismatches = Mismatches + 1
        
        NumWCLost = NumWCLost - 1  # a WC is gained here (a max of only one WC can be gained)
        
        if (Mismatches > 1) {SubsumedSch = 0; UniSchema = NULL; break}  # e.g., (100) U (111) => NULL
        
      }
      
    } else {
      
      UniSchema[i] = s1
      
    }
    
    # If more WC has been lost even after gaining (a max of) one WC, or more WC would have been lost even if a WC could
    # be gained in the future, then break.
    
    BreakProcess = ((NumWCLost[MinWCSch] > MaxWCLoss) & (Mismatches == 1)) | (NumWCLost[MinWCSch] > (MaxWCLoss + 1))
    
    if (BreakProcess) {UniSchema = NULL; SubsumedSch = 0; MoreWCLost = T; break}
    
  }
  
  if (NumWCLost[MinWCSch] > MaxWCLoss) {UniSchema = NULL; SubsumedSch = 0; MoreWCLost = T}
  
  # No mismatches means that either the resultant schema is covered by both of the uniting schemata, or one of the 
  # uniting schemata is covered by the other.
  
  if (Mismatches == 0 & !MoreWCLost) { 
    
    UniSchema = NULL
    
    if (SubsumedSch != 0) DuplSch = T   # A duplicate (or subsumed) schema could not have lost WC
    
  }
  
  return(list(UniSchema = UniSchema, SubsumedSch = SubsumedSch, DuplSch = DuplSch))
  
}

# ------------------------------------------------------------------------------------------------------------------------------------ #

# Disjunct or compress a set of schemata matrices:

# Basic idea:
# This is an optimization method where the goal is to maximally compress a set of schemata --- very much like Quine-McCluskey. 
# Unlike Quine-McCluskey, however, where smaller schemata of the same sizes always combine to become larger, in this method schemata 
# of any sizes combine to produce schemata of any other size (larger or smaller). Thus, this method may be thought of as a nonlinear 
# version of Quine-McCluskey, where the nonlinearity may be necessary for optimization, depending on the problem.

# Outline of the method:
# Two sets of schemata SchSet1 and SchSet2 in the main set of all schemata sets called 'SchemataSet' are combined to produce two sets 
# one of which shall contain schemata from SchSet1 and SchSet2 itself and the other shall contain new schemata produced from union. 
# These two sets combine further and so on until a single set remains. This set will similarly combine with other sets in 'SchemataSet' 
# until a single SchSet1 remains which is finally returned.
# Combining SchSet1 and SchSet2: Every schema in SchSet1 (main loop) is combined with every schema in SchSet2 (sub loop). 
# The resultant schemata are added to newSchSet2. Subsumed schemata in SchSet2 are marked and are prohibited from combining with 
# other schemata in SchSet1. Subsumed schemata in SchSet1 are marked as well, and are likewise prohibited from combining further. 
# Every time a schema in SchSet1 is found to be subsumed a "cleanup" operation is undertaken where all the schemata in newSchSet2
# that were produced as a result of combining that schema with those in SchSet2 are removed, and the schemata in SchSet2 that were 
# marked as subsumed during the union with that schema are reset. Once the main loop over SchSet1 is complete, all the schemata in 
# SchSet1 and SchSet2 that are not marked subsumed are put together into a single set that then becomes the new SchSet2 and those 
# in newSchSet2 become the new SchSet1. The above steps are then repeated, only this time onwards SchSet1 is also combined with itself,
# and the new schemata are added to newSchSet2 as before. See Example 2 below to see why SchSet1 needs to combine with itself. 
# The iterations are repeated until newSchSet2 is empty.

# Pruning: 
# For every schema in SchSet1, we first filter out schemata in SchSet2/SchSubSet1 that are not eligible for a union.
# (1) Conditions A and B: This condition stops schemata from combining whose result is guaranteed to be an absolutely 
# non-essential schema. Representative example: (10###) Union (###01) -> NULL; a NULL is returned because the actual resultant 
# schema (10#01) is fully subsumed by each one of the combining schemata.
# When two schemata are allowed to combine, the resultant is either an essential schema that subsumes both of the uniting 
# schemata, or a non-essential schema that bridges them; e.g., (00#) Union (01#) -> 0##; (0#0) Union (11#) -> (#10)
# respectively. The following condition further prunes which schemata can combine based on the MaxWCLoss parameter.
# (2) Condition C: This is similar to the one in Quine-McCluskey algorithm that fixes the exact difference in 
# the number of ones between the two combining schemata with the same number of WCs to be equal to 1. In our case, we allow 
# schemata with any number of WCs to combine while also control the number of WCs in the resultant schema. Thus, we need an
# appropriate method to calculate the maximum difference in NumOnes, with a default minimum equal to 0. The formula below
# calculates the max difference in NumOnes that can be allowed between the combining schemata, and is mainly controlled by
# the specified MaxWCLoss.
# Formula: MaxDiffNumOnes = MaxNumWC - MinNumWC + min(MaxWCLoss + 1, MinNumWC) + 1  
# where, Max/MinNumWC represents the max/min WC of the two schemata; the '+1' is to include cases where a 0 from one schema 
# combines with a 1 from the other to form a # (a '#' could be gained, that is).
# Purpose: Two schemata won't be allowed to combine if their NumOnes differ by an amount greater than MaxDiffOnes.
# Interpretation of the term MaxNumWC - MinNumWC:
#   (a) that many NumOnes can be absorbed by the WCs in the schema with MaxNumWC, so that much difference in NumOnes is allowed.
# Interpretation of the term min(MaxWCLoss + 1, MinNumWC) + 1:
#   (a) the first +1 means that an extra WC can be lost since a WC could be gained, so one more NumOnes on top of MaxWCLoss can be
#   absorbed by the schema with MaxNumWC; 
#   (b) the min refers to the fact that more WC can't be lost than MinNumWC, laying an upper bound on (a);
#   (c) the last +1 refers to the fact that a WC could be gained, so that's one more NumOnes.
# To further understand this method, consider the following examples specifying different schemata and parameter combinations
# (for each MaxWCLoss specification, the smaller of the two schemata could lose only up to that many WCs). 
# Examples:
# (3a1) MaxWCLoss = 0 => MaxDiffOnes = 3 - 2 + min(1, 2) + 1 = 3; 0##00# Union 1##001 -> ###001 
# (3a2) MaxWCLoss = 1 => MaxDiffOnes = 3 - 2 + min(2, 2) + 1 = 4; 0##00# Union 1##001 -> ###001 
# (3a3) MaxWCLoss = 2 => MaxDiffOnes = 3 - 2 + min(3, 2) + 1 = 4; 0##00# Union 1##001 -> ###001
# (3b1) MaxWCLoss = 0 => MaxDiffOnes = 3 - 2 + min(1, 2) + 1 = 3; 0##00# Union 1#10#1 -> ##1001
# (3b2) MaxWCLoss = 1 => MaxDiffOnes = 3 - 2 + min(2, 2) + 1 = 4; 0##00# Union 1#10#1 -> ##1001
# (3b3) MaxWCLoss = 2 => MaxDiffOnes = 3 - 2 + min(3, 2) + 1 = 4; 0##00# Union 1#10#1 -> ##1001
# (3c1) MaxWCLoss = 0 => MaxDiffOnes = 3 - 1 + min(1, 1) + 1 = 4; 0##00# Union 1110#1 -> #11001
# (3c2) MaxWCLoss = 1 => MaxDiffOnes = 3 - 1 + min(2, 1) + 1 = 4; 0##00# Union 1110#1 -> #11001
# (3c3) MaxWCLoss = 2 => MaxDiffOnes = 3 - 1 + min(3, 1) + 1 = 4; 0##00# Union 1110#1 -> #11001
# (3c4) MaxWCLoss = 0 => MaxDiffOnes = 3 - 1 + min(1, 1) + 1 = 4; 0##00# Union 1111#1 -> NULL
# (3c5) MaxWCLoss = Inf => MaxDiffOnes = 3 - 1 + min(1, 1) + 1 = 4; 0##00# Union 1111#1 -> NULL
# (3d1) MaxWCLoss = 0 => MaxDiffOnes = 3 - 2 + min(1, 2) + 1 = 3; 0##00# Union 111##1 -> NULL
# (3d2) MaxWCLoss = 1 => MaxDiffOnes = 3 - 2 + min(2, 2) + 1 = 4; 0##00# Union 111##1 -> #11001
# (3d3) MaxWCLoss = 2 => MaxDiffOnes = 3 - 2 + min(3, 2) + 1 = 4; 0##00# Union 111##1 -> #11001
# An outline of the proof for why MaxDiffOnes calculated this way is indeed the maximum:
# (Fact 1): The schema with MinNumWC could naturally have a higher NumOnes than the other schema. In other words, it is 
# always possible to construct a schema corresponding to MinNumWC with all of its literals filled with ones; such a schema 
# will have a higher NumOnes than the other schema regardless of its NumOnes. Hence, the absolute max difference in NumOnes
# between the two schemata, regardless of whether they can combine or not, is equal to (k - MinNumWC); this term is only to
# be noted and does not play a direct role in the formula.
# (Fact 2): For two schemata to combine, they can contain at most one position where one of the schemata contains a 0 and
# the other a 1. When this happens, we say that a '#' is "gained". The rest of the 0s and 1s of a schema must be positioned 
# in such a way that the corresponding position in the other schema contains a '#'. When this happens, we say that the 1s and 
# 0s involved are "absorbed", and the corresponding #s are "lost". Any number of #s can be lost. A pair of schemata that 
# satisfy these constraints is said to be "valid". 
# Our goal in this proof is to construct a valid pair of schemata A and B such that the difference in their NumOnes is maximum
# for the given MaxNumWC and MinNumWC. We assume that A contains MaxNumWC #s and B contains MinNumWC #s. First, we count the  
# number of positions in each schema of the pair that contribute to the difference. PA1 and PB1 are the number of positions where 
# there is a 0 in one schema and a 1 in the same positions in the other schema. Fact 2 implies that PA1 = PB1 = 1. Next, PA2 
# is the number of positions where there are 0s or 1s in schema A and #s in the same positions in B. Likewise, PB2 is the number 
# of positions where there are 0s or 1s in B and #s in the same positions in A. It is clear that PB2 >= PB1. All other positions 
# in A and B must necessarily contain the same symbols (0, 1 or #) in order for them to remain valid; thus, they don't contribute 
# to the difference in NumOnes. Thus, the difference in NumOnes is equal to total NumOnes in (PB1 + PB2) minus total NumOnes in 
# (PA1 + PA2). Thus, max difference can be achieved by filling PB1 and PB2 with 1s and PA1 and PA2 with 0s, resulting in the
# value equal to (PB1 + PB2) = (1 + PB2). It is clear that PB2 = MaxNumWC - MinNumWC + min(MaxWCLoss + 1, MinNumWC), for reasons
# mentioned above. Thus, MaxDiffNumOnes = MaxNumWC - MinNumWC + min(MaxWCLoss + 1, MinNumWC) + 1.
# As a special case, consider MaxNumWC = MinNumWC and MaxWCLoss = 0, as in Quine-McCluskey. In this case, MaxDiffNumOnes = 2. In the
# Quine-McCluskey procedure, the MaxDiffNumOnes is actually set at 1 (there are no lower or upper bounds). In our method, since we 
# allow non-essential schemata to form as well, we compute the MaxDiffNumOnes as 2 for schemata with the same number of WCs with 
# MaxWCLoss = 0. So, for example, (00#) and (1#1) would combine to produce the non-essential schema (#01); such combinations are not 
# allowed in Quine-McCluskey.

# Optional parameters:
# (1) Parallelization: The parameters 'Idx' and 'NumSchPerCoreList' are used for this purpose. The idea is that the list 
# 'SchemataSet' can be split into blocks and fed to parallel instances of the union function. So, 'Idx' represents the block
# id or index (an integer), and 'NumSchPerCoreList' is the list of the number of schemata per block. These parameters will be 
# prepared by the calling function and not the union function itself, which will also process the results. For an example of
# how to set these parameters, see 'ComputeIntegrateBoolNet.R'.
# (2) Randomization: The parameters 'SetMaxPI', 'MaxNumPIs' and 'RndPow' are used for this purpose. These parameters control
# whether or not there is an upper bound to the number of prime implicants (PIs) or schemata being processed, the number of
# those if there is a bound and how they are chosen respectively. The first two parameters are self-explanatory. The 'RndPow'
# parameter takes an integer greater than or equal to 0, and controls the extent to which the number of wildcards in a schema
# influences its odds of being chosen; a value of 0 implies an uniform distribution, while higher values implies a distribution
# biased more towards schemata with larger numbers of wildcards.
# (3) Compression nonlinearity: The parameters 'MaxWCLossProp' and 'MaxWCLoss' are used for this purpose. They determine the 
# extent to which schemata unions are (non)linear. Larger the values of these parameters, the more nonlinear the union tends to
# be. See examples under 'Pruning' for details on what these parameters mean.
# (4) Compressing a single set of schemata or generating the set of ALL prime implicants given a set of implicants or prime
# implicants: The parameter 'InitNewSet' is used for this purpose. Normally, when InitNewSet = F, the union method will assume
# that each schemata matrix in the set passed to it is maximally compressed. Then, during the union process, new schemata matrices
# may be generated that generally need further compression; InitNewSet is set to T to indicate that. If, on the other hand, you
# have a schemata matrix that is not maximally compressed or does not contain all possible schemata (non-essential PIs), and
# which needs maximal compression or requires to be expanded with all possible schemata, then do the following:
# Set the following parameters with the listed values: MaxWCLossProp = 1, MaxWCLoss = Inf, InitNewSet = T
# Create a list, set its first item to NULL and its second item with a matrix containing the given set of schemata.
# That is, do the following:
# SchSet = list()
# SchSet[[1]] = NULL
# SchSet[[2]] = SchemataMatrix
# NewSchemataMatrix = UnionSchemataSet(SchemataSet = SchSet, MaxWCLossProp = 1, MaxWCLoss = Inf, InitNewSet = T)
# 'NewSchemataMatrix' will compress all of the implicants in 'SchemataMatrix' and includes those that are non-essential --- all
# possible, but unnecessary, schemata that do not violate the logical conditions specified by the original set of schemata.

# Example 1: 
# {10##, 0100} U {11#0}
# 1st iter of Main loop: 10## U 11#0 => 1##0; 11#0 marked 'sub'
# 2nd iter of Main loop: SchSet2 is empty since the only schema in it is subsumed.
# New sets to combine: {1##0} U {10##, 0100}
# 1st iter of Main loop: 1##0 U 10## => NULL
# 2nd iter of Main loop: 1##0 U 0100 => #100
# New sets to combine: {#100} U {1##0, 10##, 0100}
# 1st iter of Main loop: #100 U 1##0 => NULL
# 2nd iter of Main loop: #100 U 10## => 1#00
# 3rd iter of Main loop: #100 U 0100 => NULL; 0100 marked 'sub'
# New sets to combine: {1#00} U {#100, 1##0, 10##}
# 1st iter of Main loop: 1#00 U #1#0 => NULL 
# 1st iter of Main loop: 1#00 U 1##0 => NULL; 1#00 marked 'sub'; break loop
# New sets to combine: NULL U {#100, 1##0, 10##}
# Return: {#100, 1##0, 10##}

# Example 2:
# {00, 11} U {10, 01}
# 1st iter of Main loop: 00 U 10 => #0; 00 and 10 marked 'sub'; break loop
# 2nd iter of Main loop: 11 U 01 => #1; 11 and 01 marked 'sub'; break loop
# New sets to combine: {#0, #1} U NULL
# Combine SchSet1 with itself
# 1st iter of Main loop: #0 U #1 => ##; #0 and #1 marked 'sub'; break loop
# New sets to combine: {##} U NULL
# Return: {##}

UnionSchemataSet <- function(Idx = 1, SchemataSet, NumSchPerCoreList = length(SchemataSet), 
                             MaxWCLossProp = 0, MaxWCLoss = 1, InitNewSet = F, 
                             SetMaxPI = T, MaxNumPIs = 1000, RndPow = 0) {
  
  st = ((Idx - 1) * NumSchPerCoreList[1]) + 1
  nd = st + NumSchPerCoreList[Idx] - 1
  
  SchemataSet = SchemataSet[st : nd]
  
  SchemataSet = SchemataSet[order(sapply(SchemataSet, length))]
  
  N = ncol(SchemataSet[[2]])
  
  while (length(SchemataSet) > 1) {
    
    if ( !InitNewSet ) {
      
      SchSet1 = SchemataSet[[1]]
      
      SchSet2 = SchemataSet[[2]]
      
    } else {  # length(SchemataSet) = 2 and SchemataSet[[1]] = NULL guaranteed
      
      SchSet2 = SchemataSet[[1]]
      
      SchSet1 = SchemataSet[[2]]
      
    }
    
    NewSet = InitNewSet
    
    while ( !is.null(SchSet1) ) {
      
      if (SetMaxPI) {
        
        if (nrow(SchSet1) > MaxNumPIs) {
          
          NumEnputs = apply(SchSet1, 1, function(v) sum(v < 2))
          
          randPIIdx = sample(1 : length(NumEnputs), MaxNumPIs, prob = 1 / NumEnputs ^ RndPow) 
          
          SchSet1 = SchSet1[randPIIdx, ]
          
          SchSet1 = matrix(SchSet1, ncol = N)
          
        }
        
        if ( !is.null(SchSet2) ) {
          
          if (nrow(SchSet2) > MaxNumPIs) {
            
            NumEnputs = apply(SchSet2, 1, function(v) sum(v < 2))
            
            randPIIdx = sample(1 : length(NumEnputs), MaxNumPIs, prob = 1 / NumEnputs ^ RndPow) 
            
            SchSet2 = SchSet2[randPIIdx, ]
            
            SchSet2 = matrix(SchSet2, ncol = N)
            
          }
          
        }
        
      }
      
      newSchSet2 = list()
      newSchSet2Idx = 1
      
      SchInfo1 = t(apply(SchSet1, 1, function(v) {
        enpidx = which(v < 2)
        numwc = sum(v == 2)
        numones = sum(v[enpidx])
        if (length(enpidx) > 0) mnmx = c(min(enpidx), max(enpidx))
        else mnmx = c(-Inf, Inf)  # if the schema is ####...
        c(mnmx, numones, numwc)
      }))
      
      SubSchInd1 = rep(F, nrow(SchSet1))   # subsumed schema indicator
      
      if ( !is.null(SchSet2) ) {
        
        SchInfo2 = t(apply(SchSet2, 1, function(v) {
          enpidx = which(v < 2)
          numwc = sum(v == 2)
          numones = sum(v[enpidx])
          if (length(enpidx) > 0) mnmx = c(min(enpidx), max(enpidx))
          else mnmx = c(-Inf, Inf)	  # if the schema is ####...
          c(mnmx, numones, numwc)
        }))
        
        SubSchInd2 = rep(F, nrow(SchSet2))
        
      }
      
      NumSet1Sch = nrow(SchSet1)
      
      for (i in 1 : NumSet1Sch) {  # Main loop
        
        if (SubSchInd1[i]) next  # applies when NewSet = T and SchSet1 is combining with itself
        
        sch1 = SchSet1[i, ]
        
        Info1 = SchInfo1[i, ]
        
        # Prune search
        
        newSet2IdxGenByThisSch = c()
        schSet2SubIndSetBythisSch = c()
        schSet1SubIndSetBythisSch = c()
        
        if ( !is.null(SchSet2) ) {
          
          MinMaxNumWC = sapply(SchInfo2[ , 4], function(val) c(min(Info1[4], val), max(Info1[4], val)))
          MinNumWC = MinMaxNumWC[1, ]
          MaxNumWC = MinMaxNumWC[2, ]
          MaxWCLossP = MaxWCLossProp * MinNumWC
          MaxWCLoss = sapply(MaxWCLossP, max, MaxWCLoss)  
          A = Info1[1] > SchInfo2[ , 2]
          B = (SchInfo2[ , 1] > Info1[2])
          MaxDiffNumOnes = MaxNumWC - MinNumWC + min(MaxWCLoss + 1, MinNumWC) + 1  
          ActualDiffNumOnes = abs(Info1[3] - SchInfo2[ , 3])
          C = (ActualDiffNumOnes <= MaxDiffNumOnes)
          
          SchSet2Idx = which( !SubSchInd2 & !A & !B & C )  # See Example 1 for an instance where !SubSchInd2 is put to use
          
        } else SchSet2Idx = c()
        
        if (length(SchSet2Idx) > 0) {
          
          for (j in 1 : length(SchSet2Idx)) {  # Sub loop
            
            Sch2Row = SchSet2Idx[j]
            
            sch2 = SchSet2[Sch2Row, ]
            
            UnionRes = UnionSchemata(sch1, sch2, MaxWCLossProp, MaxWCLoss)
            
            UniSch = UnionRes$UniSchema
            
            if ( !is.null(UniSch) ) {
              
              if (UnionRes$SubsumedSch == 3) {  # both sch1 and sch2 are subsumed
                
                SubSchInd1[i] = T
                SubSchInd2[Sch2Row] = T
                
                # Cleanup
                
                newSchNum = length(newSet2IdxGenByThisSch)
                
                if (newSchNum > 0) {
                  
                  newSchSet2 = newSchSet2[ -newSet2IdxGenByThisSch ]
                  
                  newSchSet2Idx = newSchSet2Idx - newSchNum
                  
                }
                
                if (length(schSet2SubIndSetBythisSch) > 0) SubSchInd2[schSet2SubIndSetBythisSch] = F
                
                newSchSet2[[newSchSet2Idx]] = UniSch
                
                newSchSet2Idx = newSchSet2Idx + 1
                
                break
                
              } else if (UnionRes$SubsumedSch == 1) {  # sch1 is subsumed
                
                SubSchInd1[i] = T
                
                # Cleanup
                
                newSchNum = length(newSet2IdxGenByThisSch)
                
                if (newSchNum > 0) {
                  
                  newSchSet2 = newSchSet2[ -newSet2IdxGenByThisSch ]
                  
                  newSchSet2Idx = newSchSet2Idx - newSchNum
                  
                }
                
                if (length(schSet2SubIndSetBythisSch) > 0) SubSchInd2[schSet2SubIndSetBythisSch] = F
                
                newSchSet2[[newSchSet2Idx]] = UniSch
                
                newSchSet2Idx = newSchSet2Idx + 1
                
                break
                
              } else if (UnionRes$SubsumedSch == 2) {  # sch2 is subsumed
                
                SubSchInd2[Sch2Row] = T
                
                newSchSet2[[newSchSet2Idx]] = UniSch
                
                newSet2IdxGenByThisSch = c(newSet2IdxGenByThisSch, newSchSet2Idx)
                
                newSchSet2Idx = newSchSet2Idx + 1
                
                schSet2SubIndSetBythisSch = c(schSet2SubIndSetBythisSch, Sch2Row)
                
              } else if (UnionRes$SubsumedSch == 0) {  # neither sch1 nor sch2 is subsumed
                
                newSchSet2[[newSchSet2Idx]] = UniSch
                
                newSet2IdxGenByThisSch = c(newSet2IdxGenByThisSch, newSchSet2Idx)
                
                newSchSet2Idx = newSchSet2Idx + 1
                
              }
              
            } else {  # UniSch is NULL
              
              if (UnionRes$DuplSch) {
                
                # (0##, 0##) are duplicates; so is 01# in the pair (0##, 01#).
                # The rationale is that no schema in Set1 will combine with a schema
                # marked 'duplicate' in Set2, since it has already maximally combined 
                # with the "master copy" in Set1.
                
                if (UnionRes$SubsumedSch == -1 | UnionRes$SubsumedSch == 1) {
                  
                  # there is a larger subsuming schema or a duplicate in SchSet2, so sch1 need not combine further
                  
                  SubSchInd1[i] = T  
                  
                  # Cleanup
                  
                  newSchNum = length(newSet2IdxGenByThisSch)
                  
                  if (newSchNum > 0) {
                    
                    newSchSet2 = newSchSet2[ -newSet2IdxGenByThisSch ]
                    
                    newSchSet2Idx = newSchSet2Idx - newSchNum
                    
                  }
                  
                  if (length(schSet2SubIndSetBythisSch) > 0) SubSchInd2[schSet2SubIndSetBythisSch] = F
                  
                  break
                  
                } else if (UnionRes$SubsumedSch == 2) {
                  
                  SubSchInd2[Sch2Row] = T
                  
                }
                
              }
              
            }  # end of processing a subsumed schema
            
          }  # end of 'for' over SchSet2
          
        }  # end of 'if (length(SchSet2Idx) > 0)'
        
        if (SubSchInd1[i] == T) next  # if True, the Sub loop would have broken and reached here
        
        ## Now combine Set1 with itself
        
        if (NewSet & (i < NumSet1Sch)) {
          
          # Prune search
          
          SchSubInfo1 = matrix(SchInfo1[-(1 : i), ], ncol = 4)
          
          MinMaxNumWC = sapply(SchSubInfo1[ , 4], function(val) c(min(Info1[4], val), max(Info1[4], val)))
          MinNumWC = MinMaxNumWC[1, ]
          MaxNumWC = MinMaxNumWC[2, ]
          MaxWCLossP = MaxWCLossProp * MinNumWC
          MaxWCLoss = sapply(MaxWCLossP, max, MaxWCLoss)
          A = Info1[1] > SchSubInfo1[ , 2]
          B = (SchSubInfo1[ , 1] > Info1[2])
          MaxDiffNumOnes = MaxNumWC - MinNumWC + min(MaxWCLoss + 1, MinNumWC) + 1  
          ActualDiffNumOnes = abs(Info1[3] - SchSubInfo1[ , 3])
          C = (ActualDiffNumOnes <= MaxDiffNumOnes)
          
          # Note the '+ i' at the end of the following line; works even if 'which' returns an empty set
          
          SchSubSet1Idx = which( !SubSchInd1[-(1 : i)] & !A & !B & C ) + i  # Note the +i at the end
          
          if (length(SchSubSet1Idx) > 0) {  # Sub loop 2
            
            for (k in 1 : length(SchSubSet1Idx)) {
              
              SchSub1Row = SchSubSet1Idx[k]
              
              schSub1 = SchSet1[SchSub1Row, ]
              
              UnionRes = UnionSchemata(sch1, schSub1, MaxWCLossProp, MaxWCLoss)
              
              UniSch = UnionRes$UniSchema
              
              if ( !is.null(UniSch) ) {
                
                if (UnionRes$SubsumedSch == 3) {  # both sch1 and sch2 are subsumed
                  
                  SubSchInd1[i] = T
                  SubSchInd1[SchSub1Row] = T
                  
                  # Cleanup
                  
                  newSchNum = length(newSet2IdxGenByThisSch)
                  
                  if (newSchNum > 0) {
                    
                    newSchSet2 = newSchSet2[ -newSet2IdxGenByThisSch ]
                    
                    newSchSet2Idx = newSchSet2Idx - newSchNum
                    
                  }
                  
                  if (length(schSet2SubIndSetBythisSch) > 0) SubSchInd2[schSet2SubIndSetBythisSch] = F
                  
                  if (length(schSet1SubIndSetBythisSch) > 0) SubSchInd1[schSet1SubIndSetBythisSch] = F
                  
                  newSchSet2[[newSchSet2Idx]] = UniSch
                  
                  newSchSet2Idx = newSchSet2Idx + 1
                  
                  break
                  
                } else if (UnionRes$SubsumedSch == 1) {  # sch1 is subsumed
                  
                  SubSchInd1[i] = T
                  
                  # Cleanup
                  
                  newSchNum = length(newSet2IdxGenByThisSch)
                  
                  if (newSchNum > 0) {
                    
                    newSchSet2 = newSchSet2[ -newSet2IdxGenByThisSch ]
                    
                    newSchSet2Idx = newSchSet2Idx - newSchNum
                    
                  }
                  
                  if (length(schSet2SubIndSetBythisSch) > 0) SubSchInd2[schSet2SubIndSetBythisSch] = F
                  
                  if (length(schSet1SubIndSetBythisSch) > 0) SubSchInd1[schSet1SubIndSetBythisSch] = F
                  
                  newSchSet2[[newSchSet2Idx]] = UniSch
                  
                  newSchSet2Idx = newSchSet2Idx + 1
                  
                  break
                  
                } else if (UnionRes$SubsumedSch == 2) {  # sch2 is subsumed
                  
                  SubSchInd1[SchSub1Row] = T
                  
                  newSchSet2[[newSchSet2Idx]] = UniSch
                  
                  newSet2IdxGenByThisSch = c(newSet2IdxGenByThisSch, newSchSet2Idx)
                  
                  newSchSet2Idx = newSchSet2Idx + 1
                  
                  schSet1SubIndSetBythisSch = c(schSet1SubIndSetBythisSch, SchSub1Row)
                  
                } else if (UnionRes$SubsumedSch == 0) {  # neither sch1 nor sch2 is subsumed
                  
                  newSchSet2[[newSchSet2Idx]] = UniSch
                  
                  newSet2IdxGenByThisSch = c(newSet2IdxGenByThisSch, newSchSet2Idx)
                  
                  newSchSet2Idx = newSchSet2Idx + 1
                  
                }
                
              } else {  # UniSch is NULL
                
                if (UnionRes$DuplSch) {
                  
                  # (0##, 0##) are duplicates; so is 01# in the pair (0##, 01#).
                  # The rationale is that no schema in set 1 will combine with a schema
                  # marked 'duplicate' in set 2, since it has already maximally combined 
                  # with the "master copy" in set 1.
                  
                  if (UnionRes$SubsumedSch == -1 | UnionRes$SubsumedSch == 1) {
                    
                    # there is a larger subsuming schema or a duplicate in SchSubSet1, so sch1 need not combine further
                    
                    SubSchInd1[i] = T
                    
                    # Cleanup
                    
                    newSchNum = length(newSet2IdxGenByThisSch)
                    
                    if (newSchNum > 0) {
                      
                      newSchSet2 = newSchSet2[ -newSet2IdxGenByThisSch ]
                      
                      newSchSet2Idx = newSchSet2Idx - newSchNum
                      
                    }
                    
                    if (length(schSet2SubIndSetBythisSch) > 0) SubSchInd2[schSet2SubIndSetBythisSch] = F
                    
                    if (length(schSet1SubIndSetBythisSch) > 0) SubSchInd1[schSet1SubIndSetBythisSch] = F
                    
                    break
                    
                  } else if (UnionRes$SubsumedSch == 2) {
                    
                    SubSchInd1[SchSub1Row] = T
                    
                  }
                  
                }
                
              }  # end of processing a subsumed schema
              
            }  # end of 'for' over SchSubSet1
            
          }  # end of 'if (length(SchSubSet1Idx) > 0)'
          
        } # end of 'if(NewSet)'
        
      }  # end of 'for' over SchSet1 (Main loop)
      
      # Form new SchSet2
      
      if (length(newSchSet2) > 0) {
        
        NewSchSet2 = matrix(unlist(newSchSet2), ncol = N, byrow = T)
        
      } else NewSchSet2 = NULL
      
      # Form new SchSet1
      
      NumSubSet1 = sum( !SubSchInd1)
      if ( !is.null(SchSet2) ) NumSubSet2 = sum( !SubSchInd2) else NumSubSet2 = 0
      NumSubSet = NumSubSet1 + NumSubSet2
      
      if (NumSubSet > 0) {
        
        NewSchSet1 = matrix(-1, nrow = NumSubSet, ncol = N)
        
        if (NumSubSet1 > 0) {
          
          NewSchSet1[1 : NumSubSet1, ] = matrix(SchSet1[ !SubSchInd1, ], ncol = N)
          
        }
        
        if (NumSubSet2 > 0) {
          
          NewSchSet1[(1 : NumSubSet2) + NumSubSet1, ] = matrix(SchSet2[ !SubSchInd2, ], ncol = N)
          
        }
        
        # Note that the assignment is reversed, which is what's intended.
        
        SchSet2 = NewSchSet1
        
        SchSet1 = NewSchSet2
        
        NewSet = T
        
      } else {
        
        SchSet1 = NewSchSet2
        SchSet2 = NULL
        
        NewSet = T
        
      }
      
    }  # end of 'while !is.null(SchSet1)'
    
    SchemataSet = SchemataSet[-2]
    
    SchemataSet[[1]] = SchSet2
    
    SchemataSet = SchemataSet[order(sapply(SchemataSet, length))]
    
    NewSet = F
    
  }  # end of 'while length(SchemataSet) > 1' 
  
  return(SchemataSet[[1]])
  
}

# ------------------------------------------------------------------------------------------------------------------------------------ #

# Conjunct or intersect a pair of schemata:

# Attempt to join two schemata via conjunction. Returns NULL if they contradict each other.
# Example: (A AND C) AND (A AND B) = (A AND B AND C); (A') AND (A AND B) = NULL.

IntersectSchemata <- function(sch1, sch2) {
  
  N = length(sch1)
  
  IntSchema = rep(-1, N)
  
  for (i in 1 : N) {
    
    s1 = sch1[i]
    s2 = sch2[i]
    
    if (s1 == s2) IntSchema[i] = s1 
    else if (s1 == 2) IntSchema[i] = s2 
    else if (s2 == 2) IntSchema[i] = s1
    else  {IntSchema = NULL; break}
    
  }
  
  return(IntSchema)
  
}

# ------------------------------------------------------------------------------------------------------------------------------------ #

# Identify superschema(ta) among a pair of schemata:

# In a pair of schemata, identify, if any, the schema that is the same as the result of the conjunction of the pair.
# That is, If Sch1 AND Sch2 => Sch2, then we say that Sch2 is a "super schema".
# Examples: (A) AND (A AND B) => (A AND B); ### AND 1## => 1##.
# It's possible that neither of the two schemata is a super of the other.
# Examples: 0# AND #0, 0# AND 1#, 0## AND #01 are pairs in which the schemata don't super each other.

# Outline of the method:
# When the two schemata are the same, they are both super schemata. 
# The function works by comparing symbols in the same position of the two schemata. If one is a literal and the other '#'
# then the one with the literal is the current superschema. If further symbol comparisons contradict it, then the procedure
# is terminated and it is concluded that neither schema is a superschema; otherwise, if the current superschema maintains its
# status throughout the procedure, it is declared as the superschema.

# Output: 
# a scalar 'super' with one of the following possible values:
# super = 3 => sch1 and sch2 are the same, that is they super each other.
# super = 1 => sch1 is super; super =2 => sch2 is super.
# super = 0 => neither sch1 nor sch2 supers the other.

# How does identifying superschemata help?
# The basic idea is that they help avoid unnecessary union operations later. For example, consider intersecting the sets
# {AB} and {A, BC, CD, BD}; read AB as (A AND B), and {A, BC} as (A OR (B AND C)). So, the full intersection calculation
# woud produce the set {AB, ABC, ABCD, ABD}, by the distributive law of Boolean logic, which upon further compression (union 
# with itself) yields {AB}. But note that {AB} is a superschema of {A}, and once this is discovered, the calculations could stop 
# right there and isolate {AB}; letting it intersect with the other schemata in the other set is unnecessary since their union 
# will simply result in {AB} itself. This way, identifying and isolating superschemata first before letting the rest of the sets 
# intersect and then compress could help circumvent many unnecessary union operations.

SuperSchema <- function(sch1, sch2) {
  
  super = 3
  
  N = length(sch1)
  
  for (i in 1 : N) {
    
    s1 = sch1[i]
    s2 = sch2[i]
    
    if (s1 != s2) {
      
      if (s1 == 2) {
        
        if (super == 3 | super == 2) super = 2 else {super = 0; break}
        
      } else if (s2 == 2) {
        
        if (super == 3 | super == 1) super = 1 else {super = 0; break}
        
      } else {super = 0; break}
      
    }
    
  }
  
  return(super)
  
}

# ------------------------------------------------------------------------------------------------------------------------------------ #

# Conjunct or intersect a set of schemata matrices:

# Outline of the method:
# Intersect (combine using the logical conjunction operation) a set of sets of schemata and compress the result. 
# The function proceeds by intersecting the first two sets in 'SchemataSet' into one set which then intersects likewise 
# with the other sets. At every step, during the intersection of SchSet1 and SchSet2, first the superschemata in both sets
# are identified and isolated; the rest of the schemata, if any, in both sets are then actually intersected. The final set
# of schemata so obtained is finally compressed by the 'UnionSchemataSet' function. This last compression step is necessary
# sometimes. For example, {##1} Intersect {11#, #11} -> {111, #11} -> further compression -> {#11}; here compression is 
# necessary. On the other hand, {##1} Intersect {1##, #11} -> {1#1, #11}; here further compression is not necessary. We 
# simply send the intersection result to the union method and let it perform any further compression needed. 
# If at any point SchSet1 and SchSet2 intersect to return a NULL, the procedure is halted and a NULL is returned. The reason
# is that a NULL result from an intersection means a logical contradiction (FALSE), and intersecting FALSE with other schemata
# is FALSE as well.

# Pruning conditions:
# We prune the search for superschemata by considering only those schemata in SchSet2, for every schema sch1 in SchSet1, 
# whose literals (enputs) are positioned within the start and end positions of the literals in sch1. 
# Condition (!A & !B & !C & !D) implements that constraint.
# Unlike schemata union, we have not implemented pruning for schemata intersection. Although, it is possible to conceive one.
# The only schemata-pairs we possibly wouldn't want to combine are those that logically contradict each other; thus a necessary
# and sufficient condition for that is for the schemata to have a 0 and 1 each at the same position. Thus, one possible pruning 
# condition would take into account the difference in NumWCs and NumOnes, just like condition C in the union function. For example,
# if k = 5 (length of a schema), sch1 has 3 WCs, sch2 has 1 WC, and the difference in NumOnes is equal to 4 (implying that sch2
# has 4 ones and sch1 has 2 zeros) then no matter how the symbols are arranged in the schemata, they would inevitably lead to a
# contradiction. A generalization of this example leads to the pruning condition: MaxDiffNumOnes = MaxNumWC.

# Optional parameters:
# The parameters 'MaxWCLossProp', 'MaxWCLoss', 'InitNewSet', 'SetMaxPI', 'MaxNumPIs' and 'RndPow' are passed to the union
# function and not made use of directly in the intersection step.

# Example:
# {1####, #1#1#, ##1##} AND {11###, 1##1#, #1#1#, ####1}
# Sets to intersect: SchSet1 = {1####, #1#1#, ##1##}; SchSet2 = {11###, 1##1#, #1#1#, ####1}
# SuperSchemata in SchSet1 = {#1#1#}; and in SchSet2 = {11###, 1##1#}
# Note that the schema #1#1# in SchSet2 is both a superschema and a duplicate.
# CombinedSchemata (isolate superschemata, ignore duplicates): {#1#1#, 11###, 1##1#}
# New sets to intersect: SchSet1 = {1####, ##1##}; SchSet2 = {####1}
# Intersection: {1###1, ##1#1}
# New CombinedSchemata: {#1#1#, 11###, 1##1#, 1###1, ##1#1}
# Compress/Self-union CombinedSchemata: {#1#1#, 11###, 1##1#, 1###1, ##1#1}
# Return: {#1#1#, 11###, 1##1#, 1###1, ##1#1}

IntersectSchemataSet <- function(SchemataSet, MaxWCLossProp = 0.0, MaxWCLoss = 1, 
                                 SetMaxPI = T, MaxNumPIs = 1000, RndPow = 0) {
  
  SchemataSet = SchemataSet[order(sapply(SchemataSet, length))]
  
  N = ncol(SchemataSet[[1]])
  
  while (length(SchemataSet) > 1) {
    
    SchSet1 = SchemataSet[[1]]
    SchSet2 = SchemataSet[[2]]
    
    # Identify superschemata in both SchSet1 and SchSet2
    
    # Collect information on the sets to be then used for pruning the search for superschemata
    
    SchMinMax1 = t(apply(SchSet1, 1, function(v) {idx = which(v < 2)
    if (length(idx) > 0) c(min(idx), max(idx))
    else c(-Inf, Inf)  # if the schema is ####...
    }))
    SchMinMax2 = t(apply(SchSet2, 1, function(v) {idx = which(v < 2)
    if (length(idx) > 0) c(min(idx), max(idx))
    else c(-Inf, Inf)  # if the schema is ####...
    }))
    
    SuperSchInd1 = rep(F, nrow(SchSet1))
    SuperSchInd2 = rep(F, nrow(SchSet2))
    
    DuplSchInd2 = rep(F, nrow(SchSet2))
    
    for (i in 1 : nrow(SchSet1)) {
      
      sch1 = SchSet1[i, ]
      
      MinMax1 = SchMinMax1[i, ]
      
      # Prune search
      
      A = MinMax1[1] > SchMinMax2[ , 2]
      B = (MinMax1[1] > SchMinMax2[ , 1]) & (MinMax1[2] > SchMinMax2[ , 2])
      C = (SchMinMax2[ , 1] > MinMax1[1]) & (SchMinMax2[ , 2] > MinMax1[2])
      D = (SchMinMax2[ , 1] > MinMax1[2])
      
      SchSet2Idx = which( !SuperSchInd2 & !A & !B & !C & !D)
      
      if (length(SchSet2Idx) == 0) next
      
      for (j in 1 : length(SchSet2Idx)) {
        
        Sch2Row = SchSet2Idx[j]
        
        sch2 = SchSet2[Sch2Row, ]
        
        super = SuperSchema(sch1, sch2)
        
        if (super == 1) {SuperSchInd1[i] = T; break}  # sch1 is a superschema
        else if (super == 2) {SuperSchInd2[SchSet2Idx[j]] = T}  # Do not place break here, else other superschemata in SchSet2 may be missed
        else if (super == 3) {  # sch1 and sch2 are exactly the same, so they are both superschemata and sch2 is duplicate
          SuperSchInd1[i] = T
          SuperSchInd2[Sch2Row] = T
          DuplSchInd2[Sch2Row] = T
          break
        }
        
      }
      
    }
    
    # Isolate superschemata and place them in the matrix 'CombinedSchemata'
    
    NumSubSet1 = sum(SuperSchInd1)
    SupNotDuplSchInd2 = (SuperSchInd2 &  ! DuplSchInd2)  # Ignore superschemata that are duplicates
    NumSubSet2 = sum(SupNotDuplSchInd2)
    NumSubSet3A = sum(! SuperSchInd1)
    NumSubSet3B = sum(! SuperSchInd2)
    NumSubSet3 =  NumSubSet3A * NumSubSet3B
    
    MaxPsblSch = NumSubSet1 + NumSubSet2 + NumSubSet3
    
    CombinedSchemata = matrix(-1, nrow = MaxPsblSch, ncol = N)
    
    if (NumSubSet1 > 0) {
      
      CombinedSchemata[1 : NumSubSet1, ] = matrix(SchemataSet[[1]][SuperSchInd1, ], ncol = N)
      
    }
    
    if (NumSubSet2 > 0) {
      
      CombinedSchemata[(1 : NumSubSet2) + NumSubSet1, ] = matrix(SchemataSet[[2]][SupNotDuplSchInd2, ], ncol = N)
      
    }
    
    # Form new SchSet1 and SchSet2 to intersect
    
    if (NumSubSet3A > 0) {
      
      SchSet1 = matrix(SchemataSet[[1]][ !SuperSchInd1, ], ncol = N)
      
    } else SchSet1 = NULL
    
    if (NumSubSet3B > 0) {
      
      SchSet2 = matrix(SchemataSet[[2]][ !SuperSchInd2, ], ncol = N)
      
    } else SchSet2 = NULL
    
    # Intersect SchSet1 and SchSet2 and append results to 'CombinedSchemata'
    
    Offset = NumSubSet1 + NumSubSet2
    
    Subset3Idx = 1
    
    if (NumSubSet3A > 0 & NumSubSet3B > 0) {
      
      for (i in 1 : nrow(SchSet1)) {
        
        sch1 = SchSet1[i, ]
        
        for (j in 1 : nrow(SchSet2)) {
          
          sch2 = SchSet2[j, ]
          
          schU = IntersectSchemata(sch1, sch2)
          
          if (! is.null(schU)) {
            
            CombinedSchemata[Subset3Idx + Offset, ] = schU
            
            Subset3Idx = Subset3Idx + 1
            
          }
          
        }
        
      }
      
    }
    
    # Remove any unoccupied rows in 'CombinedSchemata'
    
    if ((Subset3Idx + Offset) <= MaxPsblSch) {
      
      CombinedSchemata = matrix(CombinedSchemata[-((Subset3Idx + Offset) : MaxPsblSch), ], ncol = N)
      
    }
    
    # Self-union 'CombinedSchemata'
    
    if (nrow(CombinedSchemata) > 0) {
      
      # Compress the combined schemata
      
      CompSchSet = list()
      CompSchSet[[1]] = NULL
      CompSchSet[[2]] = CombinedSchemata
      
      CombinedSchemata = UnionSchemataSet(1, CompSchSet, c(2), MaxWCLossProp, MaxWCLoss, 
                                          InitNewSet = T, SetMaxPI, MaxNumPIs, RndPow)
      
      SchemataSet = SchemataSet[-2]
      
      SchemataSet[[1]] = CombinedSchemata
      
      SchemataSet = SchemataSet[order(sapply(SchemataSet, length))]
      
    } else {  # the schema being expanded here is impossible; NULL AND any schemata => NULL (FALSE or Contradiction)
      
      SchemataSet = NULL
      
      break
      
    }
    
  }
  
  if (length(SchemataSet) > 0) return(SchemataSet[[1]]) else return(NULL)
  
}

# ------------------------------------------------------------------------------------------------------------------------------------ #

# Subtract a schema from another to leave a remainder set of schemata that contain no overlap with the schema being
# subtracted.

# Basic idea:
# Subtract 'Schema' from 'FromSchema' by keeping any intersection between the two with the former, and splitting
# the latter into one or more 'Remainder' schemata in the process. In geometric terms, every remainder schema is nothing
# but a part of the 'FromSchema' that lies parallel to 'Schema' along a particular axis or dimension (corresponding to a 
# particular variable). The "parallel" nature of the remainder schemata is what that underlies the fact that the remainder
# schemata have no intersection with 'Schema' (although the remainder schemata could contain overlaps among themselves).

# Outline of the method:
# Sequentially consider each position of the schemata, then compare and process the symbols, while generating the
# remainder schemata in the process. 
# Repeat the following steps for every position from 1 to length of the schemata, until all positions are considered:
# (1) Initiate a new 'RemainderSchema' with 'FromSchema'.
# (2) If one the pair of symbols is a literal (0 or 1) and the other a wildcard (WC), or both symbols are the same literal 
# or are both WC, then do the following:
# (2a) If the WC belongs to 'FromSchema', then subtract the literal from '#' using the following rules: '#' - '1' = '0';
# '#' - '0' = '1', and place the remainder symbol at the corresponding position in the 'RemainderSchema'. Append 
# 'RemainderSchema to the list 'RemainderSchemata'. Go to the next schema position and repeat Step (1) onwards.
# (2b) If the WC belongs to 'Schema', then subtract the '#' from the literal using the following rules: '1' - '#' = '1';
# '0' - '#' = '0', and place the remainder symbol at the corresponding position in the 'RemainderSchema'. Go to the next 
# schema position and repeat Step (2) onwards.
# (2c) If both symbols are the same literal or are both WC, then place the literal/WC at the corresponding position in the 
# 'RemainderSchema'. Go to the next schema position and repeat Step (2) onwards.
# Return the list 'RemainderSchemata' if it is not empty; otherwise, return NA (since 'Schema' subsumes 'FromSchema').
# Return NULL if there is at least one position where the schemata contain opposite literals, since their intersection 
# would be NULL.

# Conceptual examples:
# SubtractSchema(1#, 11) -> {10}
# SubtractSchema(###, 111) -> {0##, #0#, ##0}
# SubtractSchema(1##, #11) -> {10#, 1#0}
# SubtractSchema(1#1#, 1011) -> {111#, 1#10}
# SubtractSchema(0##, #10) -> {00#, 0#1}
# SubtractSchema(#10, 0##) -> {110}
# SubtractSchema(##0, #11) -> NULL  # NULL intersection
# SubtractSchema(00#, 11#) -> NULL  # NULL intersection
# SubtractSchema(11#, 1##) -> NA    # latter subsumes former; NA = negative remainder
# SubtractSchema(111, ###) -> NA    # latter subsumes former; NA = negative remainder

SubtractSchema <- function(FromSchema, Schema, K = length(Schema)) {
  
  # Identify any mismatching positions in the schemata: where one has a 0 and the other 1
  
  MisMtch = sapply(1 : K, function(i) (FromSchema[i] != Schema[i]) & !(FromSchema[i] == 2 | Schema[i] == 2))
  
  if (any(MisMtch)) return(NULL)  # nothing can be subtracted (NULL intersection)
  
  # Remainder schemata could be one or more in number
  
  RemainderSchemata = list()
  Idx = 1
  
  # Step (1): Initialize the remainder
  
  RemainderSchema = FromSchema
  
  # Sequentially consider every schema position
  
  for (i in 1 : K) {
    
    s1 = FromSchema[i]
    s2 = Schema[i]
    
    if (s1 != s2) {
      
      if (s1 == 2) {  # Step (2a)
        
        RemainderSchema[i] = (s2 + 1)  %% 2
        
        # A new remainder schema is formed, so add it to the list
        
        RemainderSchemata[[Idx]] = RemainderSchema
        
        Idx = Idx + 1
        
        # Re-initialize a new remainder
        
        RemainderSchema = FromSchema
        
      } else if (s2 == 2) {  # Step (2b)
        
        RemainderSchema[i] = s1
        
      }
      
    } else {  # Step (2c)
      
      RemainderSchema[i] = s1
      
    }
    
  }
  
  if (Idx == 1) RemainderSchemata = NA  # either duplicate or FromSchema is fully subsumed
  
  return(RemainderSchemata)
  
}



















