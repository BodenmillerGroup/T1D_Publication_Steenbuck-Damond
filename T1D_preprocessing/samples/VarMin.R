# Function to minimize variance in intervention studies.
# Read paper at: https://link.springer.com/article/10.3758/s13423-021-01970-5

VarMin <- function(file, pRand = 0, save_output = F){
  
  # VarMin allocates a data point to a group, such that group 
  # differences are minimized on a number of specified variables.
  #
  # VarMin(file) uses the data specified in file to
  # match the current subject and returns the optimal group allocation and a
  # table contained the information in the CSV provided, plus the new
  # group allocation of the current subject.
  #
  # VARMIN(file, pRand = p) sets the probability of doing a random assignment, 
  # instead of VM, to p. The default value for pRand is 0, meaning that 
  # the algorithm always attempts to do VM, while keeping groups balanced. 
  # This is the recommended method, unless the experimenter has access to 
  # previous group allocations, in which case a small p_rand should be added to 
  # the algorithm, such that the current allocation cannot be inferred from
  # previous allocations (because the group balancing in VM makes every n_th 
  # participant predictable, [where n = number of groups], since participants 
  # are allocated to a shrinking set of potential groups each cycle).
  #
  # NOTE: if you set pRand > 0, the function automatically opts for random 
  # allocation of the first couple of participants, until each group has at
  # least one participant, instead of doing deterministic allocation of the
  # first participants to fill the groups as quickly as possible (since 
  # this is also predictable if previous allocations and the algorithm are known).
  #
  # VARMIN(file, SAVEFILE = T) if SAVEFILE is TRUE, the CSV in file
  # is updated and saved in the current working directory as a file called
  # 'VM_output.csv'. The file will not be saved if SAVEFILE is FALSE or not
  # specified.
  # NOTE: If you'd like to overwrite the current file, please rename your 
  # CSV as 'VM_output.csv' and make sure it is in the current working 
  # directory. In this case, please ensure that the CSV you are writing to 
  # is CLOSED, otherwise permission to overwrite the file may be denied.
  
  # Required packages
  require(pracma)
  require(utils)
  
  #Function to calculate sd of population
  sd.p=function(x){sd(x)*sqrt((length(x)-1)/length(x))}
  
  #Function to check for scalar
  is.scalar <- function(x) is.atomic(x) && length(x) == 1L
  
  # input data frame
  df <- read.csv(file = file, header = T)
  
  # removing empty columns
  df <- df[, !sapply(df, function(x)all(is.na(x))), drop=F]
  
  # Saving variables separately
  groupNum <- mean(df$NumberOfGroups, na.rm = T)
  previousAllocations <- cbind(df$Allocation[!is.na(df$Allocation)])
  varargin <- data.frame(df[, grepl( "Measure", names(df))])
  # drop NAs in varargin
  varargin <- as.matrix(varargin[rowSums(is.na(varargin)) != ncol(varargin), ])
  
  # for group number, needs to be scalar
  if(!is.scalar(groupNum)){
    stop('Specified group number needs to be scalar')
  }
  
  # if no samples have been assigned, i.e. previousAllocations is an empty
  # vector, do nothing
  # otherwise make sure that previous group allocations are a vector & numeric
  if(isempty(previousAllocations)){} else if(!is.numeric(previousAllocations) || !is.matrix(previousAllocations)){
    stop('Previous group allocations need to be numeric')
  } else {
    
    # check whether the  previous group allocations are valid 
    # (groups are numbered 1,2,3,..,number of groups):
    # 1) if they are more than the group number, input is invalid 
    # 2) if they are less than 1, input is invalid
    # 3) if they are not an integer number
    if(max(previousAllocations) > groupNum || min(previousAllocations) < 1 || sum(previousAllocations%%1)!=0){
      stop('Invalid group number in previous allocations')
    }
  }
  
  # invalid if no matching variables were provided
  if(is.null(varargin)){
    stop('Provide at least one variable on which to match')
  } else if(!is.matrix(varargin)){
    stop('All matching variables must be numeric arrays')
  } else {nRows <- nrow(varargin)
  nCols <- ncol(varargin)}
  
  # check that matching variables have more entries than previous
  # allocations (so that current samples can be allocated to a group)
  if(length(previousAllocations) >= nRows){
    stop('previous group allocations do not match the number of rows in matching variables')
  }

  # check that pRand is valid
  if(!is.scalar(pRand)){
    stop('pRand needs to be a scalar between 0 and 1')
  }
  else if (pRand < 0 || pRand > 1) {
      stop('pRand needs to be between 0 and 1')
  }
  
  # convert matching variables 
  # if the matching variables were given as separate inputs (rather than in a
  # matrix), create a matrix and standardize
  #str(scale(varargin))   
  matchingVars <- matrix(scale(varargin, center = TRUE, scale = apply(varargin, 2, sd.p)), nrow = nRows, ncol = nCols)
  
  #  allocate to group
  
  # number of samples that need to be allocated
  numAllocations = nrow(matchingVars) - nrow(previousAllocations)
  
  # list collecting group allocations
  groupAllocations = vector()
  
  for(x in 1:numAllocations) {
    
    matchingVarsCurrent = cbind(matchingVars[1:(nrow(matchingVars)-(numAllocations-x)),])
    
    # if pRand > 0, this indicates that previous allocations are known and we should randomly allocate if:
    # 1) groups don't have at least one participant
    # 2) the coin flip of pRand was successful
    
    if (pRand > 0 && (length(unique(previousAllocations)) < groupNum || runif(1) < pRand)){
      groupAllocation <- sample(1:groupNum, 1)
    }
    
    # if the number of groups to which samples were already assigned is less 
    # than the total number of groups, then  assign to one of the groups that 
    # does not have any samples in it yet
    else if(length(previousAllocations) < groupNum){
      
      # get empty groups
      emptyGroups <- setdiff(randperm(groupNum, groupNum), previousAllocations)
      
      # take first empty group as group for current participant
      groupAllocation <- emptyGroups[1]
    }   
    # otherwise, proceed with matching procedure
    else {
      # check which groups currently have the minimum number of samples, these
      # are the groups to which the current sample could be allocated
      # (i.e. possible options are:
      # (1) if every group has an equal number of samples, all groups will
      # be considered
      # (2) if some groups have an equal number samples, these groups will be
      # considered
      # (3) if only one group has the minimum amount, this is the group to
      # which the sample will be allocated)
      
      
      # sum(previousAllocations == previousAllocations') counts the
      # number of occurences of each group for the previous samples 
      numOfOccurences <- cbind(rep(as.numeric(table(previousAllocations)), length.out=length(previousAllocations)))
      
      # we know that the first n elements, where n = the number of groups, 
      # in previousAllocations are the number of occurences in groups 1 to n, 
      # because this is how we allocated in the if-statement above.
      # Therefore we can simply take the first n elements in numOfOccurences 
      # to get the number of samples in groups 1 to n
      numPerGroup <- numOfOccurences[1:groupNum]
      
      # we now find the groups that have the minimum number of samples and
      # can therefore be a potential group for the current sample
      potentialGroups <- which(numPerGroup == min(numPerGroup))
      
      # if there is only a single group that has the minimum sample size, we
      # allocate to that group directly, this is case (1) above
      if(length(potentialGroups)==1){
        groupAllocation <- potentialGroups}
      
      #% otherwise we need to find the best fit to handle cases (2) and (3)
      else{ 
        
        # expand previous allocations by 1 to make space for temporary
        # group allocation of current sample
        temporaryAllocations <- rbind(previousAllocations, NA)
        
        # create variable to store the sums of standard deviations for each
        # of the assignments to the potential groups
        tempFit <- matrix(NA, nrow = length(potentialGroups), ncol = 1)
        
        
        for(i in 1:length(potentialGroups)){
          #temporarily assign person to each potential group
          temporaryAllocations[length(temporaryAllocations)] <- potentialGroups[i]
          
          # this array will hold the mean z-score for each potential
          # group for each matching variable
          groupMeans <- matrix(NA, nrow = length(potentialGroups), ncol = ncol(matchingVarsCurrent))
          
          for(j in 1:length(potentialGroups)){
            # get mean of matching variables for the j'th group in 
            # 'potentialGroups', given this temporary assignment
            
            # these are the values of the matching variables for the
            # j'th group
            currentMatchingVars <- matrix(matchingVarsCurrent[which(temporaryAllocations == potentialGroups[j]),], ncol=ncol(matchingVarsCurrent))
            
            # enter mean z-score of j'th group given temp. assignment
            groupMeans[j,] <- colMeans(currentMatchingVars)
          }
          
          # groupMeans now contains the mean z score in each group 
          # (along the rows) for each matching variable (along columns),
          # i.e. each column contains the group means for one matching
          # variable.
          
          # now we calculate how 'bad' this assignment is by computing
          # the standard deviation between the groups for each matching
          # variable (i.e. of one column). We can then take the sum of the 
          # standard deviations along the rows as an
          # indication how well the groups are matched 
          # in terms of all matching variables taken together. 
          # The higher this sum, the worse is the fit 
          # given the current assignment to the i'th group in
          # potentialGroups.
          # tempFit stores this value for assignment to each potential
          # group.
          
          tempFit[i] <- sum(apply(X = groupMeans, FUN = sd.p, MARGIN = 2), na.rm = T) 
          
        }
        # find index for group that minimized the sum of standard
        # deivations
        bestGroupIdx <- which.min(tempFit)
        
        # index into potential groups to get final group allocation
        groupAllocation <- potentialGroups[bestGroupIdx]
        
        
      }} 
    
    # vector collecting group allocations
    groupAllocations <- cbind(groupAllocations, groupAllocation)
    
    # add current group to previous allocations
    previousAllocations <- rbind(previousAllocations, groupAllocation)
    
    # add to data frame
    df$Participant[which.max(is.na(df$Allocation))] <- length(previousAllocations)
    df$Allocation[which.max(is.na(df$Allocation))] <- groupAllocation
  }

  # remove leftover rows
  df <- df[1:length(previousAllocations),]
  
  if(save_output == T){
    write.csv(x = df, file = "VM_output.csv", row.names = FALSE)
  }
  return(list(allocations=groupAllocations, dataframe=df))
}