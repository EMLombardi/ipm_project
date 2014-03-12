### christopher.steenbock@colorado.edu
### 12 March 2014

## bin.breaks() will create bin breaks/boundaries for an arbitrary number of bins that are
## *evenly* spaced apart based upon min/max values of the state variable. The function will
## also accept user specified break points and will continue populating the bin break vector
## as described for the remaining unspecified bins. User specified break points must be in
## vector form and start from zero and ascend continuously.

bin.breaks <- function(user_spec = NULL, bins){
    if(is.null(user_spec)){                                     # if/else used to determine is user_spec'd supplied
        n.bin <- bins
        vec.bin <- rep(NA, n.bin)                               # init. container for bin breaks
        vec.bin[1] <- 0 
        vec.bin[n.bin+1] <- max(size[,1])                       # the last bin needs an upper limit specified
        bin.width <- (max(size[,1]) - min(size[,1])) / n.bin    # find bin width
        # fill vector with evenly spaced bin breaks
        for(i in 1:(n.bin-1)){
            vec.bin[i+1] <- vec.bin[i] + bin.width              # assign values to bin break vector
        }
        return(vec.bin) 
    } else {
        n.bin <- bins
        vec.bin <- rep(NA, n.bin)                               # init. container for bin breaks
        # assign user specified break points to vec.bin
        for(i in 1:length(user_spec)){
            vec.bin[i] <- user_spec[i]
        }
        vec.bin[n.bin+1] <- max(size[,1])                       # the last bin needs an upper limit specified
        bin.width <- (max(size[,1]) - min(size[,1])) / n.bin    # find bin width
        # fill vector with evenly spaced bin breaks; note that for loop begins from the last user spec'd point
        for(i in length(user_spec):(n.bin-1)){
            vec.bin[i+1] <- vec.bin[i] + bin.width              # assign values to bin break vector
        }
        return(vec.bin)    
    }
}