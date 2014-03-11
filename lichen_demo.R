## C Steenbock
## 11 March 2014

rm(list=ls())
library(dplyr)
library(popbio)

### Load and sort data --------------
dat <- read.csv("lichen_demo.csv")
dat.size <- as.data.frame(dat[, 5:6])
colnames(dat.size) <- c("sqrtszt"="t0", "sqrtszt1"="t1") # rename column header
size <- arrange(dat.size, t0) # arrange initial column in ascending order

### Functions -----------------------

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
### Estimate vital rates, populate matrix -----------------------------

bin.num <- seq(from = 4, to = 25, by = 1)
lambdas <- rep(NA, length(bin.num))
    
for(i in 1:length(bin.num)){
    
    vec.bin <- bin.breaks(user_spec = c(0, sqrt(0.05)), bins = bin.num[i])
    
    ## Initialize storage
    n.bin <- length(vec.bin)-1                  # this is a workaround
    n <- rep(NA, n.bin)                         # count of indvs per bin
    median <- rep(NA, n.bin)                    # median size per bin for F
    surv <- rep(NA, n.bin)                      # survivorship for each class
    grow <- matrix(NA, n.bin, n.bin)            # store growth probabilites for each class

    # bin, survival, growth
    for(i in 1:(length(vec.bin)-1)){
        bounds <- c(vec.bin[i], vec.bin[i+1])                           # set limits for subset according to bin breaks
        subset <- size[size[,1] > bounds[1] & size[,1] < bounds[2],]    # subset data according to bounds
        n[i] <- length(subset[,2])                                      # store number of inviduals in this bin for future ref
        median[i] <- median(subset[,1])
        surv[i] <- sum(subset[,2] != 0) / length(subset[,2])            # calculate survivorship for this class, note TRUE = 1, which is a possible return of the logical expression !=
        histo <- hist(subset$t1, breaks = vec.bin, plot = FALSE)        # store hist as object, to access counts per bin
        grow[,i] <- histo$counts/length(subset[,1])                     # $counts returns the number of individuals of a certain size class                                                                   
    }
    
    M <- matrix(NA, n.bin, n.bin)   # initiate projection matrix
    
    # populate projection martix, specific to lichen data set
    r <- 0.185                      # MLE for r from Shriver (2012)
    
    for(i in 1:length(surv)){
        M[,i] <- surv[i] * grow[,i]
        M[1,i] <- (r * median[i]) + (surv[i] * grow[1,i]) # how does growth fit into M[1,i]?
    }
    
    lambdas[i] <- lambda(M)

}

lambdas
plot(lambdas, type = 'l')



# each reproductive matrix element a[1,i] = si*ci*rc, with ci the circumference of the mid-class size for 
# class i, si the survival rate of class i, and rc is the number of recruits produced per cm length of 
# thallus circumference (estimated as a circle from thallus area)






###################
#### An attempt to parameterize the model by logistic regression, and class based vital rates
# library(ggplot2)
# 
# dat <- read.csv("lichen_demo.csv")
# t0 <- dat$sqrtszt
# t1 <- dat$sqrtszt1
# 
# par(mfrow = c(1,2)) # histogram of indv size at both times
# hist(t0)    # first census
# hist(t1)    # second census
# 
# surv <- rep(NA,length(t1)) # Initialize vector to hold 0/1's, survivorship
# # Use for loop to assign zeros and ones to surv vector based on 2nd census
# for(i in 1:length(t1)){
#     if (t1[i] != 0){
#         surv[i] <- 1
#     } else {
#         surv[i] <- 0
#         }
# }
# 
# size.surv <- as.data.frame(cbind(t0, surv)) # append size vector and surv vector
# 
# # plot surviorship v. size
# surv.plot <- ggplot(size.surv, aes(t0, surv)) + 
#     geom_point(position=position_jitter(w=0.03, h=0.03), alpha=0.2, size=5)
# surv.plot
# 
# # bin size data into classes
# size <- sort(dat$sqrtszt)
# bin.size <- split(size, ceiling(seq_along(size)/203))
# med.size <- lapply(bin.size, median)
# 
# # logistic regression model of surv as fx of size
# mod1 <- glm(surv~t0, family=binomial, data=size.surv)
# (mod1.plot <- surv.plot + geom_smooth(method="glm", family="binomial"))
# 
# # function takes median values from each size class and uses it as the ind. var
# # in the equation from the above fitted logisic regression model.  plogis is
# # necessary to transform output from log-odds to probabilites
# 
# mod1.fun <- function(x){
#     pr.surv <- plogis(mod1$coef[1]+mod1$coef[2]*x)
#     return(pr.surv)
# }
# # apply above function to median values; returns a list of probs for each class
# vital.surv <- lapply(med.size, mod1.fun)

################
################


