## C Steenbock
## 21 February 2014
## From data to a demographic model

rm(list=ls())
library(dplyr)

dat <- read.csv("lichen_demo.csv")
dat.size <- as.data.frame(dat[, 5:6])
colnames(dat.size) <- c("sqrtszt"="t0", "sqrtszt1"="t1") # rename column header
size <- arrange(dat.size, t0) # arrange initial column in ascending order

n.bin <- 5  # specify number of classes

# initialize storage
n <- rep(NA, n.bin)                         # count of indvs per bin
median <- rep(NA, n.bin)                    # median size per bin for F
surv <- rep(NA, n.bin)
grow <- matrix(NA, n.bin, n.bin)            # store growth probabilites for each size class

vec.bin <- rep(NA, n.bin)                   # init. container for bin breaks
vec.bin[1] <- 0 
vec.bin[n.bin+1] <- max(size[,1])           # the last bin needs an upper limit specified

# fill vector with evenly spaced bin breaks
for(i in 1:(n.bin-1)){
    bin.width <- (max(size[,1]) - min(size[,1])) / n.bin    # find bin width
    vec.bin[i+1] <- vec.bin[i] + bin.width                  # assign values to bin break vector
}

# bin, survival, growth
for(i in 1:(length(vec.bin)-1)){
    bounds <- c(vec.bin[i], vec.bin[i+1])                           # set limits for subset according to bin breaks
    subset <- size[size[,1] > bounds[1] & size[,1] < bounds[2],]    # subset data according to bounds
    n[i] <- length(subset[,2])                                      # store number of inviduals in this bin for future ref
    median[i] <- median(subset[,1])
    surv[i] <- sum(subset[,2] != 0) / length(subset[,2])            # calculate survivorship for this class, note TRUE = 1, which is a possible return of the logical expression !=
    histo <- hist(subset$t1, breaks = vec.bin, plot = FALSE)        # store hist as object, to access counts per bin
    grow[,i] <- histo$counts/length(subset[,1])                     # $counts returns the number of individuals of a certain size class                                                                   # at t1 from the bin set at t0, that is growth or shrinkage or neither   
}

n
surv
grow

M <- matrix(NA, n.bin, n.bin)   # initiate projection matrix
r <- 0.185                      # MLE for r from Shriver (2012)

for(i in 1:length(surv)){
    M[,i] <- surv[i] * grow[,i]
    M[1,i] <- r * median[i] * surv[i] # how does growth fit into M[1,i]?
}

M

# each reproductive matrix element a[1,i] = si*ci*rc, with ci the circumference of the mid-class size for 
# class i, si the survival rate of class i, and rc is the number of recruits produced per cm length of 
# thallus circumference (estimated as a circle from thallus area)






# ###################
# #### An attempt to parameterize the model by logistic regression, and class based vital rates
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
# 
# ################
# ################


