## C Steenbock
## 30 January 2014
## From data to a demographic model

library(ggplot2)

dat <- read.csv("lichen_demo.csv")
t0 <- dat$sqrtszt
t1 <- dat$sqrtszt1

par(mfrow = c(1,2)) # histogram of indv size at both times
hist(t0)    # first census
hist(t1)    # second census

surv <- rep(NA,length(t1)) # Initialize vector to hold 0/1's, survivorship
# Use for loop to assign zeros and ones to surv vector based on 2nd census
for(i in 1:length(t1)){
    if (t1[i] != 0){
        surv[i] <- 1
    } else {
        surv[i] <- 0
        }
}

size.surv <- as.data.frame(cbind(t0, surv)) # append size vector and surv vector

# plot surviorship v. size
surv.plot <- ggplot(size.surv, aes(t0, surv)) + 
    geom_point(position=position_jitter(w=0.03, h=0.03), alpha=0.2, size=5)
surv.plot

# bin size data into classes
size <- sort(dat$sqrtszt)
bin.size <- split(size, ceiling(seq_along(size)/136))
med.size <- lapply(bin.size, median)

# logistic regression model of surv as fx of size
mod1 <- glm(surv~t0, family=binomial, data=size.surv)
mod1.plot <- surv.plot + geom_smooth(method="glm", family="binomial")

# function takes median values from each size class and uses it as the ind. var
# in the equation from the above fitted logisic regression model.  plogis is
# necessary to transform output from log-odds to probabilites
mod1.fun <- function(x){
    pr.surv <- plogis(mod1$coef[1]+mod1$coef[2]*x)
    return(pr.surv)
}
# apply above function to median values; returns a list of probs for each class
surv.est <- lapply(med.size, mod1.fun)




