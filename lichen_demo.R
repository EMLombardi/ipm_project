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
# logistic regression model of surv as fx of size
mod1 <- glm(surv~t0, family=binomial, data=size.surv)

