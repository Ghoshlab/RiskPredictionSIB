# simulated code
# for a reproducible example
# using the methods in the paper
#
# "A weighted sample framework 
# to incorporating external calculators 
# for risk prediction" by Ghosh and Sabel

# Last updated: 7/12/21

set.seed(101)
GenderMale <- rbinom(2000,1,0.55)
Age <- rnorm(2000,50,10)
BD <- exp(rnorm(2000,0.5,0.5))
MR <- exp(rnorm(2000,0.5,1.53))
prob <- exp(0.87*log(BD)-0.019*Age-0.021*GenderMale+0.13*log(MR))/
  (1+exp(0.87*log(BD)-0.019*Age-0.021*GenderMale+0.13*log(MR)))
SLN <- rbinom(2000,1,prob)
sln.new <- data.frame(cbind(Age,BD,MR,SLN))


# Kattan 
# V  = -3.94 - 0.00034*Age-4.934*max(Age-33,0)^3 + 1.06*max(Age-56,0)^3
# - 5.67e-06*max(Age-76,0)^3 - 0.502*I(Site == "Head" or "Neck") +
# 0.8875*BD - 0.0688*max(BD-0.9,0)^3 +
# 0.0855*max(BD-1.8,0)^3
# - 0.0167*max(BD-5.5,0)^3
# + 0.772I(Clark == 4)
# 0.53*I(Clark == 5) +
# 0.535*I(Ulc == T)


#Kattan
constraints = 
  c(mean(- 0.00034*sln.new$Age),mean(-4.934*pmax(sln.new$Age-33,0)^3), 
    mean(1.06*pmax(sln.new$Age-56,0)^3),mean(-5.67e-06*pmax(sln.new$Age-76,0)^3),
    mean(0.8875*sln.new$BD),mean (-0.0688*pmax(sln.new$BD-0.9,0)^3),
    mean(0.0855*pmax(sln.new$BD-1.8,0)^3),mean(- 0.0167*pmax(sln.new$BD-5.5,0)^3))

Xmat <- cbind(sln.new$Age,pmax(sln.new$Age-33,0)^3,pmax(sln.new$Age-56,0)^3,
              pmax(sln.new$Age-76,0)^3,sln.new$BD,pmax(sln.new$BD-0.9,0)^3,
              pmax(sln.new$BD-1.8,0)^3,pmax(sln.new$BD-5.5,0)^3)


# nonnegative ridge regression
library(nnls)
lambda <- 10
newXmat <- cbind(Xmat,sqrt(lambda)*diag(1,2000,2000))
newconstraints <- c(constraints,rep(0,2000))
nnls1 <- nnls(t(newXmat),newconstraints)
wts <- coef(nnls1)
wts2 <- wts/sum(wts)

# Only 1 observations with nonzero weights!

# Now try Tejera-Vaquerizo constraints

# V = 0.5I(B \in (1,2)) + 1.9I(B \in (2,4)) + 2.45I(B > 4) + 1.5Sat 

constraints = c(0.5*mean(I(sln.new$BD >= 1 &sln.new$BD <= 2 )),
                1.9*mean(I(sln.new$BD > 2 &sln.new$BD <= 4 )),
                2.45*mean(I(sln.new$BD > 4)))

Xmat <- cbind(I(sln.new$BD >= 1 &sln.new$BD <= 2 ),
              I(sln.new$BD > 2 &sln.new$BD <= 4 ),
              I(sln.new$BD > 4))

# nonnegative ridge regression
library(nnls)
lambda <- 10
newXmat <- cbind(Xmat,sqrt(lambda)*diag(1,2000,2000))
newconstraints <- c(constraints,rep(0,2000))
nnls1 <- nnls(t(newXmat),newconstraints)
wts <- coef(nnls1)
wts2 <- wts/sum(wts)
                             

# Fit weighted regression

lr1 = glm(factor(SLN)~log(BD)+log(Age)+log(MR)+GenderMale,data=sln.new,family=binomial,weights=wts2,na.action=na.omit)
                                                
# perturbation inference


B <- 1000
coefs <- matrix(0,B,4)
for (i in 1:B) {
  set.seed(67+i)
  perturb <- rexp(2000,1)
  wts.p <- wts2*perturb
  lr.temp <- glm(factor(SLN)~log(BD)+log(Age)+log(MR)+GenderMale,data=sln.new,family=binomial,na.action=na.omit,weights=wts.p)
  coefs[i,] <- coef(lr.temp)[-1]
  cat(i,"\n")
}

se.beta <- apply(coefs,2,sd)
                            
# Tree modelling
ynew <- I(wts2 > 1/2000)

library(rpart)
sln.new$ynew <- factor(ynew)
tree1 <- rpart(ynew~log(BD)+log(Age)+GenderMale+log(MR),data=sln.new, na.action=na.omit)

