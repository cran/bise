.packageName <- "bise"

noMix <- function (obsday, count) {

aic.mix1 <- rep(NA,1)
aic.mix2 <- rep(NA,1)
aic.mix3 <- rep(NA,1)
aic.mix4 <- rep(NA,1)

# ********************************************************
#                  Function for 1 Gaussian
# ********************************************************
gauss.mix1 <- function (x) {

# definition of variables
eof <- length(obsday)
norm.a <- NULL
disc.prob.a <- NULL
part.pred.a <- NULL
pred.mix.dist <- NULL
pred.cumul.dist <- NULL
deviation <- NULL
sum.deviation <- 0
total.norm.a <- 0
total.obs <- 0
cumul <- NULL

for(i in 1:eof) {
    norm.a[i] <- dnorm(obsday[i],mean = x[1],sd = x[2])
    total.norm.a <- total.norm.a + norm.a[i]
    total.obs <- total.obs + count[i]
        if(i == 1) { cumul[i] <- count[i] }
        else {cumul[i] <- cumul[i-1] + count[i]}
}

for(j in 1:eof) {
    disc.prob.a[j] <- norm.a[j] / total.norm.a
    part.pred.a[j] <- disc.prob.a[j] * x[3] * total.obs
    pred.mix.dist[j] <- part.pred.a[j]
        if(j == 1) { pred.cumul.dist[j] <- pred.mix.dist[j] }
        else { pred.cumul.dist[j] <- pred.cumul.dist[j-1] + pred.mix.dist[j] }
    deviation[j] <- (cumul[j] - pred.cumul.dist[j])^2
    sum.deviation <- sum.deviation + deviation[j]
}
return(sum.deviation)
}
optim.mod1 <- optim(c(mean(obsday),8,0.8), fn=gauss.mix1, gr=NULL, method="L-BFGS-B", lower=c(min(obsday),0.9999,0.0001), upper=c(max(obsday),20,0.9999), control=list(maxit=5000, pgtol=0.0001, ndeps=c(0.01,0.01,0.0000000001)))

# *****************************************************************************************************
# function "gauss.mix1" modified such that the return value is now pred.mix.dist (to be drawn as lines)
# *****************************************************************************************************
gauss.mix1 <- function (x) {

# definition of variables
eof <- length(obsday)
norm.a <- NULL
disc.prob.a <- NULL
part.pred.a <- NULL
pred.mix.dist <- NULL
pred.cumul.dist <- NULL
deviation <- NULL
sum.deviation <- 0
total.norm.a <- 0
total.obs <- 0
cumul <- NULL

for(i in 1:eof) {
    norm.a[i] <- dnorm(obsday[i],mean = x[1],sd = x[2])
    total.norm.a <- total.norm.a + norm.a[i]
    total.obs <- total.obs + count[i]
        if(i == 1) { cumul[i] <- count[i] }
        else {cumul[i] <- cumul[i-1] + count[i]}
}

for(j in 1:eof) {
    disc.prob.a[j] <- norm.a[j] / total.norm.a
    part.pred.a[j] <- disc.prob.a[j] * x[3] * total.obs
    pred.mix.dist[j] <- part.pred.a[j]
        if(j == 1) { pred.cumul.dist[j] <- pred.mix.dist[j] }
        else { pred.cumul.dist[j] <- pred.cumul.dist[j-1] + pred.mix.dist[j] }
    deviation[j] <- (cumul[j] - pred.cumul.dist[j])^2
    sum.deviation <- sum.deviation + deviation[j]
}
return(pred.mix.dist)
}


# *********************************************
#   MLE for gauss.mix1
# *********************************************
pred.mix.dist1 <- gauss.mix1(optim.mod1$par)
logL <- 0
var.par <- var(count- pred.mix.dist1)
for(i in 1:length(obsday)) {
    logL <- logL + dnorm(count[i], pred.mix.dist1[i], sqrt(var.par),TRUE)
}
aic.mix1 <- -2 * logL + 3
# **********************************************
#                 end model 1
# **********************************************


# ********************************************************
#                  Function for 2 Gaussians
# ********************************************************
gauss.mix2 <- function (x) {

# definition of variables
eof <- length(obsday)
norm.a <- NULL
norm.b <- NULL
disc.prob.a <- NULL
part.pred.a <- NULL
disc.prob.b <- NULL
part.pred.b <- NULL
pred.mix.dist <- NULL
pred.cumul.dist <- NULL
deviation <- NULL
sum.deviation <- 0
total.norm.a <- 0
total.norm.b <- 0
total.obs <- 0
cumul <- NULL

for(i in 1:eof) {
    norm.a[i] <- dnorm(obsday[i],mean = x[1],sd = x[3])
    norm.b[i] <- dnorm(obsday[i],mean = x[2],sd = x[4])
    total.norm.a <- total.norm.a + norm.a[i]
    total.norm.b <- total.norm.b + norm.b[i]
    total.obs <- total.obs + count[i]
        if(i == 1) { cumul[i] <- count[i] }
        else {cumul[i] <- cumul[i-1] + count[i]}
}

for(j in 1:eof) {
    disc.prob.a[j] <- norm.a[j] / total.norm.a
    part.pred.a[j] <- disc.prob.a[j] * x[5] * total.obs
    disc.prob.b[j] <- norm.b[j] / total.norm.b
    part.pred.b[j] <- disc.prob.b[j] * x[6] * total.obs
    pred.mix.dist[j] <- part.pred.a[j] +  part.pred.b[j]
        if(j == 1) { pred.cumul.dist[j] <- pred.mix.dist[j] }
        else { pred.cumul.dist[j] <- pred.cumul.dist[j-1] + pred.mix.dist[j] }
    deviation[j] <- (cumul[j] - pred.cumul.dist[j])^2
    sum.deviation <- sum.deviation + deviation[j]
}
return(sum.deviation)
}
optim.mod2 <- optim(c(quantile(obsday, probs=0.25, names=FALSE),quantile(obsday, probs=0.75, names=FALSE),15,15,0.9,0.1), fn=gauss.mix2, gr=NULL, method="L-BFGS-B", lower=c(min(obsday),min(obsday),2,2,0.0001,0.00001), upper=c(max(obsday),max(obsday),20,20,0.9999,0.9999), control=list(maxit=5000, pgtol=0.00001, ndeps=c(0.01,0.01,0.01,0.01,0.0000000001,0.0000000001), lmm=200))

# *****************************************************************************************************
#       function "gauss.mix2" modified such that the return value is now pred.mix.dist 
# *****************************************************************************************************
gauss.mix2 <- function (x) {

# definition of variables
eof <- length(obsday)
norm.a <- NULL
norm.b <- NULL
disc.prob.a <- NULL
part.pred.a <- NULL
disc.prob.b <- NULL
part.pred.b <- NULL
pred.mix.dist <- NULL
pred.cumul.dist <- NULL
deviation <- NULL
sum.deviation <- 0
total.norm.a <- 0
total.norm.b <- 0
total.obs <- 0
cumul <- NULL

for(i in 1:eof) {
    norm.a[i] <- dnorm(obsday[i],mean = x[1],sd = x[3])
    norm.b[i] <- dnorm(obsday[i],mean = x[2],sd = x[4])
    total.norm.a <- total.norm.a + norm.a[i]
    total.norm.b <- total.norm.b + norm.b[i]
    total.obs <- total.obs + count[i]
        if(i == 1) { cumul[i] <- count[i] }
        else {cumul[i] <- cumul[i-1] + count[i]}
}

for(j in 1:eof) {
    disc.prob.a[j] <- norm.a[j] / total.norm.a
    part.pred.a[j] <- disc.prob.a[j] * x[5] * total.obs
    disc.prob.b[j] <- norm.b[j] / total.norm.b
    part.pred.b[j] <- disc.prob.b[j] * x[6] * total.obs
    pred.mix.dist[j] <- part.pred.a[j] +  part.pred.b[j]
        if(j == 1) { pred.cumul.dist[j] <- pred.mix.dist[j] }
        else { pred.cumul.dist[j] <- pred.cumul.dist[j-1] + pred.mix.dist[j] }
    deviation[j] <- (cumul[j] - pred.cumul.dist[j])^2
    sum.deviation <- sum.deviation + deviation[j]
}
return(pred.mix.dist)
}

# *********************************************
#   MLE for gauss.mix2
# *********************************************
pred.mix.dist2 <- gauss.mix2(optim.mod2$par)
logL <- 0
var.par <- var(count- pred.mix.dist2)
for(i in 1:length(obsday)) {
    logL <- logL + dnorm(count[i], pred.mix.dist2[i], sqrt(var.par),TRUE)
}
aic.mix2 <- -2 * logL + 6
# **********************************************
#                 end model 2
# **********************************************


if(aic.mix1 < aic.mix2) { 
    print("just one distribution")
    #rm(pred.mix.dist2, optim.mod2, logL, var.par, aic.mix2)
    return(as.vector(pred.mix.dist1))
    }
else { rm(pred.mix.dist1, optim.mod1, logL, var.par, aic.mix1) }


# ********************************************************
#                  Function for 3 Gaussians
# ********************************************************
gauss.mix3 <- function (x) {

# definition of variables
eof <- length(obsday)
norm.a <- NULL
norm.b <- NULL
norm.c <- NULL
disc.prob.a <- NULL
part.pred.a <- NULL
disc.prob.b <- NULL
part.pred.b <- NULL
disc.prob.c <- NULL
part.pred.c <- NULL
pred.mix.dist <- NULL
pred.cumul.dist <- NULL
deviation <- NULL
sum.deviation <- 0
total.norm.a <- 0
total.norm.b <- 0
total.norm.c <- 0
total.obs <- 0
cumul <- NULL

for(i in 1:eof) {
    norm.a[i] <- dnorm(obsday[i],mean = x[1],sd = x[4])
    norm.b[i] <- dnorm(obsday[i],mean = x[2],sd = x[5])
    norm.c[i] <- dnorm(obsday[i],mean = x[3],sd = x[6])
    total.norm.a <- total.norm.a + norm.a[i]
    total.norm.b <- total.norm.b + norm.b[i]
    total.norm.c <- total.norm.c + norm.c[i]
    total.obs <- total.obs + count[i]
        if(i == 1) { cumul[i] <- count[i] }
        else {cumul[i] <- cumul[i-1] + count[i]}
}

for(j in 1:eof) {
    disc.prob.a[j] <- norm.a[j] / total.norm.a
    part.pred.a[j] <- disc.prob.a[j] * x[7] * total.obs
    disc.prob.b[j] <- norm.b[j] / total.norm.b
    part.pred.b[j] <- disc.prob.b[j] * x[8] * total.obs
    disc.prob.c[j] <- norm.c[j] / total.norm.c
    part.pred.c[j] <- disc.prob.c[j] * x[9] * total.obs
    pred.mix.dist[j] <- part.pred.a[j] +  part.pred.b[j] +  part.pred.c[j]
        if(j == 1) { pred.cumul.dist[j] <- pred.mix.dist[j] }
        else { pred.cumul.dist[j] <- pred.cumul.dist[j-1] + pred.mix.dist[j] }
    deviation[j] <- (cumul[j] - pred.cumul.dist[j])^2
    sum.deviation <- sum.deviation + deviation[j]
}
return(sum.deviation)
}
optim.mod3 <- optim(c(quantile(obsday, probs=0.25, names=FALSE),mean(obsday),quantile(obsday, probs=0.75, names=FALSE),5,8,4,0.3,0.3,0.4), fn=gauss.mix3, gr=NULL, method="L-BFGS-B", lower=c(min(obsday),min(obsday),min(obsday),2,2,2,0.0001,0.0001,0.0001), upper=c(max(obsday),max(obsday),max(obsday),20,20,20,0.9999,0.9999,0.9999), control=list(maxit=5000, pgtol=0.00001, ndeps=c(0.01,0.01,0.01,0.01,0.01,0.01,0.0000000001,0.0000000001,0.0000000001), lmm=200))

# *****************************************************************************************************
#       function "gauss.mix3" modified such that the return value is now pred.mix.dist 
# *****************************************************************************************************
gauss.mix3 <- function (x) {

# definition of variables
eof <- length(obsday)
norm.a <- NULL
norm.b <- NULL
norm.c <- NULL
disc.prob.a <- NULL
part.pred.a <- NULL
disc.prob.b <- NULL
part.pred.b <- NULL
disc.prob.c <- NULL
part.pred.c <- NULL
pred.mix.dist <- NULL
pred.cumul.dist <- NULL
deviation <- NULL
sum.deviation <- 0
total.norm.a <- 0
total.norm.b <- 0
total.norm.c <- 0
total.obs <- 0
cumul <- NULL

for(i in 1:eof) {
    norm.a[i] <- dnorm(obsday[i],mean = x[1],sd = x[4])
    norm.b[i] <- dnorm(obsday[i],mean = x[2],sd = x[5])
    norm.c[i] <- dnorm(obsday[i],mean = x[3],sd = x[6])
    total.norm.a <- total.norm.a + norm.a[i]
    total.norm.b <- total.norm.b + norm.b[i]
    total.norm.c <- total.norm.c + norm.c[i]
    total.obs <- total.obs + count[i]
        if(i == 1) { cumul[i] <- count[i] }
        else {cumul[i] <- cumul[i-1] + count[i]}
}

for(j in 1:eof) {
    disc.prob.a[j] <- norm.a[j] / total.norm.a
    part.pred.a[j] <- disc.prob.a[j] * x[7] * total.obs
    disc.prob.b[j] <- norm.b[j] / total.norm.b
    part.pred.b[j] <- disc.prob.b[j] * x[8] * total.obs
    disc.prob.c[j] <- norm.c[j] / total.norm.c
    part.pred.c[j] <- disc.prob.c[j] * x[9] * total.obs
    pred.mix.dist[j] <- part.pred.a[j] +  part.pred.b[j] +  part.pred.c[j]
        if(j == 1) { pred.cumul.dist[j] <- pred.mix.dist[j] }
        else { pred.cumul.dist[j] <- pred.cumul.dist[j-1] + pred.mix.dist[j] }
    deviation[j] <- (cumul[j] - pred.cumul.dist[j])^2
    sum.deviation <- sum.deviation + deviation[j]
}
return(pred.mix.dist)
}

# *********************************************
#   MLE for gauss.mix3
# *********************************************
pred.mix.dist3 <- gauss.mix3(optim.mod3$par)
logL <- 0
var.par <- var(count- pred.mix.dist3)
for(i in 1:length(obsday)) {
    logL <- logL + dnorm(count[i], pred.mix.dist3[i], sqrt(var.par),TRUE)
}
aic.mix3 <- -2 * logL + 9
# **********************************************
#                 end model 3
# **********************************************


if(aic.mix2 < aic.mix3) {
    print("mixture of two distributions")
    #rm(pred.mix.dist3, optim.mod3, logL, var.par, aic.mix3)
    return(as.vector(pred.mix.dist2))
    }
else { rm(pred.mix.dist2, optim.mod2, logL, var.par, aic.mix2) }


# ********************************************************
#                  Function for 4 Gaussians
# ********************************************************
gauss.mix4 <- function (x) {

# definition of variables
eof <- length(obsday)
norm.a <- NULL
norm.b <- NULL
norm.c <- NULL
norm.d <- NULL
disc.prob.a <- NULL
part.pred.a <- NULL
disc.prob.b <- NULL
part.pred.b <- NULL
disc.prob.c <- NULL
part.pred.c <- NULL
disc.prob.d <- NULL
part.pred.d <- NULL
pred.mix.dist <- NULL
pred.cumul.dist <- NULL
deviation <- NULL
sum.deviation <- 0
total.norm.a <- 0
total.norm.b <- 0
total.norm.c <- 0
total.norm.d <- 0
total.obs <- 0
cumul <- NULL

for(i in 1:eof) {
    norm.a[i] <- dnorm(obsday[i],mean = x[1],sd = x[5])
    norm.b[i] <- dnorm(obsday[i],mean = x[2],sd = x[6])
    norm.c[i] <- dnorm(obsday[i],mean = x[3],sd = x[7])
    norm.d[i] <- dnorm(obsday[i],mean = x[4],sd = x[8])
    total.norm.a <- total.norm.a + norm.a[i]
    total.norm.b <- total.norm.b + norm.b[i]
    total.norm.c <- total.norm.c + norm.c[i]
    total.norm.d <- total.norm.d + norm.d[i]
    total.obs <- total.obs + count[i]
        if(i == 1) { cumul[i] <- count[i] }
        else {cumul[i] <- cumul[i-1] + count[i]}
}

for(j in 1:eof) {
    disc.prob.a[j] <- norm.a[j] / total.norm.a
    part.pred.a[j] <- disc.prob.a[j] * x[9] * total.obs
    disc.prob.b[j] <- norm.b[j] / total.norm.b
    part.pred.b[j] <- disc.prob.b[j] * x[10] * total.obs
    disc.prob.c[j] <- norm.c[j] / total.norm.c
    part.pred.c[j] <- disc.prob.c[j] * x[11] * total.obs
    disc.prob.d[j] <- norm.d[j] / total.norm.d
    part.pred.d[j] <- disc.prob.d[j] * x[12] * total.obs
    pred.mix.dist[j] <- part.pred.a[j] +  part.pred.b[j] +  part.pred.c[j] + part.pred.d[j]
        if(j == 1) { pred.cumul.dist[j] <- pred.mix.dist[j] }
        else { pred.cumul.dist[j] <- pred.cumul.dist[j-1] + pred.mix.dist[j] }
    deviation[j] <- (cumul[j] - pred.cumul.dist[j])^2
    sum.deviation <- sum.deviation + deviation[j]
}
return(sum.deviation)
}
optim.mod4 <- optim(c(quantile(obsday, probs=0.2, names=FALSE),quantile(obsday, probs=0.4, names=FALSE),quantile(obsday, probs=0.6, names=FALSE),quantile(obsday, probs=0.8, names=FALSE),5,8,4,9,0.3,0.3,0.3,0.1), fn=gauss.mix4, gr=NULL, method="L-BFGS-B", lower=c(min(obsday),min(obsday),min(obsday),min(obsday),2,2,2,2,0.0001,0.0001,0.0001,0.0001), upper=c(max(obsday),max(obsday),max(obsday),max(obsday),20,20,20,20,0.9999,0.9999,0.9999,0.9999), control=list(maxit=5000, pgtol=0.000001, ndeps=c(0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.0000000001,0.0000000001,0.0000000001,0.0000000001), lmm=200))

# *****************************************************************************************************
# function "gauss.mix4" modified such that the return value is now pred.mix.dist (to be drawn as lines)
# *****************************************************************************************************
gauss.mix4 <- function (x) {

# definition of variables
eof <- length(obsday)
norm.a <- NULL
norm.b <- NULL
norm.c <- NULL
norm.d <- NULL
disc.prob.a <- NULL
part.pred.a <- NULL
disc.prob.b <- NULL
part.pred.b <- NULL
disc.prob.c <- NULL
part.pred.c <- NULL
disc.prob.d <- NULL
part.pred.d <- NULL
pred.mix.dist <- NULL
pred.cumul.dist <- NULL
deviation <- NULL
sum.deviation <- 0
total.norm.a <- 0
total.norm.b <- 0
total.norm.c <- 0
total.norm.d <- 0
total.obs <- 0
cumul <- NULL

for(i in 1:eof) {
    norm.a[i] <- dnorm(obsday[i],mean = x[1],sd = x[5])
    norm.b[i] <- dnorm(obsday[i],mean = x[2],sd = x[6])
    norm.c[i] <- dnorm(obsday[i],mean = x[3],sd = x[7])
    norm.d[i] <- dnorm(obsday[i],mean = x[4],sd = x[8])
    total.norm.a <- total.norm.a + norm.a[i]
    total.norm.b <- total.norm.b + norm.b[i]
    total.norm.c <- total.norm.c + norm.c[i]
    total.norm.d <- total.norm.d + norm.d[i]
    total.obs <- total.obs + count[i]
        if(i == 1) { cumul[i] <- count[i] }
        else {cumul[i] <- cumul[i-1] + count[i]}
}

for(j in 1:eof) {
    disc.prob.a[j] <- norm.a[j] / total.norm.a
    part.pred.a[j] <- disc.prob.a[j] * x[9] * total.obs
    disc.prob.b[j] <- norm.b[j] / total.norm.b
    part.pred.b[j] <- disc.prob.b[j] * x[10] * total.obs
    disc.prob.c[j] <- norm.c[j] / total.norm.c
    part.pred.c[j] <- disc.prob.c[j] * x[11] * total.obs
    disc.prob.d[j] <- norm.d[j] / total.norm.d
    part.pred.d[j] <- disc.prob.d[j] * x[12] * total.obs
    pred.mix.dist[j] <- part.pred.a[j] +  part.pred.b[j] +  part.pred.c[j] + part.pred.d[j]
        if(j == 1) { pred.cumul.dist[j] <- pred.mix.dist[j] }
        else { pred.cumul.dist[j] <- pred.cumul.dist[j-1] + pred.mix.dist[j] }
    deviation[j] <- (cumul[j] - pred.cumul.dist[j])^2
    sum.deviation <- sum.deviation + deviation[j]
}
return(pred.mix.dist)
}

# *********************************************
#   MLE for gauss.mix4
# *********************************************
pred.mix.dist4 <- gauss.mix4(optim.mod4$par)
logL <- 0
var.par <- var(count- pred.mix.dist4)
for(i in 1:length(obsday)) {
    logL <- logL + dnorm(count[i], pred.mix.dist4[i], sqrt(var.par),TRUE)
}
aic.mix4 <- -2 * logL + 12

# **********************************************
# **********************************************
#                 end model 4
# **********************************************
# **********************************************

if(aic.mix3 < aic.mix4) {
    print("mixture of 3 distributions")
    r#m(pred.mix.dist4, optim.mod4, logL, var.par, aic.mix4)
    return(as.vector(pred.mix.dist3))
    }
else { 
    print("mixture of 4 distributions")
    #rm(pred.mix.dist3, optim.mod3, logL, var.par, aic.mix3)
    return(as.vector(pred.mix.dist4))
    }
}
