# Conditional distribution-based downweighting for robust estimation of logistic regression models

* * *

# Install pre-requisite packages

```{r,results='hide',message=FALSE,tidy=TRUE, tidy.opts=list(width.cutoff=60)}
Package.Names <- c("mvtnorm","pracma","robust","robustbase","ggplot2")
lapply(Package.Names, library, character.only = TRUE)
library(mvtnorm)
library(pracma)
library(robust)
library(robustbase)
library(ggplot2)
source("MallowsMEstRecode.R")
```

Initialize labels (balanced)
```{r}
n1 <- 100
n0 <- 100
y <- c(rep(1,n1),rep(0,n0))
```

Verify the true coefficient values
```{r}
n1big <- 100000
n0big <- 100000
ybig <- c(rep(1,n1big),rep(0,n0big))
set.seed(999999)
X1big <- matrix(rnorm( n=(2*n1big), mean=2,sd =1),nrow=n1big, ncol=2)
X0big <- matrix(rnorm( n=(2*n0big), mean=0,sd =1),nrow=n0big, ncol=2)  
TrueCoef <- coef(glm(ybig~rbind(X1big,X0big),family=binomial(link="logit")))
```

Repeated simulations (no outliers; balanced)
```{r, warnings=FALSE}
TotalReps <- 500
MLEmatrix <- matrix(0,nrow=TotalReps,ncol=3)
BYmatrix <- matrix(0,nrow=TotalReps,ncol=3)
ModMallowsmatrix <- matrix(0,nrow=TotalReps,ncol=3)
cubifmatrix <- matrix(0,nrow=TotalReps,ncol=3)
mallowsmatrix <- matrix(0,nrow=TotalReps,ncol=3)
ModMallowsCovEstArray <- array(0,dim=c(3,3,TotalReps))
CompTime.MLE <- rep(0,TotalReps)
CompTime.ModMallows <- rep(0,TotalReps)
CompTime.cubif <- rep(0,TotalReps)
CompTime.mallows <- rep(0,TotalReps)

rep.id <- 0
num.reps <- 0
rep.idStore <- rep(0,TotalReps)
while(num.reps < TotalReps){
  
  rep.id <- rep.id + 1
  set.seed(rep.id)
  X1 <- matrix(rnorm((n1*2),mean=2,sd=1),nrow=n1,ncol=2)
  X0 <- matrix(rnorm((n0*2),mean=0,sd=1),nrow=n0,ncol=2)
  X <- rbind(X1,X0)
  X1out <- X1
  Xout <- rbind(X1out,X0)
  tempGLm <- glm(y[-1]~Xout[-1,],family=binomial(link="logit"))
  tempcoef <- coef(tempGLm)
  
  if( (tempGLm$converged)&(sum((tempcoef - TrueCoef)^2) < 3) ){
    
    num.reps <- num.reps + 1
    start.MLE <- Sys.time()
    MLEmatrix[num.reps,] <- coef(glm(y~Xout,family=binomial(link="logit")))
    end.MLE <- Sys.time()
    CompTime.MLE[num.reps] <- as.numeric(difftime(end.MLE,start.MLE,units = "secs"))
    start.cubif <- Sys.time()
  cubifmatrix[num.reps,] <- coef(glmRob(y~Xout, family=binomial(link="logit"), method="cubif", bpar=0.824,cpar=2))
  end.cubif <- Sys.time()
  CompTime.cubif[num.reps] <- as.numeric(difftime(end.cubif,start.cubif,units = "secs"))
  start.mallows <- Sys.time()
  mallowsmatrix[num.reps,] <- coef(glmRob(y~Xout, family=binomial(link="logit"), method="mallows", wt.tuning=3.005))
  end.mallows <- Sys.time()
  CompTime.mallows[num.reps] <- as.numeric(difftime(end.mallows,start.mallows,units = "secs"))
  BYmatrix[num.reps,] <- coef(glmrob(y~Xout, family=binomial(link="logit"), method="BY", const =7, maxhalf=1))
  start.Modmallows <- Sys.time()
  ModMallows.init <- coef(glmRob(y~Xout, family=binomial(link="logit"), method="mallows"))
  MMobj <- CalcModMallows(y=y,X=Xout,theta.init = ModMallows.init, h=5.3)
  ModMallowsCovEstArray[,,num.reps] <- MMobj$SigmaEst
  ModMallowsmatrix[num.reps,] <- c(MMobj$theta.curr)
  end.Modmallows <- Sys.time()
  CompTime.ModMallows[num.reps] <- as.numeric(difftime(end.Modmallows,start.Modmallows,units = "secs"))
  rep.idStore[num.reps] <- rep.id
    
  }
  
}

MLEeff <- apply(MLEmatrix,1,function(s){ (s - TrueCoef)^2  })
BYeff <- apply(BYmatrix,1, function(s){ (s - TrueCoef)^2 })
ModMalloweff <- apply(ModMallowsmatrix,1,function(s){ (s - TrueCoef)^2  })
cubifeff <- apply(cubifmatrix,1,function(s){ (s - TrueCoef)^2  })
mallowseff <- apply(mallowsmatrix,1,function(s){ (s - TrueCoef)^2})

MLEeff.total <- colSums(MLEeff)
BYeff.total <- colSums(BYeff)
ModMalloweff.total <- colSums(ModMalloweff)
cubifeff.total <- colSums(cubifeff)
mallowseff.total <- colSums(mallowseff)
```

Average relative efficiency (no outliers; balanced)
```{r}
#MLE vs BY
mean(BYeff.total)/mean(ModMalloweff.total)
#MLE vs ModMallows
mean(MLEeff.total)/mean(ModMalloweff.total)
#MLE vs cubif
mean(MLEeff.total)/mean(cubifeff.total)
#MLE vs mallows
mean(MLEeff.total)/mean(mallowseff.total)
```

Assess accuracy of asymptotic SE of modified Mallows
```{r}
#True covariance matrix (based on 500 draws)
truth.Matrix <- var(ModMallowsmatrix)
ModMallows.Cov.Est.Final <- apply(ModMallowsCovEstArray,c(1,2),mean)
truth.Matrix
ModMallows.Cov.Est.Final
norm(ModMallows.Cov.Est.Final - truth.Matrix)
```


Repeated simulations (one outlier; balanced)
```{r}
TotalReps <- 500
position.out <- seq(-3.5,1,0.1)
J <- length(position.out)
MLEeff.total.Mat <- matrix(0,nrow=TotalReps,ncol=J)
BYeff.total.Mat <- matrix(0,nrow=TotalReps,ncol=J)
ModMalloweff.total.Mat <- matrix(0,nrow=TotalReps,ncol=J)
cubifeff.total.Mat <- matrix(0,nrow=TotalReps,ncol=J)
mallowseff.total.Mat <- matrix(0,nrow=TotalReps,ncol=J)
#Weight.Mat <- matrix(0,nrow=TotalReps,ncol=J)

for(j in 1:J){
  
  
MLEmatrix <- matrix(0,nrow=TotalReps,ncol=3)
BYmatrix <- matrix(0,nrow=TotalReps,ncol=3)
ModMallowsmatrix <- matrix(0,nrow=TotalReps,ncol=3)
cubifmatrix <- matrix(0,nrow=TotalReps,ncol=3)
mallowsmatrix <- matrix(0,nrow=TotalReps,ncol=3)

for(num.reps in 1:TotalReps){
  
  set.seed(rep.idStore[num.reps])
  X1 <- matrix(rnorm((n1*2),mean=2,sd=1),nrow=n1,ncol=2)
  X0 <- matrix(rnorm((n0*2),mean=0,sd=1),nrow=n0,ncol=2)
  X <- rbind(X1,X0)
  X1out <- X1
  X1out[1,] <- c(position.out[j],position.out[j])
  Xout <- rbind(X1out,X0)
  
    MLEmatrix[num.reps,] <- coef(glm(y~Xout,family=binomial(link="logit")))
    BYmatrix[num.reps,] <- coef(glmrob(y~Xout, family=binomial(link="logit"), method="BY", const =7, maxhalf=1))
    cubifmatrix[num.reps,] <- coef(glmRob(y~Xout, family=binomial(link="logit"), method="cubif", bpar=0.824,cpar=2))
  mallowsmatrix[num.reps,] <- coef(glmRob(y~Xout, family=binomial(link="logit"), method="mallows", wt.tuning=3.005))
  ModMallows.init <- coef(glmRob(y~Xout, family=binomial(link="logit"), method="mallows"))
  MMobj <- CalcModMallows(y=y,X=Xout,theta.init = ModMallows.init, h=5.3)
  ModMallowsmatrix[num.reps,] <- c(MMobj$theta.curr)
  #Weight.Mat[num.reps,j] <- MMobj$w[1]
}

MLEeff <- apply(MLEmatrix,1,function(s){ (s - TrueCoef)^2  })
BYeff <- apply(BYmatrix,1, function(s){ (s - TrueCoef)^2 })
ModMalloweff <- apply(ModMallowsmatrix,1,function(s){ (s - TrueCoef)^2  })
cubifeff <- apply(cubifmatrix,1,function(s){ (s - TrueCoef)^2  })
mallowseff <- apply(mallowsmatrix,1,function(s){ (s - TrueCoef)^2  })


MLEeff.total.Mat[,j] <- colSums(MLEeff)
BYeff.total.Mat[,j] <- colSums(BYeff)
ModMalloweff.total.Mat[,j] <- colSums(ModMalloweff)
cubifeff.total.Mat[,j] <- colSums(cubifeff)
mallowseff.total.Mat[,j] <-  colSums(mallowseff)

}
```

Plot Total MSE (averaged across 500 datasets; balanced and one outlier)
```{r}
TrimMSE.MLE.Total <- colMeans(MLEeff.total.Mat)
TrimMSE.BY.Total <- colMeans(BYeff.total.Mat)
TrimMSE.ModMallows.Total <- colMeans(ModMalloweff.total.Mat)
TrimMSE.Cubif.Total <- colMeans(cubifeff.total.Mat)
TrimMSE.Mallows.Total <- colMeans(mallowseff.total.Mat)

plot(TrimMSE.MLE.Total~position.out,pch=1,ylab="Total MSE",xlab="outlier position", main="Total MSE")
points(TrimMSE.ModMallows.Total~position.out,pch=2,col="brown")
points(TrimMSE.Cubif.Total~position.out,pch=3,col="blue")
points(TrimMSE.Mallows.Total~position.out,pch=4,col="purple")
points(TrimMSE.BY.Total~position.out,pch=5,col="red")

MSEdata <- data.frame(MSEval=c(TrimMSE.MLE.Total,TrimMSE.BY.Total,TrimMSE.ModMallows.Total,TrimMSE.Cubif.Total,TrimMSE.Mallows.Total),OutlierPos=rep(position.out,5),Method=rep(c("MLE","BY","ModMallows","CUBIF","Mallows"),rep(J,5)))
ScatterPlot <- ggplot(MSEdata, aes(x=OutlierPos, y=MSEval, shape=Method, color=Method)) +geom_point(size=4.5) + theme_bw() +theme(legend.position = "bottom") + ylab("Average Total MSE") + xlab("Outlier position")+ggtitle("Balanced dataset with one outlier")+theme(strip.text = element_text(size = 22)) + theme(text = element_text(size=22)) + theme(axis.text.x = element_text(size = 22, hjust = 1))
```

Verify the true coefficient values (for unbalanced design)
```{r}
n1big <- 40000
n0big <- 160000
ybig <- c(rep(1,n1big),rep(0,n0big))
set.seed(999999)
X1big <- matrix(rnorm( n=(2*n1big), mean=2,sd =1),nrow=n1big, ncol=2)
X0big <- matrix(rnorm( n=(2*n0big), mean=0,sd =1),nrow=n0big, ncol=2)  
TrueCoef <- coef(glm(ybig~rbind(X1big,X0big),family=binomial(link="logit")))
```

Repeated simulations (no outliers; unbalanced)
```{r}
TotalReps <- 500
MLEmatrix <- matrix(0,nrow=TotalReps,ncol=3)
BYmatrix <- matrix(0,nrow=TotalReps,ncol=3)
ModMallowsmatrix <- matrix(0,nrow=TotalReps,ncol=3)
cubifmatrix <- matrix(0,nrow=TotalReps,ncol=3)
mallowsmatrix <- matrix(0,nrow=TotalReps,ncol=3)
n1 <- 20*2
n0 <- 80*2
y <- c(rep(1,n1),rep(0,n0))

rep.id <- 0
num.reps <- 0
rep.idStore <- rep(0,TotalReps)
while(num.reps < TotalReps){
  
  rep.id <- rep.id + 1
  set.seed(rep.id)
  X1 <- matrix(rnorm((n1*2),mean=2,sd=1),nrow=n1,ncol=2)
  X0 <- matrix(rnorm((n0*2),mean=0,sd=1),nrow=n0,ncol=2)
  X <- rbind(X1,X0)
  X1out <- X1
  Xout <- rbind(X1out,X0)
  tempGLm <- glm(y[-1]~Xout[-1,],family=binomial(link="logit"))
  tempcoef <- coef(tempGLm)
  
  if( (tempGLm$converged)&(sum((tempcoef - TrueCoef)^2) < 3) ){
    
    num.reps <- num.reps + 1
    MLEmatrix[num.reps,] <- coef(glm(y~Xout,family=binomial(link="logit")))
    BYmatrix[num.reps,] <- coef(glmrob(y~Xout, family=binomial(link="logit"), method="BY", const =7, maxhalf=1))
    cubifmatrix[num.reps,] <- coef(glmRob(y~Xout, family=binomial(link="logit"), method="cubif", bpar=0.824,cpar=2))
  mallowsmatrix[num.reps,] <- coef(glmRob(y~Xout, family=binomial(link="logit"), method="mallows", wt.tuning=3.005))
  ModMallows.init <- coef(glmRob(y~Xout, family=binomial(link="logit"), method="mallows"))
  MMobj <- CalcModMallows(y=y,X=Xout,theta.init = ModMallows.init, h=5.3)
  ModMallowsmatrix[num.reps,] <- c(MMobj$theta.curr)
  rep.idStore[num.reps] <- rep.id
    
  }
  
}

MLEeff <- apply(MLEmatrix,1,function(s){ (s - TrueCoef)^2  })
BYeff <- apply(BYmatrix,1, function(s){ (s - TrueCoef)^2 })
ModMalloweff <- apply(ModMallowsmatrix,1,function(s){ (s - TrueCoef)^2  })
cubifeff <- apply(cubifmatrix,1,function(s){ (s - TrueCoef)^2  })
mallowseff <- apply(mallowsmatrix,1,function(s){ (s - TrueCoef)^2  })


MLEeff.total <- colSums(MLEeff)
BYeff.total <- colSums(BYeff)
ModMalloweff.total <- colSums(ModMalloweff)
cubifeff.total <- colSums(cubifeff)
mallowseff.total <- colSums(mallowseff)
```

Average relative efficiency (no outliers; unbalanced)
```{r}
#MLE vs BY
mean(BYeff.total)/mean(ModMalloweff.total)
#MLE vs ModMallows
mean(MLEeff.total)/mean(ModMalloweff.total)
#MLE vs cubif
mean(MLEeff.total)/mean(cubifeff.total)
#MLE vs mallows
mean(MLEeff.total)/mean(mallowseff.total)
```

Outlier and unbalanced
```{r, warning=FALSE}
position.out <- seq(-3.5,1,0.1)
J <- length(position.out)
prob.out <- pnorm(position.out,mean=2,sd=1)^2
n1 <- 20*2
n0 <- 80*2
y <- c(rep(1,n1),rep(0,n0))
TotalReps <- 500

MLEeff.total.Mat <- matrix(0,nrow=TotalReps,ncol=J)
ModMalloweff.total.Mat <- matrix(0,nrow=TotalReps,ncol=J)
cubifeff.total.Mat <- matrix(0,nrow=TotalReps,ncol=J)
mallowseff.total.Mat <- matrix(0,nrow=TotalReps,ncol=J)
for(j in 1:J){
  
  
MLEmatrix <- matrix(0,nrow=TotalReps,ncol=3)
ModMallowsmatrix <- matrix(0,nrow=TotalReps,ncol=3)
cubifmatrix <- matrix(0,nrow=TotalReps,ncol=3)
mallowsmatrix <- matrix(0,nrow=TotalReps,ncol=3)

for(num.reps in 1:TotalReps){
  
  set.seed(rep.idStore[num.reps])
  X1 <- matrix(rnorm((n1*2),mean=2,sd=1),nrow=n1,ncol=2)
  X0 <- matrix(rnorm((n0*2),mean=0,sd=1),nrow=n0,ncol=2)
  X <- rbind(X1,X0)
  X1out <- X1
  X1out[1,] <- c(position.out[j],position.out[j]) 
  Xout <- rbind(X1out,X0)
  
    MLEmatrix[num.reps,] <- coef(glm(y~Xout,family=binomial(link="logit")))
    BYmatrix[num.reps,] <- coef(glmrob(y~Xout, family=binomial(link="logit"), method="BY", const =7, maxhalf=1))
    cubifmatrix[num.reps,] <- coef(glmRob(y~Xout, family=binomial(link="logit"), method="cubif", bpar=0.824,cpar=2))
  mallowsmatrix[num.reps,] <- coef(glmRob(y~Xout, family=binomial(link="logit"), method="mallows", wt.tuning=3.005))
  ModMallows.init <- coef(glmRob(y~Xout, family=binomial(link="logit"), method="mallows"))
  MMobj <- CalcModMallows(y=y,X=Xout,theta.init = ModMallows.init, h=5.3)
  ModMallowsmatrix[num.reps,] <- c(MMobj$theta.curr)
    
  
  
}

MLEeff <- apply(MLEmatrix,1,function(s){ (s - TrueCoef)^2  })
BYeff <- apply(BYmatrix,1, function(s){ (s - TrueCoef)^2 })
ModMalloweff <- apply(ModMallowsmatrix,1,function(s){ (s - TrueCoef)^2  })
cubifeff <- apply(cubifmatrix,1,function(s){ (s - TrueCoef)^2  })
mallowseff <- apply(mallowsmatrix,1,function(s){ (s - TrueCoef)^2  })


MLEeff.total.Mat[,j] <- colSums(MLEeff)
BYeff.total.Mat[,j] <- colSums(BYeff)
ModMalloweff.total.Mat[,j] <- colSums(ModMalloweff)
cubifeff.total.Mat[,j] <- colSums(cubifeff)
mallowseff.total.Mat[,j] <-  colSums(mallowseff)
}
```


Plot Total MSE (averaged across 500 datasets)
```{r}
TrimMSE.MLE.Total <- colMeans(MLEeff.total.Mat)
TrimMSE.BY.Total <- colMeans(BYeff.total.Mat)
TrimMSE.ModMallows.Total <- colMeans(ModMalloweff.total.Mat)
TrimMSE.Cubif.Total <- colMeans(cubifeff.total.Mat)
TrimMSE.Mallows.Total <- colMeans(mallowseff.total.Mat)

plot(TrimMSE.MLE.Total~position.out,pch=1,ylab="Total MSE",xlab="outlier position", main="Total MSE")
points(TrimMSE.ModMallows.Total~position.out,pch=2,col="brown")
points(TrimMSE.Cubif.Total~position.out,pch=3,col="blue")
points(TrimMSE.Mallows.Total~position.out,pch=4,col="purple")
points(TrimMSE.BY.Total~position.out,pch=5,col="red")

MSEdata <- data.frame(MSEval=c(TrimMSE.MLE.Total,TrimMSE.BY.Total,TrimMSE.ModMallows.Total,TrimMSE.Cubif.Total,TrimMSE.Mallows.Total),OutlierPos=rep(position.out,5),Method=rep(c("MLE","BY","ModMallows","CUBIF","Mallows"),rep(J,5)))
ScatterPlot <- ggplot(MSEdata, aes(x=OutlierPos, y=MSEval, shape=Method, color=Method)) +geom_point(size=4.5) + theme_bw() +theme(legend.position = "bottom") + ylab("Average Total MSE") + xlab("Outlier position")+ggtitle("Unbalanced dataset with one outlier")+theme(strip.text = element_text(size = 22)) + theme(text = element_text(size=22)) + theme(axis.text.x = element_text(size = 22, hjust = 1))
```
