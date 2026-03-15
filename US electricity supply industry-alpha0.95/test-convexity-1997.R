#############################################
## US electriciy supply industry data
#############################################
# remove all the previous stuff in R 
rm(list = ls(all=TRUE))
require(readxl)
require(FEAR) 
if (exists(".Random.seed")) {
  save.seed=.Random.seed
  flag.seed=TRUE
} else {
  flag.seed=FALSE
}
set.seed(900001)
#######################  BEGIN Import Cleaned Data ##########################
df <- read.table("./Data/kt-data.txt", quote="\"", comment.char="")
######################## End Import Cleaned Data #############################
#

dim(df)
colnames(df)

for (Year in 1986:1998) {
  
  ii=which(df$Year==Year)
  cat(Year, ": ",length(ii), "\n")
  
}

##################################
res.alpha=matrix(NA, nrow=13, ncol=4)
res.M=matrix(NA, nrow=13, ncol=4)
years=c(1986:1998)
########################
year=1997
########################
j=year-1986+1
df=df[df$V9==year, ]
##################################

X1=df$X1=df$V4
X2=df$X2=df$V5
X3=df$X3=df$V6
Y1=df$Y1=df$V7
ID=df$ID=df$V8
Year=df$Year=df$V9

X=rbind(X1,X2,X3)
Y=rbind(Y1)
  
TX=t(X)
TY=t(Y)


N=ncol(X)
alpha = 0.95
M=floor(N*alpha)


############################# Full Frontier  ######################################
# Technical Efficiency, Nonconvex, FEAR Package
NC.full=FEAR::fdh(XOBS=X, YOBS=Y, XREF=X, YREF=Y, METRIC=2, ORIENTATION=1)
summary(NC.full)

# Technical Efficiency, Convex, FEAR Package
C.full=FEAR::dea(XOBS=X, YOBS=Y, XREF=X, YREF=Y, RTS=1, METRIC=2, ORIENTATION=1)
summary(C.full)

## Compare (Full, NConvex) Versus (Full, Convex): should be no less than 1
summary(NC.full/C.full)


############################# Order alpha ######################################
# Technical Efficiency, Order alpha, Nonconvex, FEAR Package
NC.alpha=1/FEAR::cquan(XOBS=X, YOBS=Y, XREF=X, YREF=Y, alpha=alpha, ORIENTATION=1)
summary(NC.alpha)

# Technical Efficiency, Order alpha, Convex, FEAR Package
X.eff.alpha=t(t(X)*NC.alpha)
#X.eff
C.alpha=FEAR::dea(XOBS=X, YOBS=Y, XREF=X.eff.alpha, YREF=Y, RTS=1, METRIC=2, ORIENTATION=1)
summary(C.alpha)

## Compare (Order alpha, NConvex) Versus (Order alpha, Convex): should be no less than 1
summary(NC.alpha/C.alpha)


############################# Order M ######################################
# Technical Efficiency, Order M, Nonconvex, FEAR Package
NC.M=FEAR::orderm(XOBS=X, YOBS=Y, XREF=X, YREF=Y, M=M, ORIENTATION=1)$dhat
summary(NC.M)


# Technical Efficiency, Order M, Convex, FEAR Package
X.eff.M=t(t(X)*NC.M)
#X.eff
C.M=FEAR::dea(XOBS=X, YOBS=Y, XREF=X.eff.M, YREF=Y, RTS=1, METRIC=2, ORIENTATION=1)
summary(C.M)

## Compare (Order alpha, NConvex) Versus (Order alpha, Convex): should be no less than 1
summary(NC.M/C.M)



###############  Compare Full Versus Partial Efficiency  ##############
# should be greater than or equal to 1
summary(NC.alpha/NC.full)

# should be greater than or equal to 1
summary(NC.M/NC.full)

# should be greater than or equal to 1
summary(C.alpha/C.full)

# should be greater than or equal to 1
summary(C.M/C.full)


######### Save aLl the Data For Full and Partial Efficiency ##########

res<-cbind(df,NC.full,C.full,NC.alpha,C.alpha,NC.M,C.M)
outfile=paste("./Estimates/efficiency-estimates-",year,".csv",sep="")
write.csv(res,row.names=FALSE,file=outfile)


###############    Apply KSW Test for Full Efficiency   ########

start.time <- Sys.time()
KSW.test=FEAR::test.convexity(X=X, Y=Y, ORIENTATION=1, METRIC=2, NSPLIT=100, NREP=1000,NBCR=100)
end.time <- Sys.time()


time.taken.full <- end.time - start.time
time.taken.full
time.taken.full.hours<-time.taken.full/3600
cat("Elapsed Time in Hour",time.taken.full.hours,"\n")


round(c(KSW.test$tau,KSW.test$pval),4)


######################  Save the Results into Latex  ###########################

res.full=matrix(NA, nrow=length(years), ncol=6)

res.full[j,1]=mean(C.full)
res.full[j,2]=mean(NC.full)
res.full[j,3]=KSW.test$tau[1]
res.full[j,4]=KSW.test$pval[1]
res.full[j,5]=KSW.test$tau[2]
res.full[j,6]=KSW.test$pval[2]
#res.full[j,7]=time.taken.full.hours

#### Latex for Full Efficiency
tex=formatC(years,width=6,digits=0,format="f")
for (k in 1:ncol(res.full)) {
  tex = paste(tex,"&",formatC(res.full[,k],width=6,digits = 3,format = "f"))
}
tex = paste(tex,"\\\\")

outfile=paste("./Output/test-convexity-full-",year,".tex",sep="")
write(tex,file=outfile)


###############   Apply Li Test for Partial Efficiency   ########
### USING R 4.4

require(np)

### Order Alpha
t1<-data.frame(x=NC.alpha)
t2<-data.frame(x=C.alpha)

start.time <- Sys.time()
li.alpha=np::npdeneqtest(t1,t2,boot.num=1000)
end.time <- Sys.time()

time.taken.alpha <- end.time - start.time
time.taken.alpha
time.taken.alpha.hours<-time.taken.alpha/3600
cat("Elapsed Time in Hour",time.taken.alpha.hours,"\n")


round(c(li.alpha$Tn, li.alpha$Tn.P),4)

### Order M
t1<-data.frame(x=NC.M)
t2<-data.frame(x=C.M)

start.time <- Sys.time()
li.M=np::npdeneqtest(t1,t2,boot.num=1000)
end.time <- Sys.time()

time.taken.M <- end.time - start.time
time.taken.M
time.taken.M.hours<-time.taken.M/3600
cat("Elapsed Time in Hour",time.taken.M.hours,"\n")

round(c(li.M$Tn, li.M$Tn.P),4)




######################  Save the Results into Latex  ###########################

res.alpha=matrix(NA, nrow=length(years), ncol=4)
res.M=matrix(NA, nrow=length(years), ncol=4)

res.alpha[j,1]=mean(C.alpha)
res.alpha[j,2]=mean(NC.alpha)
res.alpha[j,3]=li.alpha$Tn
res.alpha[j,4]=li.alpha$Tn.P
#res.alpha[j,5]=time.taken.alpha.hours

res.M[j,1]=mean(C.M)
res.M[j,2]=mean(NC.M)
res.M[j,3]=li.M$Tn
res.M[j,4]=li.M$Tn.P
#res.M[j,5]=time.taken.M.hours


#### Latex for Order Alpha
tex=formatC(years,width=6,digits=0,format="f")
for (k in 1:ncol(res.alpha)) {
  tex = paste(tex,"&",formatC(res.alpha[,k],width=6,digits = 3,format = "f"))
}
tex = paste(tex,"\\\\")

outfile=paste("./Output/test-convexity-alpha-",year,".tex",sep="")
write(tex,file=outfile)



#### Latex for Order M
tex=formatC(years,width=6,digits=0,format="f")
for (k in 1:ncol(res.M)) {
  tex = paste(tex,"&",formatC(res.M[,k],width=6,digits = 3,format = "f"))
}
tex = paste(tex,"\\\\")

outfile=paste("./Output/test-convexity-M-",year,".tex",sep="")
write(tex,file=outfile)