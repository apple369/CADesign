#SSD(8,13),congsider impacts of the number of adding new factors r
#fixed n2=12,beta=6,Ma=100,Mb=10,N=100

rm(list=ls())

source("SSD(8,13) Bayesian D-optimal augmented design.r")
source("SSD(8,13)  Averaged bayesian D-optimal augmented design.r")
source('select.r')
library(foreach)
library(doParallel)
###################
#run=300
detectCores( )
cl.no=detectCores(logical = F)-1 
cl=makeCluster(cl.no) 
registerDoParallel(cl)

#SSD(8,13)
X1<-c(1,1,1,1,-1,-1,-1,-1)
X2<-c(1,1,-1,-1,1,1,-1,-1)
X3<-c(1,1,-1,1,-1,-1,1,-1)
X4<-c(1,-1,-1,1,1,-1,-1,1)
X5<-c(1,-1,1,1,-1,1,-1,-1)
X6<-c(1,1,-1,-1,-1,-1,1,1)
X7<-c(1,-1,1,-1,-1,-1,1,1)
X8<-c(1,-1,-1,-1,1,1,1,-1)
X9<-c(1,-1,1,-1,1,-1,1,-1)
X10<-c(1,1,-1,-1,-1,1,-1,1)
X11<-c(1,-1,-1,1,-1,1,1,-1)
X12<-c(1,-1,1,-1,-1,1,-1,1)
X13<-c(1,1,1,-1,1,-1,-1,-1)
x14<-rep(1,8)
X_stage_1<-matrix(c(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,x14),ncol = 14)
###-------paremeters in the functions
X<-X_stage_1
tau=sqrt(5)
gamma=10
nrun<-12
r<-2
#
nfactor<-dim(X)[2]+r#add r new factors
nrun1<-dim(X)[1]
#######
AEIb<-0
IEIb<-0
AEIab<-0
IEIab<-0
a<-0

ptm=proc.time()

for (i in 1:300) {
#random setting
sam<-rnorm(6,6,1)*sample(c(-1, 1), 6, replace=TRUE)
loc1<-sample(seq(1,13),3, replace=FALSE)
if(length(seq(nfactor-r+1,nfactor))==1){loc2<-seq(nfactor-r+1,nfactor)}else{loc2<-sample(seq(nfactor-r+1,nfactor),1, replace=FALSE)}
loc3<-sample(loc1,1,replace=FALSE)
# #fixed setting
# sam<-c(8,6,9,10,4,7)
# loc1<-c(4,5,11)
# loc2<-15
# loc3<-11
Y1<-sam[1]*X_stage_1[,loc1[1]]+sam[2]*X_stage_1[,loc1[2]]+sam[3]*X_stage_1[,loc1[3]]+sam[4]*X_stage_1[,loc3]^2+sam[5]*X_stage_1[,14]+rnorm(dim(X_stage_1)[1],0,1)
#library(VGAM)
#Y1<-sam[1]*X_stage_1[,loc1[1]]+sam[2]*X_stage_1[,loc1[2]]+sam[3]*X_stage_1[,loc1[3]]+sam[4]*X_stage_1[,loc3]^2+sam[5]*X_stage_1[,14]+rlaplace(dim(X_stage_1)[1],0,1)
data<-data.frame(X_stage_1)
design_matrix<- as.matrix(data)
#################  scad,mcp,lasso
z<-select(design_matrix,Y1)

########################### pro
table1<-table(c(which(z[[1]]!=0),which(z[[2]]!=0),which(z[[3]]!=0),which(z[[4]]!=0)))

if(length(names(table1)[table1>=3])==0){primary<- 0} else{primary<-as.integer(names(table1)[table1>=3])}

# #fixed setting
# pro<-c(3,11)

#random setting
pro<-c(sample(seq(1,13)[!seq(1,13)%in%loc1],1,replace=F),loc3)

if(length(names(table1)[table1<3])==0){ secondary<- 0} else{secondary<-as.integer(names(table1)[table1<3])}
potential<-seq(1,13)[(!seq(1,13)%in%primary)&(!seq(1,13)%in%secondary)&(!seq(1,13)%in%pro)]


#############
if(primary==0){nprimary<-1+2*length(pro)+2*r}else{nprimary<-1+2*length(pro)+length(primary[!primary%in%pro])+2*r}#眒蜊
if(secondary==0){nsecondary<-0}else{nsecondary<-length(secondary)}
if(length(potential)==0){npotential<-0}else{npotential<-length(potential)}

###############
#Bayesian D-optimal:
runb=100
resultb=foreach(h=1:runb,.packages = c("ncvreg","fastclime"),.combine = rbind) %dopar%
  bayes_D_2(X,nrun,nprimary,nsecondary,npotential,tau,gamma,r) 

detb<-foreach(h=(runb+1):(2*runb),.combine = c) %dopar% resultb[[h]]
indexb<-seq(1:runb)[detb==max(detb)]
X_stage_2b<-resultb[[indexb[1]]]

Y2b<-sam[1]*X_stage_2b[,loc1[1]]+sam[2]*X_stage_2b[,loc1[2]]+sam[3]*X_stage_2b[,loc1[3]]+sam[4]*X_stage_2b[,loc3]^2+sam[5]*X_stage_2b[,14]+sam[6]*X_stage_2b[,loc2]+rnorm(dim(X_stage_2b)[1],0,1)
#Y2b<-sam[1]*X_stage_2b[,loc1[1]]+sam[2]*X_stage_2b[,loc1[2]]+sam[3]*X_stage_2b[,loc1[3]]+sam[4]*X_stage_2b[,loc3]^2+sam[5]*X_stage_2b[,14]+sam[6]*X_stage_2b[,loc2]+rlaplace(dim(X_stage_2b)[1],0,1)
Y_fullb<-c(Y1,Y2b)
X_fullb<-rbind(cbind(X_stage_1,matrix(rep(0,r*dim(X_stage_1)[1]),ncol = r)),X_stage_2b)
data_fullb<-data.frame(X_fullb,X_fullb[,pro]^2,X_fullb[,seq(nfactor-r+1,nfactor)]^2)#include of identifing r new factors
design_matrix_fullb<- as.matrix(data_fullb)

#################  scad,mcp,lasso
zb<-select(design_matrix_fullb,Y_fullb)

table2b<-table(c(which(zb[[1]]!=0),which(zb[[2]]!=0),which(zb[[3]]!=0),which(zb[[4]]!=0)))
identib<-as.integer(names(table2b)[table2b>=3])#祫屾 3 笱源楊恁堤

#######################
#Averaged Bayesian D-optimal:
runab=10
resultab=foreach(h=1:runab,.packages = c("ncvreg","fastclime","randtoolbox"),.combine = rbind) %dopar%
  quasi_integrated_bayes_D(X,nrun,nprimary,nsecondary,npotential,tau,gamma,r) 

detab<-foreach(h=(runab+1):(2*runab),.combine = c) %dopar% resultab[[h]]
indexab<-seq(1:runab)[detab==max(detab)]
X_stage_2ab<-resultab[[indexab[1]]]
#Y2ab<-sam[1]*X_stage_2ab[,loc1[1]]+sam[2]*X_stage_2ab[,loc1[2]]+sam[3]*X_stage_2ab[,loc1[3]]+sam[4]*X_stage_2ab[,loc3]^2+sam[5]*X_stage_2ab[,14]+sam[6]*X_stage_2ab[,loc2]+rlaplace(dim(X_stage_2ab)[1],0,1)
Y2ab<-sam[1]*X_stage_2ab[,loc1[1]]+sam[2]*X_stage_2ab[,loc1[2]]+sam[3]*X_stage_2ab[,loc1[3]]+sam[4]*X_stage_2ab[,loc3]^2+sam[5]*X_stage_2ab[,14]+sam[6]*X_stage_2ab[,loc2]+rnorm(dim(X_stage_2ab)[1],0,1)

Y_fullab<-c(Y1,Y2ab)

X_fullab<-rbind(cbind(X_stage_1,matrix(rep(0,r*dim(X_stage_1)[1]),ncol = r)),X_stage_2ab)
#X_full<-X_stage_2
data_fullab<-data.frame(X_fullab,X_fullab[,pro]^2,X_fullb[,seq(nfactor-r+1,nfactor)]^2)#include of identifing 2 new factors
design_matrix_fullab<- as.matrix(data_fullab)

#################  scad,mcp,lasso
zab<-select(design_matrix_fullab,Y_fullab)

table2ab<-table(c(which(zab[[1]]!=0),which(zab[[2]]!=0),which(zab[[3]]!=0),which(zab[[4]]!=0)))
identiab<-as.integer(names(table2ab)[table2ab>=3])#祫屾 3 笱源楊恁堤


#########################
true<-c(loc1,14,loc2,dim(X_stage_1)[2]+r+seq(1,length(c(pro,seq(nfactor-r+1,nfactor))))[c(pro,seq(nfactor-r+1,nfactor))==loc3])

AEIb<-all(true%in%identib)+AEIb
IEIb<-any(identib%in%true==FALSE)+IEIb

AEIab<-all(true%in%identiab)+AEIab
IEIab<-any(identiab%in%true==FALSE)+IEIab

a<-a+1


}

proc.time()-ptm
stopImplicitCluster()

#save.image("laplace_SSD(8,13)_r=2,n2=12_Bayesian & Averaged Bayesian.RData")
