#SSD(7ㄛ15),congsider impacts of the number of adding new factors r
#fixed n2=12,beta=8,Mb=100,Ma=10,N=100

rm(list=ls())

setwd("D:/OneDrive/OneDrive - mail.nankai.edu.cn/醫栠/code(1)/code")
source("SSD(7,15) Bayesian D-optimal augmented design.r")
source("SSD(7,15) Averaged bayesian D-optimal augmented design.r")
source('select.r')
library(foreach)
library(doParallel)
###################
#run=300 
detectCores(logical = F) 
cl.no=2#detectCores()-1
cl=makeCluster(cl.no)
registerDoParallel(cl)

AEIb<-0
IEIb<-0
AEIab<-0
IEIab<-0

a=0

#SSD(7,15)
X1<-c(-1,-1,-1,1,1,-1,1)
X2<-c(1,-1,1,-1,-1,1,-1)
X3<-c(-1,-1,1,-1,-1,1,1)
X4<-c(-1,-1,-1,1,-1,1,1)
X5<-c(-1,1,1,-1,1,-1,-1)
X6<-c(-1,-1,1,-1,1,-1,1)
X7<-c(1,-1,-1,-1,1,-1,1)
X8<-c(1,1,1,-1,-1,-1,1)
X9<-c(1,-1,-1,-1,1,1,1)
X10<-c(-1,-1,1,1,1,1,-1)
X11<-c(-1,-1,1,-1,1,1,1)
X12<-c(1,-1,1,1,1,-1,-1)
X13<-c(1,-1,1,1,-1,-1,1)
X14<-c(-1,1,1,1,-1,-1,1)
X15<-c(-1,1,-1,-1,1,1,1)
x16<-rep(1,7)
X_stage_1<-matrix(c(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,x16),ncol = 16,byrow = F)
#########
###-------paremeters in the functions
X<-X_stage_1
tau=sqrt(5)
gamma=10
nrun<-15
r<-2
beta<-6
#
nfactor<-dim(X)[2]+r#add r new factors
nrun1<-dim(X)[1]



ptm=proc.time()
##
for (i in 1:300) {
  
  ##fixed setting 1
  #sam <- c(8,6,11,4,7)
  #sam2<-5    #bilinear effect
  #loc1 <- c(5,10,14)
  #loc2 <- 17
  #interaction factor:x1x10(true),x5(pri)x3,x9(pri)x12
  #loc3<- c(10,5,10)
  #loc4 <-1 #interaction of other factor

  # random setting 2
  #set.seed(12)
  sam<-rnorm(6,beta,1)*sample(c(-1, 1), 6, replace=TRUE)
  sam2 <- sam[6]  #bilinear effect
  loc1 <- sample(seq(1,15),3, replace=FALSE)  #inherent factors
  if(length(seq(nfactor-r+1,nfactor))==1){loc2<-seq(nfactor-r+1,nfactor)}else{loc2<-sample(seq(nfactor-r+1,nfactor),1, replace=FALSE)}#block factor
  loc3 <-sample(loc1,3,replace=T)#loc3[1] is for true interaction, three factors in parent ones involving interaction terms
  loc4 <- sample(seq(1,15)[!seq(1,15)%in%loc1],1,replace=FALSE)#one of true interaction factors not in inherent factors
  
  
  Y1<-sam[1]*X_stage_1[,loc1[1]]+sam[2]*X_stage_1[,loc1[2]]+sam[3]*X_stage_1[,loc1[3]]+sam[4]*X_stage_1[,16]+rnorm(dim(X_stage_1)[1],0,1)+sam2*X_stage_1[,loc3[1]]*X_stage_1[,loc4]
  data<-data.frame(X_stage_1)
  design_matrix<- as.matrix(data)
  #################  scad,mcp,lasso
  z<-select(design_matrix,Y1)
  ############
 
  ########################### pro
  table1<-table(c(which(z[[1]]!=0),which(z[[2]]!=0),which(z[[3]]!=0),which(z[[4]]!=0)))
  table1 <- table1[!names(table1)==16]
  # 2: block factor in primary group
  if(length(names(table1)[table1>=3])==0){primary<- 16} else{primary<-c(as.integer(names(table1)[table1>=3]),16)}
  
  ##set1
  #pro1<-1
  #pro2<-14
  #set2
  pro1<-sample(seq(1,15)[!seq(1,15)%in%loc1],1,replace=F)
  pro2<-sample(loc1,1,replace=F)
  
  if(length(primary[primary!=16])==0){pro3 <-0}else{pro3 <- primary[primary!=16]}# x_pri: x_pri*x_other, primary in other interaction terms
  pro4 <-sample(seq(1,15)[(!seq(1,15)%in%pro3)&(!seq(1,15)%in%loc1)],length(pro3),replace=T) # x_other: x_pri*x_other, other factor in other interaction term 
  pro5 <-sample(seq(1,15)[(!seq(1,15)%in%pro3)&(!seq(1,15)%in%loc1)],2,replace=T)# x_other: x_true*x_other, other factor in other interaction term 
  #
  # 2: block factor in primary group
  if(length(names(table1)[table1<3])==0){ secondary<- 0 }else{secondary<-as.integer(names(table1)[table1<3])[as.integer(names(table1)[table1<3])!=pro1]}
  potential<-seq(1,15)[(!seq(1,15)%in%primary)&(!seq(1,15)%in%secondary)&(!seq(1,15)%in%pro1)&(!seq(1,15)%in%pro2)]
  #
  if(pro3==0){nprimary<-1+2+1+length(pro1)+2*length(pro2)+length(primary[(!primary%in%pro1)&(!primary%in%pro2)])+2*r}else{nprimary<-1+length(pro3)+2+1+length(pro1)+2*length(pro2)+length(primary[(!primary%in%pro1)&(!primary%in%pro2)])+2*r}
  
  if(length(secondary)==0){
    nsecondary<-0
  }else if(secondary==0){
    nsecondary<-0
  }else{
    nsecondary<-length(secondary)}
  if(length(potential)==0){npotential<-0}else{npotential<-length(potential)}
  ####################################
  
  ######--------------Bayesian 
  runb=100
  resultb=foreach(h=1:runb,.packages = c("ncvreg","fastclime"),.combine = rbind) %dopar%
    bayes_D_2(X,nrun,nprimary,nsecondary,npotential,tau,gamma,r) 
  
  detb<-foreach(h=(runb+1):(2*runb),.combine = c) %dopar% resultb[[h]]
  indexb<-seq(1:runb)[detb==max(detb)]
  X_stage_2b<-resultb[[indexb[1]]]
  ###############
  
  Y2b<-sam[1]*X_stage_2b[,loc1[1]]+sam[2]*X_stage_2b[,loc1[2]]+sam[3]*X_stage_2b[,loc1[3]]+sam[4]*X_stage_2b[,16]+sam[5]*X_stage_2b[,loc2]^2+rnorm(dim(X_stage_2b)[1],0,1)+sam2*X_stage_2b[,loc3[1]]*X_stage_2b[,loc4]
  Y_fullb<-c(Y1,Y2b)
  X_fullb<-rbind(cbind(X_stage_1,matrix(rep(0,r*dim(X_stage_1)[1]),ncol = r)),X_stage_2b)#imputes 0
  data_fullb<-data.frame(X_fullb,X_fullb[,pro2]^2,X_fullb[,seq(nfactor-r+1,nfactor)]^2,X_fullb[,loc3[1]]*X_fullb[,loc4],X_fullb[,pro3]*X_fullb[,pro4],X_fullb[,loc3[-1]]*X_fullb[,pro5])#include of identifing 2 new factors
  design_matrix_fullb<- as.matrix(data_fullb)
  
  #################  scad,mcp,lasso
  zb<-select(design_matrix_fullb,Y_fullb)
  
  table2b<-table(c(which(zb[[1]]!=0),which(zb[[2]]!=0),which(zb[[3]]!=0),which(zb[[4]]!=0)))
  identib<-as.integer(names(table2b)[table2b>=3])#
  ##########################
  
  ######--------------Averaged Bayesian 
  runab=10
  resultab=foreach(h=1:runab,.packages = c("ncvreg","fastclime","randtoolbox"),.combine = rbind) %dopar%
    quasi_integrated_bayes_D(X,nrun,nprimary,nsecondary,npotential,tau,gamma,r) 
  
  detab<-foreach(h=(runab+1):(2*runab),.combine = c) %dopar% resultab[[h]]
  indexab<-seq(1:runab)[detab==max(detab)]
  X_stage_2ab<-resultab[[indexab[1]]]
  ###############
  
  Y2ab<-sam[1]*X_stage_2ab[,loc1[1]]+sam[2]*X_stage_2ab[,loc1[2]]+sam[3]*X_stage_2ab[,loc1[3]]+sam[4]*X_stage_2ab[,16]+sam[5]*X_stage_2ab[,loc2]^2+rnorm(dim(X_stage_2ab)[1],0,1)+sam2*X_stage_2ab[,loc3[1]]*X_stage_2ab[,loc4]
  Y_fullab<-c(Y1,Y2ab)
  X_fullab<-rbind(cbind(X_stage_1,matrix(rep(0,r*dim(X_stage_1)[1]),ncol = r)),X_stage_2ab)#imputes 0
  data_fullab<-data.frame(X_fullab,X_fullab[,pro2]^2,X_fullab[,seq(nfactor-r+1,nfactor)]^2,X_fullab[,loc3[1]]*X_fullab[,loc4],X_fullab[,pro3]*X_fullab[,pro4],X_fullab[,loc3[-1]]*X_fullab[,pro5])#include of identifing 2 new factors
  design_matrix_fullab<- as.matrix(data_fullab)
  
  #################  scad,mcp,lasso
  zab<-select(design_matrix_fullab,Y_fullab)
  
  table2ab<-table(c(which(zab[[1]]!=0),which(zab[[2]]!=0),which(zab[[3]]!=0),which(zab[[4]]!=0)))
  identiab<-as.integer(names(table2ab)[table2ab>=3])#
  ##########################
  true<-c(loc1,16,dim(X_stage_1)[2]+r+seq(1,length(c(pro2,seq(nfactor-r+1,nfactor))))[c(pro2,seq(nfactor-r+1,nfactor))==loc2],18+2*r)#
  
  AEIb<-all(true%in%identib)+AEIb
  IEIb<-any(identib%in%true==FALSE)+IEIb
  
  AEIab<-all(true%in%identiab)+AEIab
  IEIab<-any(identiab%in%true==FALSE)+IEIab
  
  a<-a+1
  
}


proc.time()-ptm
stopImplicitCluster()

#save.image("SSD(7,15)_r=2,n2=15_Bayesian & Averaged Bayesian.RData")

