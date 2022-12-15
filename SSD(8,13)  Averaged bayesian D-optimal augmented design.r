# SSD(8,13)  Averaged bayesian D-optimal augmented design
library(randtoolbox)
quasi_integrated_bayes_D<-function(X,nrun,nprimary,nsecondary,npotential,tau,gamma,r){
  
  nfactor<-dim(X)[2]+r#add r new factors
  nrun1<-dim(X)[1]
  
  design<-matrix(runif(nrun*nfactor,min = -1,max = 1),ncol = nfactor)
  design[,14]<-rep(-1,nrun)
  model_design<-matrix(c(rep(1,nrun),design[,pro],design[,primary[!primary%in%pro]],design[,seq(nfactor-r+1,nfactor)],
                         design[,pro]^2,design[,seq(nfactor-r+1,nfactor)]^2,#primary
                         #design[,-c(primary,pro,secondary,14,15)],
                         design[,secondary],#secondary
                         design[,potential]#potential
  ),
  nrow = nrun)
  p1<-nprimary
  p2<-nsecondary
  p3<-npotential
  R<-diag(c(rep(0,p1),rep(1/(gamma^2),p2),rep(1/(tau^2),p3)))
  
  intergrated_det<-function(a)
  {
    integrated<-0
    n<-100
    x<-halton(n,nrun1*r)
    x<--1+2*(x-min(x))/(max(x)-min(x))
    # x15<-x[,1:8]
    # x16<-x[,9:16]
    for (i in 1:n) {
      model_design_1<-matrix(c(rep(1,nrun1),X[,pro],X[,primary[!primary%in%pro]],matrix(x[i,],ncol = r),
                               X[,pro]^2,matrix(x[i,],ncol = r),#primary
                               X[,secondary],#secondary
                               X[,potential]#potential
      ),
      nrow = nrun1)
      ddet<-det(t(model_design_1)%*%model_design_1+t(model_design)%*%model_design+R)
      
      integrated <-ddet+integrated 
    }
    integrated
  }
  
  level1<-c(-1,0,1)
  level2<-c(-1,1)
  
  col1<-c(pro,seq(nfactor-r+1,nfactor))
  col2<-seq(1,nfactor)[(!seq(1,nfactor)%in%col1)&(!seq(1,nfactor)==14)]
  
  ###############
  ddet_o<-intergrated_det(a)
  for (m in 1:nrun) {
    for(n in col1){
      d_ddet1<-rep(0,length(level1))
      for (i in 1:length(level1)) {
        design[m,n]<-level1[i]
        model_design<-matrix(c(rep(1,nrun),design[,pro],design[,primary[!primary%in%pro]],design[,seq(nfactor-r+1,nfactor)],
                               design[,pro]^2,design[,seq(nfactor-r+1,nfactor)]^2,#primary
                               #design[,-c(primary,pro,secondary,14,15)],
                               design[,secondary],#secondary
                               design[,potential]#potential
        ),
        nrow = nrun)
        d_ddet1[i]<-intergrated_det(a)
      }
      if(max(d_ddet1)>=ddet_o){design[m,n]<-level1[d_ddet1==max(d_ddet1)][1]}
      model_design<-matrix(c(rep(1,nrun),design[,pro],design[,primary[!primary%in%pro]],design[,seq(nfactor-r+1,nfactor)],
                             design[,pro]^2,design[,seq(nfactor-r+1,nfactor)]^2,#primary
                             #design[,-c(primary,pro,secondary,14,15)],
                             design[,secondary],#secondary
                             design[,potential]#potential
      ),
      nrow = nrun)
      ddet_o<-intergrated_det(a)
    }
    for(n in col2){
      d_ddet2<-rep(0,length(level2))
      for (i in 1:length(level2)) {
        design[m,n]<-level2[i]
        model_design<-matrix(c(rep(1,nrun),design[,pro],design[,primary[!primary%in%pro]],design[,seq(nfactor-r+1,nfactor)],
                               design[,pro]^2,design[,seq(nfactor-r+1,nfactor)]^2,#primary
                               #design[,-c(primary,pro,secondary,14,15)],
                               design[,secondary],#secondary
                               design[,potential]#potential
        ),
        nrow = nrun)
        
        d_ddet2[i]<-intergrated_det(a)
      }
      if(max(d_ddet2)>=ddet_o){design[m,n]<-level2[d_ddet2==max(d_ddet2)][1]}
      model_design<-matrix(c(rep(1,nrun),design[,pro],design[,primary[!primary%in%pro]],design[,seq(nfactor-r+1,nfactor)],
                             design[,pro]^2,design[,seq(nfactor-r+1,nfactor)]^2,#primary
                             #design[,-c(primary,pro,secondary,14,15)],
                             design[,secondary],#secondary
                             design[,potential]#potential
      ),
      nrow = nrun)
      ddet_o<-intergrated_det(a)
    }
    
  }
  
  model_design<-matrix(c(rep(1,nrun),design[,pro],design[,primary[!primary%in%pro]],design[,seq(nfactor-r+1,nfactor)],
                         design[,pro]^2,design[,seq(nfactor-r+1,nfactor)]^2,#primary
                         #design[,-c(primary,pro,secondary,14,15)],
                         design[,secondary],#secondary
                         design[,potential]#potential
  ),
  nrow = nrun)
  
  list(design,intergrated_det(a))
}

