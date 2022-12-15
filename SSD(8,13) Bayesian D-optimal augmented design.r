# SSD(8,13) Bayesian D-optimal augmented design
bayes_D_2<-function(X,nrun,nprimary,nsecondary,npotential,tau,gamma,r){
   
  nfactor<-dim(X)[2]+r#add r new factors
  nrun1<-dim(X)[1]
  
  model_design_1<-matrix(c(rep(1,nrun1),X[,pro],X[,primary[!primary%in%pro]],matrix(rep(0,r*nrun1),ncol = r),
                           X[,pro]^2,matrix(rep(0,r*nrun1),ncol = r),#primary
                           X[,secondary],#secondary
                           X[,potential]#potential
  ),
  nrow = nrun1)
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
  ddet <- det(t(model_design_1)%*%model_design_1+t(model_design)%*%model_design+R)
  
  level1<-c(-1,0,1)
  level2<-c(-1,1)
  
  col1<-c(pro,seq(nfactor-r+1,nfactor))
  col2<-seq(1,nfactor)[(!seq(1,nfactor)%in%col1)&(!seq(1,nfactor)==14)]
  
  ###############
  ddet_o<-ddet
  ddet_12<-c(0,2)
  while (ddet_12[2]>ddet_12[1]) {
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
        ddet <- det(t(model_design_1)%*%model_design_1+t(model_design)%*%model_design+R)
        d_ddet1[i]<-ddet
      }
      if(max(d_ddet1)>=ddet_o){design[m,n]<-level1[d_ddet1==max(d_ddet1)][1]}
      model_design<-matrix(c(rep(1,nrun),design[,pro],design[,primary[!primary%in%pro]],design[,seq(nfactor-r+1,nfactor)],
                             design[,pro]^2,design[,seq(nfactor-r+1,nfactor)]^2,#primary
                             #design[,-c(primary,pro,secondary,14,15)],
                             design[,secondary],#secondary
                             design[,potential]#potential
      ),
      nrow = nrun)
      ddet_o<-det(t(model_design_1)%*%model_design_1+t(model_design)%*%model_design+R)
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
        ddet <- det(t(model_design_1)%*%model_design_1+t(model_design)%*%model_design+R)
        d_ddet2[i]<-ddet
      }
      if(max(d_ddet2)>=ddet_o){design[m,n]<-level2[d_ddet2==max(d_ddet2)][1]}
      model_design<-matrix(c(rep(1,nrun),design[,pro],design[,primary[!primary%in%pro]],design[,seq(nfactor-r+1,nfactor)],
                             design[,pro]^2,design[,seq(nfactor-r+1,nfactor)]^2,#primary
                             #design[,-c(primary,pro,secondary,14,15)],
                             design[,secondary],#secondary
                             design[,potential]#potential
      ),
      nrow = nrun)
      ddet_o<-det(t(model_design_1)%*%model_design_1+t(model_design)%*%model_design+R)
    }
    
  }
    ddet_12[1]<-ddet_o
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
          ddet <- det(t(model_design_1)%*%model_design_1+t(model_design)%*%model_design+R)
          d_ddet1[i]<-ddet
        }
        if(max(d_ddet1)>=ddet_o){design[m,n]<-level1[d_ddet1==max(d_ddet1)][1]}
        model_design<-matrix(c(rep(1,nrun),design[,pro],design[,primary[!primary%in%pro]],design[,seq(nfactor-r+1,nfactor)],
                               design[,pro]^2,design[,seq(nfactor-r+1,nfactor)]^2,#primary
                               #design[,-c(primary,pro,secondary,14,15)],
                               design[,secondary],#secondary
                               design[,potential]#potential
        ),
        nrow = nrun)
        ddet_o<-det(t(model_design_1)%*%model_design_1+t(model_design)%*%model_design+R)
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
          ddet <- det(t(model_design_1)%*%model_design_1+t(model_design)%*%model_design+R)
          d_ddet2[i]<-ddet
        }
        if(max(d_ddet2)>=ddet_o){design[m,n]<-level2[d_ddet2==max(d_ddet2)][1]}
        model_design<-matrix(c(rep(1,nrun),design[,pro],design[,primary[!primary%in%pro]],design[,seq(nfactor-r+1,nfactor)],
                               design[,pro]^2,design[,seq(nfactor-r+1,nfactor)]^2,#primary
                               #design[,-c(primary,pro,secondary,14,15)],
                               design[,secondary],#secondary
                               design[,potential]#potential
        ),
        nrow = nrun)
        ddet_o<-det(t(model_design_1)%*%model_design_1+t(model_design)%*%model_design+R)
      }
      
    }
    ddet_12[2]<-ddet_o 
  }
  
  
  model_design<-matrix(c(rep(1,nrun),design[,pro],design[,primary[!primary%in%pro]],design[,seq(nfactor-r+1,nfactor)],
                         design[,pro]^2,design[,seq(nfactor-r+1,nfactor)]^2,#primary
                         #design[,-c(primary,pro,secondary,14,15)],
                         design[,secondary],#secondary
                         design[,potential]#potential
  ),
  nrow = nrun)
  
  
  ddet <- det(t(model_design_1)%*%model_design_1+t(model_design)%*%model_design+R)
  list(design,ddet)
}
