
#################  scad,mcp,lasso
select<-function(design_matrix_full,Y_full){
   # design_matrix_full<-design_matrix_fullab
   # Y_full<-Y_fullab
  library(ncvreg)
  ###
  scad_cv<-cv.ncvreg(design_matrix_full,Y_full,family = c("gaussian"),penalty = c("SCAD"))
  scad_cv$lambda.min
  mySCAD=ncvreg(design_matrix_full,Y_full,family = c("gaussian"),penalty = c("SCAD"),lambda=scad_cv$lambda.min)
  ncv_SCAD<-mySCAD$beta[-1]
  which(ncv_SCAD!=0)
  #mcp
  mcp_cv<-cv.ncvreg(design_matrix_full,Y_full,family = c("gaussian"),penalty = c("MCP"))
  mcp_cv$lambda.min
  myMCP=ncvreg(design_matrix_full,Y_full,family = c("gaussian"),penalty = c("MCP"),lambda=mcp_cv$lambda.min)
  ncv_MCP<-myMCP$beta[-1]
  which(ncv_MCP!=0)
  #lasso
  lasso_cv<-cv.ncvreg(design_matrix_full,Y_full,family = c("gaussian"),penalty = c("lasso"))
  lasso_cv$lambda.min
  mylasso=ncvreg(design_matrix_full,Y_full,family = c("gaussian"),penalty = c("lasso"),lambda=mcp_cv$lambda.min)
  ncv_lasso<-mylasso$beta[-1]
  which(ncv_lasso!=0)
  ############ DS
  library(fastclime)
  b = dantzig(design_matrix_full,Y_full, lambda = 0.01, nlambda = 100)
  # estimated regression coefficient vector
  c = dantzig.selector(b$lambdalist, b$BETA0, 25)
  which(c!=0)
  
  list(ncv_SCAD,ncv_MCP,ncv_lasso,c)
}

# z<-select(design_matrix,Y1)
# table(c(which(z[[1]]!=0),which(z[[2]]!=0),which(z[[3]]!=0),which(z[[4]]!=0)))
