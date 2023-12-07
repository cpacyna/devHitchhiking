exact.binomial=function(gender,NV,NR,cutoff=-5){
  # Function to filter out germline variants based on unmatched
  # variant calls of multiple samples from same individual (aggregate coverage
  # ideally >150 or so, but will work with less). NV is matrix of reads supporting 
  # variants and NR the matrix with total depth (samples as columns, mutations rows, 
  # with rownames as chr_pos_ref_alt or equivalent). Returns a logical vector, 
  # TRUE if mutation is likely to be germline.
  
  XY_chromosomal = grepl("X|Y",rownames(NR))
  autosomal = !XY_chromosomal
  
  if(gender=="female"){
    NV_vec = rowSums(NV)
    NR_vec = rowSums(NR)
    pval = rep(1,length(NV_vec))
    for (n in 1:length(NV_vec)){
      if(NR_vec[n]>0){
        pval[n] = binom.test(x=NV_vec[n],
                             n=NR_vec[n],
                             p=0.5,alt='less')$p.value
      }
      if (n%%1000==0){
        print(n)
      }
    }
  }
  # For male, split test in autosomal and XY chromosomal part
  if(gender=="male"){
    pval=rep(1,nrow(NV))
    NV_vec = rowSums(NV)[autosomal]
    NR_vec = rowSums(NR)[autosomal]
    pval_auto = rep(1,sum(autosomal))
    pval_XY = rep(1,sum(XY_chromosomal))
    
    for (n in 1:sum(autosomal)){
      if(NR_vec[n]>0){
        pval_auto[n] = binom.test(x=NV_vec[n],
                                  n=NR_vec[n],
                                  p=0.5,alt='less')$p.value
      }
      if (n%%1000==0){
        print(n)
      }
    }
    NV_vec = rowSums(NV)[XY_chromosomal]
    NR_vec = rowSums(NR)[XY_chromosomal]
    for (n in 1:sum(XY_chromosomal)){
      if(NR_vec[n]>0){
        pval_XY[n] = binom.test(x=NV_vec[n],
                                n=NR_vec[n],
                                p=0.95,alt='less')$p.value
      }
      if (n%%1000==0){
        print(n)
      }
    }
    pval[autosomal]=pval_auto
    pval[XY_chromosomal]=pval_XY
  }
  qval = p.adjust(pval,method="BH")
  germline = log10(qval)>cutoff
  return(germline)
}

estimateRho_gridml = function(x, mu) {
  # Estimate rho by MLE grid approach
  rhovec = 10^seq(-6,-0.5,by=0.05) # rho will be bounded within 1e-6 and 0.32
  mm = x[,2]
  #cov = c(x[,3:4])+c(x[,1:2])
  cov = c(x[,1])
  ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=mm, size=cov, rho=rhoj, prob=mu, log=T)))
  rhovec[ll==max(ll)][1]
}

logbb <- function(x, n, mu, disp) {
  lbeta(x + mu, n - x - mu + disp) - lbeta(mu, disp - mu)
}



#-------------------------------------------------
# Shearwater-like filter for WGS (post-CaVEMan)
# Tim Coorens - November 2018
#-------------------------------------------------
options(stringsAsFactors=F)

#-------------------------------------------------
# Libraries
#-------------------------------------------------

require("GenomicRanges")
library("deepSNV")
library("Rsamtools")

logbb = deepSNV:::logbb
dbetabinom = VGAM::dbetabinom

#-------------------------------------------------
# Functions
#-------------------------------------------------

estimateRho_gridml = function(x, mu) {
  # Estimate rho by MLE grid approach
  rhovec = 10^seq(-6,-0.5,by=0.05) # rho will be bounded within 1e-6 and 0.32
  mm = x[,2]
  #cov = c(x[,3:4])+c(x[,1:2])
  cov = c(x[,1])
  ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=mm, size=cov, rho=rhoj, prob=mu, log=T)))
  rhovec[ll==max(ll)][1]
}

shearwater_probability=function(patient, save=NULL, path_prefix='', rho=10^-3, mutCodes, cohortSamples, all_counts, normal_samples){
  #Function to calculcate probability of presence of mutation based on Shearwarer
  #'patient' is the name of the patient-specific subfolder
  #'path_prefix' is any prefix to the path necessary to find that subfolder
  #'rho' is the constant for the overdispersion parameter. If rho=NULL, calculate
  # it from the data (much slower)
  #'save' is path for output. If NULL, returns matrix
  # output of this function will be a matrix (muts by samples) of p values
  # 'mutCodes' is a list of all the mutations in format Chr_Ref_Pos_Alt
  # 'cohortSamples' is a list of all the samples in your cohort. It probably is a variable called "samples_xxx"
  # 'allCounts' is a matrix prepared and filled in the previous step. It has all the alleleCounter information loaded.
  # 'normalSamples' is a list of all samples in the normal panel  


  #A list of mutations in patient subdirectory (format: Chr_Ref_Pos_Alt)
  Muts_patient = mutCodes
  #A file of with the sample names belonging to this patient
  case_samples=grep(patient,cohortSamples,value = T)
  
  #Select all normal samples not belonging to this patient
  normal_panel = normal_samples[!normal_samples%in%samples&!normal_samples%in%case_samples]
  norm_all_counts = all_counts[normal_panel,,]
  
  #Set up pval matrix
  pval_mat = matrix(1,ncol=length(case_samples),nrow=length(Muts_patient))
  rownames(pval_mat)=Muts_patient
  colnames(pval_mat)=case_samples
  
  coords_proj = substr(Muts_patient,1,nchar(Muts_patient)-4)
  Alt=substr(Muts_patient,nchar(Muts_patient),nchar(Muts_patient))
  Ref=substr(Muts_patient,nchar(Muts_patient)-2,nchar(Muts_patient)-2)
  
  for (s in case_samples){
    rho_est=rep(NA,length(coords_proj))
    test_counts = all_counts[s,coords_proj,]
    for (k in 1:length(coords_proj)) {
      n = sum(test_counts[coords_proj[k],])
      x = test_counts[coords_proj[k],Alt[k]]
      
      N_indiv = rowSums(norm_all_counts[,coords_proj[k],])
      X_indiv = norm_all_counts[,coords_proj[k],c("A","C","G","T")!=Ref[k]]
      pseudo = .Machine$double.eps    
      N=sum(N_indiv)
      X=sum(X_indiv)
      
      mu = max(X,pseudo)/max(N,pseudo)
      counts = cbind(N,X)
      if(is.null(rho)) rho = estimateRho_gridml(counts,mu)
      rdisp = (1 - rho)/rho
      
      prob0 = (X + x)/(N + n); prob0[prob0==0] = pseudo
      prob1s = x/(n+pseudo); prob1s[prob1s==0] = pseudo
      prob1c = X/(N+pseudo); prob1c[prob1c==0] = pseudo
      
      prob1s = pmax(prob1s,prob1c) # Min error rate is that of the population (one-sided test)
      nu0 = prob0 * rdisp; nu1s = prob1s * rdisp; nu1c = prob1c * rdisp; 
      
      # Likelihood-Ratio Tests
      LL = logbb(x, n, nu0, rdisp) + logbb(X, N, nu0, rdisp) - logbb(x, n, nu1s, rdisp) - logbb(X, N, nu1c, rdisp)
      pvals = pchisq(-2*LL, df=1, lower.tail=F)/2 # We divide by 2 as we are performing a 1-sided test
      # Saving the result
      pval_mat[k,s] = pvals
    } 
  }
  if(is.null(save)){
    return(pval_mat)
  }else{
    write.table(pval_mat,save)
  }
}
