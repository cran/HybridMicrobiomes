FourHbootstrap <- function(x, class_grouping, core_fraction, boot_no, sample_no, replace_hosts=FALSE, reads=NULL, rarefy_each_step=TRUE) {

  rlang::englue("var: {{ core_fraction }}")
  rlang::englue("var: {{ boot_no }}")
  rlang::englue("var: {{ sample_no }}")
  rlang::englue("var: {{ replace_hosts }}")

  #Define the number of reads
  if (is.null(reads)==TRUE){
    reads<-min(sample_sums(x))
  }

  #Define the fraction of hosts a microbe must be present on to be considered 'core'
  corefactor<-core_fraction

  #Define vectors for the four models and sigma from each trial
  fraction_AND<-c()
  fraction_OR<-c()
  fraction_OUT<-c()
  fraction_MISS<-c()
  fraction_sigma<-c()

  #Define the number of hosts in each bootstrap
  mintype<-sample_no

  #If you are not going to rarefy at each bootstrap, rarefy at the beginning (this uses a single rarefaction for all bootstrap samples)
  if (rarefy_each_step==FALSE){
    rared<-phyloseq::rarefy_even_depth(x,verbose=FALSE,sample.size=reads)
  }

  #For each bootstrap sample...
  for (m in 1:boot_no){

    if (rarefy_each_step==TRUE){
    #Rarefy the OTU table to the lowest read depth across all samples, with the option of specifying a lower read depth if you want to standardize across studies
    rared<-phyloseq::rarefy_even_depth(x,verbose=FALSE,sample.size=reads)
    }

    #Convert the OTU table to a matrix
    OTU1 = as(phyloseq::otu_table(rared), "matrix")

    #Create matrices for each host class (progenitor 1 == 1, hybrid == 2, progenitor 2 == 3)
    OTU1P1<-OTU1[,which(class_grouping==1)]
    OTU1H<-OTU1[,which(class_grouping==2)]
    OTU1P2<-OTU1[,which(class_grouping==3)]

    #Subsample from the matrices for each host class
    P1sub<-sample(1:dim(OTU1P1)[2],mintype,replace=replace_hosts)
    Hsub<-sample(1:dim(OTU1H)[2],mintype,replace=replace_hosts)
    P2sub<-sample(1:dim(OTU1P2)[2],mintype,replace=replace_hosts)

    #Subsample matrices for this particular bootstrap
    OTU1P1sub<-OTU1P1[,P1sub]
    OTU1Hsub<-OTU1H[,Hsub]
    OTU1P2sub<-OTU1P2[,P2sub]

    #Find the core microbiome for each host class
    coreP1<-which(rowSums(sign(OTU1P1sub))>=round(corefactor*dim(OTU1P1sub)[2]))
    coreH<-which(rowSums(sign(OTU1Hsub))>=round(corefactor*dim(OTU1Hsub)[2]))
    coreP2<-which(rowSums(sign(OTU1P2sub))>=round(corefactor*dim(OTU1P2sub)[2]))

    #Find the names of the core microbial taxa for each host class
    namesP1<-rownames(OTU1P1sub)
    namesH<-rownames(OTU1Hsub)
    namesP2<-rownames(OTU1P2sub)

    corenamesP1<-namesP1[coreP1]
    corenamesH<-namesH[coreH]
    corenamesP2<-namesP2[coreP2]


    #Find the fraction of each 'model' represented by the hybrid microbiomes
    fraction_AND<-c(fraction_AND,length(intersect(intersect(corenamesP1,corenamesP2),corenamesH))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
    fraction_OR<-c(fraction_OR,(length(intersect(setdiff(unique(c(corenamesP1,corenamesP2)),intersect(corenamesP1,corenamesP2)),corenamesH)))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
    fraction_OUT<-c(fraction_OUT,(length(corenamesH)-length(intersect(unique(c(corenamesP1,corenamesP2)),corenamesH)))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
    fraction_MISS<-c(fraction_MISS,length(setdiff(unique(c(corenamesP1,corenamesP2)),corenamesH))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
    fraction_sigma<-c(fraction_sigma,length(intersect(corenamesP1,corenamesP2))/length(unique(c(corenamesP1,corenamesP2))))


  }


  #Make a data frame of fractions
  fractiondf<-data.frame(fraction_AND,fraction_OR,fraction_OUT,fraction_MISS,fraction_sigma)


  return(fractiondf)
}
