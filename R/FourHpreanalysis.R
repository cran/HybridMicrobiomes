FourHpreanalysis <- function(x, class_grouping, core_fraction, boot_no, sample_no, core_average_abundance = 0, core_max_abundance = 0, core_all_average_abundance = 0, core_all_max_abundance = 0, replace_hosts=FALSE, reads=NULL, rarefy_each_step=TRUE,seed=NULL) {

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
  core_fraction_P1<-c()
  core_fraction_H<-c()
  core_fraction_P2<-c()

  #Define the seed
  if (is.null(seed)==TRUE){
    seed<-sample(10000,1)
  }

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
    P1sub<-sample(1:dim(OTU1P1)[2],sample_no,replace=replace_hosts)
    Hsub<-sample(1:dim(OTU1H)[2],sample_no,replace=replace_hosts)
    P2sub<-sample(1:dim(OTU1P2)[2],sample_no,replace=replace_hosts)

    #Subsample matrices for this particular bootstrap
    OTU1P1sub<-OTU1P1[,P1sub]
    OTU1Hsub<-OTU1H[,Hsub]
    OTU1P2sub<-OTU1P2[,P2sub]
    OTU1sub<-cbind(OTU1P1,OTU1H,OTU1P2)

    #Find the core microbiome for each host class
    coreP1<-intersect(intersect(intersect(intersect(which(rowSums(sign(OTU1P1sub))>=round(corefactor*dim(OTU1P1sub)[2])), which(rowSums(OTU1P1sub)/(dim(OTU1P1sub)[2]*colSums(OTU1P1sub)[1])>=core_average_abundance)), which(apply(OTU1P1sub,1,max)/colSums(OTU1P1sub)[1]>=core_max_abundance)), which(rowSums(OTU1sub)/(dim(OTU1sub)[2]*colSums(OTU1sub)[1])>=core_all_average_abundance)), which(apply(OTU1sub,1,max)/colSums(OTU1sub)[1]>=core_all_max_abundance))
    coreH<-intersect(intersect(intersect(intersect(which(rowSums(sign(OTU1Hsub))>=round(corefactor*dim(OTU1Hsub)[2])), which(rowSums(OTU1Hsub)/(dim(OTU1Hsub)[2]*colSums(OTU1Hsub)[1])>=core_average_abundance)), which(apply(OTU1Hsub,1,max)/colSums(OTU1Hsub)[1]>=core_max_abundance)), which(rowSums(OTU1sub)/(dim(OTU1sub)[2]*colSums(OTU1sub)[1])>=core_all_average_abundance)), which(apply(OTU1sub,1,max)/colSums(OTU1sub)[1]>=core_all_max_abundance))
    coreP2<-intersect(intersect(intersect(intersect(which(rowSums(sign(OTU1P2sub))>=round(corefactor*dim(OTU1P2sub)[2])), which(rowSums(OTU1P2sub)/(dim(OTU1P2sub)[2]*colSums(OTU1P2sub)[1])>=core_average_abundance)), which(apply(OTU1P2sub,1,max)/colSums(OTU1P2sub)[1]>=core_max_abundance)), which(rowSums(OTU1sub)/(dim(OTU1sub)[2]*colSums(OTU1sub)[1])>=core_all_average_abundance)), which(apply(OTU1sub,1,max)/colSums(OTU1sub)[1]>=core_all_max_abundance))

    #Find the overall gamma richness for each host class
    allP1<-which(rowSums(sign(OTU1P1sub))>=0)
    allH<-which(rowSums(sign(OTU1Hsub))>=0)
    allP2<-which(rowSums(sign(OTU1P2sub))>=0)

    #Find the fraction of the overall gamma richness that ends up in the core of each host class
    core_fraction_P1<-c(core_fraction_P1,100*length(coreP1)/length(allP1))
    core_fraction_H<-c(core_fraction_H,100*length(coreH)/length(allH))
    core_fraction_P2<-c(core_fraction_P2,100*length(coreP2)/length(allP2))

  }


  #Make a data frame of fractions
  fractiondf<-data.frame(core_fraction_P1,core_fraction_H,core_fraction_P2)

  #If the median core fraction of the hybrid is less than 25% of the average median core fraction of the progenitors, print a warning
  if (median(core_fraction_H)<0.25*(median(core_fraction_P1)+median(core_fraction_P2))){
    print('Warning: Hybrid core microbiota represents a much smaller fraction of the hybrid y-diversity. Proceed with caution and consider dividing the hybrid population into several groups based on genetics or environments.')
  }

  return(fractiondf)
}
