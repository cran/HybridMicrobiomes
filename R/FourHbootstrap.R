FourHbootstrap <- function(x, class_grouping, core_fraction, boot_no, sample_no, core_average_abundance = 0, core_max_abundance = 0, core_all_average_abundance = 0, core_all_max_abundance = 0, replace_hosts=FALSE, reads=NULL, rarefy_each_step=TRUE,seed=NULL,dist='Jaccard') {

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
  fraction_OR_P1<-c()
  fraction_OR_P2<-c()

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

    #Find the names of the core microbial taxa for each host class
    namesP1<-rownames(OTU1P1sub)
    namesH<-rownames(OTU1Hsub)
    namesP2<-rownames(OTU1P2sub)

    corenamesP1<-namesP1[coreP1]
    corenamesH<-namesH[coreH]
    corenamesP2<-namesP2[coreP2]

    #If you are using the Jaccard inspired 4H Index...
    if (dist == 'Jaccard'){
      #Intersection is the intersection between P1, P2 and H, divided by the total number of unique core microbial taxa across all three host classes
      fraction_AND<-c(fraction_AND,length(intersect(intersect(corenamesP1,corenamesP2),corenamesH))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
      #Union is the intersection between H and the microbes on P1 and P2 but not both, divided by the total number of unique core microbial taxa. To find the microbes on P1 and P2 but not both, we find the total list of unique microbes on P1 and P2 and remove those microbes shared by P1 and P2
      fraction_OR<-c(fraction_OR,(length(intersect(setdiff(unique(c(corenamesP1,corenamesP2)),intersect(corenamesP1,corenamesP2)),corenamesH)))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
      #Gain is microbes on H that are not found on either P1 or P2, divided by the total number of unique core microbial taxa. To find the microbes found on H but not P1 or P2, we find the total number of microbes on H, and subtract those that intersect with the list of unique microbes found on P1 and P2
      fraction_OUT<-c(fraction_OUT,(length(corenamesH)-length(intersect(unique(c(corenamesP1,corenamesP2)),corenamesH)))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
      #Loss is the microbes found on P1, P2 or both that are not found on H, divided by the total number of unique core microbial taxa across all three host classes.
      fraction_MISS<-c(fraction_MISS,length(setdiff(unique(c(corenamesP1,corenamesP2)),corenamesH))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
      #sigma is the intersection of P1 and P2 divided by the total number of unique microbes found on P1 and P2
      fraction_sigma<-c(fraction_sigma,length(intersect(corenamesP1,corenamesP2))/length(unique(c(corenamesP1,corenamesP2))))
      fraction_OR_P1<-c(fraction_OR_P1,(length(intersect(setdiff(unique(c(corenamesP1)),intersect(corenamesP1,corenamesP2)),corenamesH)))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
      fraction_OR_P2<-c(fraction_OR_P2,(length(intersect(setdiff(unique(c(corenamesP2)),intersect(corenamesP1,corenamesP2)),corenamesH)))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
    }
    else if (dist == 'Sorensen'){
      #Intersection is 3x the intersection between P1, P2 and H, divided by the total number of core microbial taxa on each host, summed across all three host classes
      fraction_AND<-c(fraction_AND,3*length(intersect(intersect(corenamesP1,corenamesP2),corenamesH))/length((c(corenamesH,corenamesP1,corenamesP2))))
      #Union is 2x the intersection between H and the microbes on P1 and P2 but not both, divided by the total number core microbial taxa on each host, summed across all three host classes. To find the microbes on P1 and P2 but not both, we find the total list of unique microbes on P1 and P2 and remove those microbes shared by P1 and P2
      fraction_OR<-c(fraction_OR,(2*length(intersect(setdiff(unique(c(corenamesP1,corenamesP2)),intersect(corenamesP1,corenamesP2)),corenamesH)))/length((c(corenamesH,corenamesP1,corenamesP2))))
      #Gain is microbes on H that are not found on either P1 or P2, divided by the total number core microbial taxa on each host, summed across all three host classes. To find the microbes found on H but not P1 or P2, we find the total number of microbes on H, and subtract those that intersect with the list of unique microbes found on P1 and P2
      fraction_OUT<-c(fraction_OUT,(length(corenamesH)-length(intersect(unique(c(corenamesP1,corenamesP2)),corenamesH)))/length((c(corenamesH,corenamesP1,corenamesP2))))
      #Loss is 1x the microbes found on P1 and/or P2 that are not found on H plus another 1x the microbes found on both P1 and P2, divided by the total number core microbial taxa on each host, summed across all three host classes.
      fraction_MISS<-c(fraction_MISS,(length(setdiff(unique(c(corenamesP1,corenamesP2)),corenamesH))+length(setdiff(intersect(corenamesP1,corenamesP2),corenamesH)))/length((c(corenamesH,corenamesP1,corenamesP2))))
      #sigma is 2x the intersection of P1 and P2 divided by the total number of microbes found on P1 and P2 summed across P1 and P2
      fraction_sigma<-c(fraction_sigma,2*length(intersect(corenamesP1,corenamesP2))/length((c(corenamesP1,corenamesP2))))
      fraction_OR_P1<-c(fraction_OR_P1,(2*length(intersect(setdiff(unique(c(corenamesP1)),intersect(corenamesP1,corenamesP2)),corenamesH)))/length((c(corenamesH,corenamesP1,corenamesP2))))
      fraction_OR_P2<-c(fraction_OR_P2,(2*length(intersect(setdiff(unique(c(corenamesP2)),intersect(corenamesP1,corenamesP2)),corenamesH)))/length((c(corenamesH,corenamesP1,corenamesP2))))
    }
    else{
      print('Unrecognized model inspiration: Please use Jaccard or Sorensen')
    }

  }


  #Make a data frame of fractions
  fractiondf<-data.frame(fraction_AND,fraction_OR,fraction_OUT,fraction_MISS,fraction_sigma,fraction_OR_P1,fraction_OR_P2)


  return(fractiondf)
}
