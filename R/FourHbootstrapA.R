FourHbootstrapA <- function(x, class_grouping, core_fraction, boot_no, sample_no, core_average_abundance = 0, core_max_abundance = 0, core_all_average_abundance = 0, core_all_max_abundance = 0, replace_hosts=FALSE, reads=NULL, rarefy_each_step=TRUE,seed=NULL,dist = 'Bray-Curtis',representative = 'mean',rescale_core = FALSE, use_microViz = 'yes') {

  rlang::englue("var: {{ core_fraction }}")
  rlang::englue("var: {{ boot_no }}")
  rlang::englue("var: {{ sample_no }}")
  rlang::englue("var: {{ replace_hosts }}")


  if (use_microViz == 'yes'){
  #Check to make sure there is sample data and if not create a sample data dataframe with sample names
    x<-microViz::phyloseq_validate(x,verbose=FALSE)
  }


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

  #Put the class grouping into the phyloseq object
  sample_data(x)$classes<-class_grouping

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


    #Make a phyloseq object for each host class
    rared1<-phyloseq::prune_samples(sample_data(rared)$classes == 1,rared)
    rared2<-phyloseq::prune_samples(sample_data(rared)$classes == 2,rared)
    rared3<-phyloseq::prune_samples(sample_data(rared)$classes == 3,rared)

    #Pick the individuals from each host class for this particular bootstrap
    pick1<-phyloseq::sample_names(rared1)[sample(1:length(sample_names(rared1)),sample_no,replace=replace_hosts)]
    pick2<-phyloseq::sample_names(rared2)[sample(1:length(sample_names(rared2)),sample_no,replace=replace_hosts)]
    pick3<-phyloseq::sample_names(rared3)[sample(1:length(sample_names(rared3)),sample_no,replace=replace_hosts)]

    #Make phyloseq objects for each host class for this particular bootstrap
    otu1<-phyloseq::prune_samples(pick1,rared1)
    otu2<-phyloseq::prune_samples(pick2,rared2)
    otu3<-phyloseq::prune_samples(pick3,rared3)

    #Find the names of the core microbial taxa for each host class for this particular bootstrap
    core1<-rownames(phyloseq::otu_table(otu1))[which(rowSums(sign(phyloseq::otu_table(otu1)))>=round(corefactor*sample_no))]
    core2<-rownames(phyloseq::otu_table(otu2))[which(rowSums(sign(phyloseq::otu_table(otu2)))>=round(corefactor*sample_no))]
    core3<-rownames(phyloseq::otu_table(otu3))[which(rowSums(sign(phyloseq::otu_table(otu3)))>=round(corefactor*sample_no))]

    #Provided a host class has at least one core microbial taxon, calculate the average or median number of reads associated with each core taxon for each host class
    #Progenitor 1
    if (length(core1)>0){
      coreotu1<-phyloseq::prune_taxa(core1,otu1)
      if (representative == 'mean'){
        mean1<-rowMeans(phyloseq::otu_table(coreotu1))
      }
      else if (representative == 'median'){
        mean1<-apply(phyloseq::otu_table(coreotu1),1,median,na.rm=TRUE)/sum(apply(phyloseq::otu_table(coreotu1),1,median,na.rm=TRUE))
      }
      else{
        print('Method for finding representative microbial abundances not recognized: Please use either mean or median')
      }
    }
    #Hybrid
    if (length(core2)>0){
      coreotu2<-phyloseq::prune_taxa(core2,otu2)
      if (representative == 'mean'){
        mean2<-rowMeans(phyloseq::otu_table(coreotu2))
      }
      else if (representative == 'median'){
        mean1<-apply(phyloseq::otu_table(coreotu2),1,median,na.rm=TRUE)/sum(apply(phyloseq::otu_table(coreotu2),1,median,na.rm=TRUE))
      }
      else{
        print('Method for finding representative microbial abundances not recognized: Please use either mean or median')
      }
    }
    #Progenitor 2
    if (length(core3)>0){
      coreotu3<-phyloseq::prune_taxa(core3,otu3)
      if (representative == 'mean'){
        mean3<-rowMeans(phyloseq::otu_table(coreotu3))
      }
      else if (representative == 'median'){
        mean1<-apply(phyloseq::otu_table(coreotu3),1,median,na.rm=TRUE)/sum(apply(phyloseq::otu_table(coreotu3),1,median,na.rm=TRUE))
      }
      else{
        print('Method for finding representative microbial abundances not recognized: Please use either mean or median')
      }
    }

    #If you are going to rescale the core microbiome reads to total one, divide the average number of reads for each taxon by the total average number of reads across all taxa in the core
    if (rescale_core == TRUE){
      if (length(core1)>0){
        mean1<-mean1/sum(mean1)
      }
      if (length(core2)>0){
        mean2<-mean2/sum(mean2)
      }
      if (length(core3)>0){
        mean3<-mean3/sum(mean3)
      }

    }

    #If there is no core microbiome on any host class, set all parameters to zero
    if (length(core1) == 0 && length(core2) == 0 && length(core3) == 0){
      print('no core microbiome on any class')
      n_intersection<-0
      n_union<-0
      n_union1<-0
      n_union2<-0
      n_gain<-0
      n_loss<-0
    }

    #As long as at least one core microbial taxon on at least one host class...
    else{

      #If progenitor 1 has no core microbial taxa, arbitrarly pick a host class that DOES have a core and use those microbial taxon names... but set their abundances to zero
      if (length(core1)==0){
        if (length(core2)>0){
          mean1<-mean2
          mean1[1:length(mean1)]<-0
        }
        else{
          mean1<-mean3
          mean1[1:length(mean1)]<-0
        }
      }

      #If the hybrid has no core microbial taxa, arbitrarly pick a host class that DOES have a core and use those microbial taxon names... but set their abundances to zero
      if (length(core2)==0){
        if (length(core1)>0){
          mean2<-mean1
          mean2[1:length(mean2)]<-0
        }
        else{
          mean2<-mean3
          mean2[1:length(mean2)]<-0
        }
      }

      #If progenitor 3 has no core microbial taxa, arbitrarly pick a host class that DOES have a core and use those microbial taxon names... but set their abundances to zero
      if (length(core3)==0){
        if (length(core2)>0){
          mean3<-mean2
          mean3[1:length(mean3)]<-0
        }
        else{
          mean3<-mean1
          mean3[1:length(mean3)]<-0
        }
      }


      #make a matrix with the core microbial taxa and abundances from progenitor 1 and the hybrid
      vv<-merge(mean1,mean2,by='row.names',all=TRUE)
      vv[which(is.na(vv[,2])),2]<-0
      vv[which(is.na(vv[,3])),3]<-0
      rownames(vv)<-vv[,1]
      vv<-vv[,2:3]
      #make a matrix with the microbial taxa and abundances from progenitor 1, the hybrid and progenitor 2
      ww<-merge(vv,mean3,by='row.names',all=TRUE)
      ww[which(is.na(ww[,2])),2]<-0
      ww[which(is.na(ww[,3])),3]<-0
      ww[which(is.na(ww[,4])),4]<-0
      rownames(ww)<-ww[,1]
      ww<-ww[,2:4]
      colnames(ww)<-c('parent1','hybrid','parent2')

      #Calculate the number of reads on the progenitors to find sigma (this is regardless of the model that these reads fall into)
      #reads on both progenitors double counting the reads that they share
      sigma_total<-ww[,1]+ww[,3]
      #Only the reads that the two progenitors share
      sigma_shared<-apply(cbind(ww[,1],ww[,3]),1,FUN=min)
      #create a matrix for subtracting the shared progenitor 1/progenitor 2 reads from progenitor 1 and progenitor 2
      sigma_unique_matrix<-cbind(ww[,1],ww[,3])-sigma_shared
      #Only the reads that are on one or other of the progenitors, but not both
      sigma_unique<-sigma_unique_matrix[,1]+sigma_unique_matrix[,2]


      #for each microbial taxon, find the minimum number of reads across each host class... these are the reads shared by all three host classes for core taxa
      intersection<-apply(ww, 1, FUN = min)

      #subtract the shared reads from the matrix
      ww<-ww-intersection

      #For the reads that remain, consider progenitor 1 and the hybrid... for each microbial taxon, find the minimum number of reads across these two host classes... these are the reads shared by progenitor 1 and the hybrid for core taxa
      union1<-apply(ww[,1:2],1,FUN=min)
      #create a matrix for subtracting the progenitor 1/hybrid reads from progenitor 1 and the hybrid (subtract nothing from progenitor 2)
      union1m<-cbind(union1,union1,rep(0,length(union1)))
      ww<-ww-union1m

      #For the reads that remain, consider progenitor 2 and the hybrid... for each microbial taxon, find the minimum number of reads across these two host classes... these are the reads shared by progenitor 2 and the hybrid for core taxa
      union2<-apply(ww[,2:3],1,FUN=min)
      #create a matrix for subtracting the progenitor 2/hybrid reads from progenitor 2 and the hybrid (subtract nothing from progenitor 1)
      union2m<-cbind(rep(0,length(union2)),union2,union2)
      ww<-ww-union2m

      #The full union model are those reads shared with EITHER progenitor
      union<-union1+union2

      #Any reads remaining on the hybrid are attributed to the gain model
      gain<-ww[,2]

      #Any reads remaining on EITHER progenitor or both are attributed to the loss model
      #Loss double counting reads on both progenitors
      loss<-(ww[,1]+ww[,3])
      #Only the shared loss reads found on both progenitors
      loss_shared<-apply(cbind(ww[,1],ww[,3]),1,FUN=min)
      #create a matrix for subtracting the shared progenitor 1/progenitor 2 reads from progenitor 1 and progenitor 2
      loss_unique_matrix<-cbind(ww[,1],ww[,3])-loss_shared
      #Only the loss reads found on one or other progenitor, but not both
      loss_unique<-loss_unique_matrix[,1]+loss_unique_matrix[,2]



      #sum up the total reads for each model
      if (dist == 'Bray-Curtis'){
        #multiply the intersection reads by 3 (since there was one copy on EACH of the three host classes)
        s_intersection<-3*sum(intersection)
        #multiply the union reads by 2 (since there was one copy on a progenitor and one copy on the host)
        s_union<-2*sum(union)
        s_union1<-2*sum(union1)
        s_union2<-2*sum(union2)
        #do not multiply gain or loss by anything, since these were summed reads on individuals (or, if on both progenitors, were counted already because all remaining reads on progenitors were summed)
        s_gain<-sum(gain)
        s_loss<-sum(loss)
        #for sigma, calculate the number of shared reads relative to the number of total reads across parents
        s_sigma<-2*sum(sigma_shared)/sum(sigma_total)
      }
      else if (dist == 'Ruzicka'){
        #multiply the intersection reads by 3 (since there was one copy on EACH of the three host classes)
        s_intersection<-sum(intersection)
        #multiply the union reads by 2 (since there was one copy on a progenitor and one copy on the host)
        s_union<-sum(union)
        s_union1<-sum(union1)
        s_union2<-sum(union2)
        #do not multiply gain or loss by anything, since these were summed reads on individuals (or, if on both progenitors, were counted already because all remaining reads on progenitors were summed)
        s_gain<-sum(gain)
        s_loss<-sum(loss_shared+loss_unique)
        #for sigma, calculate the number of shared reads relative to the number of total reads across parents
        s_sigma<-sum(sigma_shared)/sum(sigma_shared+sigma_unique)
      }
      else{
        print('Unrecognized model inspiration: Please use Bray-Curtis or Ruzicka')
      }

      #turn the reads for each model into a fraction of the total reads
      n_intersection<-s_intersection/(s_intersection+s_union+s_gain+s_loss)
      n_union<-s_union/(s_intersection+s_union+s_gain+s_loss)
      n_union1<-s_union1/(s_intersection+s_union+s_gain+s_loss)
      n_union2<-s_union2/(s_intersection+s_union+s_gain+s_loss)
      n_gain<-s_gain/(s_intersection+s_union+s_gain+s_loss)
      n_loss<-s_loss/(s_intersection+s_union+s_gain+s_loss)

    }


    #Find the fraction of each 'model' represented by the hybrid microbiomes
    fraction_AND<-c(fraction_AND,n_intersection)
    fraction_OR<-c(fraction_OR,n_union)
    fraction_OUT<-c(fraction_OUT,n_gain)
    fraction_MISS<-c(fraction_MISS,n_loss)
    fraction_sigma<-c(fraction_sigma,s_sigma)
    fraction_OR_P1<-c(fraction_OR_P1,n_union1)
    fraction_OR_P2<-c(fraction_OR_P2,n_union2)


  }


  #Make a data frame of fractions
  fractiondf<-data.frame(fraction_AND,fraction_OR,fraction_OUT,fraction_MISS,fraction_sigma,fraction_OR_P1,fraction_OR_P2)


  return(fractiondf)
}
