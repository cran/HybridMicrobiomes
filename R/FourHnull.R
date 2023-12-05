FourHnull<-function(x,class_grouping,core_fraction,boot_no,sample_no,null_model=1,replace_hosts=FALSE,reads=NULL,rarefy_each_step=TRUE,seed=NULL, use_microViz = 'yes'){

  rlang::englue("var: {{ core_fraction }}")
  rlang::englue("var: {{ boot_no }}")
  rlang::englue("var: {{ sample_no }}")
  rlang::englue("var: {{ replace_hosts }}")

  #%%%%%%%%%%%%%%%%%%%%%%   READ IN PARAMETERS AND INITIALIZE VECTORS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

  if (null_model > 9){
    print('Null model does not conserve read depth and should only be used for presence/absence indices (e.g., Jaccard, UniFrac).')
  }

  if (use_microViz == 'yes'){
    #Check to make sure there is sample data and if not create a sample data dataframe with sample names
    x<-microViz::phyloseq_validate(x,verbose=FALSE)
  }

  #Define the seed
  if (is.null(seed)==TRUE){
    seed<-sample(10000,1)
  }

  if (rarefy_each_step==FALSE){
    x<-phyloseq::rarefy_even_depth(x,verbose=FALSE,sample.size=reads, rngseed = sample(1000000,1))
  }


  #Add the species code to the sample data
  phyloseq::sample_data(x)$species_code<-class_grouping

  #Set the seed
  set.seed(seed)

  #%%%%%%%%%%%%%%%%%%%%%%   NULL MODEL SIMULATIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  #For the selected number of trials...
  for (k in 1:boot_no){


    #Make a list of the individual animals you're going to randomly choose from the parent and hybrid populations
    picklist<-c()
    #Randomly select the appropriate number of animals from the progenitor 1 group
    pick<-sample(which(sample_data(x)$species_code==1),sample_no,replace=replace_hosts)
    pickmom<-pick
    #Add those animals to your list
    picklist<-c(picklist,pick)
    #Randomly select the appropriate number of animals from the progenitor 2 group
    pick<-sample(which(sample_data(x)$species_code==3),sample_no,replace=replace_hosts)
    pickdad<-pick
    #Add those animals to your list
    picklist<-c(picklist,pick)
    #Randomly select the appropriate number of animals from the hybrid group (you will write over these later)
    pick<-sample(which(sample_data(x)$species_code==2),sample_no,replace=replace_hosts)
    #Add those animals to your list
    picklist<-c(picklist,pick)
    #Record which spots are 'hybrid' (you will write your null model over these later)
    hybridspots<-pick

    #Find the columns corresponding to the P1 individuals
    moms<-which(sample_data(x)$species_code==1)

    #Find the columns corresponding to the P2 individuals
    dads<-which(sample_data(x)$species_code==3)

    #Find the columns corresponding to the hybrid individuals
    hybrids<-which(sample_data(x)$species_code==2)


    #%%%%%%%%%%%%%%%%%%%%%  CREATE NULL MODEL HYBRIDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #Select the P1 individuals for breeding (depending on null model you will use some or all of these)
    selector_mom<-sample(moms,sample_no,replace=replace_hosts)
    selector_mom2<-sample(moms,sample_no,replace=replace_hosts)

    #Select the P2 individuals for breeding (depending on null model you will use some or all of these)
    selector_dad<-sample(dads,sample_no,replace=replace_hosts)
    selector_dad2<-sample(dads,sample_no,replace=replace_hosts)

    #Select the H individuals for breeding (depending on null model you will use some or all of these)
    selector_hybrid<-sample(hybrids,sample_no,replace=replace_hosts)
    selector_hybrid2<-sample(hybrids,sample_no,replace=replace_hosts)

    if (rarefy_each_step==TRUE){
      #Rarefy the OTU table specifically for this bootstrap sample
      tempx<-phyloseq::rarefy_even_depth(x,verbose=FALSE,sample.size=reads, rngseed = sample(1000000,1))
      #Remove taxa that drop out of the rarefied OTU table
      tempx<-phyloseq::prune_taxa(row.names(tax_table(tempx)),tempx)


    }else{
      tempx<-x
    }

    #Write the phyloseq object to an OTU table that you can directly pull abundances out of
    OTU1<-otu_table(tempx)
    #Breed individuals for null models, replace the actual 'hybrid' microbial abundances with the 'null model' hybrid abundances
    for (k in 1:sample_no){
      #Add one mom and one dad microbiome readcounts and divide by 2
      if (null_model==1){
        new_abundances<-round(OTU1[,selector_mom[k]]+OTU1[,selector_dad[k]])/2
        otu_table(tempx)[,hybridspots[k]]<-new_abundances
      }
      else if (null_model==10){
        setzero1<-sample(intersect(which(OTU1[,selector_mom[k]]>0), which(OTU1[,selector_dad[k]]==0)),round(length(intersect(which(OTU1[,selector_mom[k]]>0), which(OTU1[,selector_dad[k]]==0)))/2))
        setzero2<-sample(intersect(which(OTU1[,selector_dad[k]]>0), which(OTU1[,selector_mom[k]]==0)),round(length(intersect(which(OTU1[,selector_dad[k]]>0), which(OTU1[,selector_mom[k]]==0)))/2))
        OTU1temp<-OTU1
        OTU1temp[setzero1,selector_mom[k]]<-0
        OTU1temp[setzero2,selector_dad[k]]<-0
        new_abundances<-sign(OTU1temp[,selector_mom[k]]+OTU1temp[,selector_dad[k]])
        otu_table(tempx)[,hybridspots[k]]<-new_abundances
      }
      #Add two mom and one dad microbiome readcounts and divide by 3 (useful for polyploids)
      else if (null_model==2){
        new_abundances<-round(OTU1[,selector_mom[k]]+OTU1[,selector_dad[k]]+OTU1[,selector_mom2[k]])/3
        otu_table(tempx)[,hybridspots[k]]<-new_abundances
      }
      #else if (null_model==20){
   #}
      #Add one mom and two dad microbiome readcounts and divide by 3 (useful for polyploids)
      else if (null_model==3){
        new_abundances<-round(OTU1[,selector_mom[k]]+OTU1[,selector_dad[k]]+OTU1[,selector_dad2[k]])/3
        otu_table(tempx)[,hybridspots[k]]<-new_abundances
      }
      #else if (null_model==30){
      #  otu_table(tempx)[,hybridspots[k]]<-new_abundances
      #}
      #Sample from and mix the mom population
      else if (null_model==4){
        new_abundances<-round(OTU1[,selector_mom[k]]+OTU1[,selector_mom2[k]])/2
        otu_table(tempx)[,hybridspots[k]]<-new_abundances
      }
      #else if (null_model==40){
       #}
      #Sample from and mix the dad population
      else if (null_model == 5) {
        new_abundances<-round(OTU1[,selector_dad[k]]+OTU1[,selector_dad2[k]])/2
        otu_table(tempx)[,hybridspots[k]]<-new_abundances
      }
      #else if (null_model==50){
      #}
      #Sample from the mom population without mixing
      else if (null_model == 6){
        new_abundances<-(OTU1[,selector_mom[k]])
        otu_table(tempx)[,hybridspots[k]]<-new_abundances
      }
      #Sample from the dad population without mixing
      else if (null_model == 7){
        new_abundances<-(OTU1[,selector_dad[k]])
        otu_table(tempx)[,hybridspots[k]]<-new_abundances
      }
      #No null model, sample hybrids directly
      else{
        print('Null model not recognized; randomly sampling hybrids...')
        new_abundances<-(OTU1[,selector_hybrid[k]])
        otu_table(tempx)[,hybridspots[k]]<-new_abundances
      }
    }

    #pull the selected individuals from the rarefied phyloseq object with null model replacements
    OTU1_downsample<-phyloseq::prune_samples(c(sample_names(tempx)[picklist]),tempx)
    OTU1_downsample<-phyloseq::prune_taxa(rowSums(otu_table(OTU1_downsample))>0,OTU1_downsample)
    class_grouping_downsample<-sample_data(OTU1_downsample)$species_code

    #Convert the OTU table to a matrix
    OTU1 = as(phyloseq::otu_table(OTU1_downsample), "matrix")

    #Create matrices for each host class (progenitor 1 == 1, hybrid == 2, progenitor 2 == 3)
    OTU1P1<-OTU1[,which(class_grouping_downsample==1)]
    OTU1H<-OTU1[,which(class_grouping_downsample==2)]
    OTU1P2<-OTU1[,which(class_grouping_downsample==3)]

    #Subsample from the matrices for each host class
    P1sub<-sample(1:dim(OTU1P1)[2],sample_no,replace=replace_hosts)
    Hsub<-sample(1:dim(OTU1H)[2],sample_no,replace=replace_hosts)
    P2sub<-sample(1:dim(OTU1P2)[2],sample_no,replace=replace_hosts)

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
