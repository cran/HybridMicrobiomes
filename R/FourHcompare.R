FourHcompare <- function(boots_sets, boots_types, method='ilr', permutations=1000) {

  rlang::englue("var: {{ method }}")
  rlang::englue("var: {{ permutations }}")


  #Transform four-dimensional metric based on compositional space
  if (method=='ilr'){
    input_points<-compositions::ilr(boots_sets[,1:4])
  }
  else if (method == 'clr'){
    input_points<-compositions::clr(boots_sets[,1:4])
  }
  else{
    message('unrecognized transformation method')
  }

  #Calculate distances between points
  distance_matrix<-dist(input_points)

  #Convert input data and distance matrix into a list for the PERMANOVA package
  inputpermanova<-list(Data=input_points,D=as.matrix(distance_matrix),Coefficient='Other')

  #Perform PERMANOVA
  pm<-PERMANOVA::PERMANOVA(inputpermanova,factor(boots_types),nperm=permutations)

  #Return output
  return(pm)

}
