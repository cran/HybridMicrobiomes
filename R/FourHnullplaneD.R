FourHnullplaneD <- function(boots) {

    #Collect each of the dimensions of the 4H index, and sigma
    union<-boots[,2]
    intersection<-boots[,1]
    gain<-boots[,3]
    loss<-boots[,4]
    sigma<-boots[,5]

    #For each bootstrap, calculate theta
    theta<-gain+loss
    #Use theta to calculate the expected union fraction for each bootstrap
    union_expectation<-(1-sigma)*(1-theta)
    #Use theta to calculate the expected intersection fraction for each bootstrap
    intersection_expectation<-sigma*(1-theta)

    diffI<-intersection-intersection_expectation
    diffU<-union - union_expectation

    #Make a histogram of the difference between the true intersection dimension and the expected intersection dimension
    #This will be positive when the true intersection is more than the expected intersection (hybrids are more likely to retain microbes shared by both progenitors)
    graphics::hist(intersection-intersection_expectation)
    #Mean intersection preference:
    mean_intersection_preference<-mean(intersection-intersection_expectation)
    #Standard deviation in intersection preference:
    sd_intersection_preference<-stats::sd(intersection-intersection_expectation)
    #p-value (fraction of bootstraps where the true intersection was LESS than the expected intersection):
    p_value_intersection_preference<-length(which(intersection-intersection_expectation<0))/length(union)

    expectationdf<-data.frame(sigma,theta,union,union_expectation,intersection,intersection_expectation,diffI,diffU)
    return(list("trials"=expectationdf,"mean_intersection_preference"=mean_intersection_preference,"standard_deviation_intersection_preference"=sd_intersection_preference,"p_value_intersection_preference"=p_value_intersection_preference))
}
