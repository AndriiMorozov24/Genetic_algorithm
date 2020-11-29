nmain <- function(.name,.k,.POP){
  (WD <- getwd())
  if (!is.null(WD)) setwd(WD)
  name <- .name
  df <- read.csv(name, sep = ";", row.names = "Cities")
  
  library(ecr) # package for PMX and OX crossovering
  library(rlist)
  
  N = nrow(df) # number of tasks
  POP = .POP # number of "persons" in population
  k = .k # number of iterations
  p.mut <<- 0.05 # probability of the mutation
  p.recomb <<- 0.8 # probability of the recombination
  
  Fitness <- function(.vector){
    sum <- 0
    for(i in 1:N){
      if(i != N){
        sum <- sum + df[[.vector[i]]][.vector[i+1]]
      }else{
        sum <- sum + df[[.vector[i]]][.vector[1]]
      }
    }
    return(sum)
  } 
  
  FirstGenerationOFPopulation <- function(){
    temp <- list()
    for(i in 1:POP){
      temp[[i]] <- sample(1:N,N,replace = F) # generate a "pearson"
    }
    return(temp) # return population
  }
  
  TournamentSelection <- function(Population){
    final <- list()
    for(i in 1:(length(Population)/2)-1){
      tempR <- sample(1:length(Population), 2, replace = FALSE) # generate two random "pearson" from population
      candidate1 <- Fitness(Population[[tempR[1]]])
      candidate2 <- Fitness(Population[[tempR[2]]])
      if (candidate1 < candidate2){
        final <- list.append(final,Population[[tempR[1]]])
      }else {
        final <- list.append(final,Population[[tempR[2]]])
      }
      Population <- list.remove(Population,c(tempR[1],tempR[2]))
    }
    return(final)
  } # selection function
  #===============================================
  # SumOfFitness <- function(.listOfParents){
  #   sum = 0
  #   for(i in 1:length(.listOfParents)){
  #     temp <- Fitness(.listOfParents[[i]])
  #     sum <- sum + temp
  #   }
  #   return(sum)
  # }
  # 
  # RouletteP <- function(.listOfParents){
  #   sum <- SumOfFitness(.listOfParents)
  #   temp <- rep(0,POP*2)
  #   for(i in 1:POP*2){
  #     fi <- Fitness(.listOfParents[[i]])
  #     temp[i] <- fi / sum
  #   }
  #   return(temp)
  # }
  # 
  # RouletteSelection <- function(Population){
  #   final <- list()
  #   choi <- RouletteP(Population)
  #   for(i in 1:POP*2){
  #     if(choi[i] >= median(choi)){
  #       final <- list.append(final,Population[[i]])
  #     }
  #   }
  #   return(final)
  # }
  #==================================================
  # package "ecr" for PMX and OX recombinator
  PMX <- function(Population){
    .crossPMX <- list()
    temp1 <- list.subset(Population, c(1:(POP/2))) 
    temp2 <- list.subset(Population, c((1+(POP/2)):POP)) 
    for(i in 1:POP){
      tempR <- sample(1:(POP/2),2,replace = FALSE) # generate random parents
      temp <- list(temp1[[tempR[1]]],temp2[[tempR[2]]]) 
      .crossPMX[[i]] <- recPMX(temp) # recombination that giving two childs
    }
    return(.crossPMX) # return POP*2 kids
  } # Partially-Mapped-Crossover (PMX) recombinator. 
  
  OX <- function(Population){
    .crossOX <- list()
    temp1 <- list.subset(Population, c(1:(POP/2)))
    temp2 <- list.subset(Population, c((1+(POP/2)):POP))
    for(i in 1:POP){
      tempR <- sample(1:(POP/2),2,replace = FALSE)
      temp <- list(temp1[[tempR[1]]],temp2[[tempR[2]]])
      .crossOX[[i]] <- recOX(temp)
    }
    return(.crossOX)
  } # Ordered-Crossover (OX) recombinator.
  #==========================================
  Swap <- function(vector) {
    tempR <- sample(1:length(vector), 2, replace = FALSE)
    index <- which(vector %in% tempR)
    if((index[1]==1) && (index[2]==length(vector))){
      final <- rev(vector)
    }else if(index[1]==1){
      temp2 <- list.subset(vector, c((index[2]+1):length(vector)))
      .rev <- rev(list.subset(vector, c(index[1]:index[2])))
      final <- c(.rev,temp2)
    }else if(index[2]==length(vector)){
      temp1 <- list.subset(vector, c(1:(index[1]-1)))
      .rev <- rev(list.subset(vector, c(index[1]:index[2])))
      final <- c(temp1,.rev)
    }else{
      temp1 <- list.subset(vector, c(1:(index[1]-1)))
      temp2 <- list.subset(vector, c((index[2]+1):length(vector)))
      .rev <- rev(list.subset(vector, c(index[1]:index[2])))
      final <- c(temp1, .rev, temp2)
    }
    return(final)
  }
  
  Mutation <- function(Population){
    for(i in 1:length(Population)){
      .rand <- runif(1,0,1)
      if(.rand < p.mut){ # check probability of the mutation
        Population[[i]] <- Swap(Population[[i]]) # swap them (mutate)
      }else{
        next
      }
    }
    return(Population) # population after mutation
  }
  
  Genetic <- function(Population){
    .rand <- runif(1,0,1)
    final <- list()
    tempPOP <- Population 
    if(.rand < p.recomb){ # check probability of the recombination
      .crossPMX <- PMX(tempPOP) 
      .crossPMX <- list.flatten(.crossPMX)
      temp <- Mutation(.crossPMX)
      # .crossOX <- OX(Population) # rekombinacja typu OX
      # .crossOX <- list.flatten(.crossOX)
      # temp <- Mutation(.crossOX) #mutacja
      tempF <- TournamentSelection(temp) # select the best kids
      test <- c(Population,tempF) # add them to parents
      final <- TournamentSelection(test) # select size = POP best "persons"
    }else{ # if not recombination
      temp <- Mutation(tempPOP)
      test <- c(Population,temp)
      final <- TournamentSelection(test)
    }
    return(final) # return final population
  }
  
  BestInPopulation <- function(Population){
    temp <- rep(0,POP)
    for(i in 1:POP){
      temp[i] <- Fitness(Population[[i]])
    }
    return(min(temp))
  } # return the best "pearson" from given population
  
  BestScheduling <- function(Population){
    temp <- rep(0,POP)
    for(i in 1:POP){
      temp[i] <- Fitness(Population[[i]])
    }
    ind <- which.min(temp)
    return(Population[[ind]])
  } # return the best scheduling
  
  main <- function(.FirstPopulation){
    FinalSchedule <- NULL
    for(i in 1:k){ # generate generations = k, k = stop condition
      if(i == 1){ # for the first generation
        .new <- Genetic(.FirstPopulation) # run the genetic algorithm
        tempMIN <- BestInPopulation(.new) # the best schedule time from the returned population
        dE <- tempMIN - curMIN 
        if(dE < 0){ # if we found better scheduling time
          curMIN <<- tempMIN # override with the new optimal scheduling time
          cat("Currently the best = ", curMIN,"\n", file = "output_TSP.txt", append = TRUE)
          FinalSchedule <- BestScheduling(.new) # save schedule
        }else{
          cat("No better solution for the generation = ", i, "\n", file = "output_TSP.txt", append = TRUE)
        }
      }else{
        if(i >= (k*0.5)){ # the point at which we began to change the probability of mutation and recombination
          p.mut <<- p.mut + (0.75/(k*0.5))
          p.recomb <<- p.recomb - (0.65/(k*0.5))
        }
        .new <- Genetic(.new)
        tempMIN <- BestInPopulation(.new)
        dE <- tempMIN - curMIN
        if(dE < 0){
          curMIN <<- tempMIN
          cat("Currently the best = ", curMIN,"\n", file = "output_TSP.txt", append = TRUE)
          FinalSchedule <- BestScheduling(.new)
        }else{
          cat("No better solution for the generation = ", i, "\n", file = "output_TSP.txt", append = TRUE)
        }
      }
    }
    return(FinalSchedule)
  }
  
  cat("Current file = ", name, "\n", file = "output_TSP.txt", append = TRUE)
  cat("Output for the number of iterations = ",k, "\n", file = "output_TSP.txt", append = TRUE)
  cat("Output for the number of persons in population = ",POP, "\n", file = "output_TSP.txt", append = TRUE)
  FirstPopulation <- FirstGenerationOFPopulation() # generate first population
  curMIN <<- BestInPopulation(FirstPopulation) # the best schedule time for the first population
  cat("Optimal schedule time for the start population = ",curMIN, "\n", file = "output_TSP.txt", append = TRUE)
  .go <- main(FirstPopulation) # run program
  cat("Optimal schedule", "\n", file = "output_TSP.txt", append = TRUE)
  capture.output(.go, file = "output_TSP.txt", append = TRUE)
  print("end")
}

Names <- c("TSP_127.csv")
niter <- c(10000)
POP = 100

for(i in 1:length(Names)){
  for(j in 1:length(niter)){
    nmain(Names[i],niter[j],POP)
  }
}

