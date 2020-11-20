(WD <- getwd())
if (!is.null(WD)) setwd(WD)
name <- "Dane_S2_200_20.csv"
df <- read.csv(name, sep = ";", row.names = "Zadanie")

library(ecr) # package for PMX i OX crossovering
library(rlist)

N = nrow(df) # number of tasks
NM = ncol(df) # number of machines 
POP = 10 # number of "persons" in population
k = 10 # number of iterations
p.mut <<- 0.05 # probability of the mutation
p.recomb <<- 0.75 # probability of the recombination

ListFromListOfList <- function(Population){
  temp <- list()
  for(i in 1:length(Population)){
    temp <- list.append(temp,Population[[i]][[1]],Population[[i]][[2]])
  }
  return(temp)
} # transform list of lists to single list

Fitness <- function(.vector){
  temp <- rep(0,NM) # how much time each machine will work with given schedule
  for (i in 1:N){
    ind <- .vector[i] # number of current task
    for(j in 1:NM){
      if (j!=1) { # for other machines
        if (temp[j] < temp[j-1]){ # calculating "makespan"
          temp[j] <- temp[j-1] + df[[j]][[ind]] # starting time of the next task on the machine [j]
        }else {
          temp[j] <- temp[j] + df[[j]][[ind]] # starting time of the next task on the machine [j]
        }
      }else{ # for the first machine
        temp[j] <- temp[j] + df[[j]][[ind]] # starting time of the next task on the first machine
      }
    }
  }
  return(max(temp)) # we are interested to accomplish all tasks, so we choose max of time 
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
Swap <- function(vectorM.,vector.){ 
  for (j in 1:N){
    if (vectorM.[j] == vector.[1]){
      Find <- j
    }
    if (vectorM.[j] == vector.[2]){
      Sind <- j
    }
  }
  temp <- vectorM.[Find]
  vectorM.[Find] <- vectorM.[Sind]
  vectorM.[Sind] <- temp
  return(vectorM.)
} # function that swap 2 random tasks

Mutation <- function(Population){
  for(i in 1:length(Population)){
    .rand <- runif(1,0,1)
    if(.rand < p.mut){ # check probability of the mutation
      tempR <- sample(1:N,2,replace = F) # generate two random chromosomes
      Population[[i]] <- Swap(Population[[i]],tempR) # swap them (mutate)
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
    .crossPMX <- ListFromListOfList(.crossPMX)
    temp <- Mutation(.crossPMX)
    # .crossOX <- OX(Population) # rekombinacja typu OX
    # .crossOX <- ListFromListOfList(.crossOX)
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
        cat("Currently the best = ", curMIN,"\n", file = "output.txt", append = TRUE)
        FinalSchedule <- BestScheduling(.new) # save schedule
      }else{
        cat("No better solution for the generation = ", i, "\n", file = "output.txt", append = TRUE)
      }
    }else{
      if(i >= (k*0.2)){ # the point at which we began to change the probability of mutation and recombination
        p.mut <<- p.mut + 0.005
        p.recomb <<- p.recomb - 0.005
      }
      .new <- Genetic(.new)
      tempMIN <- BestInPopulation(.new)
      dE <- tempMIN - curMIN
      if(dE < 0){
        curMIN <<- tempMIN
        cat("Currently the best = ", curMIN,"\n", file = "output.txt", append = TRUE)
        FinalSchedule <- BestScheduling(.new)
      }else{
        cat("No better solution for the generation = ", i, "\n", file = "output.txt", append = TRUE)
      }
    }
  }
  return(FinalSchedule)
}

cat("Current file = ", name, "\n", file = "output.txt", append = TRUE)
cat("Output for the number of iterations = ",k, "\n", file = "output.txt", append = TRUE)
cat("Output for the number of persons in population = ",POP, "\n", file = "output.txt", append = TRUE)
FirstPopulation <- FirstGenerationOFPopulation() # generate first population
curMIN <<- BestInPopulation(FirstPopulation) # the best schedule time for the first population
cat("Optimal schedule time for the start population = ",curMIN, "\n", file = "output.txt", append = TRUE)
.go <- main(FirstPopulation) # run program
cat("Optimal schedule", "\n", file = "output.txt", append = TRUE)
capture.output(.go, file = "output.txt", append = TRUE)
