nmain <- function(.name,.k,.POP){
  (WD <- getwd())
  if (!is.null(WD)) setwd(WD)
  name <- .name
  df <- read.csv(name, sep = ";", row.names = "Zadanie")
  
  library(ecr) # package for PMX and OX crossovering
  library(rlist)
  
  N = nrow(df) # number of tasks
  NM = ncol(df) # number of machines
  POP = .POP # number of chromosomes in population
  k = .k # number of iterations
  p.mut <<- 0.1 # probability of the mutation
  p.recomb <<- 0.85 # probability of the recombination
  b <- (p.recomb-0.6)/(k*0.66)
  
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
      temp[[i]] <- sample(1:N,N,replace = F) # generate a chromosome
    }
    return(temp) # return population
  }
  
  TournamentSelection <- function(Population){
    final <- list()
    for(i in 1:(length(Population)/2)){
      tempR <- sample(1:length(Population), 2, replace = FALSE) # generate two random indexes of "chromosome" from population
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
  #=================================================
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
  
  hOX <- function(.parent1,.parent2,.off,.cut){
    count <- 0
    i <- .cut+1
    j <- 1
    choi <- T
    temp <- T
    w <- length(.parent1) - length(.off)
    repeat{
      if(.parent2[i] %in% .off){
        i <- i+1
        if(i>length(.parent1)){
          i <- 1
        }
      }else{
        if (choi == T){
          .off <- list.append(.off,.parent2[i])
          i <- i+1
          count <- count + 1
          if ((i>length(.parent1)) && (count != length(.parent1[(max(.cut)+1):length(.parent1)]))){
            i <- 1
            temp <- F
          }
          if(count == length(.parent1[(max(.cut)+1):length(.parent1)])){
            choi <- F
            if(temp == T){
              i <- 1
            }
          }
        }else{
          .off <- list.insert(.off, j, .parent2[i])
          j <- j+1
          i <- i+1
          count <- count + 1
        }
      }
      if (count == w){
        break
      }
    }
    return(.off)
  }
  
  hPMX <- function(.ind,.off,.cut1,.cut2){
    i <- 1
    repeat{
      if(i < min(.ind) || i > max(.ind)){
        repeat{
          if(.off[i] %in% .cut1){
            .off[i] <- .cut2[match(.off[i],.cut1)]
          }else{
            i <- i+1
            break
          }
        }
      }else{
        i <- i+1
      }
      if(i>length(.off)){
        break
      }
    }
    return(.off)
  }
  
  recOwn <- function(.population, choi){
    final <- list()
    i <- 1
    b <- length(.population)
    repeat{
      ind <- sample(1:length(.population),2,replace = F)
      parent1 <- .population[[ind[1]]]
      parent2 <- .population[[ind[2]]]
      .cut <- sample(2:(length(parent1)-1),2,replace = F)
      .off1 <- parent1[min(.cut):max(.cut)]
      .off2 <- parent2[min(.cut):max(.cut)]
      if(choi=="OX"){
        final[[i]] <- hOX(parent1,parent2,.off1,max(.cut))
        i<-i+1
        final[[i]] <- hOX(parent2,parent1,.off2,max(.cut))
        i<-i+1
        .population <- list.remove(.population,c(ind[1],ind[2]))
      }else{
        .toff1 <- c(parent2[1:(min(.cut)-1)],.off1,parent2[(max(.cut)+1):length(parent2)])
        .toff2 <- c(parent1[1:(min(.cut)-1)],.off2,parent1[(max(.cut)+1):length(parent2)])
        final[[i]] <- hPMX(.cut,.toff1,.off1,.off2)
        i<-i+1
        final[[i]] <- hPMX(.cut,.toff2,.off2,.off1)
        i<-i+1
        .population <- list.remove(.population,c(ind[1],ind[2]))
      }
      if(i>=b){
        break
      }
    }
    return(final)
  }
  # package "ecr" for PMX and OX recombinator
  recLib <- function(Population, choi){
    .cross <- list()
    temp1 <- list.subset(Population, c(1:(POP/4)))
    temp2 <- list.subset(Population, c((1+(POP/4)):(POP/2)))
    for(i in 1:(POP/4)){
      if(choi=="OX"){
        if(length(temp1)==1){
          temp <- list(temp1[[1]],temp2[[1]])
          .cross[[i]] <- recOX(temp)
        }else{
          tempR <- sample(1:length(temp1),2,replace = FALSE)
          temp <- list(temp1[[tempR[1]]],temp2[[tempR[2]]])
          .cross[[i]] <- recOX(temp)
          temp1 <- list.remove(temp1,tempR[1])
          temp2 <- list.remove(temp2,tempR[2])
        }
      }else{
        if(length(temp1)==1){
          temp <- list(temp1[[1]],temp2[[1]])
          .cross[[i]] <- recPMX(temp)
        }else{
          tempR <- sample(1:length(temp1),2,replace = FALSE)
          temp <- list(temp1[[tempR[1]]],temp2[[tempR[2]]])
          .cross[[i]] <- recPMX(temp)
          temp1 <- list.remove(temp1,tempR[1])
          temp2 <- list.remove(temp2,tempR[2])
        }
      }
    }
    return(.cross)
  }
  
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
  } # function that swap 2 random genes
  
  Mutation <- function(Population){
    for(i in 1:length(Population)){
      .rand <- runif(1,0,1)
      if(.rand < p.mut){ # check probability of the mutation
        tempR <- sample(1:N,2,replace = F) # generate two random indexes of gene to be swaped
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
    bestP <- TournamentSelection(Population) # selection best parents
    if(.rand < p.recomb){ # check probability of the recombination
      # kids <- recOwn(bestP,"OX")
      # kids <- Mutation(kids)
      
      kids <- recOwn(bestP,"PMX")
      kids <- Mutation(kids)
      
      # kids <- recLib(bestP,"PMX")
      # kids <- list.flatten(kids)
      # kids <- Mutation(kids)
      
      # kids <- recLib(bestP,"OX")
      # kids <- list.flatten(kids)
      # kids <- Mutation(kids)
      
      final <- c(bestP,kids) # add kids to parents
    }else{ # if not recombination
      kids <- Mutation(bestP)
      final <- c(bestP,kids)
    }
    return(final) # return final population
  }
  
  BestInPopulation <- function(Population){
    temp <- rep(0,POP)
    for(i in 1:POP){
      temp[i] <- Fitness(Population[[i]])
    }
    ind <- which.min(temp)
    tempFinalSchedule <<- Population[[ind]]
    return(min(temp))
  } # return the best "chromosome" from given population and save schedule for that chromosome
  
  main <- function(.FirstPopulation){
    FinalSchedule <<- NULL
    for(i in 1:k){ # generate generations = k, k = stop condition
      if(i == 1){ # for the first generation
        .new <- Genetic(.FirstPopulation) # run the genetic algorithm
        tempMIN <- BestInPopulation(.new) # the best schedule time from the returned population
        dE <- tempMIN - curMIN 
        if(dE < 0){ # if we found better scheduling time
          curMIN <<- tempMIN # override with the new optimal scheduling time
          cat("Currently the best = ", curMIN,"\n", file = "output_PMX.txt", append = TRUE)
          FinalSchedule <<- tempFinalSchedule
        }else{
          cat("No better solution for the generation = ", i, "\n", file = "output_PMX.txt", append = TRUE)
        }
      }else{
        if(i >= (k*0.33)){ # the point at which we began to change the probability of mutation and recombination
          #p.mut <<- p.mut + (0.09/(k*0.5))
          p.recomb <<- p.recomb - b
        }
        .new <- Genetic(.new)
        tempMIN <- BestInPopulation(.new)
        dE <- tempMIN - curMIN
        if(dE < 0){
          curMIN <<- tempMIN
          cat("Currently the best = ", curMIN,"\n", file = "output_PMX.txt", append = TRUE)
          FinalSchedule <<- tempFinalSchedule
        }else{
          cat("No better solution for the generation = ", i, "\n", file = "output_PMX.txt", append = TRUE)
        }
      }
    }
    return(FinalSchedule)
  }
  
  cat("Current file = ", name, "\n", file = "output_PMX.txt", append = TRUE)
  cat("Output for the own PMX recombinator", "\n", file = "output_PMX.txt", append = TRUE)
  cat("Output for the number of iterations = ",k, "\n", file = "output_PMX.txt", append = TRUE)
  cat("Output for the number of persons in population = ",POP, "\n", file = "output_PMX.txt", append = TRUE)
  FirstPopulation <- FirstGenerationOFPopulation() # generate first population
  curMIN <<- BestInPopulation(FirstPopulation) # the best schedule time for the first population
  cat("Optimal schedule time for the start population = ",curMIN, "\n", file = "output_PMX.txt", append = TRUE)
  .go <- main(FirstPopulation) # run program
  cat("Optimal schedule", "\n", file = "output_PMX.txt", append = TRUE)
  capture.output(.go, file = "output_PMX.txt", append = TRUE)
  print("end")
}

Names <- c("Dane_S2_50_10.csv","Dane_S2_100_20.csv","Dane_S2_200_20.csv")
#niter <- c(1000,2000,3000)
niter <- c(1000)
POP <- c(500,1000,2000)

for(i in 1:length(Names)){
  for(j in 1:length(niter)){
    for(x in 1:length(POP)){
      nmain(Names[i],niter[j],POP[x])
    }
  }
}
