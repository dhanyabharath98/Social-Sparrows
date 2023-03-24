#No plasticity
#Case1: No environmental fluctuations
#Case2: Env_fluct within a generation (a,f)
#Case3: Env_fluct only between generations (a,f)
#Case4: Env_fluct both within and between generations (a,f)

#a: environmental autocorrelation = 0
#b: environmental autocorrelation = 0.5
#c: environmental autocorrelation = 0.99

#clear the environment
rm(list=ls())
#load packages ggplot2
library(ggplot2)
library(MASS)
library(vioplot)
library(RColorBrewer)

## Structure of Population matrix (describing each row - columns are individuals) ##

# 1. LRN Elevation gene
# 2. LRN Slope gene
# 3. LRN Noise gene
# 4. Memory factor gene (alpha) {0,1}
# 5. Type of information gene [0,1], 0 = proportion of scroungers, 1 = payoffs
# 6. Strategy expressed current time step [0,1], 0 = scrounger, 1 = producer
# 7. Fraction of producing played

#############
############################################
################ Simulation ################
#HELP PAGE
############################################
### Function running the evolutionary simutation
## Arguments: ##
### Demographic and scenario setup parameters:
# Rounds      How many replicate simulations to run. Default 10.
# N           Adult population size. Default 300.
# S           Time steps per generation. Default 30.
# G           Number of generations. Default 300.
# mrate       Vector of mutation rates for the 5 genes. Default c(0.1,0.1,0.1,0,0.01).
#  -The order of the genes are: 1:Elev 2:Slope 3:Noise 4:Memory 5:InfoType.
# pcost       Cost of being plastic (fraction of intake). Scales with LRN Slope gene. Default 0.
# scost       Cost of switching tactic (fraction of intake)
# mcost       Cost of using memory  (fraction of intake). Scales with memory gene. Default 0.
# msize       Size of the mutational effect
# dem_stoch   Demographic stochasticity / background mortality (proportion of individuals randomly killed before selection). Default 0.
# sel_trunc   Truncation selection (the individuals below the sel_trunc quantile of food items get no fitness). Default 0.25.

### Game theoretical and environmental parameters:
# f           Food items per patch. Default 20.
# a           Producer's advantage. Default 6.
# p           Probability of finding food when producing. Default 0.25.
# env_fluct   Environmental fluctuations. Default NULL. Can contain:
#  -Parameters to be fluctuating: The above three ("f", "a", "p")
#  -Type of fluctuations: "random" or "autocorr": env_fluct_type
#  -Time scale of fluctuations: Among generations ("g") or among time steps ("s"): env_fluct_scale can be "g" or "s" or both g and s
# rw_par      Rescorla-Wagner learning rule shape parameter. 0 => flat (Info <- 0.5 for all Vp/(Vp+Vs)); 1 => linear (Info=Vp/(Vp+Vs)). Default NULL (not currently implemented).
## Outputs: ##
# Totrec, an array with dimensions [8, N, G, Rounds]
#  First dimension are the 7 loci, and the tactic expressed in the final time step. 1:Strategy type. 2:Elev. 3:Diff. 4:InflPt 5:Steepness 6:Memory 7:InfoType 8:CurrentTactic.
#  These are printed for each individual (2nd dimension) in each generation (3rd dimension) and each replicate (4th dimension)
# 

###############################
#Function autocorr: Create a random sequence of n values with an autocorrelation of rho#
autocorr <- function(rho,n){
  tmp.r <- matrix(rho,n,n) # Create an nxn matrix with values rho
  tmp.r <- tmp.r^abs(row(tmp.r)-col(tmp.r))
  x <- mvrnorm(1,rep(0,n),tmp.r) # Yes the first argument should be 1 and not n, weirdly.
  return(x)
}

#Simulation function#
sim1.1 <- 
  function(Rounds=10,N=300,S=30,G=300,mrate=c(0.1,0,0.1,0,0),
           msize=0.1,pcost=0,mcost=0,scost=0,tac=0.99,
           dem_stoch=0,sel_trunc=0.25,f=20,a=6,
           p=0.25,env_fluct=NULL,env_fluct_type="autocorr",env_fluct_scale=c("s","g"),rw_par=NULL){
    #initialising the final matrix 
    Totrec <- array(0,dim=c(7,N,G,Rounds))
    
    ### START REPLICATE ###
    
    for (z in 1:Rounds){
      
      ## Variables used for recording information about the population
      stepfood <- rep(0,N) #Vector for recording each individual's food gain each time step. Reset each time step, used to calculate evaluation of information
      Stratrec <- array(0,dim=c(7,N,G)) #An array (3D-matrix) that will contain all population matrices from all generations
      Offspring <- matrix(0,7,N)
      
      rownames(Stratrec) <- c("Elev","Slope","Noise","alpha","InfoType","CurrStrat","TotalProd")
      Stratrec[1,,1] <- rnorm(N,0,1) # Initiate elevation genes ~Normal(0,1)
      Stratrec[2,,1] <- 0 #rnorm(N,0,1) # Initiate slope genes ~Normal(0,1)
      Stratrec[3,,1] <- rexp(N,1) # Initiate noise genes ~Exp(1)
      Stratrec[4,,1] <- 0 #runif(N,0,1) # Initiate memory genes ~Uniform(0,1)
      Stratrec[5,,1] <- 0 #Dhanya: prob = relative probability? Should be okay as long as they are equal?
      
      # Reset storage vectors
      vs <- rep(4,N)#Dhanya: Scrounger payoff
      vp <- rep(4,N)#Dhanya: Producer payoff
      prs <- rep(0.5,N)#Dhanya: Proportion of scroungers in the population
      Info <- rbind(prs,vp,vs) #rowbind vector containing the two types of information (Proportion scrounger, Payoff produce, Payoff scrounge)
      Intake <- rep(0,N)
      
      ### START GENERATION ###
      gen_autocorr <- autocorr(tac,G)
      
      for(i in 1:G){
        
        #if(i > (G/2)){
          #mrate[2] <- 0.1
          #mrate[4] <- 0.1
          #mrate[5] <- 0.01
          #if(i == (G/2)){
            #Initialize plasticity and memory
            #Stratrec[2,,i] <- rnorm(N,0,1) # Initiate slope genes ~Normal(0,1)
            #Stratrec[4,,i] <- runif(N,0,1) # Initiate memory genes ~Uniform(0,1)
            #Stratrec[5,,i] <- sample(c(0,1),size=N,replace=T,prob=c(0.5,0.5)) #Initiate info genes
          #}
        #}else{
          #mrate[2] <- 0
          #mrate[4] <- 0
          #mrate[5] <- 0
        #}
        
        #For each generation:
        Intake[] <- 0 #Reset the food intake counter and patch counter
        Stratrec[7,,i] <- 0
        
        #Calculate food abundances this generation (if env_fluct!=NULL)
        if(any(env_fluct_scale=="g")){
          
          if(any(env_fluct=="f")){
            f_G <- gen_autocorr[i] + f
          }else{f_G <- f}
          
          if(any(env_fluct=="a")){
            a_G <- gen_autocorr[i] + a
            a_G <- min(a_G,f_G) # Finder's share can't be larger than patch richness.
          }else{a_G <- a}
          
          if(any(env_fluct=="p")){
            if(env_fluct_type=="random"){
            p_G <- rbeta(1,p+1,2-p) # Probability of finding patches this generation - use beta distribution so stays positive
            }
          }else{p_G <- p}
        }else{ # If environmental parameters don't vary across generations, let each generation have the same environment.
          f_G <- f
          a_G <- a
          p_G <- p
          
        }
        
        ### START TIME STEP ###
        step_autocorr <- autocorr(tac,S)
        
        for(j in 1:S){ 
          #For each time step:
          stepfood <- rep(0,N) #Reset food intake vector used for information
          
          if(any(env_fluct_scale=="s")){
            
            if(any(env_fluct=="f")){
              f_S <- step_autocorr[j] + f_G
            }else{f_S <- f_G}
            
            if(any(env_fluct=="a")){
              a_S <- step_autocorr[j] + a_G
              a_S <- min(a_S,f_S) # Finder's share can't be larger than patch richness.
            }else{a_S <- a_G}
            
            if(any(env_fluct=="p")){
              p_S <- 0.1*sin(0.1*j)+p_G
            }else{p_S <- p_G}
            
          }else{ # If environments don't fluctuate across time steps, give each time step the environment of that generation.
            f_S <- f_G
            a_S <- a_G
            p_S <- p_G
          }
          
          #Calculate probability of playing each tactic based on individual personality
          
          #If using proportion of scroungers as information:
          IT1 <- which(Stratrec[5,,i]==0) # Indices of individuals using Info Type 1
          liabp <- Stratrec[1,IT1,i]+Stratrec[2,IT1,i]*Info[1,IT1]+rnorm(length(IT1),0,Stratrec[3,IT1,i]) #Calculate liability of producing
          Stratrec[6,IT1,i] <- ifelse(liabp>0,1,0)
          
          #If using payoffs as information
          IT2 <- which(Stratrec[5,,i]==1) # Indices of individuals using info type 2
          payoff <- Info[2,IT2]/(Info[3,IT2]+Info[2,IT2]) # Add a nonlinear component here as well? (rwpar)
          liabp <-  Stratrec[1,IT2,i] + Stratrec[2,IT2,i]*payoff+rnorm(length(IT2),0,Stratrec[3,IT2,i]) # Calculate liability of producing
          Stratrec[6,IT2,i] <- ifelse(liabp>0,1,0)
          
          Prod <- which(Stratrec[6,,i] == 1) #Figuring out indices of all individuals playing the producer tactic this time step
          Scro <- which(Stratrec[6,,i] == 0) #Figuring out indices of all individuals playing the scrounger tactic this time step
          Stratrec[7,,i] <- Stratrec[7,,i]+Stratrec[6,,i] # Keep track of how many times each individual has produced.
          
          found <- Prod[which(runif(length(Prod))<p_S)] # Indices of all producers that found food
          
          stepfood[found] <- a_S+(f_S-a_S)/(1+length(Scro)) # These producers get food this time step.
          
          stepfood[Scro] <- length(found)*(f_S-a_S)/(1+length(Scro)) # Scroungers share all produced food among themselves
          
          Intake <- Intake+stepfood
          
          #Update information about %scrounger/payoffs to be used as input to strategise next time step
          # Those who use %scrounger:
          Info[1,IT1] <- Info[1,IT1]*Stratrec[4,IT1,i] + (1-Stratrec[4,IT1,i])*length(Scro)/N #Record %scrounger
          # Those who use payoffs
          updateP <- IT2[IT2 %in% Prod] #first those who produced this time
          Info[2,updateP] <- Info[2,updateP]*Stratrec[4,updateP,i] + (1-Stratrec[4,updateP,i])*stepfood[updateP]
          updateS <- IT2[IT2 %in% Scro] #then those who scrounged this time
          Info[3,updateS] <- Info[3,updateS]*Stratrec[4,updateS,i] + (1-Stratrec[4,updateS,i])*stepfood[updateS]
          
          Info[2,][which(Info[2,]==0 & Info[3,]==0)] <- Info[3,][which(Info[2,]==0 & Info[3,]==0)] <- 0.001
          
        }
        
        ### END TIME STEP ###
        
        Stratrec[7, ,i] <- Stratrec[7,,i]/S ### Calculate fraction of producing events out of all time steps
        
        ### SELECTION & REPRODUCTION ###
        
        ## Subtracting costs of plasticity/memory
        Intake[which(Stratrec[7,,i]!=1 & Stratrec[7,,i]!=0)] <- Intake[which(Stratrec[7,,i]!=1 & Stratrec[7,,i]!=0)] * (1-scost) # Cost of tactic switching
        #Intake <- ifelse(Intake<Intake*(1-pcost*Stratrec[2,,i]),Intake,Intake*(1-pcost*Stratrec[2,,i])) # Plasticity cost
        #Intake <- apply(cbind(Intake,Intake*(1-pcost*Stratrec[2,,i])),MARGIN=1,FUN=min,na.rm=T) # Plasticity cost
        Intake <- Intake*(1-pcost*abs(Stratrec[2,,i]))
        Intake <- Intake*(1-mcost*abs(Stratrec[4,,i]))
        Intake[Intake < 0] <- 0 #Dhanya: To remove negative values of intake
        ## Selection and (asexual) reproduction mechanism.
        if(sum(Intake)>0) {
          Intake[sample(1:N,N*dem_stoch)] <- 0 #Randomly kill a proportion of the individuals.
          Intake[which(Intake < quantile(Intake,sel_trunc,na.rm = T))] <- 0 #Truncated selection, choosing only the 1-quantile best.
          sel <- Intake/sum(Intake) # Determine relative fitness based on food gain
          ix <- sample(1:N,N,replace=TRUE,prob=sel) # Actual reproduction. New individuals are sampled from old population, where individuals with high fitness have higher prob. of being picked/sampled.
          Offspring[1:7,] <- Stratrec[1:7,ix,i] # Give new individuals the strategy of their parents
        } else {
          print(paste("No survivors in generation",i))
          break
        }
        
        ### END SELECTION AND REPRODUCTION ###
        
        ### MUTATIONS ###
        if(i<G){
          ## If changing elev 
          mut1 <- sample(N,size=round(N*mrate[1]),replace=F)
          #Offspring[1,mut1] <- rnorm(length(mut1),Offspring[1,mut1],rexp(length(mut1),1/msize)) # Or perhaps mutations should look like exp(log(old gene) + mutation~N(0,1))
          Offspring[1,mut1] <- rnorm(length(mut1),Offspring[1,mut1],msize)
          
          ## If changing slope
          mut2 <- sample(N,size=round(N*mrate[2]),replace=F)
          #Offspring[2,mut2] <- rnorm(length(mut2),Offspring[2,mut2],rexp(length(mut2),1/msize))
          Offspring[2,mut2] <- rnorm(length(mut2),Offspring[2,mut2],msize)
          
          ## If changing noise
          mut3 <- sample(N,size=round(N*mrate[3]),replace=F)
          #Offspring[3,mut3] <- rnorm(length(mut3),Offspring[3,mut3],rexp(length(mut3),1/msize))
          Offspring[3,mut3] <- rnorm(length(mut3),Offspring[3,mut3],msize)
          Offspring[3,which(Offspring[3,]<0)] <- 0
          
          #If changing memory
          mut4 <- sample(N,size=round(N*mrate[4]),replace=F)
          #Offspring[4,mut4] <- rnorm(length(mut4),Offspring[4,mut4],rexp(length(mut4),1/msize))
          Offspring[4,mut4] <- rnorm(length(mut4),Offspring[4,mut4],msize)
          Offspring[4,which(Offspring[4,]<0)] <- 0
          Offspring[4,which(Offspring[4,]>1)] <- 1
          
          ##If changing infotype
          mut5 <- sample(N,size=round(N*mrate[5]),replace=F)
          #Offspring[5,mut5] <- as.numeric(!(Offspring[5,mut5]))
          Offspring[5,mut5] <- ifelse(Offspring[5,mut5]==0,1,0)
          
          ### END MUTATIONS ###
          Stratrec[1:7,,i+1] <- Offspring
        }
      }
      
      ### END GENERATION ###
      #Save output
      Totrec[, , ,z] <- Stratrec
      rownames(Totrec) <- rownames(Stratrec)
      print(z)
    }
    
    ### END REPLICATE ###
    params <- list(Rounds,N,S,G,mrate,msize,pcost,mcost,scost,tac,dem_stoch,sel_trunc,f,a,p,env_fluct)
    names(params) <- c("Rounds","N","S","G","mrate","msize","pcost","mcost","scost","tac","dem_stoch","sel_trunc","f","a","p","env_fluct")
    return(list(Totrec,params))
    
    
  }

##### Running and saving simulations #####
###Case1a###
start.time <- Sys.time()

data1a <- sim1.1(Rounds=10,N=1000,S=100,G=3000,mrate=c(0.1,0,0.1,0,0),
                     msize=0.1,pcost=0.01,mcost=0.01,scost=0.01,tac=0,
                     dem_stoch=0.1,sel_trunc=0.2,f=20,a=6,
                     p=0.5,env_fluct=NULL,env_fluct_type=NULL,env_fluct_scale=NULL,rw_par=NULL)

Sys.time()-start.time
saveRDS(data1a,file=paste0("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/1a/data1a.Rdata"))
###

###Case1b###
start.time <- Sys.time()

data1b <- sim1.1(Rounds=10,N=1000,S=100,G=3000,mrate=c(0.1,0,0.1,0,0),
                msize=0.1,pcost=0.01,mcost=0.01,scost=0.01,tac=0.5,
                dem_stoch=0.1,sel_trunc=0.2,f=20,a=6,
                p=0.5,env_fluct=NULL,env_fluct_type=NULL,env_fluct_scale=NULL,rw_par=NULL)

Sys.time()-start.time
saveRDS(data1b,file=paste0("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/1b/data1b.Rdata"))
###

###Case1c###
start.time <- Sys.time()

data1c <- sim1.1(Rounds=10,N=1000,S=100,G=3000,mrate=c(0.1,0,0.1,0,0),
                 msize=0.1,pcost=0.01,mcost=0.01,scost=0.01,tac=0.99,
                 dem_stoch=0.1,sel_trunc=0.2,f=20,a=6,
                 p=0.5,env_fluct=NULL,env_fluct_type=NULL,env_fluct_scale=NULL,rw_par=NULL)

Sys.time()-start.time
saveRDS(data1c,file=paste0("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/1c/data1c.Rdata"))
###

###Case2a###
start.time <- Sys.time()

data2a <- sim1.1(Rounds=10,N=1000,S=100,G=3000,mrate=c(0.1,0,0.1,0,0),
                     msize=0.1,pcost=0.01,mcost=0.01,scost=0.01,tac=0,
                     dem_stoch=0.1,sel_trunc=0.2,f=20,a=6,
                     p=0.5,env_fluct=c("a", "f"),env_fluct_type="autocorr",env_fluct_scale="s",rw_par=NULL)

Sys.time()-start.time
saveRDS(data2a,file=paste0("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/2a/data2a.Rdata"))
###

###Case2b###
start.time <- Sys.time()

data2b <- sim1.1(Rounds=10,N=1000,S=100,G=3000,mrate=c(0.1,0,0.1,0,0),
                 msize=0.1,pcost=0.01,mcost=0.01,scost=0.01,tac=0.5,
                 dem_stoch=0.1,sel_trunc=0.2,f=20,a=6,
                 p=0.5,env_fluct=c("a", "f"),env_fluct_type="autocorr",env_fluct_scale="s",rw_par=NULL)

Sys.time()-start.time
saveRDS(data2b,file=paste0("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/2b/data2b.Rdata"))
###

###Case2c###
start.time <- Sys.time()

data2c <- sim1.1(Rounds=10,N=1000,S=100,G=3000,mrate=c(0.1,0,0.1,0,0),
                 msize=0.1,pcost=0.01,mcost=0.01,scost=0.01,tac=0.99,
                 dem_stoch=0.1,sel_trunc=0.2,f=20,a=6,
                 p=0.5,env_fluct=c("a", "f"),env_fluct_type="autocorr",env_fluct_scale="s",rw_par=NULL)

Sys.time()-start.time
saveRDS(data2c,file=paste0("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/2c/data2c.Rdata"))
###

###Case3a###
start.time <- Sys.time()

data3a <- sim1.1(Rounds=10,N=1000,S=100,G=3000,mrate=c(0.1,0,0.1,0,0),
                     msize=0.1,pcost=0.01,mcost=0.01,scost=0.01,tac=0,
                     dem_stoch=0.1,sel_trunc=0.2,f=20,a=6,
                     p=0.5,env_fluct=c("a", "f"),env_fluct_type="autocorr",env_fluct_scale="g",rw_par=NULL)

Sys.time()-start.time
saveRDS(data3a,file=paste0("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/3a/data3a.Rdata"))
###

###Case3b###
start.time <- Sys.time()

data3b <- sim1.1(Rounds=10,N=1000,S=100,G=3000,mrate=c(0.1,0,0.1,0,0),
                 msize=0.1,pcost=0.01,mcost=0.01,scost=0.01,tac=0.5,
                 dem_stoch=0.1,sel_trunc=0.2,f=20,a=6,
                 p=0.5,env_fluct=c("a", "f"),env_fluct_type="autocorr",env_fluct_scale="g",rw_par=NULL)

Sys.time()-start.time
saveRDS(data3b,file=paste0("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/3b/data3b.Rdata"))
###

###Case3c###
start.time <- Sys.time()

data3c <- sim1.1(Rounds=10,N=1000,S=100,G=3000,mrate=c(0.1,0,0.1,0,0),
                 msize=0.1,pcost=0.01,mcost=0.01,scost=0.01,tac=0.99,
                 dem_stoch=0.1,sel_trunc=0.2,f=20,a=6,
                 p=0.5,env_fluct=c("a", "f"),env_fluct_type="autocorr",env_fluct_scale="g",rw_par=NULL)

Sys.time()-start.time
saveRDS(data3c,file=paste0("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/3c/data3c.Rdata"))
###

###Case4a###
start.time <- Sys.time()

data4a <- sim1.1(Rounds=10,N=1000,S=100,G=3000,mrate=c(0.1,0,0.1,0,0),
                     msize=0.1,pcost=0.01,mcost=0.01,scost=0.01,tac=0,
                     dem_stoch=0.1,sel_trunc=0.2,f=20,a=6,
                     p=0.5,env_fluct=c("a", "f"),env_fluct_type="autocorr",env_fluct_scale=c("s","g"),rw_par=NULL)

Sys.time()-start.time
saveRDS(data4a,file=paste0("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/4a/data4a.Rdata"))
###

###Case4b###
start.time <- Sys.time()

data4b <- sim1.1(Rounds=10,N=1000,S=100,G=3000,mrate=c(0.1,0,0.1,0,0),
                 msize=0.1,pcost=0.01,mcost=0.01,scost=0.01,tac=0.5,
                 dem_stoch=0.1,sel_trunc=0.2,f=20,a=6,
                 p=0.5,env_fluct=c("a", "f"),env_fluct_type="autocorr",env_fluct_scale=c("s","g"),rw_par=NULL)

Sys.time()-start.time
saveRDS(data4b,file=paste0("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/4b/data4b.Rdata"))
###

###Case4c###
start.time <- Sys.time()

data4c <- sim1.1(Rounds=10,N=1000,S=100,G=3000,mrate=c(0.1,0,0.1,0,0),
                 msize=0.1,pcost=0.01,mcost=0.01,scost=0.01,tac=0.99,
                 dem_stoch=0.1,sel_trunc=0.2,f=20,a=6,
                 p=0.5,env_fluct=c("a", "f"),env_fluct_type="autocorr",env_fluct_scale=c("s","g"),rw_par=NULL)

Sys.time()-start.time
saveRDS(data4c,file=paste0("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/4c/data4c.Rdata"))
###
#####

#####################################
## Generating data from saved file ##
#####################################
data <- readRDS("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/4c/data4c.Rdata")
Totrec <- data[[1]]
#Totrec format:
#First dimension: Gene values (1: Elev, 2: Slope, 3: Noise, 4: Memory, 5: Info Type, 6: CurrStrat, 7: Fraction Prod played)
#Second dimension: Individuals
#Third dimension: Generations
#Fourth dimension: Rounds
Rounds <- length(Totrec[1,1,1,]) 
G <- length(Totrec[1,1,,1]) #No of Generations in simulation
N <- length(Totrec[1,,1,1]) #Population size in simulation
loci <- length(Totrec[,1,1,1])-2
f <- data[[2]]$"f"
a <- data[[2]]$"a"
p <- data[[2]]$"p"
Elevation<-Slope<-Noise<-Inflection<-Memory<-Tactic<-rbind(rep(0,G),rep(0,G))
Information<-rep(0,G)

propstrats<-matrix(data=0,ncol=9,nrow=G) 
colnames(propstrats) <- c("Non-plast", "Step %p", "Step Vp/Vs", "Logistic %p", "Logistic Vp/Vs", "Learning %p", "Learning Vp/Vs", "Fixed","%Prod in pop")

Proprec<-array(0,dim=c(G,9,Rounds))
Meanrec<-array(0,dim=c(11,G,Rounds))

for (z in 1:Rounds) {
  for (i in 1:G) {
    Elevation[1,i] <- mean(Totrec[1, ,i,z])
    Elevation[2,i] <- sd(Totrec[1, ,i,z])
    
    Slope[1,i] <- mean(Totrec[2, ,i,z])
    Slope[2,i] <- sd(Totrec[2, ,i,z])
    
    Noise[1,i] <- mean(Totrec[3, ,i,z]) #Only calculate for the ones that are plastic?
    Noise[2,i] <- sd(Totrec[3, ,i,z])
    
    Memory[1,i] <- mean(Totrec[4, ,i,z])
    Memory[2,i] <- sd(Totrec[4, ,i,z])
    
    Tactic[1,i] <- mean(Totrec[7, ,i,z])
    Tactic[2,i] <- sd(Totrec[7, ,i,z])
    
    Information[i] <- mean(Totrec[5, ,i,z])
    
    #Could still get something like this by the following:
    propstrats[i,] <- c(
      length(which(abs(Totrec[2,,i,z]) < 0.01 & abs(Totrec[1,,i,z])<5*Totrec[3,,i,z]))/N, # Non-plast, mixed
      length(which(Totrec[3,,i,z] < 0.01 & abs(Totrec[2,,i,z])>0.01 & Totrec[4,,i,z] < 0.01 & Totrec[5,,i,z]==0))/N, #Step %scr
      length(which(Totrec[3,,i,z] < 0.01 & abs(Totrec[2,,i,z])>0.01 & Totrec[4,,i,z] < 0.01 & Totrec[5,,i,z]==1))/N, #Step payoff
      length(which(Totrec[3,,i,z] > 0.01 & abs(Totrec[2,,i,z])>0.01 & Totrec[4,,i,z] < 0.01 & Totrec[5,,i,z]==0))/N, #Log %scr
      length(which(Totrec[3,,i,z] > 0.01 & abs(Totrec[2,,i,z])>0.01 & Totrec[4,,i,z] < 0.01 & Totrec[5,,i,z]==1))/N, #Log payoff)
      length(which(Totrec[4,,i,z] > 0.01 & abs(Totrec[2,,i,z])>0.01 & Totrec[5,,i,z]==0))/N, #Learn %scr)
      length(which(Totrec[4,,i,z] > 0.01 & abs(Totrec[2,,i,z])>0.01 & Totrec[5,,i,z]==1))/N, #Learn payoff)
      length(which(abs(Totrec[2,,i,z]) < 0.01 & abs(Totrec[1,,i,z])>5*Totrec[3,,i,z]))/N, # Fixed
      mean(Totrec[7,,i,z]) # Mean proportion producer (across all time steps)
    )
  }
  
  Final <-
    rbind(Elevation, Slope, Noise, Memory, Tactic, Information)
  
  
  Proprec[, , z] <- propstrats
  Meanrec[, , z] <- Final
  ## Meanrec format:
  # Dim 1: 
  # 1) Mean elev
  # 2) Sd elev
  # 3) Mean slope
  # 4) Sd slope
  # 5) Mean noise
  # 6) Sd noise
  # 7) Mean memory
  # 8) Sd memory
  # 9) Mean tactic use (% producing)
  # 10) Sd tactic use
  # 11) Mean info type
  
}

#Plot of strategies from each round
#par(mfrow=c(2,5),mar=c(4,4,1,1))
#for (z in 1:Rounds) {
#  matplot(propstrats, type="l",col=c(1:8,"Gold"),lwd=2,ylab="Proportion of each strategy in population",xlab="Generations",ylim=c(0,1),xlim=c(0,G),main=z)
#}
#legend("topleft",legend=c("Non-plast mix","Step %scr","Step payoff","Log %scr","Log payoff","Learn %scr","Learn payoff","Fixed","Fraction Prod"),col=c(1:8,"Gold"),pch=1,cex=0.9)
#lines(1:G,rep(a/f+(1/N),G),lty=2) # If you want to draw the line for the ESS you need to define the parameters a, f and N.

par(mfrow=c(1,1),mar=c(4,4,2,1))
#Summary plot of all gene values across all reps
plot(1:G,apply(Meanrec[1,,],1,mean),type="l",ylim=c(-3,3),ylab="Mean gene values",xlab="Generation",cex.lab=1) #Elevation
for(r in 1:10){
  lines(1:G,Meanrec[1,,r],col=alpha(col=1,0.1))#Elevation
  lines(1:G,Meanrec[3,,r],col=alpha(col=2,0.1))#Slope
  lines(1:G,Meanrec[5,,r],col=alpha(col=3,0.1))#Noise
  lines(1:G,Meanrec[7,,r],col=alpha(col=4,0.1))#Memory
  lines(1:G,Meanrec[11,,r],col=alpha(col=7,0.1))#Info type
}
lines(1:G,apply(Meanrec[3,,],1,mean,na.rm=T),col=2) #Slope
lines(1:G,apply(Meanrec[5,,],1,mean,na.rm=T),col=3) #Noise
lines(1:G,apply(Meanrec[7,,],1,mean),col=4) #Memory
lines(1:G,apply(Meanrec[11,,],1,mean,na.rm=T),col=7) #Info type
#lines(1:G,apply(Meanrec[9,,],1,mean),col=6) #Proportion producing
legend("topleft",legend=c("Elevation","Noise"),col=c(1,3),pch=1,cex=1)
#lines(1:G,rep(a/f+(1/N),G),lty=2) # If you want to draw the line for the ESS you need to define the parameters a, f and N.
abline(h=0,col="gray",lty = "solid",lwd = 2)
#abline(v=1500,col="gray",lty = "dashed",lwd = 2)
#lines(1:G,apply(Meanrec[12,,],1,mean),col=7) #Proportion using %scrounger as information

#Confidence intervals example - can add for the others too but it might get messy...
#lowerbounds <- apply(Meanrec[1,,],1,mean)-apply(Meanrec[2,,],1,mean,na.rm=T)
#upperbounds <- apply(Meanrec[1,,],1,mean)+apply(Meanrec[2,,],1,mean,na.rm=T)
#polygon(c(1:G,rev(1:G)),c(lowerbounds,rev(upperbounds)),col=alpha(1,alpha=0.1),border=NA)

#lowerbounds <- apply(Meanrec[3,,],1,mean)-apply(Meanrec[4,,],1,mean,na.rm=T)
#upperbounds <- apply(Meanrec[3,,],1,mean)+apply(Meanrec[4,,],1,mean,na.rm=T)
#polygon(c(1:G,rev(1:G)),c(lowerbounds,rev(upperbounds)),col=alpha(2,alpha=0.1),border=NA)

#lowerbounds <- apply(Meanrec[5,,],1,mean)-apply(Meanrec[6,,],1,mean,na.rm=T)
#upperbounds <- apply(Meanrec[5,,],1,mean)+apply(Meanrec[6,,],1,mean,na.rm=T)
#polygon(c(1:G,rev(1:G)),c(lowerbounds,rev(upperbounds)),col=alpha(3,alpha=0.1),border=NA)

#lowerbounds <- apply(Meanrec[7,,],1,mean)-apply(Meanrec[8,,],1,mean,na.rm=T)
#upperbounds <- apply(Meanrec[7,,],1,mean)+apply(Meanrec[8,,],1,mean,na.rm=T)
#polygon(c(1:G,rev(1:G)),c(lowerbounds,rev(upperbounds)),col=alpha(4,alpha=0.1),border=NA)

#lowerbounds <- apply(Meanrec[11,,],1,mean)-apply(Meanrec[6,,],1,mean,na.rm=T)
#upperbounds <- apply(Meanrec[11,,],1,mean)+apply(Meanrec[6,,],1,mean,na.rm=T)
#polygon(c(1:G,rev(1:G)),c(lowerbounds,rev(upperbounds)),col=alpha(7,alpha=0.1),border=NA)

#Quantiles?
#matplot(t(apply(Meanrec[1,,],1,quantile,na.rm=T)),pch=1,type="l",xlab="Generation",ylab="Elevation (population means)",
#col=c("Light grey","Dark grey","Black","Dark grey","Light grey"),lty=1,ylim=c(-1,1)) #0, 25, 50, 75, 100 percentiles

#Producing vs scrounging
producers <- apply(Proprec[,9,],1,mean)
scroungers <- 1-producers
prod_sd <- apply(Proprec[,9,],1,sd)
scro_sd <- apply(Proprec[,9,],1,sd)
plot(1:G,producers,type="l",ylim=c(0,1),col="Blue",lwd=2,xlab="Generation",ylab="Percent tactic use",cex.lab=1)
lines(1:G,scroungers,col="Red",lwd=2)
polygon(c(1:G,rev(1:G)),c(producers-prod_sd,producers+prod_sd),col=alpha("Blue",alpha=0.1),border=NA)
polygon(c(1:G,rev(1:G)),c(scroungers-scro_sd,scroungers+scro_sd),col=alpha("Red",alpha=0.1),border=NA)
legend("topleft", legend=c("Producers","Scroungers"),col=c("Blue","Red"),pch=1,cex=1)
lines(1:G,rep(a/f+(1/N),G),lty=2) #ESS proportion producer

#Summary plot of all strategies
Meanprops <- matrix(0,nrow=G,ncol=9)
Sdprops <- matrix(0,nrow=G,ncol=9) #Not sure if standard deviations here are useful...

for (i in 1:9){
  for (j in 1:G){
    Meanprops[j,i] <- mean(Proprec[j,i,])
    Sdprops[j,i] <- sd(Proprec[j,i,]) 
  }
}

#Average strategy use per generation across rounds
matplot(Meanprops[,1:8], type="l",pch=1,col=c(1:8),lwd=2,ylab="Proportion of each strategy in population",
        xlab="Generations",cex.lab=1,ylim=c(0,1),xlim=c(0,G))
for(i in 1:8){
  polygon(c(1:G,rev(1:G)),c(Meanprops[,i]-Sdprops[,i],rev(Meanprops[,i]+Sdprops[,i])),col=alpha(i,alpha=0.1),border=NA)
}
#polygon(c(1:G,rev(1:G)),c(Meanprops[,7]-Sdprops[,7],Meanprops[,7]+Sdprops[,7]),col=alpha("Orange",alpha=0.1),border=NA)
#polygon(c(1:G,rev(1:G)),c(Meanprops[,8]-Sdprops[,8],Meanprops[,8]+Sdprops[,8]),col=alpha("Dark grey",alpha=0.1),border=NA)
legend("topleft", legend = c("Mixed","Pure"), col=c(1,8),pch=1,cex=1)
abline(h=0,col="gray",lty = "solid",lwd = 2)

## 20 randomly sampled reaction norms from the final time steps of 10 replicate populations
par(mfrow=c(2,5),mar=c(1,1,1,1),oma=c(2,2,0,0))
etas <- rnorm(10000)
infos <- seq(0,1,by=0.01)
for(z in sample(1:Rounds,min(10,Rounds))){
  plot(infos,infos,type="n",xlab="",ylab="")# xaxt="n",yaxt="n"
  #axis(side=1,)
  for (i in sample(1:N,50,replace=F)){
    isfixed <- logical(1)
    ys <- numeric(101)
    for(j in 1:length(infos)){
      ys[j] <- mean(Totrec[1,i,G,z]+infos[j]*Totrec[2,i,G,z]+Totrec[3,i,G,z]*etas>0)
    }
    isfixed <- ifelse(abs(Totrec[1,i,G,z])<5*Totrec[3,i,G,z],T,F)
    if(as.numeric(isfixed) == 0){
      lines(seq(0,1,0.01),ys,lty=1,col=2)#pure strategies = red
    }
    else{
      lines(seq(0,1,0.01),ys,lty=1,col=12)#mixed strategies = blue
    }
  }
}
#mtext("Probability of producing",side=2,line=0.2,outer=T)
#mtext("Information",side=1,line=0.2,outer=T)

###Summary plots for environmental scales###
pal <- brewer.pal(10,"Paired")

##Noise##
noise_scale <- array(0, dim=c(10,4))
colnames(noise_scale) <- c("Case 1","Case 2","Case 3","Case 4")
noise_scale[,1] <- apply(Meanrec[5,(G-100):G,], FUN = mean, MARGIN = 2)
noise_scale[,2] <- apply(Meanrec[5,(G-100):G,], FUN = mean, MARGIN = 2)
noise_scale[,3] <- apply(Meanrec[5,(G-100):G,], FUN = mean, MARGIN = 2)
noise_scale[,4] <- apply(Meanrec[5,(G-100):G,], FUN = mean, MARGIN = 2)
saveRDS(noise_scale,file=paste0("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/Scale/noise_scale.Rdata"))
noise_scale <- readRDS("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/Scale/noise_scale.Rdata")
vioplot(noise_scale, col = alpha(2:5,alpha=0.2),ylim=c(0,1.25),xlab="Scale of environmental variation",ylab = "Noise",colMed = "black",rectCol = "grey",cex.lab=2,border=NA)
points(jitter(c(rep(1,10),rep(2,10),rep(3,10),rep(4,10))),noise_scale,add=T, col = pal[1:10], pch=20, cex=1.5)
abline(h=0,col="gray",lty = "solid",lwd = 2)
##

##Elevation##
elev_scale <- array(0, dim=c(10,4))
colnames(elev_scale) <- c("Case 1","Case 2","Case 3","Case 4")
elev_scale[,1] <- apply(Meanrec[1,(G-100):G,], FUN = mean, MARGIN = 2)
elev_scale[,2] <- apply(Meanrec[1,(G-100):G,], FUN = mean, MARGIN = 2)
elev_scale[,3] <- apply(Meanrec[1,(G-100):G,], FUN = mean, MARGIN = 2)
elev_scale[,4] <- apply(Meanrec[1,(G-100):G,], FUN = mean, MARGIN = 2)
saveRDS(elev_scale,file=paste0("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/Scale/elev_scale.Rdata"))
elev_scale <- readRDS("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/Scale/elev_scale.Rdata")
vioplot(elev_scale, col = alpha(2:5,alpha=0.2),ylim=c(-4,0),xlab="Scale of environmental variation",ylab = "Elevation",colMed = "black",rectCol = "grey",cex.lab=2,border=NA)
points(jitter(c(rep(1,10),rep(2,10),rep(3,10),rep(4,10))),elev_scale,add=T, col = pal[1:10], pch=20, cex=1.5)
abline(h=0,col="gray",lty = "solid",lwd = 2)
##

##Pure strategies##
pure_scale <- array(0, dim=c(10,4))
colnames(pure_scale) <- c("None","Within","Between","Both")
pure_scale[,1] <- apply(Proprec[(G-100):G,8,], FUN = mean, MARGIN = 2)
pure_scale[,2] <- apply(Proprec[(G-100):G,8,], FUN = mean, MARGIN = 2)
pure_scale[,3] <- apply(Proprec[(G-100):G,8,], FUN = mean, MARGIN = 2)
pure_scale[,4] <- apply(Proprec[(G-100):G,8,], FUN = mean, MARGIN = 2)
saveRDS(pure_scale,file=paste0("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/Scale/pure_scale.Rdata"))
pure_scale <- readRDS("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/Scale/pure_scale.Rdata")
vioplot(pure_scale, col = alpha(2:5,alpha=0.2),xlab="Scale of environmental variation",ylab = "Proportion of pure strategies",colMed = "black",rectCol = "grey",cex.lab=2,border=NA)
points(jitter(c(rep(1,10),rep(2,10),rep(3,10),rep(4,10))),pure_scale,add=T, col = pal[1:10], pch=20, cex=1.5)
abline(h=0,col="gray",lty = "solid",lwd = 2)
##

##Mixed strategies##
mixed_scale <- array(0, dim=c(10,4))
colnames(mixed_scale) <- c("None","Within","Between","Both")
mixed_scale[,1] <- apply(Proprec[(G-100):G,1,], FUN = mean, MARGIN = 2)
mixed_scale[,2] <- apply(Proprec[(G-100):G,1,], FUN = mean, MARGIN = 2)
mixed_scale[,3] <- apply(Proprec[(G-100):G,1,], FUN = mean, MARGIN = 2)
mixed_scale[,4] <- apply(Proprec[(G-100):G,1,], FUN = mean, MARGIN = 2)
saveRDS(mixed_scale,file=paste0("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/Scale/mixed_scale.Rdata"))
mixed_scale <- readRDS("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/Scale/mixed_scale.Rdata")
vioplot(mixed_scale, col = alpha(2:5,alpha=0.2),ylim=c(0,0.5),xlab="Scale of environmental variation",ylab = "Proportion of mixed strategies",colMed = "black",rectCol = "grey",cex.lab=2,border=NA)
points(jitter(c(rep(1,10),rep(2,10),rep(3,10),rep(4,10))),mixed_scale,add=T, col = pal[1:10], pch=20, cex=1.5)
abline(h=0,col="gray",lty = "solid",lwd = 2)
##

###Summary plots for environmental autocorrelation###
pal <- brewer.pal(10,"Paired")

##Noise##
noise_ac <- array(0, dim=c(10,3))
colnames(noise_ac) <- c("0","0.5","0.99")
noise_ac[,1] <- apply(Meanrec[5,(G-100):G,], FUN = mean, MARGIN = 2)
noise_ac[,2] <- apply(Meanrec[5,(G-100):G,], FUN = mean, MARGIN = 2)
noise_ac[,3] <- apply(Meanrec[5,(G-100):G,], FUN = mean, MARGIN = 2)
saveRDS(noise_ac,file=paste0("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/Autocorrelation/noise_ac.Rdata"))
noise_ac <- readRDS("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/Autocorrelation/noise_ac.Rdata")
vioplot(noise_ac, col = alpha(2:5,alpha=0.2),ylim=c(0,1.25),xlab="Environmental autocorrelation",ylab = "Noise",colMed = "black",rectCol = "grey",cex.lab=2,border=NA)
points(jitter(c(rep(1,10),rep(2,10),rep(3,10))),noise_ac,add=T, col = pal[1:10], pch=20, cex=1.5)
abline(h=0,col="gray",lty = "solid",lwd = 2)
##

##Elevation##
elev_ac <- array(0, dim=c(10,3))
colnames(elev_ac) <- c("0","0.5","0.99")
elev_ac[,1] <- apply(Meanrec[1,(G-100):G,], FUN = mean, MARGIN = 2)
elev_ac[,2] <- apply(Meanrec[1,(G-100):G,], FUN = mean, MARGIN = 2)
elev_ac[,3] <- apply(Meanrec[1,(G-100):G,], FUN = mean, MARGIN = 2)
saveRDS(elev_ac,file=paste0("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/Autocorrelation/elev_ac.Rdata"))
elev_ac <- readRDS("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/Autocorrelation/elev_ac.Rdata")
vioplot(elev_ac, col = alpha(2:5,alpha=0.2),ylim=c(-4,0),xlab="Environmental autocorrelation",ylab = "Elevation",colMed = "black",rectCol = "grey",cex.lab=2,border=NA)
points(jitter(c(rep(1,10),rep(2,10),rep(3,10))),elev_ac,add=T, col = pal[1:10], pch=20, cex=1.5)
abline(h=0,col="gray",lty = "solid",lwd = 2)
##

##Pure strategies##
pure_ac <- array(0, dim=c(10,3))
colnames(pure_ac) <- c("0","0.5","0.99")
pure_ac[,1] <- apply(Proprec[(G-100):G,8,], FUN = mean, MARGIN = 2)
pure_ac[,2] <- apply(Proprec[(G-100):G,8,], FUN = mean, MARGIN = 2)
pure_ac[,3] <- apply(Proprec[(G-100):G,8,], FUN = mean, MARGIN = 2)
saveRDS(pure_ac,file=paste0("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/Autocorrelation/pure_ac.Rdata"))
pure_ac <- readRDS("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/Autocorrelation/pure_ac.Rdata")
vioplot(pure_ac, col = alpha(2:5,alpha=0.2),xlab="Environmental autocorrelation",ylab = "Proportion of pure strategies",colMed = "black",rectCol = "grey",cex.lab=2,border=NA)
points(jitter(c(rep(1,10),rep(2,10),rep(3,10))),pure_ac,add=T, col = pal[1:10], pch=20, cex=1.5)
abline(h=0,col="gray",lty = "solid",lwd = 2)
##

##Mixed strategies##
mixed_ac <- array(0, dim=c(10,3))
colnames(mixed_ac) <- c("0","0.5","0.99")
mixed_ac[,1] <- apply(Proprec[(G-100):G,1,], FUN = mean, MARGIN = 2)
mixed_ac[,2] <- apply(Proprec[(G-100):G,1,], FUN = mean, MARGIN = 2)
mixed_ac[,3] <- apply(Proprec[(G-100):G,1,], FUN = mean, MARGIN = 2)
saveRDS(mixed_ac,file=paste0("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/Autocorrelation/mixed_ac.Rdata"))
mixed_ac <- readRDS("/Users/MyPC/Desktop/Social sparrows/Code/NEW/Non-plastic/Autocorrelation/mixed_ac.Rdata")
vioplot(mixed_ac, col = alpha(2:5,alpha=0.2),ylim=c(0,0.5),xlab="Environmental autocorrelation",ylab = "Proportion of mixed strategies",colMed = "black",rectCol = "grey",cex.lab=2,border=NA)
points(jitter(c(rep(1,10),rep(2,10),rep(3,10))),mixed_ac,add=T, col = pal[1:10], pch=20, cex=1.5)
abline(h=0,col="gray",lty = "solid",lwd = 2)
##

cor.test(noise_scale[,1], elev_scale[,1], method=c("pearson"))
cor.test(noise_scale[,2], elev_scale[,2], method=c("pearson"))
cor.test(noise_scale[,3], elev_scale[,3], method=c("pearson"))
cor.test(noise_scale[,4], elev_scale[,4], method=c("pearson"))

cor.test(noise_ac[,1], elev_ac[,1], method=c("pearson"))
cor.test(noise_ac[,2], elev_ac[,2], method=c("pearson"))
cor.test(noise_ac[,3], elev_ac[,3], method=c("pearson"))
