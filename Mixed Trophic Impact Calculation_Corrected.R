#**************************************************************************
#**** File: Calculate Mixed Trophic Impact From Ecosystem Model Output ****
#****   Project: Ecosystem Modelling of the Oceanic Gulf of Mexico     ****
#****                    Developer: Matt Woodstock                     ****
#**************************************************************************

## Citation ##
#* Woodstock, M.S., T.T. Sutton, T. Frank, Y. Zhang. (2021). 
#*   An early warning sign: trophic structure changes in the 
#*   oceanic Gulf of Mexico from 2011-2018. Ecological Modelling.

## Contact Email ##
#* fishesofthedeep@gmail.com - Matt Woodstock


## Notes ##
#* Most of this has to do with taking Ecopath with Ecosim model output and calculating MTI (Part 1)
#* This methodology could be used for any ecosystem where you have both diet and predation mortality matrices (Part 2)

## Clear Working Directory
rm(list=ls())

## Load Packages ##
packages = c("MASS","xlsx")

package.check <- lapply(packages,FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

#Set working directory#
#* For EwE output, @indir is your mc_Ecosimscenario directory
indir<-"LOCATION OF YOUR MONTE CARLO RESULTS"
outdir <-"YOUR OUTPUT DIRECTORY" #*N files created = N iterations:: Best to make a separate directory

## Load Species names ##
## You need to create these files ##
#* I provide an example. Take note of the formatting issues with importing special characters in R
load(file="Species Information.RData")

#species<-as.matrix(read.csv("FILE NAME.csv",header=F)) #List of functional groups in the model
#species_adj<-as.matrix(read.csv("FILE NAME.csv",header=F)) #Functional groups spelled as R will input them (i.e., special characters become periods)

## Global Parameters ##
n_sim<-  #Number of iterations
n_spec<-length(species[,2]) #Number of functional groups
n_fisher<- #Number of fisheries in the model
year_start<- #Model Start Year
year_end<-  #Model End Year
n_year<- #Number of years of the simulation
n_year_concern<-  #Could be different/implemented if you include a spin up to the years of concern


#### Part 1:  Gather Ecopath with Ecosim Output and Calculate MTI For Each Year and Iteration####
#* If you already have a diet matrix, skip this
mti_array<-array(0,dim=c(n_spec+n_fisher,n_spec+n_fisher,n_year_concern,n_sim)) #* Placeholder for future calculation in loop
MTI<-array(0,dim=c(n_spec+n_fisher,n_spec+n_fisher,n_year_concern,n_sim)) #* Actual result that includes indirect impacts
n_count<-0 #Use a counter to control simulations

##Load Information##
for (file_count in 1:n_sim){
  if (file_count < 10){
    foldername<-paste0(indir,"/mc_output_trial000",file_count,sep="") #Working Directory where output is located
    if ((dir.exists(foldername))){
      setwd(foldername)
      n_count<-n_count+1
    }
  }
  if (file_count < 100){
    foldername<-paste0(indir,"/mc_output_trial00",file_count,sep="")
    if ((dir.exists(foldername))){
      setwd(foldername)
      n_count<-n_count+1
    }
  }
  if (file_count < 1000){
    foldername<-paste0(indir,"/mc_output_trial0",file_count,sep="")
    if ((dir.exists(foldername))){
      setwd(foldername)
      n_count<-n_count+1
    }
  }
  if (file_count < 10000){
    foldername<-paste0(indir,"/mc_output_trial",file_count,sep="")
    if ((dir.exists(foldername))){
      setwd(foldername)
      n_count<-n_count+1
    }
  }
  
  #Input diet composition#
  diet_array<-array(0,dim=c(n_spec+n_fisher,n_spec+n_fisher,n_year_concern))
  for (a in 1:length(species[,2])){
    filename<-paste("prey_",species[a,2],"_annual.csv",sep="")
    if (file.exists(filename)){
      diet<-read.csv(filename,skip=9,header=T) #Necessary data always starts on row 10
      prey_taxa<-colnames(diet) #Collected prey items of predator
      for (b in 1:length(diet[,1])){
        for (c in 1:length(species[,2])){
          for (d in 2:length(diet[1,])){
            if (species_adj[c,2]==prey_taxa[d]){
              diet_array[c,a,b]<-diet[b,d] #You are creating a diet matrix by skipping all non-preys
            }
          }
        }
      }
    }
  }
  
  if (n_fisher>0){
    filename<-"catch-fleet-group_annual.csv" #Fishery Catches
    fish_contrib<-array(0,dim=c(n_spec,n_fisher,n_year_concern))
    
    #** If you have alot of fisheries or a long model, creat a subset to reduce run time **
    if (file.exists(filename)){
      fishery<-read.csv(filename,skip=9,header=T)
      for (z in 1:length(fishery[,1])){
        for (y in year_start:year_end){ #Input years your model starts and ends
          for (x in 1:n_fisher){
            for (w in 1:n_spec){
              if (fishery[z,1]==y && fishery[z,2]==x && fishery[z,3]==w){
                fish_contrib[w,x,y-(year_start-1)]<-fishery[z,4]
              }
            }
          }
        }
      }
    }
    
    ## Calculate "Diet Composition" for fishery and append to diet_array
    for (a in ((n_spec+1):(n_spec+n_fisher))){
      for (y in year_start:year_end){
        for (x in 1:n_fisher){
          for (w in 1:n_spec){
            if (mean(fish_contrib[,x,-(year_start-1)])>0)
              diet_array[w,a,y-(year_start-1)]<-fish_contrib[w,x,y-(year_start-1)]/sum(fish_contrib[,x,y-(year_start-1)]) #Calculate proporitional contribution of a species to fishery
          }
        }
      }
    }
  }
  #************************Update 3/12/2021**************************
  #*This is the process to calculate Fij in RPath (https://github.com/NOAA-EDAB/Rpath) put
  #*in the context of this exercise.
  
  # Need biomass and Q/B
  biomass <- read.csv("biomass_annual.csv",skip=9)
  qb <- read.csv("consumption-biomass_annual.csv",skip=9)
  
  #* First columns are Year. The dimensions should be = [n_year,n_spec]
  biomass <- biomass[,-1] 
  qb <- qb[,-1]
  
  #* Calculate consumption
  bqb <- biomass * qb
  
  Prodij <- array(0,dim=dim(diet_array))
  Prodi <- matrix(0,nrow = dim(diet_array)[1],ncol=dim(diet_array)[3])
  Predcontij <- array(0,dim=dim(diet_array))
  Predcontji <- array(0,dim=dim(diet_array))
  direct <- array(0,dim=dim(diet_array))
  
  
  for (year in 1:n_year_concern){
    Prodij[,,year] <- diet_array[,,year] * bqb[col(as.matrix(diet_array[,,year]))]
    
    #* Replace final diet row with total fishery contribution
    #* This step is required because we previously converted the landings+discards into a percent
    if (n_fisher > 0){
      Prodij[(1:n_spec),((n_spec+1):(n_spec+n_fisher)),year]<-fish_contrib[,,year]
    }
    #* Net Production
    Prodi[,year] <- rowSums(Prodij[,,year])
    
    for(spec in 1:dim(Prodij)[1]){
      Predcontij[spec,,year] <- Prodij[spec,,year] / Prodi[spec,year]
    }
    Predcontij[,,year][is.na(Predcontij[,,year])] <- 0
  
    for (prey in 1:dim(Predcontij)[1]){
      for (pred in 1:dim(Predcontij)[2]){
        Predcontji[prey,pred,year] <- Predcontij[pred,prey,year] 
      }
    }
    Predcontji[is.na(Predcontji)] <- 0
    
  #* Direct MTI
    direct[,,year] <- diet_array[,,year] - Predcontji[,,year]
  }
  #************************************************************************

  setwd(outdir)
  #* If you are not interested in having excel files, comment those lines out
  #* Will significantly increase run time, but may be wise to have in case the loop cannot finish for any reason
  for (year in 1:n_year_concern){
    filename<-paste("MTI simulation ", n_count,".xlsx",sep="")
    identity.matrix <- diag(ncol(direct[,,year]))
    MTI[,,year,n_count] <- MASS::ginv(identity.matrix - direct[,,year]) - identity.matrix #Calculation of total MTI for each interaction
    if (year == 1){ #Year one creates the file
      sheet_name<-paste("MTI Year ",year,sep="")
      write.xlsx(MTI[,,year,n_count],file=filename,
                 sheetName = sheet_name,append=F)
    }
    if (year > 1){ #Years 2+ add sheets onto existing file
      sheet_name<-paste("MTI Year ",year,sep="")
      write.xlsx(MTI[,,year,n_count],file=filename,
                 sheetName = sheet_name,append=T)
    }
  }
  print(c(file_count)) #Will want to keep track of your loop#
  Sys.sleep(0.01)
}

#* Save all data to reload later
#* Will save you time instead of reloading all files
save(MTI,file="Mixed Trophic Impact Input.RData")

#****************************Done****************************************



#### Part 2: If you already have diet and predation matrices and just want to calculate MTI ####

  setwd(outdir)

  #* If you are not interested in having excel files, comment those lines out
  #* Will significantly increase run time, but may be wise to have in case the loop cannot finish for any reason
  for (year in 1:n_year_concern){
    filename<-paste("MTI simulation ", n_count,".xlsx",sep="")
    identity.matrix <- diag(ncol(direct[,,year]))
    MTI[,,year,n_count] <- MASS::ginv(identity.matrix - direct[,,year]) - identity.matrix #Calculation of total MTI for each interaction
    if (year == 1){ #Year one creates the file
      sheet_name<-paste("MTI Year ",year,sep="")
      write.xlsx(MTI[,,year,n_count],file=filename,
               sheetName = sheet_name,append=F)
    } else { #Years 2+ add sheets onto existing file
      sheet_name<-paste("MTI Year ",year,sep="")
      write.xlsx(MTI[,,year,n_count],file=filename,
               sheetName = sheet_name,append=T)
    }
  }
  print(c(file_count)) #Will want to keep track of your loop#
  Sys.sleep(0.01)
}