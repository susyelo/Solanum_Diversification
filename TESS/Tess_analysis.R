###################################################################################################
### Using the package TESS to assess diversification changes throught time in Solanum ---
###################################################################################################
library(TESS)
library(ape)
library(tools)

### Using the Sarkinen phylogeny
### NO idea why Ape is not recognizing the tree as Ultrametric. Before everything was working but after updgrading R the trees are recongnize as not ultrametric
Solanum_tree_sarkinen<-chronoMPL(read.tree("./Data/Solanum_Sarkinen_non_negative.tre"))
samplingFraction <- (Solanum_tree_sarkinen$Nnode + 1) / 1200

sarkinen_analysis<-tess.analysis(Solanum_tree_sarkinen,
                                 empiricalHyperPriors = TRUE,
                                 samplingProbability = samplingFraction,
                                 estimateNumberMassExtinctions = FALSE,
                                 MAX_ITERATIONS = 100000000,
                                 MAX_TIME = 24*60*60,
                                 MIN_ESS = 500,
                                 dir = "Output/Sarkinen_tree/")
numExpectedRateChanges=1
output_sarkinen <- tess.process.output("./Output/Sarkinen_tree/", Solanum_tree_sarkinen, numExpectedRateChanges = numExpectedRateChanges)


### ESS
ESS_lambda<-effectiveSize(output_sarkinen$numSpeciationCategories)
ESS_mu<-effectiveSize(output_sarkinen$numExtinctionCategories)

pdf("./Figures/TESS_Sarkinen_tree.pdf")
layout.mat <- matrix(1:4,nrow=2,ncol=2,byrow=TRUE)
layout(layout.mat)
tess.plot.output(output_sarkinen,
                 fig.types = c("speciation rates",
                               "speciation shift times",
                               "extinction rates",
                               "extinction shift times"),
                 las=2)
dev.off()



### Using the PASTIS trees
runs<-file_path_sans_ext(list.files("./Data/Pastis_trees/"))


for(i in 1:length(runs)){
  Solanum_tree_fully<-chronoMPL(read.tree(paste("./Data/Pastis_trees/",runs[i], ".tre", sep="")))
  samplingFraction <- (Solanum_tree_fully$Nnode + 1) / 1200
  tess.analysis(Solanum_tree_fully,
                empiricalHyperPriors = TRUE,
                samplingProbability = samplingFraction,
                estimateNumberMassExtinctions = FALSE,
                MAX_ITERATIONS = 100000000,
                MAX_TIME = 24*60*60,
                MIN_ESS = 500,
                dir = paste("./Output/Pastis_trees/",runs[i], sep=""))
  
}

##################################################################
## Reading output and check effective size for each run ---
##################################################################

E_size<-data.frame(Run=NULL,numSpeciationCategories=NULL,numExtinctionCategories=NULL)
numExpectedRateChanges=1

output_runs<-list()

for (i in 1:length(runs)){
  output_runs[[i]] <- tess.process.output(paste("./Output/Pastis_trees/",runs[i],sep="_"),
                                          numExpectedRateChanges = numExpectedRateChanges)
}

names(output_runs)<-paste("Run",runs,sep="_")

for (i in 1:length(runs)){
  
  TMP_DF<-data.frame(names(output_runs[i]))
  TMP_DF$numSpeciationCategories<-effectiveSize(output_runs[[i]]$numSpeciationCategories)
  TMP_DF$numExtinctionCategories<-effectiveSize(output_runs[[i]]$numExtinctionCategories)
  
  E_size<-rbind(E_size,TMP_DF)
}

##############################################################################
### Joining all the mcmc runs into a single one with pseudoposterior samples ---
##############################################################################

### 1. Reading the files of the first run and append the files of the other runs

files<-list.files("./Output/Pastis_trees/run_1/", pattern="*.txt")

file_tmp<-lapply(list.files("./Output/Pastis_trees/run_1/", pattern="*.txt", full.names = TRUE), readLines)
names(file_tmp)<-list.files("./Output/Pastis_trees/run_1/", pattern="*.txt")

for (dir in 2:length(runs)){
  for (file in files){
    file2<-readLines(paste("./Output/Pastis_trees/run_",runs[dir],"/",file,sep=""))[-1]
    file_tmp[[file]]<-append(file_tmp[[file]],file2)
    
  }
}

## 2. Writing the final output files ---
for (file in files){
  write(file_tmp[[file]], file=(paste("./Output/Pastis_trees/all_runs/", file, sep="")))
}

## 3. Reading back the joining mcmc samples and plot the summary
Solanum_tree<-read.tree("./Data/Solanum_Sarkinen_non_negative.tre")
#Solanum_tree_fully<-read.tree("./Data/Input_BAMM2.5/1/tree.tre")

## 4. Reading the fusion of all the files ----
output_all <- tess.process.output("./Output/Pastis_trees/all_runs/", tree = Solanum_tree,
                                  numExpectedRateChanges = numExpectedRateChanges)


pdf("./Figures/TESS_Solanum.pdf")
layout.mat <- matrix(1:4,nrow=2,ncol=2,byrow=TRUE)
layout(layout.mat)
tess.plot.output(output_all,
                 fig.types = c("speciation rates",
                               "speciation shift times",
                               "extinction rates",
                               "extinction shift times"),
                 las=2)
dev.off()