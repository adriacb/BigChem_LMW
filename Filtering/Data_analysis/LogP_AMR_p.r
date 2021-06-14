library(BiocManager)
library(ggplot2)
library(rJava)
library(rcdklibs)
library(rcdk)
library(Rcpi)
library(readr)
library(foreach)
library(doParallel)
library(plyr)
inputdir = ''
workingdir = ''

setwd(workingdir)
start_time <- Sys.time()
X100smi <- read_csv(inputdir, 
                    col_names = FALSE)

#setup parallel backend to use many processors
cl <- makeCluster(6)
registerDoParallel(cl)

results<- foreach(i = 1:(length(X100smi$X1)), .packages = "tcltk", .combine='rbind') %dopar% 
  {
    require(rcdk)
    require(Rcpi)
    if(!exists("pb")) pb <- tkProgressBar("P task", min=1, max = (length(X100smi$X1)))
    setTkProgressBar(pb, i)
    log2(i)
    mol = parse.smiles(as.character(X100smi$X1[i]))
    tempDF <- extractDrugALOGP(mol)
    ALogP <- tempDF[,1]
    AMR <- tempDF[,3]
    return (data.frame(ALogP, AMR))
  }

stopCluster(cl)

#results
# Histogram with density plot

p1<-ggplot(results, aes(x=ALogP)) + 
  geom_histogram(aes(y=..density..), fill ="#69b3a2", bins = 20)+
  geom_density(alpha=.1, fill="#FF6666")+ggtitle("A. LogP") +theme_classic()
png("ALogP.png")
print(p1)
dev.off()

p2<-ggplot(results, aes(x=AMR)) + 
  geom_histogram(aes(y=..density..), fill = "#404080", bins = 20)+
  geom_density(alpha=.1, fill="#FF6666")+ggtitle("Molar Refractivity") + theme_classic()
png("AMR.png")
print(p2)
dev.off()


end_time <- Sys.time()
end_time - start_time