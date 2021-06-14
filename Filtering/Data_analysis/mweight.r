library(BiocManager)
library(ggplot2)
library(rJava)
library(rcdklibs)
library(rcdk)
library(Rcpi)
library(readr)
library(foreach)
library(doParallel)

inputdir = ''
workingdir = ''

setwd(workingdir)
start_time <- Sys.time()
X100smi <- read_csv(inputdir, 
                    col_names = FALSE)
# MW DF

#setup parallel backend to use many processors
cl <- makeCluster(6)
registerDoParallel(cl)
total = (length(X100smi$X1))
results<- foreach(i = 1:total, .packages = "tcltk", .combine='rbind') %dopar% 
  {
    require(rcdk)
    require(Rcpi)
    if(!exists("pb")) pb <- tkProgressBar("P task", min=1, max = total)
    setTkProgressBar(pb, i)
    log2(i)
    mol = parse.smiles(as.character(X100smi$X1[i]))
    tempDF <- extractDrugWeight(mol)
    MW <- tempDF[,1]
    return (data.frame(MW))
  }

stopCluster(cl)

print("Clusters stoped.")
print("Removing DF.")
print("Starting plots.")

p<-ggplot(results, aes(x=MW)) + 
  geom_histogram(aes(y=..count..), fill ="#69b3a2", bins = 20)+
  ggtitle("Molecular Weight")+ theme_classic()
png("weight.png")
print(p)
dev.off()

end_time <- Sys.time()
end_time - start_time