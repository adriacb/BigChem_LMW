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

print("DF created.")
#setup parallel backend to use many processors
cl <- makeCluster(6)
registerDoParallel(cl)
total = (length(X100smi$X1))
results<- foreach(i = 1:total, .packages = "tcltk", .combine='rbind') %dopar% 
  {
    require(rcdk)
    require(Rcpi)
    if(!exists("pb")) pb <- tkProgressBar("P task", min=1, max = total/6)
    setTkProgressBar(pb, i)
    log2(i)
    mol = parse.smiles(as.character(X100smi$X1[i]))
    tempDF <- extractDrugHBondAcceptorCount(mol)
    tempDF1 <- extractDrugHBondDonorCount(mol)
    HBA <- tempDF[,1]
    HBD <- tempDF1[,1]
    return (data.frame(HBA, HBD))
  }

stopCluster(cl)
print("Clusters stoped.")
print("Removing DF.")
remove(X100smi)
print("Starting plots.")


p1<-ggplot(results, aes(x=HBA)) + 
  geom_histogram(aes(y=..count..), fill ="#69b3a2", bins = 6)+
  ggtitle("Hydrogen Bond Acceptors")+ theme_classic()
png("hba.PNG")
print(p1)
dev.off()

p2<-ggplot(results, aes(x=HBD)) + 
  geom_histogram(aes(y=..count..), fill = "#404080", bins = 5)+
  ggtitle("Hydrogen Bond Donors") + theme_classic()
png("hbd.png")
print(p2)
dev.off()

end_time <- Sys.time()
end_time - start_time
