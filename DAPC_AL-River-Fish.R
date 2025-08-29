library(adegenet)
library(ape)
library(poppr)
library(na.tools)
rm(list=ls())
setwd("~/Fish/AL-River_migratory-fish/white-crappie/white_crappie_m3M1n1_p4r7_mac3_multi/")

#rename any .genepop to .gen to make packages happy
genepops<-list.files(pattern="*.genepop")
gens<-gsub(".genepop$",".gen", genepops)
file.rename(genepops,gens)


data <- read.genepop(file="populations.haps.gen")

#Determine number of clusters
grp <- find.clusters(data, max.n.clust=5, n.pca = 250, n.clust = 2)

dapc <- dapc(data, grp$grp)
pdf("DAPC_clustered-by-population.pdf", width = 5, height = 5)
scatter.dapc(dapc,
             scree.da = TRUE,
             scree.pca = TRUE,
             posi.pca = "bottomright", 
             posi.da = "top", 
             cstar = 0,
             clab = 0,
             cell =0,
             col=c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","grey","darkblue","green"), 
             grp=data$pop,
             posi.leg="topright",
             leg=TRUE,
             
)
dev.off()
pdf
scatter.dapc(dapc,
             scree.da = TRUE,
             scree.pca = TRUE,
             posi.pca = "bottom", 
             posi.da = "top", 
             cstar = 0,
             clab = 0,
             cell =0,
             col=c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","grey","darkblue","green"), 
             grp=data$pop,
             posi.leg="bottomright",
             leg=TRUE,
             
)
pdf("DAPC_plotted-by-cluster.pdf")
scatter.dapc(dapc,
             scree.da = TRUE,
             scree.pca = TRUE,
             posi.pca = "topright", 
             posi.da = "top", 
             cstar = 0,
             clab = 0,
             cell =0,
             col=c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","grey","darkblue","green"), 
             grp=grp$grp,
             posi.leg="bottomright",
             leg=TRUE,
             
)
dev.off()


####Freshwater Drum#############
rm(list=ls())
setwd("~/Fish/AL-River_migratory-fish/freshwater-drum/freshwater-drum_m3M2n2_p4r7_mac3/")
data <- read.genepop(file="freshwater-drum_m3M2n2_p4r7_mac3.haps.gen")

#Determine number of clusters
grp <- find.clusters(data, max.n.clust=5, n.pca = 250)
#grp <- find.clusters(data, max.n.clust=5, n.pca = 100, n.clust = 3) ##values based on perlim run
#dapc
#dapc <- dapc(data, grp$grp, n.pca = 2, n.da = 2)
dapc <- dapc(data, grp$grp)
pdf("DAPC_clustered-by-population.pdf", width = 5, height = 5)
scatter.dapc(dapc,
             scree.da = TRUE,
             scree.pca = TRUE,
             posi.pca = "bottomright", 
             posi.da = "top", 
             cstar = 0,
             clab = 0,
             cell =0,
             col=c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","grey","darkblue","green"), 
             grp=data$pop,
             posi.leg="topright",
             leg=TRUE,
             
)
dev.off()
pdf
scatter.dapc(dapc,
             scree.da = TRUE,
             scree.pca = TRUE,
             posi.pca = "bottom", 
             posi.da = "top", 
             cstar = 0,
             clab = 0,
             cell =0,
             col=c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","grey","darkblue","green"), 
             grp=data$pop,
             posi.leg="bottomright",
             leg=TRUE,
             
)
pdf("DAPC_plotted-by-cluster.pdf")
scatter.dapc(dapc,
             scree.da = TRUE,
             scree.pca = TRUE,
             posi.pca = "topright", 
             posi.da = "top", 
             cstar = 0,
             clab = 0,
             cell =0,
             col=c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","grey","darkblue","green"), 
             grp=grp$grp,
             posi.leg="bottomright",
             leg=TRUE,
             
)
dev.off()
