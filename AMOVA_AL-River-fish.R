library(poppr)
library(ade4)



####Freshwater Drum multiple SNP####

rm(list=ls())
setwd("~/Fish/AL-River_migratory-fish/freshwater-drum/freshwater-drum_m3M2n2_p4r7_mac3/")
#rename any .genepop to .gen to make packages happy
genepops<-list.files(pattern="*.genepop")
gens<-gsub(".genepop$",".gen", genepops)
file.rename(genepops,gens)
drum.geneid = read.genepop("freshwater-drum_m3M2n2_p4r7_mac3.haps.gen")
drum.genclone<-as.genclone(drum.geneid)
drum.genclone$pop
strata(drum.genclone)<-(as.data.frame(drum.geneid$pop))
drum.genclone$strata

##AMOVA
drum.genclone
drum.site.amova.pop = poppr.amova(drum.genclone, ~drum.geneid.pop, cutoff = 0.5, method = "ade4")


####Print Results
drum.site.amova.pop

##Randomization Test
drum.site.amova.pop.rtest<-randtest(drum.site.amova.pop,nrepet = 1000)
drum.site.amova.pop.rtest
plot(drum.site.amova.pop.rtest)



####White crappie multiple SNP####

rm(list=ls())
setwd("~/Fish/AL-River_migratory-fish/freshwater-drum/freshwater-drum_m3M2n2_p4r7_mac3/")
#rename any .genepop to .gen to make packages happy
genepops<-list.files(pattern="*.genepop")
gens<-gsub(".genepop$",".gen", genepops)
file.rename(genepops,gens)
drum.geneid = read.genepop("freshwater-drum_m3M2n2_p4r7_mac3.haps.gen")
drum.genclone<-as.genclone(drum.geneid)
drum.genclone$pop
strata(drum.genclone)<-(as.data.frame(drum.geneid$pop))
drum.genclone$strata

##AMOVA
drum.genclone
drum.site.amova.pop = poppr.amova(drum.genclone, ~drum.geneid.pop, cutoff = 0.5, method = "ade4")


####Print Results
drum.site.amova.pop

##Randomization Test
drum.site.amova.pop.rtest<-randtest(drum.site.amova.pop,nrepet = 1000)
drum.site.amova.pop.rtest
plot(drum.site.amova.pop.rtest)
