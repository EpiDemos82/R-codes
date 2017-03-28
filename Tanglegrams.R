#Tanglegrams (Need to clean up!)
library(dendextendRcpp)
library(dendextend)
library(ape)
library(phytools)
library(phangorn)

#Treescape comparison of tree space
library(treescape)
require(RColorBrewer)
require(scatterD3)
require(ggplot2)
require(rgl)
require(sandwich)
require(msm)
require(reshape2)
library(ggrepel)

setwd("/Trees")

#Importing core genome and accessory (binary) trees inferred using RAxML
SC.core <- read.tree(file=file.path("CoreML.tre"))
SC.acc <- read.tree(file=file.path("AccML.tre"))


####**Manipulating trees for KC comparisons and randomization**####

#laderizing core and access trees (sorting by branch length)
SC.core <- midpoint(SC.core)
SC.acc <- midpoint(SC.acc)
SC.core <- ladderize(SC.core)
SC.acc <- ladderize(SC.acc)

#Plotting the ladderized trees - just for a quick check
plot(SC.core, show.tip.label = TRUE, cex=.3)
plot(SC.acc, show.tip.label = TRUE, cex=.3)

#Generating a random-binary-ultrametric tree with the same number of tips and the same tip lables
tiplabs <- SC.core$tip.label #tip labels 
n.tips <- length(tiplabs) #number of taxa 
ran.utree.binary <- rcoal(n.tips, tip.label = tiplabs, br = "coalescent") #generating nuetral coalescent tree with n taxa
SC.ran <- rtree(n.tips, tip.label = tiplabs, rooted = TRUE) 
SC.ran <- midpoint(SC.ran) #midpoint rooting
SC.ran <- ladderize(SC.ran) #ladderizing 

#TESTING SIGNIFICANCE BY RANDOMLY GENERATING 1000 TREES AND COMPARING KC METRICS 
KC.core.acc<-treeDist(SC.core, SC.acc) #KC of core and accessory
KC.core.ran<-treeDist(SC.core, SC.ran) #KC of core and random tree
#setting up the matrix to store results and looping the 1000 comparisons of random trees to the core tree
randomization <-matrix(nrow=1000,ncol=1)
for (n in 1:1000){
  SC.ran <- rtree(n.tips, tip.label = tiplabs, rooted = TRUE)
  SC.ran <- midpoint(SC.ran)
  SC.ran <- ladderize(SC.ran)
  randomization[n,] <- treeDist(SC.core, SC.ran)
  cat(".")
}
mean <- round(mean(randomization), digits=1) #mean KC statistics from randomizations
hist(randomization) #distribution of KC statistics from tree randomization comparisons 
error <- qt(0.975,df=length(randomization)-1)*sd(randomization)/sqrt(length(randomization))
left <- round(mean(randomization)-error, digits = 1) #Standard errors 
right <- round(mean(randomization)+error, digits = 1)
#p-value
p.value <- length(randomization[randomization<KC.core.acc])/1000 #calculating the p-value
p.value <- ifelse(p.value == 0,0.001,p.value) #this is for instances when p <<0.001
####-----------------------####

####**TRANSFORMING TREES AND CREATING TANGLEGRAMS**####

#Converting binary and access trees to ultrametric - requirement for tanglgram
core.utree = chronos(SC.core, lambda = 0, model = "correlated")
acc.utree = chronos(SC.acc, lambda = 0, model = "correlated")

#randomly resolve polytomies of binary and access trees
core.utree.binary <- multi2di(core.utree, random = FALSE)
acc.utree.binary <- multi2di(acc.utree, random = FALSE)

#Confriming that trees are ultrametric and binary (no polytomies)
is.ultrametric(core.utree.binary)
is.ultrametric(core.utree.binary)
is.ultrametric(acc.utree.binary)
is.binary.tree(acc.utree.binary)

#Plotting final untrametric, binary core/access/random trees for a quick look
#plot(core.utree.binary, show.tip.label = TRUE, cex=.3)
#plot(acc.utree.binary, show.tip.label = TRUE, cex=.3)
#plot(ran.utree.binary, show.tip.label = TRUE, cex=.3)


pdf(file='Core_Accessory_Comp_Tangle.pdf') #PDF of tanglegrams
###**TANGLEGRAMS**####
#Tanglegrams of Core-Access-Random
tanglegram(core.utree.binary, acc.utree.binary,
           main_left = "Core",
           main_right = "Accessory",
           common_subtrees_color_branches=TRUE,
           axes=FALSE,
           remove_nodePar=TRUE,
           sort = TRUE, k_branches = as.integer(n.tips),
           lab.cex = 0.001, edge.lwd = 1,  
           margin_inner= 3.5)
mtext(paste("Kendall-Colijn = ",round(KC.core.acc,1),"(p<",p.value,")"),side=1,line=2)

tanglegram(core.utree.binary, ran.utree.binary,
           main_left = "Core",
           main_right = "Neutral Coalescent",
           common_subtrees_color_branches=TRUE,
           cex_main_left = 1.5,
           cex_main_right = 1.5,
           axes=FALSE,
           remove_nodePar=TRUE,
           sort = TRUE, k_branches = as.integer(n.tips),
           lab.cex = 0.001, edge.lwd = 1,  
           margin_inner= 3.5)
mtext(paste("Kendall-Colijn = ",mean," (95% CI:",left,"-",right,")",sep=""),side=1,line=2)

dev.off()
###**PRINTING KC COMPARISON STATS TO SCREEN FOR ASSESSMENT**####
round(KC.core.acc,1)
p.value
paste(mean," (",left,"-",right,")",sep="")
