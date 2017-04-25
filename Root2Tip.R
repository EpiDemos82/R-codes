library(ape)
library(phytools)

setwd("/RAxML/") #Specify your working directory where all of your ML trees are located

#Getting tree files
tree.files <- list.files(".", pattern=".tre") #obtaining names of all trees in folder
tree.names<-gsub('.tre',"", tree.files) ##Remove all of the extensions off the files

#Setting up a matrix to store correlation data
num.trees <- length(tree.files) #obtaining the number of trees
corr.results <- matrix(nrow = num.trees, ncol = 7) #creating matrix to store R2T results

##Setting up functions for script ----
## Find the best location on that edge.
f <- function(x) {
  dist <- x * dist[best.edge.parent, ] + (1 - x) * dist[best.edge.child, ]
  objective(tip.dates, dist)
}

## Do root-to-tip regressions for every possible choice of root,
## that is, compute the objective function for every choice of root.
choice.f <- function(row) {
  if(row %in% missing.indices) # don't compute the objective function for tips that're missing data
    return (-Inf)
  else
    return (objective(tip.dates[valid.indices], dist[row, valid.indices])) # Only do the regression over data that exist
}

#Pick function to use (use the correlation for most)
##You can change the method --->>> This could be added to a funcition
#objective <- function(x, y) cor.test(y, x)$estimate #correlation
objective <- function(x, y) summary(lm(y ~ x))$r.squared #rsquared
#objective <- function(x, y) -summary(lm(y ~ x))$sigma^2 #rms

count=1 #counter
pdf(file='RootToTip2.pdf')
####looping over all trees----
for(i in tree.names){
  treepath <- file.path(".",paste(i,".tre",sep=""))  
  #input file
  #tree <- read.newick(file="/Users/tazarian/Documents/BrianColab/R2T/full_set.revised.final_tree_renamed.tre") #This is for testing a tree 
  tree <- read.newick(file=treepath) # I updated this to "read.newick" which is more robust then read tree
  tips <- tree$tip.label #obtain tip labels from phylogeny and abstract year from last time point
  tip.info <- strsplit(tips, "_") # split string by '_'
  tip.dates <- as.numeric(sapply(tip.info, "[", 3))  #THIS NEEDS TO BE CHANGED BASED ON NAMING
  rooted.tre <- rtt(tree, tip.dates, ncpu = 1, objective = "correlation", opt.tol = .Machine$double.eps^0.25)  
  
  missing.indices <- which(is.na(tip.dates))
  valid.indices <- which(!is.na(tip.dates))
  tip.lengths <- node.depth.edgelength(rooted.tre)
  dist <- dist.nodes(rooted.tre)[, 1:(rooted.tre$Nnode + 1)] 
  
  # Apply the objective function
  fits <- unlist(lapply(1:nrow(dist), choice.f))
  
  ## Find the best one (highest value of objective function).
  fit.edge <- apply(rooted.tre$edge, 2, function(e) fits[e])
  obj.edge <- apply(fit.edge, 1, mean)
  
  ## Compatibility with Path-O-Gen: find the last maximum, not the first.
  best.edge <- length(obj.edge) - which.max(rev(obj.edge)) + 1
  best.edge.parent <- rooted.tre$edge[best.edge, 1]
  best.edge.child <- rooted.tre$edge[best.edge, 2]
  best.edge.length <- rooted.tre$edge.length[best.edge]
  
  best.pos <- optimize(f, c(0, 1), maximum = TRUE, tol = .Machine$double.eps^0.25)$maximum
  
  ## Reroot the tree at the optimal location
  new.root <- list(edge = matrix(c(2L, 1L), 1, 2), tip.label = "new.root", edge.length = 1, Nnode = 1L, root.edge = 1)
  class(new.root) <- "phylo"
  t <- bind.tree(rooted.tre, new.root, where = best.edge.child, position = best.pos * best.edge.length)
  t <- collapse.singles(t)
  t <- root(t, "new.root")
  t <- drop.tip(t, "new.root")
  
  # Construct the linear model of the data
  tip.lengths <- node.depth.edgelength(t)
  
  distances <- tip.lengths[valid.indices]
  times <- tip.dates[valid.indices]
  model <- lm(times ~ distances)
  
  #Plot diagnostics from regression
  #layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page 
  #plot(model)
  #dev.off()
  
  ##Checking fit for linear regression of Root-to-Tip----
  #coefficients(model) # model coefficients
  #confint(model, level=0.95) # CIs for model parameters 
  #anova(model) # anova table 
  AnovaTable <- anova(model) #creating the ANOVA table to get the p-value for the reg
  
  #Variables to print on tree - p-value, correlation coeff, r2, and RMS
  model<-lm(times ~ distances)
  pval <- AnovaTable[1,5] #pvalue for root-to-tip correlation
  corl <- cor.test(times,distances)$estimate #correlation
  rsqaured <- summary(lm(times ~ distances))$r.squared #rsquared
  rms <- summary(lm(times ~ distances))$sigma^2 #rms
  slope <- summary(lm(times ~ distances))$coefficients[2]
  CI1 <- round(confint(model, level = 0.95)[2], 2) #Upper CI bound
  CI2 <- round(confint(model, level = 0.95)[4], 2) #Lower CI bound
  
  ##separate figures ----
  #dev.off()
  #plot(rooted.tre, cex=.5, main="Sequence Cluster SC02_3")
  #mtext(paste("p-value",round(pval,3)),side=1,line=1)
  #mtext(paste("correlation",round(corl,3)),side=1,line=2)
  #mtext(paste("r-sqaured",round(rsqaured,3)),side=1,line=3)
  #mtext(paste("RMS",round(rms,3)),side=1,line=4)
  
  #dev.off()
  #plot(times,distances, ylab="Distance", xlab="Year")
  #abline(lm(distances~times), col="red") # regression line (y~x) 
  #lines(lowess(times,distances), col="blue") # lowess line (x,y)
  
  rooted.tre <- ladderize(rooted.tre)
  ##Layout of bothfigures - phylogeny and root to tip----
  par(mfrow=c(1,2))
  plot(rooted.tre, cex=.3)
  mtext(paste("p-value",round(pval,3)),side=1,line=0)
  mtext(paste("correlation",round(corl,3)),side=1,line=1)
  mtext(paste("r-sqaured",round(rsqaured,3)),side=1,line=2)
  mtext(paste("RMS",round(rms,3)),side=1,line=3)
  mtext(paste("rate",round(slope,3)," (95% CI ",CI1,"-",CI2,")"),side=1,line=4)
  plot(times,distances, ylab="Distance", xlab="Year")
  abline(lm(distances~times), col="red") # regression line (y~x) 
  lines(lowess(times,distances), col="blue") # lowess line (x,y)
  mtext(paste("Population: ",i), side=3, outer=TRUE, line=-3)
  
  #Filling the matrix
  corr.results[count,1] <- i
  corr.results[count,2] <- round(corl,3)
  corr.results[count,3] <- round(pval,3)
  corr.results[count,4] <- round(rsqaured,3)
  corr.results[count,5] <- round(rms,3)
  corr.results[count,6] <- round(slope,3)
  corr.results[count,7] <- paste("(",CI1,"-",CI2,")", sep = "")
   
  count=count+1
  cat(".")
  
  #remove('treepath','tree','tips','tip.info','tip.dates','rooted.tre','missing.indices','valid.indices','tip.lengths','dist','fits','fit.edge','obj.edge','new.root','t','best.edge','best.edge.parent','best.edge.child','best.edge.length','best.pos','distances','times','model','AnovaTable','pval','corl','rsqaured','rms')
}
dev.off()

colnames(corr.results) <- c("Clade","corl","pval","rsquared","rms","slope","CI")
corr.results
write.table(corr.results, 'RootToTip2.txt', sep="\t")
