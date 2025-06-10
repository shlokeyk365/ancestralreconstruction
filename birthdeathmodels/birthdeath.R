#loading the needed libraries for the phylo analysis
library(phytools)
library(ape)

#generate random phylogenetic tree w/ 12 species 
tree<-pbtree(n=12,scale=100)
par(mfrow=c(2,1))
plotTree(tree,ftype="off",mar=c(4.1,4.1,2.1,1.1))

#compute the Lineage-Through-Time (LTT) plot data for the tree without directly plotting it. 
obj<-ltt(tree,plot="False")
abline(v=obj$times,lty="dotted",
       col=make.transparent("blue",0.5))
axis(1,cex.axis=0.8)
mtext("a)",line=1,at=-10)

#Plot the LTT graph based on the computer data (obj)
plot(obj,mar=c(5.1,4.1,2.1,1.1),bty="n",
     log.lineages=FALSE,las=1,cex.axis=0.8)
abline(v=obj$times,lty="dotted",
       col=make.transparent("blue",0.5))
mtext("b)", line=1,at=-10)

#load phylogeny tree for darters from local file 
darter.tree<-read.tree("etheostoma_percina_chrono.tre")
#plot the darter tree using the ape packages plot function, and not phytools due to many compile errors and not working with my version of R studio
plotTree(darter.tree,ftype="i",
         fsize=0.4,type="fan",lwd=1,part=0.88)

getwd()
setwd("/Users/shlokkulkarni/Downloads")
darter.tree<-read.tree("etheostoma_percina_chrono.tre")

class(darter.tree)
str(darter.tree)
plot.phylo(darter.tree, type="fan", cex=0.6)

#compute heights of darter tree/annotating temporal axis
h<-max(nodeHeights(darter.tree))
obj<-axis(1,pos=-2,at=h-c(0,5,10,15,20),
          cex.axis=0.5,labels=FALSE)
text(obj,rep(-5,length(obj)),h-obj,
     cex=0.6)
text(mean(obj),-8,"time (mybp)",cex=0.8)

#compute and plot the LTT graph for the darter tree on a logarithmic scale
darter.ltt<-ltt(darter.tree, plot=FALSE)
par(mar=c(5.1,4.1,2.1,2.1))
plot(darter.ltt,log.lineages=FALSE,log="y",
     col="blue",lwd=2,bty="n",las=1,
     cex.axis=.8)

#extra steps taken by paper to convert darter tree to fully bifurcating tree by resolving polytomies 
darter.tree<-multi2di(darter.tree)
darter.ltt<-ltt(darter.tree,plot=FALSE)
darter.ltt


#estimating speciation/extinction rates from the simulated trees with and without extinction

tree.noExtinction<-pbtree(b=0.039,n=100,t=100,
                          method="direct")
tree.withExtinction<-pbtree(b=0.195,d=0.156,
                            n=100,t=100,method="direct")
tree.reconstructed<-drop.tip(tree.withExtinction,
                             getExtinct(tree.withExtinction))
tree.reconstructed$root.edge<100-
  max(nodeHeights(tree.reconstructed))
tree.reconstructed<-rootedge.to.singleton(
  tree.reconstructed
)

#create the ltt objects from each tree: no extinction, with extinction, and reconstructed tree
ltt.noE<-ltt(tree.noExtinction,plot=FALSE)
ltt.wE<-ltt(tree.withExtinction,plot=FALSE)
ltt.recon<-ltt(tree.reconstructed,plot=FALSE)

#graphing the LTTs - for the different scenariors to compare patterns of speciation and extinction
par(lend=1,mar=c(5.1,4.1,2.1,2.1))
plot(ltt.noE,bty="n",log.lineages=FALSE,log="y",
     lwd=2,xlim=c(0,110),las=1,cex.axis=0.8)
plot(ltt.wE,log.lineages=FALSE,lty="dotted",
     col="black",lwd=2,add=TRUE)
plot(ltt.recon,log.lineages=FALSE,
     lwd=2,add=TRUE,col="darkgrey")

#adding a legend to the plot describing each scenario
legend(x="topleft",lty=c("solid", "dotted", "solid"),
       lwd=2,col=c("black","black","darkgrey"),
       legend=c("no extinction",
                "extinction(w/ extinct lineages",
                "extinction (reconstructed)"),
       bty="n",cex=0.7)


#fitting the birth death model to the darter tree and adjusting for incomplete sampling. 
bd.model<-fit.bd(darter.tree)
bd.model

sampling.f<-201/216
sampling.f
bd.model<-fit.bd(darter.tree,
                 rho=sampling.f)
bd.model

tree.reconstructed<-collapse.singles(tree.reconstructed) #simulate missing taxa by randomly dropping tips from reconstructed tree
tree.missing<-drop.tip(tree.reconstructed,
                       sample(tree.reconstructed$tip.label,50))
#create ltt objects for the reconstruted tree and the tree with missing taxa
ltt.recon<-ltt(tree.reconstructed,plot=FALSE)
ltt.missing<-ltt(tree.missing,plot=FALSE)
#plot the LTT graphs for the reconstructed tree and the tree with missing taxa
par(mar=c(5.1,4.1,2.1,2.1),lend=2)
plot(ltt.recon,bty="n",log.lineages=FALSE,log="y",
     lwd=2,col="darkgrey",las=1,cex.axis=0.8)
plot(ltt.missing,log.lineages=FALSE,lty="dotted",
     lwd=2,add=TRUE)

#adding a legend to show two scenarios: reconstructed vs missing taxa
legend(x="topleft",lty=c("solid","dotted"),
       lwd=c(2,2),col=c("darkgrey","black"),
       legend=c("reconstructed phylogeny",
                "phylogeny with missing taxa"),
       bty="n",cex=0.8)

