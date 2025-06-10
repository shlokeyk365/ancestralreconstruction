library(phytools)
library(ape)

dir.create("birthdeathmodels/plots", showWarnings = FALSE)

# simulate tree with 20 species
tree <- pbtree(n=20, scale=150)

# initial tree visualization
png("birthdeathmodels/plots/initial_tree.png", width=800, height=600)
par(mfrow=c(2,1))
plotTree(tree, ftype="i", mar=c(4.1,4.1,2.1,1.1), show.tip.label = TRUE, 
         cex = 0.8, label.offset = 0.5)

obj <- ltt(tree, plot="False")
abline(v=obj$times, lty="dotted", col=make.transparent("red",0.5))
axis(1, cex.axis=0.8)
mtext("a)", line=1, at=-10)

plot(obj, mar=c(5.1,4.1,2.1,1.1), bty="n", log.lineages=FALSE, 
     las=1, cex.axis=0.8)
abline(v=obj$times, lty="dotted", col=make.transparent("red",0.5))
mtext("b)", line=1, at=-10)
dev.off()

darter.tree <- tree

# fan plot
png("birthdeathmodels/plots/fan_tree.png", width=800, height=600)
plotTree(darter.tree, ftype="i", fsize=0.4, type="fan", lwd=1, part=0.88)

h <- max(nodeHeights(darter.tree))
obj <- axis(1, pos=-2, at=h-c(0,5,10,15,20,25), cex.axis=0.5, labels=FALSE)
text(obj, rep(-5,length(obj)), h-obj, cex=0.6)
text(mean(obj), -8, "time (mybp)", cex=0.8)
dev.off()

# ltt plot
darter.ltt <- ltt(darter.tree, plot=FALSE)

png("birthdeathmodels/plots/ltt_plot.png", width=800, height=600)
par(mar=c(5.1,4.1,2.1,2.1))
plot(darter.ltt, log.lineages=FALSE, log="y", col="red", lwd=2, 
     bty="n", las=1, cex.axis=.8)
dev.off()

# resolve polytomies
darter.tree <- multi2di(darter.tree)
darter.ltt <- ltt(darter.tree, plot=FALSE)

# simulate different scenarios
tree.noExtinction <- pbtree(b=0.05, n=150, t=150, method="direct")
tree.withExtinction <- pbtree(b=0.25, d=0.20, n=150, t=150, method="direct")
tree.reconstructed <- drop.tip(tree.withExtinction, getExtinct(tree.withExtinction))
tree.reconstructed$root.edge <- 150 - max(nodeHeights(tree.reconstructed))
tree.reconstructed <- rootedge.to.singleton(tree.reconstructed)

ltt.noE <- ltt(tree.noExtinction, plot=FALSE)
ltt.wE <- ltt(tree.withExtinction, plot=FALSE)
ltt.recon <- ltt(tree.reconstructed, plot=FALSE)

# comparison of scenarios
png("birthdeathmodels/plots/comparison_plot.png", width=800, height=600)
par(lend=1, mar=c(5.1,4.1,2.1,2.1))
plot(ltt.noE, bty="n", log.lineages=FALSE, log="y", lwd=2, 
     xlim=c(0,160), las=1, cex.axis=0.8)
plot(ltt.wE, log.lineages=FALSE, lty="dotted", col="red", lwd=2, add=TRUE)
plot(ltt.recon, log.lineages=FALSE, lwd=2, add=TRUE, col="darkred")

legend(x="topleft", lty=c("solid", "dotted", "solid"), lwd=2,
       col=c("black","red","darkred"), 
       legend=c("no extinction", "extinction (w/ extinct lineages)", 
                "extinction (reconstructed)"), bty="n", cex=0.7)
dev.off()

# fit birth-death model
bd.model <- fit.bd(darter.tree)
print("birth-death model results:")
print(bd.model)

# adjust for sampling
sampling.f <- 0.85
bd.model <- fit.bd(darter.tree, rho=sampling.f)
print("birth-death model with sampling:")
print(bd.model)

tree.reconstructed <- collapse.singles(tree.reconstructed)
tree.missing <- drop.tip(tree.reconstructed, 
                        sample(tree.reconstructed$tip.label, 75))

ltt.recon <- ltt(tree.reconstructed, plot=FALSE)
ltt.missing <- ltt(tree.missing, plot=FALSE)

# final comparison
png("birthdeathmodels/plots/final_comparison.png", width=800, height=600)
par(mar=c(5.1,4.1,2.1,2.1), lend=2)
plot(ltt.recon, bty="n", log.lineages=FALSE, log="y", lwd=2, 
     col="darkred", las=1, cex.axis=0.8)
plot(ltt.missing, log.lineages=FALSE, lty="dotted", lwd=2, add=TRUE)

legend(x="topleft", lty=c("solid","dotted"), lwd=c(2,2), 
       col=c("darkred","black"), 
       legend=c("reconstructed phylogeny", "phylogeny with missing taxa"),
       bty="n", cex=0.8)
dev.off()

cat("\nplots saved to birthdeathmodels/plots/\n")
