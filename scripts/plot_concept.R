png("/home/himanshu/.gemini/antigravity-ide/brain/d6303ffc-29ea-4fc1-9499-09d29b4a7ebc/signal_to_noise_concept.png", width=2200, height=1500, res=300)
par(mar=c(5, 5, 4, 2), bg="#ffffff")

set.seed(42)
igg <- rnorm(500, mean=19.0, sd=0.4)
trf2 <- rnorm(500, mean=28.0, sd=0.4)

h1 <- hist(igg, plot=FALSE, breaks=30)
h2 <- hist(trf2, plot=FALSE, breaks=30)

plot(h1, col=rgb(0.19, 0.51, 0.74, 0.6), xlim=c(17, 30), ylim=c(0, max(h1$counts, h2$counts)*1.15),
     main="Welch's t-Test Concept: High Signal vs Low Noise (p < 0.0001)",
     xlab="Log2 LFQ Intensity", ylab="Frequency (Protein Replicates)",
     cex.main=1.1, cex.lab=1.0)
plot(h2, col=rgb(0.87, 0.18, 0.15, 0.6), add=TRUE)

abline(v=mean(igg), col="#08519c", lty=2, lwd=2)
abline(v=mean(trf2), col="#a50f15", lty=2, lwd=2)

arrows(x0=mean(igg), y0=55, x1=mean(trf2), y1=55, code=3, length=0.1, lwd=2, col="#252525")
text(x=(mean(igg)+mean(trf2))/2, y=58, labels="SIGNAL (Log2FC Difference = +9.0)", font=2, cex=0.8)

legend("topright", legend=c("IgG Control (Low/Noise)", "TERF2 Target (High/Signal)"),
       fill=c(rgb(0.19, 0.51, 0.74, 0.6), rgb(0.87, 0.18, 0.15, 0.6)), bty="n")

dev.off()
cat("Concept diagram created cleanly in R!\n")
