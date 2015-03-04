args <- commandArgs(trailingOnly=T)

counts <- read.table(args[1], row.names=1)

weights <- read.table(args[2], row.names=1)

variances <- read.table(args[3], row.names=1)

num_individuals <- as.numeric(args[4])

pdf(args[5])

par(mfrow=c(2,2))

hist(counts[,1]/num_individuals, breaks="FD", xlim=c(0,4), xlab="average normalized counts", main="Per k-mer average normalized counts")

hist(weights[,1], xlab="weights", main="Per k-mer inferred weights", breaks="FD", xlim=c(0,4))

hist(variances[,1]/counts[,1], breaks="FD", xlim=c(0,0.005), xlab="variance/count", main="Per k-mer ratio of variance\n to average count")

hist(variances[,1], xlab="variance", main="Per k-mer count variance", breaks="FD", xlim=c(0,4))

dev.off()