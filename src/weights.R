args <- commandArgs(trailingOnly=T)

counts <- read.table(args[1], row.names=1)

weights <- read.table(args[2], row.names=1)

variances <- read.table(args[3], row.names=1)

num_individuals <- as.numeric(args[4])

pdf(args[5])

par(mfrow=c(2,2))

hist(2 * counts[,1]/num_individuals, breaks="FD", xlim=c(0,25), xlab="average normalized counts", main="Per k-mer average normalized counts")

hist(weights[,1], xlab="weights", main="Per k-mer inferred weights", breaks="FD", xlim=c(0,6.0))

hist(variances[,1]/(2 * counts[,1]), breaks="FD", xlim=c(0,0.002), xlab="variance/count", main="Per k-mer ratio of variance\n to average count")

hist(variances[,1], xlab="variance", main="Per k-mer count variance", breaks="FD", xlim=c(0,3))

dev.off()