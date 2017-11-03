library(ggplot2)
library(plyr)

args <- commandArgs(trailingOnly = TRUE)

data=read.table(args[1], sep="\t", header=FALSE, check.names=FALSE, stringsAsFactors=F)
colnames(data) <- c("chr", "position", "count")
head(data)
mpileup_data = read.table(args[2], sep="\t", header=FALSE, check.names=FALSE, stringsAsFactors=F, quote="", comment.char="")
colnames(mpileup_data) <- c("chr", "position", "nt", "count", "seq", "seq_q")
head(mpileup_data)
coverage <- data[, c("position", "count")]
coverage$nt <- mpileup_data$nt

## ATCG average coverage
tmp <- data.frame(table(coverage$nt))
colnames(tmp) <- c("nt", "count")
tmp2 <- aggregate(coverage$count, by=list(nt=coverage$nt), FUN=sum)
colnames(tmp2) <- c("nt", "total_cov")
tmp2 <- merge(tmp2, tmp, by="nt")
tmp2$average_cov <- tmp2$total_cov/tmp2$count
tmp2$average_cov <- tmp2$average_cov/tmp2$average_cov[tmp2$nt=="A"]

write.table(tmp2, file="coverage_nt.tsv", sep="\t", quote=F, row.names=F)

# plot
set_size <- 15
p <- ggplot(tmp2, aes(x=nt, y=average_cov, fill = tmp2$nt)) +
    geom_bar(stat="identity", color = "white") +
 	scale_fill_manual(values=c("A"="#3cba54", "T"="#f4c20d", "C"="#db3236", "G"="#4885ed")) +
    ggtitle("Average Coverage") +
    xlab("") +
	ylab("Normalized average coverage") +
    theme(text=element_text(colour="black", size=set_size)) +
    theme(plot.title = element_text(size=set_size)) +
    theme(legend.text=element_text(size=set_size), legend.title=element_text(size=set_size, face="bold")) +
	theme(legend.position="none")

ggsave(plot=p, height=5, width=5, dpi=300, filename="avg_cov.png")

