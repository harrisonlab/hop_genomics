#!/usr/bin/env Rscript
#install.packages("tidyverse")
#install.packages("ggpubr")
#library('tidyverse')
#library("ggpubr")

library('ggplot2')
library('tools')
library('tidyverse')
library('ggpubr')


args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = 16
  args[3] = 8
}
output <- paste(file_path_sans_ext(args[1]), "_snp_quality_stats.pdf", sep="")
df1 = data.frame(read_tsv(args[1]))
head(df1)


df1$PASS <- df1$FILTER
df1$PASS[df1$PASS!="PASS"] <- "FAIL"
df1$PASS <- as.factor(df1$PASS)
# Histogram with density plot
#    geom_histogram(aes(y=..density..), alpha=0.2, bins=50)

p0 <- ggplot(df1, aes(x=QD)) +
  geom_density(alpha=.2) +
  theme(legend.position="none")
ggsave("../QD_all.pdf", p0)
p1 <- ggplot(df1, aes(x=QD, colour=PASS, fill=PASS)) +
  geom_density(alpha=.2) +
  theme(legend.position="none")
p2 <- ggplot(df1, aes(x=FS, colour=PASS, fill=PASS)) +
  geom_density(alpha=.2) +
  theme(legend.position="none") +
  scale_x_continuous(trans='log10')
p3 <- ggplot(df1, aes(x=SOR, colour=PASS, fill=PASS)) +
  geom_density(alpha=.2) +
  theme(legend.position="none")
p4 <- ggplot(df1, aes(x=MQ, colour=PASS, fill=PASS)) +
  geom_density(alpha=.2) +
  theme(legend.position="none")
p5 <- ggplot(df1, aes(x=MQRankSum, colour=PASS, fill=PASS)) +
  geom_density(alpha=.2) +
  theme(legend.position="none")
p6 <- ggplot(df1, aes(x=ReadPosRankSum, colour=PASS, fill=PASS)) +
  geom_density(alpha=.2) +
  theme(legend.position="none")
p <- ggarrange(p1, p2, p3, p4,p5, p6, ncol = 2, nrow = 3)
          #labels = c("A", "B", "C", "D", "E", "F"),
ggsave(output, p)