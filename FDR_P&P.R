#FDR
library(data.table)
library(dplyr)
library(CMplot)
gwasresult <- fread("GAPIT.Association.GWAS_Results.MLM.CIR10d.csv",header = T)
p <- gwasresult$P.value
p.plot <-  gwasresult[,c(1,2,3,4)]
FDR_P <- p.adjust(p,method = "fdr",n=length(p))
new_result <- cbind(gwasresult[,c(1:6)],FDR_P)
fwrite(new_result,file="new.fdr.MLMresult.csv")

#lambda
install.packages("qqman")
library(qqman)
p_value=gwasresult$P.value
z = qnorm(p_value/ 2)
lambda = round(median(z^2, na.rm = TRUE) / 0.454, 3)
lambda


#Manhanttan plot & QQplot
keySNPs <- read.table("highlightsnps.txt",header = F)
highlightSNPs <- keySNPs$V2
QTL <- keySNPs$V1

CMplot(p.plot,plot.type="m",
       col=c("grey30","grey60"),
       highlight=highlightSNPs,
       highlight.text = QTL,
       highlight.col="red",
       highlight.cex=1.5,
       highlight.pch=19,
       chr.labels=paste("Chr",c(1:12),sep=""),
       r=0.4,
       #chr.den.col=c("darkgreen","yellow","red"),
       file.name = "p.manhattanplot",
       file="pdf",
       dpi=600,
       file.output=TRUE,
       verbose=TRUE)

CMplot(p.plot,plot.type="q",
       file.name = "p.qqplot",
       file="pdf",
       dpi=600,file.output=TRUE,verbose=TRUE)

CMplot(fdr_p.plot,plot.type="m",
       col=c("grey30","grey60"),
       highlight=highlightSNPs,
       highlight.text = QTL,
       highlight.col="red",
       highlight.cex=1.5,
       highlight.pch=19,
       chr.labels=paste("Chr",c(1:12),sep=""),
       r=0.4,
       #chr.den.col=c("darkgreen","yellow","red"),
       file.name = "FDR.manhattan.plot",
       file="pdf",
       dpi=600,file.output=TRUE,verbose=TRUE)

CMplot(fdr_p.plot,plot.type="q",
       file.name = "FDR.qqplot",
       file="pdf",
       dpi=600,file.output=TRUE,verbose=TRUE)