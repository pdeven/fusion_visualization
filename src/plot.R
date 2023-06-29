#!/usr/bin/env Rscript
library(circlize)
library(optparse)

option_list = list(
  make_option(c("-l", "--left_bed"),help="regions showing breakpoint of left gene"),
  make_option(c("-r", "--right_bed"),help="regions showing breakpoint of right gene"),
  make_option(c("-o", "--pdf"),help="output pdf path"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
left <- opt$left_bed
right <- opt$right_bed
pdf_output <- opt$pdf
df1=read.csv(left, sep = "\t")
df2=read.csv(right, sep = "\t")
df=rbind(df1, df2)

pdfWidth=list("pdfWidth", "numeric", 200)
pdfHeight=list("pdfHeight", "numeric", 100)

pdf(pdf_output , width = pdfWidth, height = pdfHeight, title = "Fusions")

circos.initializeWithIdeogram()
#circos.initializeWithIdeogram(plotType = c("labels", "axis"))

circos.genomicLabels(df, 
                     labels.column = 4, 
                     side = "outside",
                     col = as.numeric(factor(df[[1]])), 
                     cex = 0.1,
                     line_col = as.numeric(factor(df[[1]])))
circos.genomicLink(df1, df2, col = rand_color(nrow(df1), transparency = 0), 
                   border = NA, lwd=1)
mtext("    GRCh37\n", outer=TRUE,  cex=0.8, line=-2,side=3,adj=0,font = 4)
mtext("     Fusion Genes", outer=TRUE,  cex=0.7, line=-2.5,side=3,adj=0,font = 3)
dev.off()

