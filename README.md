# HormoneSignalling

[DESEQ2](https://lashlock.github.io/compbio/R_presentation.html)

[DESEQ_Modeldata](https://bioconnector.github.io/workshops/data.html)

```{r}

#Install packages and load libraries

#install.packages("htmltools")
#library(htmltools)
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

library( "DESeq2" )
library(ggplot2)


#Download data

countsName <- "http://bioconnector.org/workshops/data/airway_scaledcounts.csv"
download.file(countsName, destfile = "airway_scaledcounts.csv", method = "auto")

countData <- read.csv('airway_scaledcounts.csv', header = TRUE, sep = ",")
head(countData)


metaDataName <- "http://bioconnector.org/workshops/data/airway_metadata.csv"
download.file(metaDataName, destfile = "airway_metadata.csv", method = "auto")

metaData <- read.csv('airway_metadata.csv', header = TRUE, sep = ",")
metaData

Construct DESEQDataSet Object
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~dex, tidy = TRUE)
							  
							  
## converting counts to integer mode

dds

#Now we’re ready to run DESEQ function
dds <- DESeq(dds)

res <- results(dds)
head(results(dds, tidy=TRUE))

#Summary of differential gene expression
summary(res) 

#Sort summary list by p-value
res <- res[order(res$padj),]
head(res)

#plotCounts
par(mfrow=c(2,3))

plotCounts(dds, gene="ENSG00000152583", intgroup="dex")
plotCounts(dds, gene="ENSG00000179094", intgroup="dex")
plotCounts(dds, gene="ENSG00000116584", intgroup="dex")
plotCounts(dds, gene="ENSG00000189221", intgroup="dex")
plotCounts(dds, gene="ENSG00000120129", intgroup="dex")
plotCounts(dds, gene="ENSG00000148175", intgroup="dex")

#Volcano Plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

PCA

vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="dex") #using the DESEQ2 plotPCA fxn we can

```



1. [Molecular Mechanisms Underlying Abscisic Acid/Gibberellin Balance in the Control of Seed Dormancy and Germination in Cereals](https://www.frontiersin.org/articles/10.3389/fpls.2018.00668/full)
![img](https://www.frontiersin.org/files/Articles/362906/fpls-09-00668-HTML/image_m/fpls-09-00668-g001.jpg)
![img](https://www.frontiersin.org/files/Articles/362906/fpls-09-00668-HTML/image_m/fpls-09-00668-g002.jpg)

2. [Abscisic Acid and Abiotic Stress Tolerance in Crop Plants](https://www.frontiersin.org/articles/10.3389/fpls.2016.00571/full)
![img](https://www.frontiersin.org/files/Articles/190245/fpls-07-00571-HTML-r1/image_m/fpls-07-00571-g001.jpg)

3. [Hormone-regulated defense and stress response networks contribute to heterosis in Arabidopsis F1 hybrids](https://www.pnas.org/doi/full/10.1073/pnas.1519926112)
![img](https://www.pnas.org/cms/10.1073/pnas.1519926112/asset/a45b2234-12d9-4b07-9c19-5ef9f0ede0a1/assets/graphic/pnas.1519926112fig02.jpeg)

4. [Abscisic acid (ABA) sensitivity regulates desiccation tolerance in germinated Arabidopsis seeds](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.12785)
![img](https://nph.onlinelibrary.wiley.com/cms/asset/d76c8c9a-22b4-4bb4-9694-5b5ac86a93b9/nph12785-fig-0003-m.jpg)

5. [Transcriptional and physiological data revealed cold tolerance in a photo-thermo sensitive genic male sterile line Yu17S](https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-022-03437-8)
![img](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs12870-022-03437-8/MediaObjects/12870_2022_3437_Fig3_HTML.png?as=webp)

6. [Identification of major candidate genes for multiple abiotic stress tolerance at seedling stage by network analysis and their validation by expression profiling in rice (Oryza sativa L.)](https://link.springer.com/article/10.1007/s13205-022-03182-7)
![img](https://media.springernature.com/full/springer-static/image/art%3A10.1007%2Fs13205-022-03182-7/MediaObjects/13205_2022_3182_Fig1_HTML.png?as=webp)
![img](https://media.springernature.com/full/springer-static/image/art%3A10.1007%2Fs13205-022-03182-7/MediaObjects/13205_2022_3182_Fig2_HTML.png?as=webp)

7.[Transcriptomic and metabolomic landscape of quinoa during seed germination](https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-022-03621-w#Sec10)
![img](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs12870-022-03621-w/MediaObjects/12870_2022_3621_Fig2_HTML.png?as=webp)
