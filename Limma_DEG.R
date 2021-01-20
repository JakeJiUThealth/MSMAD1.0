library(GEOquery)
library(limma)
gse <- getGEO('GSE65067')  #GSE104775; GSE75431
myMatrix <- gse[[1]]@assayData$exprs  #extract gene expression matrix
Ftrd<-gse[[1]]@featureData
ID_labels<-Ftrd@data    #annotations for row-name
myeset<-log2(myMatrix)  #log2 normalization
rownum<-nrow(myeset)    #gene numbers
colnum<-10                #5 control : 5 treatment
design <- cbind(WT=1,TrvsWT=c(0,0,0,0,0,1,1,1,1,1))
colnames(design) <- c("WT","ADvsWT")
myeset1<-myeset[1:rownum,c(1:5,10:14)]   # 5 for WT and 5 for AD
fit <- lmFit(myeset1, design)
fit <- eBayes(fit)
cal_DEGs<-topTable(fit, number = nrow(fit), coef="ADvsWT", adjust="BH", confint = TRUE)  #DEGs 
write.csv(cal_DEGs,'DEG_GSE65067.csv')

