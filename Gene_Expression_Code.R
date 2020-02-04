library(Biobase)

dat<-dat[rowSums(is.na(dat))<ncol(dat),]

#impute missing values with KNN imputation method 
library(DMwR)
dat<-knnImputation(dat,k=10,scale=T,meth="weighAvg",disData=NULL)


#need to obtain sample classes "dasatinib-sensitive pancreatic cells", "dasatinib-resistant pancreatice cancer cells"

dat3.samples<-dat3			#subset data with sample cell type 
colnames(dat3.samples)<-cell.type

sensitive<-dat.samples[,1:8]		#subset data to give two data frames; one for each 
resistant<-dat.samples[,9:16]			cell type

sensitive.m<-apply(sensitive,1,mean)
resistant.m<-apply(resistant,1,mean)
fold<-sensitive.m-resistant.m

# Volcano Plot

p.trans<- -1*log10(original.pvalues)

plot(range(p.trans),range(fold),type='n',xlab='-1*log10(p-value)',ylab='fold change',main='Volcano Plot\n Dasatinib-Sensitive vs Dasatinib-Resistant')
points(p.trans,fold,col='black',pch=21,bg=1)
points(p.trans[(p.trans> -log10(.05)&fold>log2(2))],fold[(p.trans> -log10(.05)&fold>log2(2))],col=1,bg=2,pch=21)
points(p.trans[(p.trans> -log10(.05)&fold< -log2(2))],fold[(p.trans> -log10(.05)&fold< -log2(2))],col=1,bg=3,pch=21)
abline(v= -log10(.05))
abline(h= -log2(2))
abline(h=log2(2))





test<-dat.loading[c(5:8,13:16),]			#Using 4 sensitive and 4 resistant cell types for testing set



#print a table for prediction results
pred.table<-table(cell.class,pred)
plot(dat.nn2$fitted.values,type="n",ylab="Fitted Values",xlab="",axes=F,main="Visualization of Neural Network\n Predicting Pancreatic Cancer Cell Type")


points(dat.nn2$fitted.values,pch=23,cex=1.5,bg="grey")
5.522324e-06 1.611060e-05 1.678626e-05 1.816544e-05 2.140999e-05 

#find the p-value and fold change for these five genes
probes<-c("ILMN_1783337","ILMN_1682343","ILMN_3304086","ILMN_1740772","ILMN_2403006")
fold[probes]
dat.genes2[probes]



ILMN_1682343    NONO	 0.3381725	1.611060e-05
ILMN_3304086 ANAPC10	 0.4805388	1.678626e-05
ILMN_1740772   APBB3	-1.1376750	1.816544e-05
ILMN_2403006    TJP1	 1.5622788	2.140999e-05

ensembl ids: 
ENSG00000242612
