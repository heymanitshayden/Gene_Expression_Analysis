library(Biobase) library(GEOquery) #to open GEO .SOFT file  gds5627<-getGEO('GDS5627',destdir="/Users/haydenthomas/downloads") dat<-Table(gds5627) # gather data.frame from file # the first column of data.frame is probe IDs, set row names as probe IDs row.names(dat)<-dat[,1]  colnames(dat) [1] "ID_REF"     "IDENTIFIER" "GSM1435684" "GSM1435685" [5] "GSM1435686" "GSM1435687" "GSM1435688" "GSM1435689" [9] "GSM1435690" "GSM1435691" "GSM1435692" "GSM1435693"[13] "GSM1435694" "GSM1435695" "GSMa1435696" "GSM1435697"[17] "GSM1435698" "GSM1435699" "GSM1435700" "GSM1435701"# need to remove first two columns  dat<-dat[,3:20] colnames(dat) [1] "GSM1435684" "GSM1435685" "GSM1435686" "GSM1435687" [5] "GSM1435688" "GSM1435689" "GSM1435690" "GSM1435691" [9] "GSM1435692" "GSM1435693" "GSM1435694" "GSM1435695"[13] "GSM1435696" "GSM1435697" "GSM1435698" "GSM1435699"[17] "GSM1435700" "GSM1435701"#remove rows where all columns values are NA

dat<-dat[rowSums(is.na(dat))<ncol(dat),]

#impute missing values with KNN imputation method 
library(DMwR)
dat<-knnImputation(dat,k=10,scale=T,meth="weighAvg",disData=NULL)
# look for outliers - Pearson's correlation matrixlibrary(gplots)dat.cor<-cor(dat)layout(matrix(c(1,1,1,1,1,1,1,1,2,2),5,2,byrow=TRUE))par(oma=c(5,7,1,1))cx<-rev(colorpanel(25,"yellow","blue","black"))leg<-seq(min(dat.cor,na.rm=T),max(dat.cor,na.rm=T),length=10)image(dat.cor,main="Correlation Plot of Dasatinib-Resistant/Sensitive \n Pancreatic Cells",axes=F,col=cx)axis(1,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]],cex.axis=0.9,las=2)axis(2,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]],cex.axis=0.9,las=2)image(as.matrix(leg),col=cx,axes=F)tmp<-round(leg,2)axis(1,at=seq(0,1,length=length(leg)),labels=tmp,cex.axis=2)plot(dat.cor)axis(1,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]],cex.axis=0.9,las=2)axis(2,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]],cex.axis=0.9,las=2)#cv vs. mean plot for outlier analysis dat.mean<-apply(dat,2,mean)dat.sd<-sqrt(apply(dat,2,var))dat.cv<-dat.sd/dat.meanplot(dat.mean,dat.cv,main="Dasatinib Resistant/Dasatinib Sensitive Cells \n CV vs. Mean Plot",xlab="Mean",ylab="CV",col="blue",cex=1.5,type="n")points(dat.mean,dat.cv,bg="green",col=1,pch=21)text(dat.mean,dat.cv,label=dimnames(dat)[[2]],pos=1,cex=0.5)#Avg. correaltion plot for outlier analysis dat.avg<-apply(dat.cor,1,mean)par(oma=c(3,0.1,0.1,0.1))plot(c(1,length(dat.avg)),range(dat.avg),type="n",xlab="",ylab="Average R",main="Average Correlation of Dasatinib Risistant/Sensitive \n Pancreatic Cancer Cells",axes=F)points(dat.avg,bg="blue",col=1,pch=21,cex=1.25)axis(1,at=c(1:length(dat.avg)),labels=dimnames(dat)[[2]],las=2,cex.lab=0.4,cex.axis=0.6)axis(2)abline(v=seq(0.5,62.5,1),col="black")#Cluster Dendrogram of Samples for oulier analysis trans.dat<-t(dat)dat.dist<-dist(trans.dat,method="euclidean")dat.clust<-hclust(dat.dist,method="single")plot(dat.clust,labels=names(dat2),cex=0.75,main="Cluster Dendrogram of Pancreatic Cancer Cell Samples",xlab="Sample Name")#Remove the two outlier samples "GSM1435691", "GSM1435694"dat2<-subset(dat,select=-c(GSM1435691,GSM1435694))#look for avg. signal of probesthreshold<-mean(dat[[1]])-sd(dat[[1]])threshold[1] 6.066021# subset data with values above thresholddat3<-subset(dat2[,],rowMeans(dat2)>threshold)dim(dat3)[1] 46889    16

#need to obtain sample classes "dasatinib-sensitive pancreatic cells", "dasatinib-resistant pancreatice cancer cells" samples<-Columns(gds5627)samples<-samples[-c(8,11),]			#remove samples determined as outlierscell.type<-as.matrix(samples[,3])# find the fold difference per gene for sensitive vs resistant cell lines

dat3.samples<-dat3			#subset data with sample cell type 
colnames(dat3.samples)<-cell.type

sensitive<-dat.samples[,1:8]		#subset data to give two data frames; one for each 
resistant<-dat.samples[,9:16]			cell type

sensitive.m<-apply(sensitive,1,mean)
resistant.m<-apply(resistant,1,mean)
fold<-sensitive.m-resistant.m #t-test for feature selection functiont.test.pvalues<-function(x,s1,s2){	x1<-x[s1]	x2<-x[s2]	x1<-as.numeric(x1)	x2<-as.numeric(x2)	t.out<-t.test(x1,x2,alternative="two.sided",var.equal=T)	out<-as.numeric(t.out$p.value)	return(out)}original.pvalues<-apply(dat3,1,t.test.pvalues,s1=cell.type=="dasatinib-sensitive pancreatic cancer cells",s2=cell.type=="dasatinib-resistant pancreatic cancer cells")ttest.adjust<-p.adjust(original.pvalues,method="bonferroni",n=length(original.pvalues))new.threshold<-0.001dat.genes<-original.pvalues[original.pvalues<=new.threshold]length(dat.genes)[1] 17980[1] 10148 # when p<=0.01[1] 3766 # when p<=0.001

# Volcano Plot

p.trans<- -1*log10(original.pvalues)

plot(range(p.trans),range(fold),type='n',xlab='-1*log10(p-value)',ylab='fold change',main='Volcano Plot\n Dasatinib-Sensitive vs Dasatinib-Resistant')
points(p.trans,fold,col='black',pch=21,bg=1)
points(p.trans[(p.trans> -log10(.05)&fold>log2(2))],fold[(p.trans> -log10(.05)&fold>log2(2))],col=1,bg=2,pch=21)
points(p.trans[(p.trans> -log10(.05)&fold< -log2(2))],fold[(p.trans> -log10(.05)&fold< -log2(2))],col=1,bg=3,pch=21)
abline(v= -log10(.05))
abline(h= -log2(2))
abline(h=log2(2))




# A more conservative approach I used for feature selectiont.test.statistic<-function(x,s1,s2){	x1<-x[s1]	x2<-x[s2]	x1<-as.numeric(x1)	x2<-as.numeric(x2)	t2.out<-t.test(x1,x2,alternative="two.sided",var.equal=T)	out2<-as.numeric(t2.out$statistic)	return(out2)}t.statistic.list2<-NULLoriginal.t.test<-apply(dat3,1,t.test.statistic,s1=cell.type=="dasatinib-sensitive pancreatic cancer cells",s2=cell.type=="dasatinib-resistant pancreatic cancer cells")set.seed(1)for(i in 1:100){	print(i)	ran<-sample(c(1:ncol(dat3)),replace=F)	dat.random<-dat3[,ran]	iterated.test<-apply(dat.random,1,t.test.statistic,s1=cell.type=="dasatinib-sensitive pancreatic cancer cells",s2=cell.type=="dasatinib-resistant pancreatic cancer cells")	max.t<-abs(max(iterated.test))	print(max.t)	t.statistic.list<-c(t.statistic.list,max.t)}new.threshold<-as.numeric(quantile(t.statistic.list,0.95))dat.genes<-original.t.test[abs(original.t.test)>threshold]length(dat.genes)[1] 510dat.genes.mat<-as.matrix(dat.genes)dat4<-dat3[c(rownames(dat.genes.mat)),]#now find p-values genes.pvalues<-apply(dat4,1,t.test.pvalues,s1=cell.type=="dasatinib-sensitive pancreatic cancer cells",s2=cell.type=="dasatinib-resistant pancreatic cancer cells")pvalues.adjust<-p.adjust(genes.pvalues,method="bonferroni",n=length(genes.pvalues))threshold2<-0.01dat.genes2<-pvalues.adjust[pvalues.adjust<=threshold2]> length(dat.genes2)[1] 399hist(dat.genes2,main="Histogram of P-values",xlab="p-value",border="black",col="light green",ylim=NULL,axes=TRUE)#subset data.frame with retained genes only dat.genes2.mat<-as.matrix(dat.genes2)use.dat<-dat4[c(rownames(dat.genes2.mat)),]#PCA plot pcadata<-use.datcolnames(pcadata)<-cell.typedat.pca<-prcomp(t(pcadata))dat.loading<-dat.pca$x[,1:3]plot(dat.loading[,1],dat.loading[,2],xlab="p1",ylab="p2",main="PCA Plot of Pancreatic Cancer Cells")points(dat.loading[,1][cell.type=="dasatinib-resistant pancreatic cancer cells"],dat.loading[,2][cell.type=="dasatinib-resistant pancreatic cancer cells"],col="green",pch=16,cex=1.5)points(dat.loading[,1][cell.type=="dasatinib-sensitive pancreatic cancer cells"],dat.loading[,2][cell.type=="dasatinib-sensitive pancreatic cancer cells"],col="blue",pch=16,cex=1.5)leg.names<-c("Resistant","Sensitive")leg.col=c("green","blue")legend(-2,-2,leg.names,leg.col,horiz=F)#classification by neural netlibrary(nnet)samp<-sample(c(1:ncol(dat.loading)),8,replace=F)	#Take 8 random samples to train ANN
test<-dat.loading[c(5:8,13:16),]			#Using 4 sensitive and 4 resistant cell types for testing setdat.n<-data.frame(dat.loading,cell.type.factor)		#The data frame containing the top 3 PCAs for classification
dat.nn2<-nnet(cell.type.factor~.,data=dat.n,subset=sample,size=4,rang=0.1,decay=5e-4,maxit=200) #Build of the NN using train
pred<-predict(dat.nn2,dat.n,type="class")		#Prediction of cell type (class) for the test setcell.class<-as.character(cell.type.factor[c(1:4,9:12)])	#need to gather 8 class identifiers to match dimensions of test set

#print a table for prediction resultslibrary(gridExtra)
pred.table<-table(cell.class,pred)    #ANN Visualization of Fitted Values: 
plot(dat.nn2$fitted.values,type="n",ylab="Fitted Values",xlab="",axes=F,main="Visualization of Neural Network\n Predicting Pancreatic Cancer Cell Type")axis(2)axis(1,at=c(1:8),labels=c("sensitive","resistant","sensitive","sensitive","sensitive","sensitive","resistant","resistant"),cex.axis=0.75,las=2)


points(dat.nn2$fitted.values,pch=23,cex=1.5,bg="grey")#find the top 5 differentially expressed genes: dat.genes2.sorted<-sort(dat.genes2,decreasing=FALSE)genes.ordered[1:5][1]  46 356 277 320 204 # These are the row numbers for gene IDs #print the top five genes: print(dat.genes2.sorted[1:5])ILMN_1783337 ILMN_1682343 ILMN_3304086 ILMN_1740772 ILMN_2403006 
5.522324e-06 1.611060e-05 1.678626e-05 1.816544e-05 2.140999e-05 #get gene IDs library(illuminaHumanv4.db)probeID<-c("ILMN_1682343","ILMN_3304086","ILMN_2403006","ILMN_3242673","ILMN_1842394")data.frame(Gene=unlist(mget(x=probeID,envir=illuminaHumanv4SYMBOL)))

#find the p-value and fold change for these five genes
probes<-c("ILMN_1783337","ILMN_1682343","ILMN_3304086","ILMN_1740772","ILMN_2403006")
fold[probes]
dat.genes2[probes]


                Gene	Fold Change	P-Value(adj.)ILMN_1783337   DECR2	-0.9098388 	5.522324e-06
ILMN_1682343    NONO	 0.3381725	1.611060e-05
ILMN_3304086 ANAPC10	 0.4805388	1.678626e-05
ILMN_1740772   APBB3	-1.1376750	1.816544e-05
ILMN_2403006    TJP1	 1.5622788	2.140999e-05

ensembl ids: 
ENSG00000242612
