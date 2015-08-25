require(gdata) #paczka z read.xls
require(ggplot2)# wykresy
require(clv)
require(reshape)
require(chemometrics)
require(gplots)
require(haplo.stats)
require(MASS)
require(glmnet)
require(plyr)
require(Hmisc)
source("DataAnalysis_BK.R")



raw<-read.xls("polsenior_komplet_ver6_BK1.xls")
#Sprawdzenie czy identyfikatory nie powtarzaja sie.
max(table(raw[,"ID"]))
max(table(raw[,"ID1"]))
max(table(raw[,"ID2"]))
multipleEntries<-which(table(raw[,"KOD"])!=1)#Pacjenci wielokrotnie w bazie
numEntries<-table(raw[,"KOD"]) #Liczba wystapien pacjenta w bazie
numEntries[multipleEntries] #Ile razy wystepuja wielokrotnie
multipleEntrSUM<-sum(numEntries[multipleEntries])#suma wielokrotnych wpisow

#Przykladowa zawartosc wielokrotenie badanych pacjentow
raw[which(raw[,"KOD"]=="02/02/242"),]
raw[which(raw[,"KOD"]=="02/09/138"),]
raw[which(raw[,"KOD"]=="30/26/252"),]
as.Date(raw[which(raw[,"KOD"]=="30/26/252"),"DATA"],origin="1899-12-30") # data liczona w dniach od 1899-12-30
raw[which(raw[,"KOD"]=="30/26/252"),c(1:13,23:34,36,42:43,52:55)] # sprawdzenie wybranych pol

#Wybranie wylacznie pojedynczych wpisow
singleEntries<-which(is.na(pmatch(raw[,"KOD"],names(multipleEntries),duplicates.ok=T)))
singleRaw<-raw[singleEntries,]
dim(singleRaw) #Sprawdzenie rozmiarow proby
names(singleRaw)<-toupper(names(singleRaw)) # nazwy kolumn tylko z duzych liter

singleRaw[which(singleRaw[,"WIELKOSC.MIEJSCOWOSCI"]=="Szprotawa"),"WIELKOSC.MIEJSCOWOSCI"]="miasto do 20tys"
singleRaw[which(singleRaw[,"GR.WIEKOWA"]=="90 i wiecej"),"GR.WIEKOWA"]="90 lat i wiecej"
df_NoNA<-removeNA(singleRaw)

single_NoNA<-df_NoNA$data

print(typeof(single_NoNA))
print(table(single_NoNA[,"PLEC"])) # liczba kobiet i mezczyzn
mInd<-which(single_NoNA[,"PLEC"]=="mezczyzna") #indeksy mezczyzn
kInd<-which(single_NoNA[,"PLEC"]=="kobieta") #indeksy kobiet

#singleRaw_df<-data.frame(singleRaw)
#Okreslenie typu danych w poszczegolnych polach
dataTypes<-sapply(single_NoNA[1,],class)
numericData_Ind<-which(dataTypes=="numeric" | dataTypes=="integer") # wybor danych liczbowych
numericData_Ind<-numericData_Ind[-c(1,2,3)] # usuniecie pol z indeksami
numericData_Ind<-numericData_Ind[-32] # usuniecie pola z data
numericData_Ind<-numericData_Ind[-c(28,29)] #usuniecie pol z polimorfizmami Bsm
#usuniecie transformacji
attr_trans<-which(!is.na(pmatch(names(numericData_Ind),c("LN.TEST","LOG.TEST","TEST_NMOL_L","E2_PMOL_L","LN.WITD","LOG.WITD"))))
numericData_Ind<-numericData_Ind[-attr_trans]
numericData<-single_NoNA[,numericData_Ind]

#Histogramy wszystkich atrybutow z podzialem na plcie
numericDataFlat<-melt(numericData[mInd,])
gg_hist_m<-ggplot(numericDataFlat,aes(x = value)) + facet_wrap(~variable,scales = "free_x") + geom_histogram()
print(gg_hist_m)
dev.print(device=png, file="HistogramyAtrybutow_Mezczyzni.png", width=1240, height=720)
dev.new()
numericDataFlat_k<-melt(numericData[kInd,])
gg_hist_k<-ggplot(numericDataFlat_k,aes(x = value)) + facet_wrap(~variable,scales = "free_x") + geom_histogram()
print(gg_hist_k)
dev.print(device=png, file="HistogramyAtrybutow_Kobiet.png", width=1240, height=720)

numericData_m<-numericData[mInd,]
numericData_k<-numericData[kInd,]

#Usuniecie jednowymiarowych wartosci odstajacych
numericNoOut1D_m<-lapply(numericData_m,outlierDet_Hampel)
numericNoOut1D_k<-lapply(numericData_k,outlierDet_Hampel)
allOutID_m<-lapply(numericNoOut1D_m,f<-function(x){x[[2]]})

#Histogramy wszystkich atrybutow i ich logarytmow
#plotHists(numericData_m)

#dodanie logarytmow naturalnych wybranych atrybutow
numericDataWithLogs_m<-addLog(numericData_m,c("FSH","DHEA.S","FEI","TGC","INS","WITAMINA.D","ICTP"))

#Wyznaczenie warto?ci odstaj?cych
AllOutliers<-lapply(numericDataWithLogs_m,outlierDet_Hampel)

#Progi dla wartosci odstawjacych
OutlierRanges_1D_m<-lapply(AllOutliers,function(x)x[[3]])

#Wyciagniecie indeksow
AllOutliersID<-lapply(AllOutliers,function(x)x[[2]])
t(names(AllOutliersID))
NotUsedForOutliers<-c(1,3,6,8,18,19,21,23) # kt?re atrybuty zosta?y zast?pione przez logarytmy
outliersStackIDs<-stack(AllOutliersID)$ind %in% names(AllOutliersID)[-NotUsedForOutliers]
outliersID_1D_m<-unique(stack(AllOutliersID)$value[outliersStackIDs])

numericDataWithLogsNorm_m<-data.frame(sapply(numericDataWithLogs_m,normalizeZRobust))
source("VitD_PCA.R")
numericDataWLogsNormPCA_m<-PCA_bplot_sdev(numericDataWithLogsNorm_m)
numericDataWLogsNormPCAdf_m<-data.frame(numericDataWLogsNormPCA_m$x)
plotHists(numericDataWLogsNormPCAdf_m)

AllOutliers_PCA_m<-lapply(numericDataWLogsNormPCAdf_m,outlierDet_Hampel)
AllOutliersRanges_PCA_m<-lapply(AllOutliers_PCA_m, function(x)x[[3]])
AllOutliersID_PCA_m<-lapply(AllOutliers_PCA_m, function(x)x[[2]])
outliersStackIDs_PCA<-stack(AllOutliersID_PCA_m)$ind %in% names(AllOutliersID_PCA_m)
outliersID_PCA_m<-unique(stack(AllOutliersID_PCA_m)$value[outliersStackIDs_PCA])
sum(outliersID_PCA_m %in% outliersID_1D_m)# 53 probki, wyrzucone zarowno przez PCA jak i 1D

#53 pokrywajace sie , PCA wyznaczylo 104 wartosci odstajace, 1D 81 wartosci odstajacych
#> length(outliersID_1D_m)
#[1] 81
#> length(outliersID_PCA_m)
#[1] 104


#Probka 145 musi byc odrzucona przed sama analiza PCA - zbyt mocno wplywa na dane.

numericDataPCA_no145_m<-PCA_bplot_sdev(numericDataWithLogsNorm_m[-145,])
numericDataPCA_no145df_m<-data.frame(numericDataPCA_no145_m$x)
plotHists(numericDataPCA_no145df_m)

AllOutliers_PCA_no145_m<-lapply(numericDataPCA_no145df_m[,-c(28,29,30)],outlierDet_Hampel) # nie sa brane pod uwage 3 ostatnie skladowe
AllOutliersRanges_PCA_no145_m<-lapply(AllOutliers_PCA_no145_m, function(x)x[[3]])
AllOutliersID_PCA_no145_m<-lapply(AllOutliers_PCA_no145_m, function(x)x[[2]])
outliersStackIDs_PCA_no145<-stack(AllOutliersID_PCA_no145_m)$ind %in% names(AllOutliersID_PCA_no145_m)[-c(28,29,30)]
outliersID_PCA_no145_m<-unique(stack(AllOutliersID_PCA_no145_m)$value[outliersStackIDs_PCA_no145])

outliersPCA_no145_IDcorrected<-strtoi(outlierNames_PCA_no145)
outliersPCA_no145_IDcorrected %in% outliersID_1D_m


#Wykorzystanie Mahalanobis distance!!!
md_NumData<-Moutlier(numericData_m,quantile=0.975,plot=FALSE) # standardowy kwantyl 97.5%
plot(md_NumData_m$md,md_NumData_m$rd)
lines(c(md_NumData_m$cutoff,md_NumData_m$cutoff),c(0,150), col="red")
lines(c(0,17),c(md_NumData_m$cutoff,md_NumData_m$cutoff), col="red")
dev.new()
par(mfcol=c(2,1))
hist(md_NumData$md)
hist(md_NumData$rd)
outliers_NumData_m_md0975<-which(md_NumData$md>md_NumData$cutoff)
outliers_NumData_m_rd0975<-which(md_NumData$rd>md_NumData$cutoff)

#Wykorzystanie Mahalanobis distance dla danych z logarytmami!!
md_NumDataWLogs_m<-Moutlier(numericDataWithLogs_m,quantile=0.975,plot=TRUE) # standardowy kwantyl 97.5%
plot(md_NumDataWLogs_m$md,md_NumDataWLogs_m$rd)
lines(c(md_NumDataWLogs_m$cutoff,md_NumDataWLogs_m$cutoff),c(0,150), col="red")
lines(c(0,17),c(md_NumDataWLogs_m$cutoff,md_NumDataWLogs_m$cutoff), col="red")
dev.new()
par(mfcol=c(2,1))
hist(md_NumDataWLogs_m$md)
hist(md_NumDataWLogs_m$rd)
outliers_NumDataWLogs_m_md0975<-which(md_NumDataWLogs_m$md>md_NumDataWLogs_m$cutoff)
outliers_NumDataWLogs_m_rd0975<-which(md_NumDataWLogs_m$rd>md_NumDataWLogs_m$cutoff)

#Klastrowanie hierarchiczne metoda Warda
numericDataNorm_m<-data.frame(lapply(numericData_m,normalizeZRobust))
numDataNorm_m_HClust<-hclust(dist(numericDataNorm_m),method="ward.D")
plot(numDataNorm_m_HClust)

#Inne klastrowania
numDataNorm_m_HClust_UPGMA<-hclust(dist(numericDataNorm_m),method="average")
plot(numDataNorm_m_HClust_UPGMA)
numDataNorm_m_HClust_UPGMA<-hclust(dist(numericDataNorm_m[-outliers_NumData_m_md0975,]),method="average")
plot(numDataNorm_m_HClust_UPGMA)
numDataNorm_m_HClust_WPGMC<-hclust(dist(numericDataNorm_m),method="median")
plot(numDataNorm_m_HClust_WPGMC)
numDataNorm_m_HClust_WPGMC<-hclust(dist(numericDataNorm_m[-outliers_NumData_m_md0975,]),method="median")
plot(numDataNorm_m_HClust_WPGMC)
numDataNorm_m_HClust_UPGMC<-hclust(dist(numericDataNorm_m[-outliers_NumData_m_md0975,]),method="centroid")
plot(numDataNorm_m_HClust_WPGMC)
plot(numDataNorm_m_HClust_UPGMC)
numDataNorm_m_HClust_cmpl<-hclust(dist(numericDataNorm_m[-outliers_NumData_m_md0975,]),method="complete")
plot(numDataNorm_m_HClust_cmpl)
numDataNorm_m_HClust_cmpl<-hclust(dist(numericDataNorm_m),method="complete")
plot(numDataNorm_m_HClust_cmpl)
numDataNorm_m_HClust_cmpl_noOut<-hclust(dist(numericDataNorm_m[-outliers_NumData_m_md0975,]),method="complete")
plot(numDataNorm_m_HClust_cmpl_noOut)
numDataNorm_m_HClust_single<-hclust(dist(numericDataNorm_m),method="single")
plot(numDataNorm_m_HClust_single)
numDataNorm_m_HClust_single_noOut<-hclust(dist(numericDataNorm_m[-outliers_NumData_m_md0975]),method="single")
plot(numDataNorm_m_HClust_single_noOut)
numDataNorm_m_HClust_single_noOut<-hclust(dist(numericDataNorm_m[-outliers_NumData_m_md0975,]),method="single")
plot(numDataNorm_m_HClust_single_noOut)
numDataNorm_m_HClust_wardD2<-hclust(dist(numericDataNorm_m),method="ward.D2")
plot(numDataNorm_m_HClust_wardD2)
numDataNorm_m_HClust_wardD2_noOut<-hclust(dist(numericDataNorm_m[-outliers_NumData_m_md0975,]),method="ward.D2")
plot(numDataNorm_m_HClust_wardD2_noOut)
numDataNorm_m_noOut_HClust_wardD<-hclust(dist(numericDataNorm_m[-outliers_NumData_m_md0975,]),method="ward.D2")
plot(numDataNorm_m_noOut_HClust_wardD)

#Heat map - mapa cieplna po klastrowaniu
rowv<-as.dendrogram(hclust(dist(numericDataNorm_m[-outliers_NumData_m_md0975,]),method="ward.D2"))
colv<-as.dendrogram(hclust(dist(t(numericDataNorm_m[-outliers_NumData_m_md0975,])),method="ward.D2"))
dev.new()
heatmap.2(as.matrix(numericDataNorm_m[-outliers_NumData_m_md0975,]),Rowv=rowv,Colv=colv,trace="none",symm=F,symkey=F,symbreaks=T, scale="none")
dev.new()
heatmap(as.matrix(numericDataNorm_m[-outliers_NumData_m_md0975,]),Rowv=rowv,Colv=colv,trace="none")

#Ocena klastrowania - zbior danych razem z wartosciami odstajacymi
numDataNorm_m_HCL_ward_cut<-cutree(numDataNorm_m_HClust,k=2:10)
numDataNorm_m_HCL_ward_clv<-BK_clv(data=numericDataNorm_m,clust=numDataNorm_m_HCL_ward_cut,dis='euclidean')
dev.new()
par(mfcol=c(1,2))
plot(2:10,numDataNorm_m_HCL_ward_clv$Davies.Bouldin,type='b')
plot(2:10,numDataNorm_m_HCL_ward_clv$Dunn,type='b')

#Ocena klastrowania - zbior danych bez wartosci odstajacych
numDataNorm_m_noOut_HClust_wardD<-hclust(dist(numericDataNorm_m[-outliers_NumData_m_md0975,]),method="ward.D2")
numDataNorm_m_noOut_HCL_ward_clust<-cutree(numDataNorm_m_noOut_HClust_wardD,k=2:10)
numDataNorm_m_noOut_HCL_ward_clv<-BK_clv(data=numericDataNorm_m[-outliers_NumData_m_md0975,],clust=numDataNorm_m_noOut_HCL_ward_clust,dis='euclidean',ValDistIntra='average',ValDistInter='hausdorff')
dev.new()
par(mfcol=c(1,2))
plot(2:10,numDataNorm_m_noOut_HCL_ward_clv$Davies.Bouldin,type='b')
plot(2:10,numDataNorm_m_noOut_HCL_ward_clv$Dunn,type='b')

#Ocena klastrowania - zbior danych bez wartosci odstajacych - klastrowanie metoda UPGMA
numDataNorm_m_noOut_HClust_UPGMA<-hclust(dist(numericDataNorm_m[-outliers_NumData_m_md0975,]),method="average")
numDataNorm_m_noOut_HClust_UPGMA_clust<-cutree(numDataNorm_m_noOut_HClust_UPGMA,k=2:10)
numDataNorm_m_noOut_HClust_UPGMA_clv<-BK_clv(data=numericDataNorm_m[-outliers_NumData_am_md0975,],clust=numDataNorm_m_noOut_HClust_UPGMA_clust,dis='euclidean',ValDistIntra='centroid',ValDistInter='centroid')
dev.new()
par(mfcol=c(1,2))
plot(2:10,numDataNorm_m_noOut_HClust_UPGMA_clv$Davies.Bouldin,type='b')
plot(2:10,numDataNorm_m_noOut_HClust_UPGMA_clv$Dunn,type='b')

#Wizualizacja klastrow
d<-data.frame(md_NumData_m$md[-outliers_NumData_m_md0975], md_NumData_m$rd[-outliers_NumData_m_md0975],numDataNorm_m_noOut_HCL_ward_clust)
names(d)<-c("md","rd",paste("cl",colnames(numDataNorm_m_noOut_HCL_ward_clust),sep=""))
gg<-ggplot(data=d)
gg+geom_point(aes(y=rd,x=md,colour=factor(cl4),shape=factor(cl4),size=5))

#ANALIZA PCA na znormalizowanych danych, bez warto?ci odstaj?cych

numDataNorm_m_noOut_PCA<-PCA_bplot_sdev(numericDataNorm_m[-outliers_NumData_m_md0975,])
biplot(numDataNorm_m_noOut_PCA,xlabs=d$cl4) # obserwacja z etykietami nadanymi z klastrowania.
biplot(numDataNorm_m_noOut_PCA,xlabs=d$cl6)

#Analiza korelacji pomi?dzy atrybutami

numAttributes_corMat<-cor(numData_m_noOut) # korelacja pearona - liniowa
rowv<-as.dendrogram(hclust(dist(numAttributes_corMat),method="ward.D2"))
colv<-as.dendrogram(hclust(dist(t(numAttributes_corMat)),method="ward.D2"))
dev.new()
heatmap.2(as.matrix(numAttributes_corMat,Rowv=rowv,Colv=colv))


numAttributes_corMat_Spearman<-cor(numData_m_noOut, method="spearman") # korelacja rangowa - spearmana
rowv<-as.dendrogram(hclust(dist(numAttributes_corMat_Spearman),method="ward.D2"))
colv<-as.dendrogram(hclust(dist(t(numAttributes_corMat_Spearman)),method="ward.D2"))
dev.new()
heatmap.2(as.matrix(numAttributes_corMat_Spearman,Rowv=rowv,Colv=colv))

biplot(numDataNorm_m_noOut_PCA,xlabs=all_d$NADCISNIENIE.TETNICZE)

#UTWORZENIE CALEJ MACIERZY DANYCH - tylko mezczyzni i bez wartosci odstajacych
all_d<-data.frame(single_NoNA[row.names(numData_m_noOut),])
#Poprawienie pola WIELKOSC.MIEJSCOWOSCI
all_d$WIELKOSC.MIEJSCOWOSCI<-factor(all_d$WIELKOSC.MIEJSCOWOSCI,levels(all_d$WIELKOSC.MIEJSCOWOSCI)[c(6,4,1,3,2,5)])

#Poprawienie pola GR.WIEKOWA
all_d[which(all_d[,"GR.WIEKOWA"]=="90 i wiecej"),"GR.WIEKOWA"]="90 lat i wiecej"
all_d$GR.WIEKOWA<-factor(all_d$GR.WIEKOWA)

#Uzupelnienie brakujacej wartosci w CG1.ROZPOZNANA.WCZESNIEJ.CUKRZYCA
#Pacjent ma poziom glukozy w normie - brak cukrzycy
all_d["177","CG1.ROZPOZNANA.WCZESNIEJ.CUKRZYCA"]="nie"
all_d$CG1.ROZPOZNANA.WCZESNIEJ.CUKRZYCA<-factor(all_d$CG1.ROZPOZNANA.WCZESNIEJ.CUKRZYCA)
#Uzupelnienie NASLONECZNIENIE - trzy sprawdzone probki sa z wrzesnia, wiec naslonecznienie "TAK"
all_d[which(all_d$NASLONECZNIENIE==""),"NASLONECZNIENIE"]="TAK"
all_d$NASLONECZNIENIE<-factor(all_d$NASLONECZNIENIE)
#Uzupelnienie nadcisnienia - pacjent "177" - grupa 4 - odsetek z nadcisnieniem takisam jak 
#w calym zbiorze, 11/37 i 76/254 ok 0.29 - wieksze prawdopodobienstwo na "TAK"
all_d[which(all_d$NADCISNIENIE.TETNICZE==""),"NADCISNIENIE.TETNICZE"]="tak"
all_d["17","NADCISNIENIE.TETNICZE"]="tak"
all_d$NADCISNIENIE.TETNICZE<-factor(all_d$NADCISNIENIE.TETNICZE)


#Poprawienie factorow _ PROBA
try<-lapply(all_d[,nonGenData_Ind],factor)

#SPRAWDZENIE DANYCH GENETYCZNYCH
lapply(all_d[seq(41,71)],table)

#ZAMIANA WARTOSCI NOMINALNYCH - NIEGENETYCZNYCH - NA CIAGI BINARNE
nominalData_Ind<-which(dataTypes=="factor")
#Usuniecie atrybutow:
#1 - KOD,2 - PLEC, 3- GR.WIEKOWA (zamiana na kategorie stopniowane) 5 - WOJEWODZTWO, 6 - WIELKOSC.MIEJSCOWOSCI (zamiana na kategorie stopniowane) 
#10 - ZAKRES.WITAMINY.D 11 - VDR..CDX., 16-27 - MODELE, 28 - MIESIAC,31 - WITD - zle ustalony prog, 32-34 - OBZM, OZZM, OMPMC - brak zmiennosci
nominalData_Ind<-nominalData_Ind[-c(1,2,3,5,6,10,11,16:28,31,32:34)]

#Dane niegenetyczne

nonGenData_Ind<-nominalData_Ind[c(1:4,9:12)] #wyciete dane genetyczne
lapply(all_d[,nonGenData_Ind],table)

nonGenData_Factors<-data.frame(lapply(all_d[,nonGenData_Ind],factor))
row.names(nonGenData_Factors)<-row.names(all_d)
lapply(nonGenData_Factors,table)
nonGenData_m_noOut<-model.matrix(~ .+0, nonGenData_Factors)

#Dolaczenie pol z GR.WIEKOWA i WIELKOSC.MIEJSCOWOSCI
nonGenData_m_noOut<-data.frame(nonGenData_m_noOut,as.numeric(all_d[,"WIELKOSC.MIEJSCOWOSCI"]),as.numeric(all_d[,"GR.WIEKOWA"]))
row.names(nomData_m_noOut)<-row.names(all_d)
names(nonGenData_m_noOut)[dim(nonGenData_m_noOut)[2]-1]<-"WIELKOSC.MIEJSCOWOSCI"
names(nonGenData_m_noOut)[dim(nonGenData_m_noOut)[2]]<-"GR.WIEKOWA"
nonGenData_m_noOut$GR.WIEKOWA<-normalizeZRobust(nonGenData_m_noOut$GR.WIEKOWA)
nonGenData_m_noOut$WIELKOSC.MIEJSCOWOSCI<-normalizeZRobust(nonGenData_m_noOut$WIELKOSC.MIEJSCOWOSCI)

#Skalowanie 0 i 1 do kwantylu 5% i 95% danych numerycznych
scalingRange_5_95<-quantile(unlist(numDataNorm_m_noOut),probs=c(0.05,0.95))#probs=c(0.25,0.75))
numNonGenDataNorm_m_noOut_5_95<-nonGenData_m_noOut*diff(scalingRange_5_95)+scalingRange_5_95[1]

#Skalowanie 0 i 1 do kwantyli 25% i 75%
scalingRange_25_75<-quantile(unlist(numDataNorm_m_noOut),probs=c(0.25,0.75))
numNonGenDataNorm_m_noOut_25_75<-nonGenData_m_noOut*diff(scalingRange_25_75)+scalingRange_25_75[1]


#Dane genetyczne
genData_Ind<-which(regexpr(pattern = "*VDR*",text = names(all_d))==1)
genData<-data.frame(lapply(all_d[,genData_Ind],factor))
genData<-genData[,-1]
row.names(genData)<-row.names(all_d)

#Poprawienie wartosci 
levels(genData[,"VDR..BSM."])[1]<-"BB"
levels(genData[,"VDR..BSM."])[2]<-"Bb"
levels(genData[,"VDR..BSM."])[3]<-"bb"

levels(genData[,"VDR..FOK."])[1]<-"FF"
levels(genData[,"VDR..FOK."])[2]<-"Ff"
levels(genData[,"VDR..FOK."])[3]<-"ff"

genData.Miss<-(genData$VDR..APA.=="" | genData$VDR..TAQ.=="")
genData.c<-genData[!genData.Miss,] #complete genData
genData.c$VDR..TAQ.<-factor(genData.c$VDR..TAQ.)
genData.c$VDR..APA.<-factor(genData.c$VDR..APA.)

levels(genData.c$VDR..TAQ.)<-c(0,1,2)
levels(genData.c$VDR..APA.)<-c(0,1,2)
levels(genData.c$VDR..FOK.)<-c(0,1,2)
levels(genData.c$VDR..BSM.)<-c(0,1,2)

genData.c.bin<-BK_binGenData(genData.c)

genData.c2<-geno1to2(genData.c,locus.label = names(genData))
genData.haplo.em<-haplo.em(genData.c2,locus.label = names(genData))


genData_freq<-lapply(genData,table)


genData_m_noOut<-model.matrix(~ .+0, genData)



#Binaryzacja danych nominalnych
#Jezeli nie ma brakujacych wartosci to mozna uzyc - nomData_m_noOut<-model.matrix(~ .+0, all_d[nominalData_Ind])
#nomData_m_noOut <-data.frame(lapply(all_d[nominalData_Ind],function(x){model.matrix(~x+0, x)}))



#POLCZENIE DANYCH NUMERYCZNYCH i NOMINALNYCH NIEGENETYCZNYCH

numNonGen_Numerized_m_noOut_5_95<-data.frame(numDataNorm_m_noOut,numNonGenDataNorm_m_noOut_5_95)
numNonGen_Numerized_m_noOut_25_75<-data.frame(numDataNorm_m_noOut,numNonGenDataNorm_m_noOut_25_75)


#ANALIZY NA PELNYCH DANYCH


numNonGen_5_95_HClust_wardD2<-hclust(dist(numNonGen_Numerized_m_noOut_5_95),method="ward.D2")
plot(numNonGen_5_95_HClust_wardD2)

numNonGen_25_75_HClust_wardD2<-hclust(dist(numNonGen_Numerized_m_noOut_25_75),method="ward.D2")
plot(numNonGen_25_75_HClust_wardD2)



#Heat map - mapa cieplna po klastrowaniu
rowv<-as.dendrogram(hclust(dist(numNonGen_Numerized_m_noOut_5_95),method="ward.D2"))
colv<-as.dendrogram(hclust(dist(t(numNonGen_Numerized_m_noOut_5_95)),method="ward.D2"))
dev.new()
heatmap.2(as.matrix(numNonGen_Numerized_m_noOut_5_95),Rowv=rowv,Colv=colv,trace="none",symm=F,symkey=F,symbreaks=T, scale="none",margins=c(10,10))
dev.new()
heatmap.2(as.matrix(cor(numNonGen_Numerized_m_noOut_5_95)),margins=c(10,10))

#Heat map - mapa cieplna po klastrowaniu
rowv<-as.dendrogram(hclust(dist(numNonGen_Numerized_m_noOut_25_75),method="ward.D2"))
colv<-as.dendrogram(hclust(dist(t(numNonGen_Numerized_m_noOut_25_75)),method="ward.D2"))
dev.new()
heatmap.2(as.matrix(numNonGen_Numerized_m_noOut_25_75),Rowv=rowv,Colv=colv,trace="none",symm=F,symkey=F,symbreaks=T, scale="none",margins=c(10,10))
dev.new()
heatmap.2(as.matrix(cor(numNonGen_Numerized_m_noOut_25_75)),margins=c(10,10))


#Ocena klastrowania
numNonGen_5_95_HClust_wardD2_cut<-cutree(numNonGen_5_95_HClust_wardD2,k=2:10)
numNonGen_5_95_HClust_wardD2_clv<-BK_clv(data=numNonGen_Numerized_m_noOut_5_95,clust=numNonGen_5_95_HClust_wardD2_cut,dist='euclidean')
dev.new()
par(mfcol=c(1,2))
plot(2:10,numNonGen_5_95_HClust_wardD2_clv$Davies.Bouldin,type='b')
plot(2:10,numNonGen_5_95_HClust_wardD2_clv$Dunn,type='b')

#Ocena klastrowania
numNonGen_25_75_HClust_wardD2_cut<-cutree(numNonGen_25_75_HClust_wardD2,k=2:10)
numNonGen_25_75_HClust_wardD2_clv<-BK_clv(data=numNonGen_Numerized_m_noOut_25_75,clust=numNonGen_25_75_HClust_wardD2_cut,dist='euclidean')
dev.new()
par(mfcol=c(1,2))
plot(2:10,numNonGen_25_75_HClust_wardD2_clv$Davies.Bouldin,type='b')
plot(2:10,numNonGen_25_75_HClust_wardD2_clv$Dunn,type='b')
c75<-data.frame(numNonGen_25_75_HClust_wardD2_cut)
names(c75)<-names(cN)


#Wykorzystanie Mahalanobis distance - nie jest mozliwe - macierz kowariancji jest singular
md_numNonGenData_5_95<-Moutlier(numNonGen_Numerized_m_noOut_5_95,quantile=0.975,plot=FALSE) #- standardowy kwantyl 97.5%
plot(md_numNonGenData_5_95$md,md_numNonGenData_5_95$rd)
lines(c(md_numNonGenData_5_95$cutoff,md_numNonGenData_5_95$cutoff),c(0,150), col="red")
lines(c(0,17),c(md_numNonGenData_5_95$cutoff,md_numNonGenData_5_95$cutoff), col="red")
dev.new()
par(mfcol=c(2,1))
hist(md_numNonGenData_5_95$md)
hist(md_numNonGenData_5_95$rd)
outliers_numNonGenData_5_95_md0975<-which(md_numNonGenData_5_95$md>md_numNonGenData_5_95$cutoff)
outliers_numNonGenData_5_95_rd0975<-which(md_numNonGenData_5_95$rd>md_numNonGenData_5_95$cutoff)

#Analiza z wykorzystaniem PCA
numNonGen_Numerized_m_noOut_PCA<-PCA_bplot_sdev(numNonGen_Numerized_m_noOut)
dev.new()
biplot(numNonGen_Numerized_m_noOut_PCA,xlabs=d$cl4) # obserwacja z etykietami nadanymi z klastrowania.
biplot(numNonGen_Numerized_m_noOut_PCA,xlabs=d$cl6)

#budowa modelu liniowego
T.lm_data<-data.frame(numNonGen_Numerized_m_noOut_25_75[!genData.Miss,c(4,1:3,5:dim(numNonGen_Numerized_m_noOut_25_75)[2])],genData.c.bin)
row.names(T.lm_data)<-row.names(genData.c)
T.lm_formula<-formula(T.lm_data)
trainID<-sample(1:dim(T.lm_data)[1],120)




T.lm<-lm(T.lm_formula,data=T.lm_data,subset=trainID)
T.predict<-predict(T.lm,T.lm_data[-trainID,-1])
plot(T.lm_data$TESTOSTERON[-trainID],T.predict)
lines(x = c(-2,2),y=c(-2,2),col = 'red')

#Bez nadmiarowych kolumn
sing<-c(7,29,46,49,52,55,35,41,30)#wynikaja z tego, ze 0 we wczesniejszych kolumnach juz definiuja co bedzie w tych + hyperandrogenism + FAI
T.rlm_data<-T.lm_data[,-sing]
T.rlm_formula<-formula(T.rlm_data)


#T.rlm<-rlm(T.lm_formula,data=T.lm_data,subset=trainID,method="MM")
shufleID<-sample(1:193,193)

#Wybor danych do cross-validacji k=4
cv.1<-shufleID[1:48]
cv.2<-shufleID[49:(2*48)]
cv.3<-shufleID[(2*48+1):(3*48)]
cv.4<-shufleID[(3*48+1):193]
cv<-list(cv.1=cv.1,cv.2=cv.2,cv.3=cv.3,cv.4=cv.4)

T.lm_1<-lm(T.rlm_formula,data=T.lm_data,subset=unlist(cv[-1]))
T.lm_2<-lm(T.rlm_formula,data=T.lm_data,subset=unlist(cv[-2]))
T.lm_3<-lm(T.rlm_formula,data=T.lm_data,subset=unlist(cv[-3]))
T.lm_4<-lm(T.rlm_formula,data=T.lm_data,subset=unlist(cv[-4]))

T.rlm_1<-rlm(T.rlm_formula,data=T.lm_data,subset=unlist(cv[-1]),init = T.lm_1$coefficients)
T.rlm_2<-rlm(T.rlm_formula,data=T.lm_data,subset=unlist(cv[-2]),init = T.lm_1$coefficients)
T.rlm_3<-rlm(T.rlm_formula,data=T.lm_data,subset=unlist(cv[-3]),init = T.lm_1$coefficients)
T.rlm_4<-rlm(T.rlm_formula,data=T.lm_data,subset=unlist(cv[-4]),init = T.lm_1$coefficients)


#skalowanie danych genetycznych
genData.c.num<-genData.c.bin*diff(scalingRange_25_75)+scalingRange_25_75[1]
T.lm_data.N<-data.frame(numNonGen_Numerized_m_noOut_25_75[!genData.Miss,c(4,1:3,5:dim(numNonGen_Numerized_m_noOut_25_75)[2])],genData.c.num)
T.rlm_data.N<-T.lm_data.N[,-sing]


#Bez skalowania danych genetycznych i niegenetycznych danych binarnych.
T.final_data<-data.frame(numDataNorm_m_noOut[!genData.Miss,],nonGenData_m_noOut[!genData.Miss,],genData.c.bin)
T.final_data<-T.final_data[,c(4,1:3,5:55)] # testosteron na pierwszej pozycji
names(T.final_data)
rm_col<-c(7,29,35,41,30,46,49,52,55)#wynikaja z tego, ze 0 we wczesniejszych kolumnach juz definiuja co bedzie w tych + hyperandrogenism + FAI
#rm_col<-c(7,30)#usuniecie FAI i hyperandrogenizmu - 
T.final_data<-T.final_data[,-rm_col]
T.final_data<-T.final_data[!(cN.c$cl4==4 | cN.c$cl4==3 ),]#Usuniecie pacjentow z podywzszonym FSH - zaborzone funkcjonowanie jader
#T.final_data<-T.final_data[!(cN.c$cl4==4),]
T.final_formula<-formula(T.final_data)

T.init.lm<-lm(formula=T.final_formula,data=T.final_data,subset=sample(1:(dim(T.final_data)[1]),round(3/4*(dim(T.final_data)[1]))))#inicjalizacja parametrów dla rlm
T.init.lm<-coef(T.init.lm)			
T.final_rlm<-BK_lm(formula=T.final_formula,data=T.final_data,K=10,init=T.init.lm,method="rlm")
T.final_lm<-BK_lm(formula=T.final_formula,data=T.final_data,K=10,method="lm")
#Rysowanie
par(mar=c(9,4,2,2))
# ord<-order(T.final_lm$coef.mean$Mean.coef,decreasing=T)
#ord<-1:length(T.final_lm$coef.mean$Mean.coef)
ord<-order(T.final_rlm$coef.mean$p.value)
plot.bar<-data.frame(T.final_lm$coef.mean$Mean.coef[ord],T.final_rlm$coef.mean$Mean.coef[ord])
plot.CI_low<-data.frame(T.final_lm$coef.mean$CI_low[ord],T.final_rlm$coef.mean$CI_low[ord])
plot.CI_high<-data.frame(T.final_lm$coef.mean$CI_high[ord],T.final_rlm$coef.mean$CI_high[ord])
barp<-barplot(as.matrix(t(plot.bar)),names.arg=row.names(T.final_lm$coef.mean)[ord],col=c("red","blue"),beside=T,las=2,cex.names=0.65,ylim=c(-1,1),legend.text=c("Ordinary Least Squares","Tukey's Linear Model"))
errbar(barp,t(plot.bar),yplus=t(plot.CI_high),yminus=t(plot.CI_low),add=T,type="n",cap=0.01)
#FORWARD FEATURE SELECTION
T.lm.step.forward<-BK_cv_step(formula=TESTOSTERON~.0,data=T.final_data,direction="forward",K = 10,scope=T.final_formula)


#BACKWARD FEATURE SELECTION
T.lm_step.backward<-BK_cv_step(formula=T.final_formula,data=T.final_data,direction="backward",K = 10,scope=T.final_formula)


#LASSO REGRESSION

T.final_data_lasso<-data.frame(numDataNorm_m_noOut[!genData.Miss,],nonGenData_m_noOut[!genData.Miss,],genData.c.bin)
T.final_data_lasso<-T.final_data_lasso[,c(4,1:3,5:55)] # testosteron na pierwszej pozycji
rm_col<-c(7,30)#usuniecie FAI i hyperandrogenizmu
T.final_data_lasso<-T.final_data_lasso[,-rm_col]
T.final_data_lasso<-T.final_data_lasso[!(cN.c$cl4==4 | cN.c$cl4==3 ),]#Usuniecie pacjentow z podywzszonym FSH - zaborzone funkcjonowanie jad
#T.final_data_lasso<-T.final_data_lasso[!(cN.c$cl4==4 ),]
T.final_lasso<- cv.glmnet(as.matrix(T.final_data_lasso[,-1]), T.final_data_lasso$TESTOSTERON,standardize=FALSE)

#Podsumowanie korelacji Testosteronu:
T.final_lasso_coef<-data.frame(as.matrix(coefficients(T.final_lasso,s="lambda.min")),as.matrix(coefficients(T.final_lasso,s="lambda.1se")))


#*************
VitD.final_data<-data.frame(numDataNorm_m_noOut[!genData.Miss,],nonGenData_m_noOut[!genData.Miss,],genData.c.bin)
VitD.final_data<-VitD.final_data[,c(21,1:20,22:55)]#Witamina D jako pierwsza kolumna
VitD.rm_col<-c(29,35,41,46,49,52,55) #usuniecie nadmiarowych kolumn - MAKROREGIONwschodni,PORY.ROKUzima,FENOTYP.OTYLOSCI.FELAOMWD,+ genetyczne
VitD.final_data<-VitD.final_data[,-VitD.rm_col]
VitD.final_data<-VitD.final_data[!(cN.c$cl4==4 | cN.c$cl4==3 ),]
VitD.final_formula<-formula(VitD.final_data)

VitD.final_lm<-BK_lm(formula=VitD.final_formula,data=VitD.final_data,K = 10)
VitD.final_rlm<-BK_lm(formula=VitD.final_formula,data=VitD.final_data,K = 10,init=coef(VitD.final_lm$model.list$m2),method="rlm")


#Rysowanie
par(mar=c(9,4,2,2))
# ord<-order(T.final_lm$coef.mean$Mean.coef,decreasing=T)
#ord<-1:length(T.final_lm$coef.mean$Mean.coef)
ord<-order(VitD.final_rlm$coef.mean$p.value)
plot.bar<-data.frame(VitD.final_lm$coef.mean$Mean.coef[ord],VitD.final_rlm$coef.mean$Mean.coef[ord])
plot.CI_low<-data.frame(VitD.final_lm$coef.mean$CI_low[ord],VitD.final_rlm$coef.mean$CI_low[ord])
plot.CI_high<-data.frame(VitD.final_lm$coef.mean$CI_high[ord],VitD.final_rlm$coef.mean$CI_high[ord])
barp<-barplot(as.matrix(t(plot.bar)),names.arg=row.names(VitD.final_lm$coef.mean)[ord],col=c("red","blue"),beside=T,las=2,cex.names=0.65,ylim=c(-1,1),legend.text=c("Ordinary Least Squares","Tukey's Linear Model"))
errbar(barp,t(plot.bar),yplus=t(plot.CI_high),yminus=t(plot.CI_low),add=T,type="n",cap=0.01)
#FORWARD FEATURE SELECTION
VitD.lm_step.forward<-BK_cv_step(formula=WITAMINA.D~.0,data=VitD.final_data,direction="forward",K = 10,scope=VitD.final_formula)


#BACKWARD FEATURE SELECTION
VitD.lm_step.backward<-BK_cv_step(formula=VitD.final_formula,data=VitD.final_data,direction="backward",K = 10,scope=VitD.final_formula)

#LASSO REGRESSION

VitD.final_data_lasso<-data.frame(numDataNorm_m_noOut[!genData.Miss,],nonGenData_m_noOut[!genData.Miss,],genData.c.bin)
VitD.final_data_lasso<-VitD.final_data_lasso[,c(4,1:3,5:55)] # testosteron na pierwszej pozycji
#rm_col<-c(7,30)#usuniecie FAI i hyperandrogenizmu
#T.final_data_lasso<-T.final_data_lasso[,-rm_col]
VitD.final_data_lasso<-VitD.final_data_lasso[!(cN.c$cl4==4 | cN.c$cl4==3 ),]#Usuniecie pacjentow z podywzszonym FSH - zaborzone funkcjonowanie jad
#T.final_data_lasso<-T.final_data_lasso[!(cN.c$cl4==4 ),]
VitD.final_lasso<- cv.glmnet(as.matrix(VitD.final_data_lasso[,-1]), VitD.final_data_lasso$TESTOSTERON,standardize=FALSE)

#Podsumowanie korelacji Testosteronu:
VitD.final_lasso_coef<-data.frame(as.matrix(coefficients(VitD.final_lasso,s="lambda.min")),as.matrix(coefficients(VitD.final_lasso,s="lambda.1se")))
















#W danych nie ma FAI
T.lm_formula_noFAI<-formula(T.lm_data[trainID,-7]) #FAI - ind 7
T.lm_noFAI<-lm(T.lm_formula_noFAI,data=T.lm_data,subset=trainID)
T.predict_noFAI<-predict(T.lm_noFAI,T.lm_data[-trainID,-c(1,7)])
plot(T.lm_data$TESTOSTERON[-trainID],T.predict_noFAI)
lines(x = c(-2,2),y=c(-2,2),col = 'red')
points(T.lm_data$TESTOSTERON[trainID],predict(T.lm_noFAI,T.lm_data[trainID,-c(1,7)]),col='green')



#
TrnInd=sample(1:255,200)
testosteronNN<-nnet(all_d_numNorm[TrnInd,-c(4,20:99)],all_d_numNorm[TrnInd,4], size=20) #train neural network
#Making predictions on the test set
outputs<-predict(testosteronNN,all_d_numNorm[-TrnInd,-c(4,20:99)])

#Evaluation
##############################################

#####ANALIZA ZBIORU KOBIET
numericData_k<-numericData[kInd,]
rownames(numericData_k)<-rownames(numericData[kInd,])
numericData_k$INS<-log(numericData_k$INS)
numericData_k$TESTOSTERON<-log(numericData_k$TESTOSTERON)
numericData_k$TGC<-log(numericData_k$TGC)
numericData_k$ESTRADIOL<-log(numericData_k$ESTRADIOL)
numericData_k$FAI<-log(numericData_k$FAI) #transformacja powodowala, ze macierz danych byla zle uwarunkowana do obliczen Mahalanobisa, wiec transformacje tego elementu przeprowadzono pozniej.
numericData_k$FEI<-log(numericData_k$FEI)
numericData_k$FAT<-log(numericData_k$FAT)
numericData_k$GLUKOZA<-log(numericData_k$GLUKOZA)


numDataNorm_k<-data.frame(sapply(numericData_k,normalizeZRobust)) #Normalizacja
#BK_plotHists(numDataNorm_k)


#Identyfikacja wartości odstających metodą Mahalanobisa
md_NumData_k<-Moutlier(numericData_k,quantile=0.975,plot=FALSE)
plot(md_NumData_k$md,md_NumData_k$rd)
lines(c(md_NumData_k$cutoff,md_NumData_k$cutoff),c(0,150), col="red")
lines(c(0,17),c(md_NumData_k$cutoff,md_NumData_k$cutoff), col="red")
par(mfcol=c(2,1))
hist(md_NumData_k$md)
hist(md_NumData_k$rd)
outliers_NumData_k_md0975<-which(md_NumData_k$md>md_NumData_k$cutoff)
outliers_NumData_k_rd0975<-which(md_NumData_k$rd>md_NumData_k$cutoff)

numDataNorm_k_noOut<-numDataNorm_k[-outliers_NumData_k_md0975,]
#PCA bez wartosci odstajacych
numDataNorm_k_noOut_PCA<-BK_PCA_bplot_sdev(numDataNorm_k_noOut)

#Klasteryzacja
numDataNorm_k_HClust_wardD2_noOut<-hclust(dist(numDataNorm_k_noOut),method="ward.D2")


#Heat map - mapa cieplna po klastrowaniu
rowv<-as.dendrogram(hclust(dist(numDataNorm_k_noOut),method="ward.D2"))
colv<-as.dendrogram(hclust(dist(t(numDataNorm_k_noOut)),method="ward.D2"))
dev.new()
heatmap.2(as.matrix(numDataNorm_k_noOut),Rowv=rowv,Colv=colv,trace="none",symm=F,symkey=F,symbreaks=T, scale="none",margins=c(10,10))
heatmap.2(as.matrix(cor(numDataNorm_k_noOut)),trace="none",symm=F,symkey=F,symbreaks=T, scale="none",margins=c(10,10))
#Ocena klastrowania

numDataNorm_k_noOut_cut<-cutree(numDataNorm_k_HClust_wardD2_noOut,k=2:10)
numDataNorm_k_noOut_clv<-BK_clv(data=numDataNorm_k_noOut,clust=numDataNorm_k_noOut_cut,dist='euclidean')
dev.new()
par(mfcol=c(1,2))
plot(2:10,numDataNorm_k_noOut_clv$Davies.Bouldin,type='b')
plot(2:10,numDataNorm_k_noOut_clv$Dunn,type='b')
par(mfcol=c(1,1))
biplot(numDataNorm_k_noOut_PCA,xlabs = numDataNorm_k_noOut_cut[,"4"])


#UTWORZENIE CALEJ MACIERZY DANYCH - tylko kobiety i bez wartosci odstajacych
all_d_k<-data.frame(single_NoNA[row.names(numDataNorm_k_noOut),])
#Poprawienie pola WIELKOSC.MIEJSCOWOSCI
all_d_k[c("136","137"),"WIELKOSC.MIEJSCOWOSCI"]="wies"
all_d_k$WIELKOSC.MIEJSCOWOSCI<-factor(all_d_k$WIELKOSC.MIEJSCOWOSCI,levels(all_d$WIELKOSC.MIEJSCOWOSCI))

#Poprawienie pola GR.WIEKOWA
all_d_k[which(all_d[,"GR.WIEKOWA"]=="90 i wiecej"),"GR.WIEKOWA"]="90 lat i wiecej"
all_d_k$GR.WIEKOWA<-factor(all_d_k$GR.WIEKOWA)

#Uzupelnienie danych o cukrzycy
all_d_k[40,"CG1.ROZPOZNANA.WCZESNIEJ.CUKRZYCA"]="nie"
all_d_k[41,"CG1.ROZPOZNANA.WCZESNIEJ.CUKRZYCA"]="tak"
#Uzupelnienie danych o naslonecznieniu
all_d_k[which(all_d_k$NASLONECZNIENIE==""),"NASLONECZNIENIE"]="TAK" #Wszystkie probki brakujace byly z wrzesnia - dla ktorego NASLONECZNIENIE=TAK
#uzupelnienie danych o nadcisnieniu
#CZtery pacjentki dwie z gr 1 po jednej z gr1 i 3. We wszystkch grupach prawdopodobienstwo nadcisnienia jest wieksze wiec wszystkie brakujace nadcisnienia ustawione na "TAK"
table(all_d_k$NADCISNIENIE.TETNICZE,numDataNorm_k_noOut_cut[,"4"])
all_d_k[all_d_k$NADCISNIENIE.TETNICZE=="","NADCISNIENIE.TETNICZE"]="tak"

#ZAMIANA WARTOSCI NOMINALNYCH - NIEGENETYCZNYCH - NA CIAGI BINARNE
nominalData_Ind<-which(dataTypes=="factor")
#Usuniecie atrybutow:
#1 - KOD,2 - PLEC, 3- GR.WIEKOWA (zamiana na kategorie stopniowane) 5 - WOJEWODZTWO, 6 - WIELKOSC.MIEJSCOWOSCI (zamiana na kategorie stopniowane) 
#10 - ZAKRES.WITAMINY.D 11 - VDR..CDX., 16-27 - MODELE, 28 - MIESIAC,31 - WITD - zle ustalony prog, 32-34 - OBZM, OZZM, OMPMC - brak zmiennosci
nominalData_Ind<-nominalData_Ind[-c(1,2,3,5,6,10,11,16:28,31,32:34)]
#Dane niegenetyczne

nonGenData_Ind<-nominalData_Ind[c(1:4,9:12)] #wyciete dane genetyczne
lapply(all_d_k[,nonGenData_Ind],table)

nonGenData_k_Factors<-data.frame(lapply(all_d_k[,nonGenData_Ind],factor))
rownames(nonGenData_k_Factors)<-rownames(all_d_k)
lapply(nonGenData_k_Factors,table)
nonGenData_k_noOut<-model.matrix(~ .+0, nonGenData_k_Factors) 

#Dolaczenie danych o grupie wiekowej i rozmiarach miejscowosci
nonGenData_k_noOut<-data.frame(nonGenData_k_noOut,as.numeric(all_d_k[,"WIELKOSC.MIEJSCOWOSCI"]),as.numeric(all_d_k[,"GR.WIEKOWA"]))
names(nonGenData_k_noOut)[dim(nonGenData_k_noOut)[2]-1]<-"WIELKOSC.MIEJSCOWOSCI"
names(nonGenData_k_noOut)[dim(nonGenData_k_noOut)[2]]<-"GR.WIEKOWA"
nonGenData_k_noOut$GR.WIEKOWA<-normalizeZRobust(nonGenData_k_noOut$GR.WIEKOWA)
nonGenData_k_noOut$WIELKOSC.MIEJSCOWOSCI<-normalizeZRobust(nonGenData_k_noOut$WIELKOSC.MIEJSCOWOSCI)

#Skalowanie 0 i 1 do kwantylu 5% i 95% danych numerycznych
scalingRange_k_5_95<-quantile(unlist(numDataNorm_k_noOut),probs=c(0.05,0.95))#probs=c(0.25,0.75))
numNonGenDataNorm_k_noOut_5_95<-nonGenData_k_noOut*diff(scalingRange_k_5_95)+scalingRange_k_5_95[1]

#Skalowanie 0 i 1 do kwantyli 25% i 75%
scalingRange_k_25_75<-quantile(unlist(numDataNorm_k_noOut),probs=c(0.25,0.75))
numNonGenDataNorm_k_noOut_25_75<-nonGenData_k_noOut*diff(scalingRange_k_25_75)+scalingRange_k_25_75[1]


#POLCZENIE DANYCH NUMERYCZNYCH i NOMINALNYCH NIEGENETYCZNYCH

numNonGen_Numerized_k_noOut_5_95<-data.frame(numDataNorm_k_noOut,numNonGenDataNorm_k_noOut_5_95)
numNonGen_Numerized_k_noOut_25_75<-data.frame(numDataNorm_k_noOut,numNonGenDataNorm_k_noOut_25_75)


#ANALIZY NA PELNYCH DANYCH


numNonGen_k_5_95_HClust_wardD2<-hclust(dist(numNonGen_Numerized_k_noOut_5_95),method="ward.D2")
plot(numNonGen_k_5_95_HClust_wardD2)

numNonGen_k_25_75_HClust_wardD2<-hclust(dist(numNonGen_Numerized_k_noOut_25_75),method="ward.D2")
plot(numNonGen_k_25_75_HClust_wardD2)



#Heat map - mapa cieplna po klastrowaniu
rowv<-as.dendrogram(hclust(dist(numNonGen_Numerized_k_noOut_5_95),method="ward.D2"))
colv<-as.dendrogram(hclust(dist(t(numNonGen_Numerized_k_noOut_5_95)),method="ward.D2"))
dev.new()
heatmap.2(as.matrix(numNonGen_Numerized_k_noOut_5_95),Rowv=rowv,Colv=colv,trace="none",symm=F,symkey=F,symbreaks=T, scale="none",margins=c(10,10))
dev.new()
heatmap.2(as.matrix(cor(numNonGen_Numerized_k_noOut_5_95)),margins=c(10,10))

numNonGen_Numerized_k_noOut_5_95_HClust_wardD2<-hclust(dist(numNonGen_Numerized_k_noOut_5_95),method="ward.D2")
numNonGen_Numerized_k_noOut_5_95_cut<-cutree(numNonGen_Numerized_k_noOut_5_95_HClust_wardD2,k=2:10)
numNonGen_Numerized_k_noOut_5_95_clv<-BK_clv(data=numNonGen_Numerized_k_noOut_5_95,clust=numNonGen_Numerized_k_noOut_5_95_cut,dist='euclidean')

#Heat map - mapa cieplna po klastrowaniu
rowv<-as.dendrogram(hclust(dist(numNonGen_Numerized_k_noOut_25_75),method="ward.D2"))
colv<-as.dendrogram(hclust(dist(t(numNonGen_Numerized_k_noOut_25_75)),method="ward.D2"))
dev.new()
heatmap.2(as.matrix(numNonGen_Numerized_k_noOut_25_75),Rowv=rowv,Colv=colv,trace="none",symm=F,symkey=F,symbreaks=T, scale="none",margins=c(10,10))
dev.new()
heatmap.2(as.matrix(cor(numNonGen_Numerized_k_noOut_25_75)),margins=c(10,10))

numNonGen_Numerized_k_noOut_25_75_HClust_wardD2<-hclust(dist(numNonGen_Numerized_k_noOut_25_75),method="ward.D2")
numNonGen_Numerized_k_noOut_25_75_cut<-cutree(numNonGen_Numerized_k_noOut_25_75_HClust_wardD2,k=2:10)
numNonGen_Numerized_k_noOut_25_75_clv<-BK_clv(data=numNonGen_Numerized_k_noOut_25_75,clust=numNonGen_Numerized_k_noOut_25_75_cut,dist='euclidean')

#Wizualizacja wartosci odstajacych i klastrow
d_k<-data.frame(md_NumData_k$md[-outliers_NumData_k_md0975], md_NumData_k$rd[-outliers_NumData_k_md0975],numDataNorm_k_noOut_cut)
names(d_k)<-c("md","rd",paste("cl",colnames(numDataNorm_k_noOut_cut),sep=""))
gg<-ggplot(data=d_k)
gg+geom_point(aes(y=rd,x=md,colour=factor(cl4),shape=factor(cl4),size=5))

plot(2:10,numNonGen_Numerized_k_noOut_25_75_clv$Davies.Bouldin,type="b")
plot(2:10,numNonGen_Numerized_k_noOut_25_75_clv$Dunn,type="b")
plot(2:10,numNonGen_Numerized_k_noOut_5_95_clv$Dunn,type="b")
plot(2:10,numNonGen_Numerized_k_noOut_5_95_clv$Davies.Bouldin,type="b")









#Dane genetyczne

genData_k_Ind<-which(regexpr(pattern = "*VDR*",text = names(all_d_k))==1)
genData_k<-data.frame(lapply(all_d_k[,genData_k_Ind],factor))
genData_k<-genData_k[,-1]
row.names(genData_k)<-row.names(all_d_k)

#Poprawienie wartosci 
levels(genData_k[,"VDR..BSM."])[1]<-"BB"
levels(genData_k[,"VDR..BSM."])[2]<-"Bb"
levels(genData_k[,"VDR..BSM."])[3]<-"bb"

levels(genData_k[,"VDR..FOK."])[1]<-"FF"
levels(genData_k[,"VDR..FOK."])[2]<-"Ff"
levels(genData_k[,"VDR..FOK."])[3]<-"ff"

genData_k.Miss<-(genData_k$VDR..APA.=="" | genData_k$VDR..TAQ.=="")
genData_k.c<-genData_k[!genData_k.Miss,] #complete genData
genData_k.c$VDR..TAQ.<-factor(genData_k.c$VDR..TAQ.)
genData_k.c$VDR..APA.<-factor(genData_k.c$VDR..APA.)

levels(genData_k.c$VDR..TAQ.)<-c(0,1,2)
levels(genData_k.c$VDR..APA.)<-c(0,1,2)
levels(genData_k.c$VDR..FOK.)<-c(0,1,2)
levels(genData_k.c$VDR..BSM.)<-c(0,1,2)

genData_k.c.bin<-BK_binGenData(genData_k.c)


genData_k_freq<-lapply(genData_k,table)
genData_k.c_freq<-lapply(genData_k.c,table)


genData_k_noOut<-model.matrix(~ .+0, genData_k.c)




#Budowa modelu liniowego
T.final_data_k<-data.frame(numDataNorm_k_noOut[!genData_k.Miss,],nonGenData_k_noOut[!genData_k.Miss,],genData_k.c.bin)
T.final_data_k<-T.final_data_k[,c(4,1:3,5:56)]
rm_col_k<-c(7,29,38,42,47,50,53,56,35,30) #usunieta kolumna FAI (7), hyperandrogenizm, oraz kolumny, ktore powoduja osobliwosc macierzy
T.final_data_k<-T.final_data_k[,-rm_col_k]
T.final_formula_k<-formula(T.final_data_k)


T_k.init.lm<-lm(formula=T.final_formula_k,data=T.final_data_k,subset=sample(1:193,120))#inicjalizacja parametrów dla rlm
T_k.init.lm<-coef(T_k.init.lm)  		
T_k.final_rlm<-BK_lm(formula=T.final_formula_k,data=T.final_data_k,K=10,init=T_k.init.lm,method="rlm")
T_k.final_lm<-BK_lm(formula=T.final_formula_k,data=T.final_data_k,K=10,method="lm")

cvfit <- glmnet::cv.glmnet(as.matrix(T.final_data_k[,-1]), T.final_data_k$TESTOSTERON)
coef(cvfit, s = "lambda.1se")


#Rysowanie
par(mar=c(9,4,2,2))
# ord<-order(T.final_lm$coef.mean$Mean.coef,decreasing=T)
#ord<-1:length(T.final_lm$coef.mean$Mean.coef)
ord<-order(T_k.final_rlm$coef.mean$p.value)
plot.bar<-data.frame(T_k.final_lm$coef.mean$Mean.coef[ord],T_k.final_rlm$coef.mean$Mean.coef[ord])
plot.CI_low<-data.frame(T_k.final_lm$coef.mean$CI_low[ord],T_k.final_rlm$coef.mean$CI_low[ord])
plot.CI_high<-data.frame(T_k.final_lm$coef.mean$CI_high[ord],T_k.final_rlm$coef.mean$CI_high[ord])
barp<-barplot(as.matrix(t(plot.bar)),names.arg=row.names(T_k.final_lm$coef.mean)[ord],col=c("red","blue"),beside=T,las=2,cex.names=0.65,ylim=c(-1,1),legend.text=c("Ordinary Least Squares","Tukey's Linear Model"))
errbar(barp,t(plot.bar),yplus=t(plot.CI_high),yminus=t(plot.CI_low),add=T,type="n",cap=0.01)
#FORWARD FEATURE SELECTION
T_k.lm.step.forward<-BK_cv_step(formula=TESTOSTERON~.0,data=T_k.final_data,direction="forward",K = 10,scope=T_k.final_formula)


#BACKWARD FEATURE SELECTION
T_k.lm_step.backward<-BK_cv_step(formula=T.final_formula_k,data=T.final_data_k,direction="backward",K = 10,scope=T.final_formula_k)


#LASSO REGRESSION

T_k.final_data_lasso<-data.frame(numDataNorm_k_noOut[!genData_k.Miss,],nonGenData_k_noOut[!genData_k.Miss,],genData_k.c.bin)
T_k.final_data_lasso<-T_k.final_data_lasso[,c(4,1:3,5:dim(T_k.final_data_lasso)[2])] # testosteron na pierwszej pozycji
rm_col<-c(7,30)#usuniecie FAI i hyperandrogenizmu
T_k.final_data_lasso<-T_k.final_data_lasso[,-rm_col]
T_k.final_lasso<- cv.glmnet(as.matrix(T_k.final_data_lasso[,-1]), T_k.final_data_lasso$TESTOSTERON,standardize=FALSE)

#Podsumowanie korelacji Testosteronu:
T_k.final_lasso_coef<-data.frame(as.matrix(coefficients(T_k.final_lasso,s="lambda.min")),as.matrix(coefficients(T_k.final_lasso,s="lambda.1se")))



