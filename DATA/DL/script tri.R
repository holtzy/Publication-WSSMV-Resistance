# script analyse QTLS
rm(list=ls())


# le fichier DL
setwd("C:/Users/david/Dropbox/Publi_Mosaique/DATA/DATA/DS")
setwd("C:/Users/david/Dropbox/VIRUS/COMPTE_RENDU_PRIVES/FD/FICHIERS")

# ouverture du fichier de phenotypes

pheno_DL<-read.table("phenotypage.csv",
                     sep=";", dec=".", header=TRUE, na.strings = "NA")

names(pheno_DL)

# ouverture du ficher de genotypes

geno_DL<-read.table("genotypage.csv",
                    sep=";", dec=".", header=TRUE, na.strings = "NA")

geno_DL[1:2, 1:2]

# les p values
data_DL<-read.table("bilan_simple_marker",
                      sep=",", dec=".", header=TRUE, na.strings = "NA")
                    
                    
levels(data_DL$variable)

data_DL_V <-data_DL[which(data_DL$variable=="P2016_ElisaPondPos_Blup"),]

CHR7A<-data_DL_V[which(data_DL_V$LG=="7A"),]
CHR7A<-CHR7A[which(!(is.na(CHR7A$LOD))),]

head(CHR7A)

CHR7B<-data_DL_V[which(data_DL_V$LG=="7B"),]
CHR7B<-CHR7B[which(!(is.na(CHR7B$LOD))),]

CHR2A<-data_DL_V[which(data_DL_V$LG=="2A"),]
CHR2A<-CHR2A[which(!(is.na(CHR2A$LOD))),]

plot(CHR7A$Distance, CHR7A$LOD  )
plot(CHR7B$Distance, CHR7B$LOD  )
plot(CHR2A$Distance, CHR2A$LOD  )


MAXI_7A<-max(CHR7A$LOD, na.rm=TRUE)
MAXI_7B<-max(CHR7B$LOD, na.rm=TRUE)
MAXI_2A<-max(CHR2A$LOD, na.rm=TRUE)

# la liste des marqueurs 
as.character(CHR7A[CHR7A$LOD==MAXI_7A,"marqueur"])
as.character(CHR7B[CHR7B$LOD==MAXI_7B,"marqueur"])
as.character(CHR2A[CHR2A$LOD==MAXI_2A,"marqueur"])

NomMArQ_7A<-as.character(CHR7A[CHR7A$LOD>=MAXI_7A-1,"marqueur"])
NomMArQ_7B<-as.character(CHR7B[CHR7B$LOD>=MAXI_7B-0.5,"marqueur"])
NomMArQ_2A<-as.character(CHR2A[CHR2A$LOD>=MAXI_2A-0.1,"marqueur"])

NomMArQ_7A<-gsub("@",".",NomMArQ_7A)
NomMArQ_7A<-gsub("\\|",".",NomMArQ_7A)

NomMArQ_7B<-gsub("@",".",NomMArQ_7B)
NomMArQ_7B<-gsub("\\|",".",NomMArQ_7B)

NomMArQ_2A<-gsub("@",".",NomMArQ_2A)
NomMArQ_2A<-gsub("\\|",".",NomMArQ_2A)

# pour avoir le génotype
row.names(geno_DL)<-as.character(geno_DL[,1])
  
head(colnames(geno_DL))
haplo<-geno_DL[, c(NomMArQ_7A, NomMArQ_7B, NomMArQ_2A) ]  
haplo<-geno_DL[, NomMArQ_7B]

dim(haplo)

GENO_PHENO<-merge(haplo,pheno_DL[,c("Geno","P2016_ElisaPondPos_Blup")], by.x=0, by.y=1 )
GENO_PHENO<-GENO_PHENO[order(GENO_PHENO$"P2016_ElisaPondPos_Blup"),]


write.csv2(GENO_PHENO, file = "GENO_PHENO_DB_7B.csv", row.names = TRUE, sep=";" ,col.names = TRUE)
