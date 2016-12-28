

			#--------------------------------------------------------------------------------------------------
			#   
			#		SCRIPT R : CALCULATE INTERACTION BETWEEN ALL PAIRS OF MARKERS
			#
			#					author: Holtz Yan
			#---------------------------------------------------------------------------------------------------



##-- GOAL
	# To answer the request of a reviewer, we decided to calculate the interaction between all pairs of markers. To avoid computation time, we used one marker per slice of 5cM.
 
##-- INPUT (in the right order)
	# Genotyping matrix
	# Phenotyping matrix
	# Genetic map
	# Chromosome of reference

## -- Récupération des Arguments
args <- commandArgs(trailingOnly = TRUE)
fic_geno=args[1]
fic_pheno=args[2] 
fic_map=args[3]
ref_chromo=args[4]
ref_pheno=args[5]




#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# -------------------------------------------------
# PARTIE 1 : READ INPUT.
# -------------------------------------------------


# --- Genotyping matrix
# Le format du génotpage doit respecter les règles suivantes : 1/ Nom du fichier = fic_genotypage.csv 2/ en tete des colonnes = nom des geno 3/ données manquantes = "-" 4/ Séparateur=";" 
genotype<-read.table(fic_geno, sep = ";" , header = F, na.strings = "-")
genotype=as.matrix(genotype)
colnames(genotype)=genotype[1,]
genotype=as.data.frame(genotype[-1 , ])
names(genotype)[1]<-"geno"
print("--- Your genotyping matrix looks correct. Dimension of the matrix are :")
print(dim(genotype))


# --- Map
# Format de la carte : 3 colonnes : LG, nom du marqueur, position dans le LG
map <- read.table(fic_map , header=T , dec = ".", na.strings = "-" , check.names=F)
colnames(map) <- c("LG", "marqueur", "Distance","group_physique","Posi_physique")
rownames(map) <- map$marqueur
map$LG <- as.factor(map$LG)
map=map[ which(map$marqueur%in%colnames(genotype)) , ]
print("--- Your genetic map looks correct. Dimension of the map are :")
print(dim(map))
# I remove the rare singlet markers that make a bug:
map=map[-grep("singlet", map$marqueur) , ]



# --- Phenotyping matrix
BLUP<-read.table(fic_pheno, header = TRUE, sep=";")
colnames(BLUP)[1]="geno"
print("--- Your Phenotyping matrix looks correct. Dimension of the matrix are :")
print(dim(BLUP))
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#






#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# -------------------------------------------------
# PARTIE 2 : MERGE GENOTYPE AND PHENOTYPE
# -------------------------------------------------

BLUP[,1]<-as.character(BLUP[,1])
genotype[,1]<-as.character(genotype[,1])
don<-merge(BLUP,genotype, by.x=1 , by.y=1, all=T)
print("--- Nombre d'individu communs entre pheno et géno :")
dim(don)[1]

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#





#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# -------------------------------------------------
# PARTIE 2 : CALCULATE INTERACTIONS.
# -------------------------------------------------


# Create an empty summary file!
bilan=data.frame(matrix(0,1000000,13))
colnames(bilan)=c("pheno","chromo_ref","min_Kref","max_Kref","mark_ref","chromo_other","min_Kother","max_Kother","mark_Kother","R2","pval_mark1", "pval_mark2", "pval_inter")

# features of the chromosome of reference. We are going to slice it by 5cM! --> print the number of slices:
max_map=max(map$Distance[map$LG==ref_chromo] , na.rm=T)
print("the number of slices for this chromosome is:")
print(max_map/5)

# With which chromosome do I have to do the comparison? (only chromosome higher than the reference chromosome)
vec_chromo=levels(map$LG)
K_for_comparison=vec_chromo[c(match(ref_chromo , vec_chromo):length(vec_chromo))]
print("lets calculate interaction with chromosomes:")
print(K_for_comparison)

# Let's start the slicing!
num=0
for(min_Kref in seq(0,max_map,5)){
	
	# Find a random marker in the slice of the reference chromosome:
	max_Kref=min_Kref+5
	markers_in_slice=map$marqueur[which(map$LG==ref_chromo & map$Distance>=min_Kref & map$Distance<max_Kref) ]
	marker1=markers_in_slice[ sample(1:length(markers_in_slice),1)]
	print(paste("you are on slice: ",min_Kref," ----> ",max_Kref," -- size of chromo is: ",max_map,sep=""))
	
	# We need to compare it with slices of all other chromosomes!
	for(chromo in K_for_comparison ){
		print("chromo en cours:")
		print(chromo)
		my_max=max(map$Distance[map$LG==chromo] , na.rm=T)

		# Slice this new chromosome
		for(min_Kother in seq(0,my_max,5)){
			
			# find a marker for the slice:
			max_Kother=min_Kother+5
			markers_in_slice=map$marqueur[which(map$LG==chromo & map$Distance>=min_Kother & map$Distance<max_Kother) ]
			if(length(markers_in_slice)==0){next}
			marker2=markers_in_slice[ sample(1:length(markers_in_slice),1)]

			# calculate interaction
			if( length(which(c(as.character(marker1),as.character(marker2))%in%colnames(don)))!=2){next}
			data=don[ , c(as.character(ref_pheno), as.character(marker1), as.character(marker2))]
			model=lm(data[,1] ~ data[,2] * data[,3])
			sum=summary(model)
			my_r2=round(sum$r.squared,3)
			my_pval=if(nrow(sum$coefficients)>3)  as.vector(round(sum$coefficients[c(2:4),4],6))  else rep(NA,3)
			
			# fill the summary table
			num=num+1
			bilan[num , c(1:13)]=c(as.character(ref_pheno), ref_chromo , min_Kref, max_Kref, as.character(marker1), chromo, min_Kother, max_Kother, as.character(marker2), my_r2, my_pval )
		
	}}}
bilan=bilan[bilan$chromo_ref==ref_chromo , ]



# Save the result!
write.table(bilan, file=paste("result_epistasie_K_",ref_chromo,"_",ref_pheno,".txt",sep=""), col.names=T, row.names=T , quote=F)











