#R script for BRIC LED RNA transcript
#Adapted from Wyseq by Gbolaga Olanrewaju
#6/5/2023

suppressMessages(require(edgeR))
suppressMessages(require(ggplot2))
suppressMessages(require(grid))
suppressMessages(require(limma))
library(KEGGREST)
library(ggrepel)  #to avoid text overlap in plot
library(ape)   #For PCA
library(biomaRt)
library(DESeq2)
library(ggforce)
suppressMessages(require(NMF))
suppressMessages(require(plyr))
suppressMessages(require(reshape2))
suppressMessages(require(scales))
suppressMessages(require(gridExtra))
suppressMessages(require(stringr))
suppressMessages(require(dplyr))
suppressMessages(require(pheatmap))
suppressMessages(require(tidyr))
suppressMessages(require(cowplot))
suppressMessages(require(gplots))
suppressMessages(require(WGCNA))
suppressMessages(require(data.table))
suppressMessages(require(genefilter))
suppressMessages(require(ParamHelpers))
suppressMessages(require(S4Vectors))
suppressMessages(require(IRanges))
suppressMessages(require(AnnotationDbi))
options(stringsAsFactors = F)
allowWGCNAThreads()

plant_species <- "arabidopsis"
#*********************************************************************
#to average read counts between technical replicates
#***********************************************************************
#
#To remove the first 3 columns which contained file description
# List all the txt files in the directory

#Let your files be in your working directory
files <- list.files(pattern = "*.txt$")

# Loop over each file
for(file in files) {
  # Read the file, skipping the first three rows
  data <- read.delim(file, header = FALSE)
  
  # Write the data back to the original txt file with tab-separated columns
  write.table(data, file, row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
}

#Specify the column names again

files <- list.files(pattern = "*.txt$")

# Define your custom column names
column_names <- c("GeneID", "Chr", "Start", "End", "Strand", "Length", "Count")  # replace with your column names

# Loop over each file
for(file in files) {
  # Read the file, skipping the first three rows
  data <- read.delim(file, skip = 3, header = FALSE)
  
  # Assign the custom column names to the data
  names(data) <- column_names
  
  # Write the data back to the original txt file with tab-separated columns
  write.table(data, file, row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#THIS CREATES AN AVERAGE COUNT FOR EACH TECHNICAL REPLICATES IN THE FOLDER
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
files <- list.files(pattern = "*.txt$")

# Extract base names from the file names (remove "_R1" or "_R2" and the extension ".txt")
base_names <- unique(sub("(R1|R2).*", "", files))

# Create an empty list to hold the data from each pair of replicates
data_list <- list()

for (name in base_names) {
  # Get replicates
  replicate_files <- grep(name, files, value = TRUE)
  
  # Read R1 and R2 count files using read.delim, skipping the first three rows
  R1 <- read.delim(replicate_files[grep("R1", replicate_files)], header = TRUE)
  R2 <- read.delim(replicate_files[grep("R2", replicate_files)], header = TRUE)
  
  # Calculate the average count
  average_count <- (R1$Count + R2$Count) / 2
  
  # Create a data frame with R1 counts, R2 counts, and the average count
  data <- cbind(R1, R2_Counts = R2$Count, Average_Counts = average_count)
  
  # Add it to the list
  data_list[[name]] <- data
}

# Write each data frame to a separate txt file
for (name in base_names) {
  write.table(data_list[[name]], paste0(name, "_combined_counts.txt"), row.names = FALSE, sep = "\t")
}

#......................................................
# Create a data frame with your sample information
#.......................................................

sample_info <- data.frame(Sample = c("FT-A2.txt", "FT-B2.txt", "FT-C2.txt", "FT-D2.txt", "FT-E2.txt", "FT-F2.txt","GC-A2.txt", "GC-B2.txt", "GC-C2.txt", "GC-D2.txt", "GC-E2.txt", "GC-F2.txt" ),
                          condition = c("Flight", "Flight", "Flight", "Flight", "Flight", "Flight", "Ground", "Ground", "Ground", "Ground", "Ground", "Ground"))

# Write the data frame to a tab-delimited text file
write.table(sample_info, file = "sample_info.txt", quote = FALSE, sep = "\t", row.names = FALSE)



#..................................................
#Concatenate spaceflight and Ground control
#................................................
# Creating vector of suffixes
suffixes <- c("A2", "B2", "C2", "D2", "E2", "F2")
#Create vector of file names
files <- c(paste0("FT-", suffixes, ".txt"), paste0("GC-", suffixes, ".txt"))

# Read all files
data_list <- lapply(files, function(x) read.table(x, header = TRUE))

# Concatenate all data frames

# Read all files, and keep only the 1st and 7th columns which are GeneID and Average count
data_list <- lapply(files, function(x) {
  df <- read.table(x, header = TRUE)
  df <- df[, c(1, 7)]
  colnames(df) <- c("GeneID", x)
  return(df)
})

# Bind all data frames together by columns
library(dplyr)
counts <- bind_cols(data_list)

# Set the GeneID as row names
rownames(counts) <- counts$GeneID

# Remove the redundant GeneID columns
counts <- counts[, -seq(3, ncol(counts), by = 2)]

#Create a lib_size list
b <- as.matrix(counts [,-1])
b <- ceiling(b)   #specify b to upper integer
rownames(b)<-counts[,1]

b<-b[rowSums(!as.matrix(b)) < ncol(b), ]  #Remove genes with no count
colnames(b) <- gsub("\\.txt$", "", colnames(b))   #Remove .txt from column names
##......................................
# merge exons of same gene as one. Find the aggregrate
#...................................................................
gene_keys <- sapply(strsplit(rownames(b), "\\."), function(x) {
  if (length(x) >= 1) {
    return(x[1])  # return the gene part
  } else {
    return(NA)  # return NA for rows where extraction fails
  }
})

f <- as.data.frame(b)
f$Gene <- gene_keys
b<- rowsum(f[, -ncol(f)], group = f$Gene, na.rm = TRUE)

#Create your lib_size

d<- colSums(b)
lib_size <- data.frame(
  Sample = names(d),
  TotalSum = d)
lib_sizes<-as.data.frame(d)



#.....................
#Normalize counts with Voom
#....................

y <- DGEList(counts = b[rowSums(cpm(as.matrix(b)) > 3) >= 3, ])

y <- calcNormFactors(y, method = "TMM")
design<- model.matrix(~ condition, data = sample_info)

v <- voomWithQualityWeights(y, design = design, method = "reml", plot = TRUE)


###############
#Run DESeq
###############

dds <- DESeqDataSetFromMatrix(countData = b,
                              colData = sample_info,
                              design = ~ condition)
dds <- DESeq(dds)  #run DESeq

res <- results(dds)

# Normalize counts for PCoA explained variance calculation
rld <- rlog(dds)
dist_matrix <- dist(t(assay(rld))) # Calculate distance matrix
pcoa <- pcoa(dist_matrix)

#Plot distance matrix heat map

output_file <- "heatmap.png"
# Plot heatmap of the distance matrix
heatmap(as.matrix(dist_matrix), main = "Distance Matrix Heatmap")
png(output_file)

#Hierachical cluster
output_file <- "Cluster.png"
hc <- hclust(dist_matrix)
plot(hc, main = "Dendrogram")
png(output_file)

#.....................
# MAKE PCoA
#....................

nsamples <- ncol(v)
nprobes <- nrow(v)
v.plot <- as.matrix(v)

dd <- matrix(0,nrow=nsamples,ncol=nsamples)
topindex <- nrow(v)-500 +1L
for (i in 2L:(nsamples))
  for (j in 1L:(i-1L))
    dd[i,j]=sqrt(mean(sort.int((v.plot[,i]-v.plot[,j])^2,partial=topindex)[topindex:nprobes]))

mds.out <- as.data.frame(suppressWarnings(cmdscale(as.dist(dd),k=3)))
mds.out$Fraction <- strtrim(colnames(b),1)
mds.out$Sample <- factor(rownames(lib_sizes), levels = rownames(lib_sizes))
eigen_values <- pcoa$values$Eigenvalues  #Calculate the Explained variance from eigen values


ggplot(mds.out, aes(x=V1, y=V2, color=Fraction)) + 
  scale_color_manual(values = c("#f79646", "#5b9bd5"))+
  ylim(-1.5, 1.5) +
  xlim(-1.5, 1.5) +
  ylab(paste0("PCo2 (Explained Variance: ", round((eigen_values[2] / sum(eigen_values)) * 100, 2), "%)")) +
  xlab(paste0("PCo1 (Explained Variance: ", round((eigen_values[1] / sum(eigen_values)) * 100, 2), "%)")) +
  #ylab(bquote('Leading' ~Log[2]~ 'Fold Change Dim 2')) +
  #xlab(bquote('Leading' ~Log[2]~ 'Fold Change Dim 1')) +
  theme_bw() +
  theme(
    text = element_text(size = 16, color = "black"),
    axis.text = element_text(color = "black")
  ) +
  ggforce::geom_mark_ellipse(aes(color = Fraction)) +
  geom_text(aes(label = sample_info$Sample))
fit<-lmFit(v,design)
ggsave("PCA.png", width = 8, height = 6)


#.....................
# Create Contrast Matrix
#....................
design.pairs<-function(levels){
  n <- length(levels) 
  design <- matrix(0,n,choose(n,2)) 
  rownames(design) <- levels 
  colnames(design) <- 1:choose(n,2) 
  k <- 0 
  for(i in 1:(n-1)){ 
    for(j in (i+1):n){ 
      k <- k+1 
      design[i,k] <- 1 
      design[j,k] <- -1 
      colnames(design)[k] <- paste(levels[i],"_v_",levels[j],sep="")
    }
  }
  design
}

levels=colnames(design)
contMatrix<-design.pairs(levels = levels)
DE<-contrasts.fit(fit,contMatrix)
DE<-eBayes(DE)



##############################################################################
#Calculate the CPM, Annotate to TAIR and Annotation Dbi and merge sample data
##############################################################################
# Calculate CPM values
cpm_counts <- cpm(b)

# Split the counts based on the condition
ground_counts <- cpm_counts[, sample_info$condition == "Ground"]
flight_counts <- cpm_counts[, sample_info$condition == "Flight"]

# Calculate the average CPM for each gene in each condition
avg_cpm_G <- rowMeans(ground_counts)
avg_cpm_F <- rowMeans(flight_counts)

# Add these to your result data frame
res$avg_cpm_G <- avg_cpm_G
res$avg_cpm_F <- avg_cpm_F

#.....................
#Annotate to ENSEMBLE
#....................

# Specify the Ensembl dataset and organism
ensembl_dataset <- "athaliana_eg_gene"
organism <- "athaliana"

# Connect to the Ensembl database
ensembl <- useMart(biomart = "plants_mart",
                   dataset = "athaliana_eg_gene",
                   host = "https://plants.ensembl.org")

# Convert row names to gene symbols
symbols <- getBM(attributes = c('ensembl_gene_id','description','external_gene_name','entrezgene_id'),
                 filters = "ensembl_gene_id",
                 values = rownames(res),
                 mart = ensembl)

Annotations<-symbols

# Merge annotations with DESeq2 result

# Convert DESeqResults object to data frame
res_df <- as.data.frame(res)
#Make genID the first column of res_df
res2_df <- data.frame(RowHeader = rownames(res_df), res_df)

# Remove row names from the data frame
rownames(res2_df) <- NULL
print(res2_df)

#Rename column 1 to match ensemble dataframe
names(res2_df)[1] <- "ensembl_gene_id"

w<-merged_df <- merge(res2_df, symbols, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")


counts_df <- as.data.frame(b)
counts_df <- data.frame(RowHeader = rownames(counts_df), counts_df)
rownames(counts_df) <- NULL
names(counts_df)[1] <- "ensembl_gene_id"
w<-merged_df <- merge(w, counts_df, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")


##Tidy data column
q<- w[,10:12]
w<-w[,-c(10:12)]
w <- cbind(w[, 1:2], q, w[, 2:ncol(w)])

w<- w[!is.na(w$padj), , drop = FALSE]

write.csv (w, file ="BRIC.csv", row.names = FALSE)



###############
#write outputs
###############


print("Writing counts data to file.")
Counts_CPM<-cpm(as.matrix(b))
NormalizedCountData<-v$E
write.csv(b,file = "AllCounts_RawReadCounts.csv")
write.csv(NormalizedCountData,file = "AllCounts_VoomNormalized.csv")
write.csv(Counts_CPM,file = "AllCounts_asCPM.csv")
write.csv(quantitysummaries, "Star_summaries.csv")


###############
#GO and KEGG annotation using the existing pipeline
###############

  
source("C:/Users/Gbolaga/Desktop/BRIC LED RNA RAW DATA/WySeq/GO_BP_functions_Functional (1).R")

decisions<-as.data.frame(strsplit(colnames(contMatrix),split = "_"))
decisions<-as.data.frame(t(decisions))

print("Let's get it started in ha")
for(contrasts in colnames(contMatrix)){
  print(paste("Beginning DE, GO, and Pathway Analysis for ", contrasts,".",sep = ""))
  assign("temp",topTable(DE,coef =paste(contrasts),adjust.method = "BH",p.value = 0.05,number=length(DE$coefficients)))
  if(exists("Annotations")){
    temp<-merge(Annotations,temp,by.x="ensembl_gene_id",by.y="row.names")
    temp$Row.names <- NULL
    if(length(temp$logFC)>0){
      assign(paste(contrasts,"_DE",sep = ""),temp)
      assign(paste(contrasts,"_GO",sep = ""),GO_Wrapper(temp,out.prefix=paste(contrasts),SpeciesInput=plant_species))
      assign(paste(contrasts,"_KEGG",sep = ""),PathwayD(temp,out.suffix = paste(contrasts),SpeciesInput = plant_species))
      
    }else{
      print("Bad exception case.")
    }
  }
  else{
    print("No annotation file available. Pathway and GO analysis will not be performed.")
    assign(paste(contrasts,"_DE",sep = ""),temp)
  }
  print(paste("Finished DE, GO, and Pathway Analysis for ", contrasts,".",sep = ""))
  print("Writing differential expression data to file.")
  write.csv(temp,file=paste(contrasts,"DE.csv",sep = ""))
  if(exists("dframe")){
    print("and again")
    deframe=data.frame(contrasts,sum(topTable(DE,coef =paste(contrasts),adjust.method = "BH",p.value = 1,number=length(DE$coefficients))[5]<0.05),sum(topTable(DE,coef =paste(contrasts),adjust.method = "BH",p.value = 1,number=length(DE$coefficients))[4]<0.05))
    names(deframe)=c("contrast","SigDEGenes_adj","SigDEGenes")
    dframe<-rbind(dframe,deframe)
  }
  if(!exists("dframe")){
    print("First time through")
    dframe=data.frame(contrasts,sum(topTable(DE,coef =paste(contrasts),adjust.method = "BH",p.value = 1,number=length(DE$coefficients))[5]<0.05),sum(topTable(DE,coef =paste(contrasts),adjust.method = "BH",p.value = 1,number=length(DE$coefficients))[4]<0.05))
    names(dframe)<-c("contrast","SigDEGenes_adj","SigDEGenes")
  }
  gc()
}  
  
  
  
  
