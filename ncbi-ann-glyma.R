# Reading in the CSV file
csv <- read.csv('data/FPKM_full-table.csv', header=T)
df <- data.frame()

# Function
geneSummary <- function(x){
  library(rentrez)
  
  # searches through NCBI gene db
  r_search <- entrez_search(db="gene", term=x) 
  
  # takes the UID of the first search input
  summ <- entrez_summary(db='gene', id=r_search$ids[1])
}

# Create a list of gene names from gene_name column
gene_list <- csv$gene_name

# Looping through each name in the gene list and creating a new dataframe
for(gene in gene_list){
  if(gene == '.'){
    temp <- data.frame('gene_name' = gene, 'glyma' = NA, 'desc' = NA)
  }else{
    summ <- geneSummary(gene)
    temp <- data.frame('gene_name' = gene, 'glyma' = summ$otheraliases, 'desc' = summ$description)
  }
  df <- rbind(df, temp)
}

# Create CSV file of annotations
write.csv(df, 'data/Full_GlymaAnnotations.csv')