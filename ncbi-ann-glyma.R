# Reading in the CSV file
# csv <- read.csv('data/FPKM_full-table.csv', header=T)

# Function
time <- function(){
  cat(format(Sys.time(), "%a %b %d %Y   %X"), "\n")
}

geneSummary <- function(x){
  library(rentrez)
  
  # searches through NCBI gene db
  r_search <- entrez_search(db="gene", term=x) 
  
  # takes the UID of the first search input
  summ <- entrez_summary(db='gene', id=r_search$ids[1])
}

annotate <- function(csv, output.name){
  # Create a list of gene names from gene_name column
  gene_list <- csv$gene_name
  
  # for large dataset
  n <- 1000
  df <- data.frame()
  full.df <- data.frame()
  
  for(j in 1:ceiling(nrow(csv)/n)) {
    nrange <- (n*(j-1)+1):(min(n*j, nrow(csv)))
    time()
    cat(paste0(j, "/", ceiling(nrow(csv)/n)))
    cat("\nrange: ", range(nrange))
    
    # Looping through each name in the gene list and creating a new dataframe    
    for(i in nrange){
      gene <- gene_list[i]
      
      # Progress Check
      if(i/max(nrange)==.25){
        time()
        print("25%")
      }else if(i/max(nrange)==.5){
        time()
        print("50%")
      }else if(i/max(nrange)==.75){
        time()
        print("75%")
      }else if(i/max(nrange)==1){
        time()
        print("100%")
      }
      
      if(gene == '.'){
        temp <- data.frame('uid' = NA, 'gene_name' = gene, 'alias' = NA, 'desc' = NA)
      }else{
        summ <- geneSummary(gene)
        temp <- data.frame('uid' = summ$uid, 'gene_name' = gene, 'alias' = summ$otheraliases, 'desc' = summ$description)
      }
      df <- rbind(df, temp)
    }
    full.df <- rbind(full.df, df)
    
    # Create CSV file of annotations
    write.csv(full.df, output.name, row.names=F)
  }
  return(full.df)
}

# Rename Glyma Terms
rename.glyma <- function(x){
  library(stringr)
  glyma.terms <- strsplit(na.omit(x['alias']),split=",")
  glyma.terms <- unlist(sapply(glyma.terms, head, 1))
  glyma.terms <- str_replace(glyma.terms, "GLYMA", "Glyma")
  glyma.terms <- glyma.terms[grepl("^Glyma", glyma.terms)] 
  return(glyma.terms)
}

# Annotate GO terms
glycine.db <- function(){
  # Database for Glycine Max annotations
  library(AnnotationHub)
  hub <- AnnotationHub()
  # query(hub, c("orgdb","glycine")) # to view the databases
  orgdb <- hub[['AH100839']]
  
  return(orgdb)
}

clean.list <- function(x){
  x[which(x == ".")] <- NA
  clean.list <- na.omit(x)
  return(clean.list)
}
annoGO <- function(x, orgdb){
  # get gene_name --> GO terms
  search.list <- x$gene_name
  search.list <- clean.list(search.list)
  
  # Annotate list of glyma terms
  # gene2cat <- select(orgdb, search.list, "ALIAS")
  
  # Collapse list of GO-terms to each ID
  df <- AnnotationDbi::select(orgdb, search.list, columns(orgdb), "ALIAS")
  collapsed.df <- df %>% group_by(ENTREZID, GENENAME) %>% 
    summarise(GO = paste(GO, collapse=','))
  colnames(collapsed.df)[1] <- 'uid'
  
  return(collapsed.df)
}

GO.profile <- function(x, orgdb){
  library(clusterProfiler)
  
  # get gene_name --> GO terms
  search.list <- x$gene_name
  search.list <- clean.list(search.list)
  
  # creating dataframe of each gene with list of annotations
  df <- AnnotationDbi::select(orgdb, search.list, columns(orgdb), "ALIAS")
  
  # Entrez gene ID
  yy <- groupGO(unique(df$ENTREZID), orgdb, ont="MF", level=2)
  yy <- as.data.frame(yy[order(-yy$Count),])
  yy <- subset(yy, Count>5)
  
  ego <- enrichGO(df$ENTREZID, orgdb, ont="MF", pAdjustMethod = "none")
  # ego.df <- as.data.frame(ego)
  
  return(list(yy = yy, ego = ego))
}
