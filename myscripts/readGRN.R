setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../GeneOntologyData/dataToUse")

library(tidyverse)
require(R.matlab)
geneNetFile <- "athal-network_seidr.txt"
net <- read_delim(geneNetFile,quote="",delim ="\t") %>% 
  mutate(IRP=as.double(str_split_fixed(`irp_score;irp_rank`, ";", 2)[,1])) %>%
  select(Source, Target, Type, IRP)

genes <- unique(c(net$Source,net$Target))

genesAnnot <- read_delim("dataToUse/geneHash.txt",delim="\t")

genesnoAnnot <- genes[!genes%in%genesAnnot$gene]

maxIdx = max(genesAnnot$idx)

contiSeq <- seq(1+maxIdx,length(genesnoAnnot)+maxIdx) 

AllgenesHash <- tibble("gene"=genesnoAnnot,"idx"=contiSeq) %>%
  bind_rows(genesAnnot) %>% arrange(idx)


geneLinks.r <- net %>% inner_join(AllgenesHash, by=c("Source"="gene")) %>%
  inner_join(AllgenesHash,by=c("Target"="gene")) %>% select(idx.x,idx.y,IRP) %>%
  rename("Source"=idx.x,"Target"=idx.y)


# Matlab code to transform edge list to adj Matrix
# Example with weights just being 1
# edgelist = [1 2;2 3;2 4];
# edgelist = unique(edgelist, 'rows');
# sz = max(edgelist(:));
# A = sparse(edgelist(:,1), edgelist(:,2), 1, sz, sz);
# Example with weights in the third column
#edgelist = [1 2 0.1;2 3 0.2;2 4 0.3];
#edgelist = unique(edgelist, 'rows');
#sz = max(max(edgelist(:, 1:2)));
#A = sparse(edgelist(:,1), edgelist(:,2), edgelist(:,3), sz, sz);


paths = "dataToUse"
write_tsv(AllgenesHash,paste(paths,"AllgeneHash.txt",sep="/"))
write_tsv(geneLinks.r,paste(paths,"geneInt_3col.txt",sep="/"))
write_tsv(net%>%select(Source,Target,IRP),paste(paths,"geneStr_3col.txt",sep="/"))

