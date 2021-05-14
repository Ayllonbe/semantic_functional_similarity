
# This file work for the GOA analysis in NewGOA, HPHash, PILL, and NMFGO
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../GeneOntologyData/dataToUse")



library(stringr)
read_gaf <- function(filepath,
                     database = NULL, accession = "UNIPROT", 
                     filter.NOT = T, 
                     filter.evidence = c("ND", "IPI", "IEA", "NAS"),
                     ontology = NULL, propagate = T) {
  message("reading GOA file ", filepath, "...")
  # read GOA file 
  gaf.colnames <- c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", 
                    "GO_ID", "DB_Reference", "Evidence_Code", "With_From",
                    "Aspect", "DB_Object_Name", "DB_Object_Synonym",
                    "DB_Object_Type", "Taxon", "Date", "Assigned_By", 
                    "Annotation_Extension", "Gene_Product_Form_ID")
  goa <- suppressMessages(
    readr::read_tsv(filepath,comment = "!", col_names = gaf.colnames,
    col_types=readr::cols(
      DB = col_character(),
      DB_Object_ID = col_character(),
      DB_Object_Symbol = col_character(),
      Qualifier = col_character(),
      GO_ID = col_character(),
      DB_Reference = col_character(),
      Evidence_Code = col_character(),
      With_From = col_character(),
      Aspect = col_character(),
      DB_Object_Name = col_character(),
      DB_Object_Synonym = col_character(),
      DB_Object_Type = col_character(),
      Taxon = col_character(),
      Date = col_double(),
      Assigned_By = col_character(),
      Annotation_Extension = col_logical(),
      Gene_Product_Form_ID = col_logical()
    )))
  # optionally, filter annotations with the qualifier NOT
  if (filter.NOT) 
    goa <- goa[!grepl("NOT", goa$Qualifier),]
  # optionally, filter out evidence 
  if (!is.null(filter.evidence) & length(filter.evidence) > 0) 
    goa <- goa[!goa$Evidence_Code %in% filter.evidence,]
  # map to accession 
  if (!is.null(accession) & accession != "UNIPROT") {
    map <- suppressMessages(AnnotationDbi::select(
      database, keys = goa$DB_Object_ID, keytype = "UNIPROT", 
      columns = accession))
    goa[[accession]] <- map[[accession]][match(goa$DB_Object_ID, map$UNIPROT)]
    # remove NAs
    goa <- goa[!is.na(goa[[accession]]),]
  } else {
    goa$UNIPROT <- goa$DB_Object_ID
  }
  # read the ontology
  if (!is.null(ontology) & propagate) {
    if (!"ontology_index" %in% class(ontology))
      stop("Ontology must be of class ontology_index")
    goa$ancestors <- ontology$ancestors[goa$GO_ID]
    # filter out terms missing ancestors (deprecated)
    goa <- goa[lengths(goa$ancestors) > 0,]
    goa <- tidyr::unnest(goa, goa$ancestors)
    # replace column
    goa[["GO_ID"]] <- goa[["ancestors"]]
    goa <- goa[, -ncol(goa)]
  }
  return(goa)
}

transformToMat<-function(gaf){

  require(reshape2)
  require(tidyverse) 
  require(GO.db)
  require(org.At.tair.db)
  require(UniprotR)
  
GOA <- read_gaf(gaf, filter.evidence = c("ND", "IPI", "NAS"))

# Get only BP

GOA_bp <- GOA %>% filter(Aspect=="P")



GOA_2col <- unique(GOA_bp %>% dplyr::select(DB_Object_ID,GO_ID))
# Some genes (not much) are lost
mappingFile="ARATH_3702_idmapping_20-06.dat"
cat(paste("Mapping file used:",mappingFile))
mapping <- suppressMessages(read_tsv(mappingFile,comment = "!", col_names = FALSE) %>% dplyr::rename(ACC=X1, MAPPING=X2, ID=X3))

found <- mapping %>% filter(ID %in% GOA_2col$DB_Object_ID) %>% dplyr::select(ACC,ID) %>% dplyr::rename(Found = ID)

found.loc <- mapping %>% filter(MAPPING %in% "Araport") %>% filter(ACC %in% found$ACC) %>% dplyr::select(ACC,ID) %>% dplyr::rename(AGIp = ID) 

found.loc <- inner_join(found,found.loc,by=c("ACC"="ACC")) %>% dplyr::select(Found,AGIp)

found.uni <- mapping %>% filter(MAPPING %in% "Araport") %>% filter(ACC %in% GOA_2col$DB_Object_ID) %>% dplyr::select(ACC,ID) %>% dplyr::rename(Found=ACC,AGIp = ID)

found.ara <- bind_rows(found.loc,found.uni)

cat(paste(length(unique(GOA_2col$DB_Object_ID)) - length(unique(found.ara$Found)), "genes are lost since they have not Araport IDs..."))

GOA_2col_new <- inner_join(GOA_2col ,found.ara, by=c("DB_Object_ID"="Found")) %>% dplyr::select(AGIp,GO_ID)


GOA_cols <- GOA_2col_new %>% mutate("namespace"=Ontology(GO_ID))
#obsolete <- GOA_cols %>% filter(is.na(namespace))%>% dplyr::select(AGIp,GO_ID)
GOA_cols <- GOA_cols %>% filter(!is.na(namespace))
#' Preparing the relation between the ancestors and the GO terms annotating your genes
bp <- unique(GOA_cols[GOA_cols$namespace=="BP",2])$GO_ID
anc.bp <- as.list(GOBPANCESTOR)[bp]
anc.bp [sapply(anc.bp , is.null)] <- NULL
anc.bp <- tibble(melt(anc.bp, value.name = "ANC")) %>% filter(ANC!="all") %>% dplyr::rename(GO=L1)

nonDirect <- anc.bp %>% inner_join(GOA_cols,by=c("GO"="GO_ID")) %>% dplyr::select(AGIp,ANC) %>% dplyr::rename(GO_ID=ANC)

GOA_2col_new <-unique(GOA_cols %>% dplyr::select(AGIp,GO_ID)%>% bind_rows(nonDirect))


gos <- unique(GOA_2col_new$GO_ID)
genes <- unique(GOA_2col_new$AGIp) 

resT <- matrix(data=0,nrow = length(genes), ncol=length(gos))
rownames(resT) = genes
colnames(resT) = gos
m <- as.matrix(GOA_2col_new)
resT[m] = 1

dim(resT)

res <- list("mat"=resT,"col"=GOA_2col_new)
return(res)
}



gaf <- "tair.gaf"
folder1 <- "10.5281-zenodo.1442457_1st_Oct_2018/"
folder2 <- "10.5281-zenodo.4041891_10th_Sep_2020/"

his <- transformToMat(paste(folder1,gaf,sep = ""))
rec <- transformToMat(paste(folder2,gaf,sep = ""))


his.f <- his$mat[rownames(his$mat)%in%rownames(rec$mat),]
his.fo <- his.f[order(rownames(his.f)),]
rec.o <- rec$mat[order(rownames(rec$mat)),]

his.fo.f <-  his.fo[,colnames(his$mat)%in%colnames(rec$mat)]
his.fo.fo <-  his.fo.f[,order(colnames(his.fo.f))]
rec.o.o <- rec.o[,order(colnames(rec.o))]

ncol(his.fo.fo)
ncol(rec.o.o)



hBPidx  <- which(colnames(his.fo.fo)%in%colnames(rec.o.o))
rBPidx  <- which(colnames(rec.o.o)%in%colnames(his.fo.fo))

GOBPNames.his <- colnames(his.fo.fo)[hBPidx]
GOBPNames <- colnames(rec.o.o)
hproteinNames <- rownames(his.fo.fo)
rproteinNames <- rownames(rec.o.o)



####



parentLists <- as.list(GOBPPARENTS)[GOBPNames]

parentLists <-lapply(parentLists,function(x){return(tibble("parent"=x,"type"=names(x)))})

parent.df <- suppressMessages(tibble(melt(parentLists))) %>% filter(type=="isa") %>% 
  filter(parent!="all")  %>% dplyr::select(parent,L1,type)

geneHash <-  tibble("gene"=hproteinNames, "idx"=seq(1,length(hproteinNames))) 

geneHash.recent <- rproteinNames[!rproteinNames%in%hproteinNames]
geneHash.recent <- tibble("gene"=geneHash.recent, "idx"=seq(length(hproteinNames)+1,length(hproteinNames)+length(geneHash.recent)))
geneHash <- bind_rows(geneHash,geneHash.recent)


gohash <- tibble("go"=GOBPNames.his, "idx"=seq(1,length(GOBPNames.his))) 
gohash.recent <- GOBPNames[!GOBPNames%in%GOBPNames.his]
gohash.recent <- tibble("go"=gohash.recent, "idx"=seq(length(GOBPNames.his)+1,length(GOBPNames.his)+length(gohash.recent)))
gohash <- bind_rows(gohash,gohash.recent)
InterBP <- gohash$go %>% str_replace("GO:", "") %>% as.double()
goidxLinks.r <- parent.df %>% inner_join(gohash, by=c("parent"="go")) %>%
  inner_join(gohash,by=c("L1"="go")) %>% dplyr::select(idx.x,idx.y) %>%
  dplyr::rename("From"=idx.x,"To"=idx.y)

annotationIdx.h <- his$col %>% inner_join(gohash, by=c("GO_ID"="go")) %>%
  inner_join(geneHash, by=c("AGIp"="gene")) %>% dplyr::select(idx.y,idx.x) %>%
  dplyr::rename("gene"=idx.y,"go"=idx.x)

annotationIdx.r <- rec$col %>% inner_join(gohash, by=c("GO_ID"="go")) %>%
  inner_join(geneHash, by=c("AGIp"="gene")) %>% dplyr::select(idx.y,idx.x) %>%
  dplyr::rename("gene"=idx.y,"go"=idx.x)



####################### For Python #####################

annot <- rec$col %>% inner_join(gohash, by=c("GO_ID"="go"))  %>% dplyr::select(AGIp, GO_ID)


GOvec <- gohash$go

go2desc <- as.list(GOBPOFFSPRING)[GOvec]

go2desc <- tibble(GO=names(go2desc),DESC=go2desc) %>% unnest_longer(DESC) 

go2desc2prot <- unique(go2desc %>% inner_join(annot, by=c("DESC"="GO_ID")) %>%filter(!is.na(AGIp)) %>% dplyr::select(AGIp,GO) %>% rename(GO="GO_ID"))

moreSpecificAnnot <- dplyr::setdiff(annot,go2desc2prot) 

moreSpecificAnnot %>% count(GO_ID)

################## Write ###############

inputForFunSim <- rec$col %>% group_by(AGIp) %>%
  summarise(data =paste(cur_data()$GO_ID,collapse = ",")) %>% 
  rowwise()
inputForSemSim <- rec$col %>% group_by(GO_ID) %>%
  summarise(data =paste(cur_data()$AGIp,collapse = ",")) %>% 
  rowwise()


paths = "dataToUse"
if(!dir.exists(paths)){
  dir.create(paths)
}
write_tsv(inputForSemSim,paste(paths,"GOannotationForPython.txt",sep="/"))
write_tsv(inputForFunSim,paste(paths,"annotationForPython.txt",sep="/"))
write_tsv(parent.df,paste(paths,"go.txt",sep="/"))
write.table(InterBP,paste(paths,"InterBP.txt",sep="/"),
            row.names = FALSE, col.names = FALSE)
write_tsv(goidxLinks.r,paste(paths,"goIdx.txt",sep="/"))
write_tsv(gohash,paste(paths,"goHash.txt",sep="/"))
write_tsv(geneHash,paste(paths,"geneHash.txt",sep="/"))
write_tsv(annotationIdx.h,paste(paths,"annotationHistory.txt",sep="/"))
write_tsv(annotationIdx.r,paste(paths,"annotationRecent.txt",sep="/"))




#writeMat(paste(paths,"/",gaf,".mat",sep=''),rGOBP=rec.o.o,
#         hGOBP=his.fo.fo,hBPidx=hBPidx,rBPidx=rBPidx,
#         InterBP=InterBP,GOBPNames=GOBPNames,hproteinNames=hproteinNames,
#         rproteinNames=rproteinNames)



