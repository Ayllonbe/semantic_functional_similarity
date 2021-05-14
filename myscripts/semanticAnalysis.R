library(ontologyIndex)
library(tidyverse)
library(Rcpp)

read_gaf <- function(filepath,
                     filter.NOT = T,
                     filter.evidence = c("ND", "IPI", "IEA", "NAS")) {
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
                      .default =  col_character()
                    )))
  # optionally, filter annotations with the qualifier NOT
  if (filter.NOT)
    goa <- goa %>% filter(!Qualifier%in% "NOT")
  # optionally, filter out evidence
  if (!is.null(filter.evidence) & length(filter.evidence) > 0)
    goa <- goa%>% filter(!Evidence_Code %in% filter.evidence)
  return(goa)
}



getGOTermsDepths <- function(go.obo, vec1=c("GO:0008150","GO:0003674","GO:0005575"), removeObs=TRUE){
  go2depth <- tibble("id" = go.obo$id, "depth" = rep(0, length(go.obo$id)))
  flag = 1
  level = 1
  while(flag == 1){
    vec2 <- c()
    go2depth <- go2depth %>% mutate(depth =replace(depth,id%in%vec1 ,level))
    vec2 <- c(vec2,unique(unlist(go.obo$children[vec1])))
    
    if(length(vec2) == 0){
      flag=0
    }
    vec1 <- vec2
    level=level+1
  }
  if(removeObs){
    go2depth <- go2depth %>% filter(depth!=0)
  }
  return(go2depth)
  
}


getGOTermsICMaz <- function(go.obo,go2depth,vec1=c("GO:0008150","GO:0003674","GO:0005575")){
  
  go2icMaz<- tibble("id" = vec1,
                    "alpha" = rep(1., length(vec1)),
                    "beta"  = rep(0., length(vec1)) )
  for(x in 2:max(go2depth$depth)){
    vec1 <- (go2depth %>% filter(depth==x))$id
    id2par <- go.obo$parents[vec1] %>% tibble("id"=vec1, "par" = .) %>% unnest(par)
    par.chlength <- lapply(go.obo$children[unique(id2par$par)],length) %>% tibble("par"=unique(id2par$par), "ch"=.) %>% unnest(ch)
    id2par <- id2par %>% inner_join(par.chlength, by=c("par"="par")) %>% inner_join(go2icMaz,by=c("par"="id"))
    id2par <- id2par %>% mutate("Division" = alpha/ch)
    id2par <- unique(id2par %>% mutate("betaC"=floor(log10(Division))) %>%
                       mutate("alpha" = Division/10^betaC)%>%
                       select(id, alpha, beta, betaC) %>%
                       group_by(id) %>%
                       mutate(beta=sum(beta)+sum(betaC)) %>%
                       mutate(alpha=prod(alpha)) %>% ungroup(id) %>%
                       select(id, alpha, beta))
    
    go2icMaz <- go2icMaz %>% bind_rows(id2par)
  }
  
  go2icMaz <- go2icMaz %>% mutate("ICMaz" = -log(alpha)-beta*log(10))
  
  return(go2icMaz %>% select(id, ICMaz) %>% rename("GO_ID"=id))
}



#' ### Compute the IC



getGOTermsICRes <- function(GOA, go.obo, subont=c("GO:0008150","GO:0003674","GO:0005575")){
  
  go.freq <- table(GOA$GO_ID)
  go.ic <- tibble()
  for(x in 1:length(subont)){
    go.ic <- bind_rows(go.ic, tibble("GO_ID"= names(go.freq[get_descendants(go.obo,subont[x])]),
                                     "ICRes"=-log(go.freq[get_descendants(go.obo,subont[x])]/go.freq[subont[x]])) %>%
                         filter(is.na(ICRes)==FALSE)) 
  }
  return(go.ic%>% mutate(ICRes=as.double(ICRes)))
}

getGOTermsICAgg <- function(go.ic, go.t2anc){
  go.sw <-  go.ic %>% mutate(SW=1./(1.+exp(1./IC)))%>% select(GO_ID, SW)
  go.aggIC <- unique(go.sw %>% inner_join(go.t2anc,by=c("GO_ID"="ANC")) %>%
                       rename( ANC = GO_ID, GO_ID = "GO_ID.y") %>%
                       filter(GO_ID %in% go.ic$GO_ID) %>%  arrange(GO_ID)%>% 
                       group_by(GO_ID) %>% mutate(aggIC = sum(SW)) %>% select(GO_ID,aggIC))
  
  return(go.aggIC %>% inner_join(go.sw, by=c("GO_ID"="GO_ID")))
  
}


ontologyFile <- "http://release.geneontology.org/2021-02-01/ontology/go-basic.obo"
go.obo <-get_OBO(ontologyFile, propagate_relationships = "is_a", extract_tags = "minimal")


versionGOA <- "http://release.geneontology.org/2021-02-01/annotations/tair.gaf.gz"

GOA <- read_gaf(versionGOA,
                filter.NOT = FALSE,
                filter.evidence = c("ND","NR"))
go.t2anc <-  dplyr::tibble("GO_ID"= go.obo$id, "ANC" = go.obo$ancestors) %>% unnest(ANC)
GOA <- GOA %>% inner_join(go.t2anc, by=c("GO_ID"="GO_ID")) %>%
  dplyr::select(DB_Object_ID, DB_Object_Symbol,Qualifier, ANC, Evidence_Code, DB_Object_Synonym) %>%
  dplyr::rename(GO_ID=ANC)

go2depth <- getGOTermsDepths(go.obo)
go2ICMaz <- getGOTermsICMaz(go.obo, go2depth)
go2ICRes <- getGOTermsICRes(GOA, go.obo)
go2ICagg <- getGOTermsICAgg(go2ICRes %>% rename(IC=ICRes),go.t2anc)

go2info <- go2ICRes %>% inner_join(go2depth, by=c("GO_ID"="id"))  %>%
  inner_join(go2ICMaz, by=c("GO_ID"="GO_ID"))  %>%
  inner_join(go2ICagg, by=c("GO_ID"="GO_ID")) %>%
  arrange(depth) %>%
  mutate(id=seq(1, nrow(go2ICRes))) %>%
  select(id, GO_ID, depth, ICRes, ICMaz, SW, aggIC)




# Until here is perfect

go.bp <- go2info %>% filter(GO_ID %in%  get_descendants(go.obo ,"GO:0008150"))

# CPP Function to get the Most Informative Common Ancestor (MICA)
# here type is the column of the feature that you want use to get MICA (depth or ICMaz)
cppFunction('DataFrame CA(CharacterVector& query, DataFrame& go2info, String& type, List& t2anc) {
            CharacterVector   goID = go2info["GO_ID"];
            NumericVector   score = go2info[type];
            std::map<String, int> mapGO2id;
            std::map<String, std::vector<int>> mapGO2parentID;
            for(int i=0;i<goID.size();i++){
              mapGO2id.insert(std::pair<String, int>(goID.at(i),i));
              mapGO2parentID.insert(std::pair<String, std::vector<int>>(goID.at(i),std::vector<int>()));
            }
            for(int i=0;i<goID.size();i++){
              String go1 = goID.at(i);
              CharacterVector chv = t2anc[go1];
              for(String anc: chv){
                mapGO2parentID.at(goID[i]).push_back(mapGO2id.at(anc));
              }
            }
            std::vector<std::string> ch1;
            std::vector<std::string> ch2;
            std::vector<std::string> ch3;
            std::vector<int> pivot;
            for(int i=0; i<query.size()-1; i++){
              std::vector<int> v1 = mapGO2parentID.at(query[i]);
              for(int j=i+1 ;j<query.size(); j++){
                std::vector<int> v2 = mapGO2parentID.at(query[j]);
                ch1.push_back(as<std::string>(query[i]));
                ch2.push_back(as<std::string>(query[j]));
                std::vector<int> v_intersection;
                std::set_intersection(v1.begin(), v1.end(),
                          v2.begin(), v2.end(),
                          std::back_inserter(v_intersection));
                int z = -1;
                double d = -1;
                for(int x : v_intersection){
                  double dep = score[x];
                  if(dep>d){
                  //std::cout<<x<<std::endl;
                  d = dep;
                  z = x;
                  }
                }
            if(z!=-1){
            ch3.push_back(as<std::string>(goID[z]));
            }else{
             ch3.push_back("root");
            }
              }
            }

            DataFrame df = DataFrame::create(Named("vec1") = ch1,
                                             Named("vec2") = ch2,
                                             Named("Common") = ch3);
            return df;
            }')


res.bp <- tibble(CA(go.bp$GO_ID, go.bp, "ICMaz", go.obo$ancestors))



cppFunction('DataFrame songSim(CharacterVector& query, DataFrame& go2info, List& t2anc) {
            CharacterVector   goID = go2info["GO_ID"];
            NumericVector   score = go2info["aggIC"];
            NumericVector   sw = go2info["SW"];
            std::map<String, int> mapGO2id;
            std::map<String, std::vector<int>> mapGO2parentID;
            for(int i=0;i<goID.size();i++){
              mapGO2id.insert(std::pair<String, int>(goID.at(i),i));
              mapGO2parentID.insert(std::pair<String, std::vector<int>>(goID.at(i),std::vector<int>()));
            }
            for(int i=0;i<goID.size();i++){
              String go1 = goID.at(i);
              CharacterVector chv = t2anc[go1];
              for(String anc: chv){
                mapGO2parentID.at(goID[i]).push_back(mapGO2id.at(anc));
              }
            }
            std::vector<std::string> ch1;
            std::vector<std::string> ch2;
            std::vector<double> semsim;
            std::vector<int> pivot;
            for(int i=0; i<query.size()-1; i++){
              std::vector<int> v1 = mapGO2parentID.at(query[i]);
              for(int j=i+1 ;j<query.size(); j++){
                std::vector<int> v2 = mapGO2parentID.at(query[j]);
                ch1.push_back(as<std::string>(query[i]));
                ch2.push_back(as<std::string>(query[j]));
                std::vector<int> v_intersection;
                std::set_intersection(v1.begin(), v1.end(),
                          v2.begin(), v2.end(),
                          std::back_inserter(v_intersection));
                double d = 0.;
                for(int x : v_intersection){
                  d = d+2.*sw[x];
                  }
                semsim.push_back(d/(score[i]+score[j]));
              }
            }

            DataFrame df = DataFrame::create(Named("vec1") = ch1,
                                             Named("vec2") = ch2,
                                             Named("AIC") = semsim);
            return df;
            }')

aic.bp <- tibble(songSim(go.bp$GO_ID, go.bp, go.obo$ancestors))
cppFunction('DataFrame DisjunctiveCA(CharacterVector& query, DataFrame& go2info, List& t2anc) {
            CharacterVector   goID = go2info["GO_ID"];
            std::map<String, int> mapGO2id;
            std::map<String, std::vector<int>> mapGO2parentID;
            for(int i=0;i<goID.size();i++){
              mapGO2id.insert(std::pair<String, int>(goID.at(i),i));
              mapGO2parentID.insert(std::pair<String, std::vector<int>>(goID.at(i),std::vector<int>()));
            }
            for(int i=0;i<goID.size();i++){
              String go1 = goID.at(i);
              CharacterVector chv = t2anc[go1];
              for(String anc: chv){
                mapGO2parentID.at(goID[i]).push_back(mapGO2id.at(anc));
              }
            }
            std::vector<std::string> ch1;
            std::vector<std::string> ch2;
            std::vector<std::string> ch3;
            for(int i=0; i<query.size()-1; i++){
              std::vector<int> v1 = mapGO2parentID.at(query[i]);
              for(int j=i+1 ;j<query.size(); j++){
                std::vector<int> v2 = mapGO2parentID.at(query[j]);
                
                std::vector<int> v_intersection;
                std::set_intersection(v1.begin(), v1.end(),
                          v2.begin(), v2.end(),
                          std::back_inserter(v_intersection));
        
                int n = v_intersection.size();
                for(auto i = 0; i < n-1; i++){
                 int ca1 = v_intersection.at(i);
                 for(auto j = i+1; j < n; j++){
                 int ca2 = v_intersection.at(j);
                  if(std::find(mapGO2parentID.at(goID[ca1]).begin(), mapGO2parentID.at(goID[ca1]).end(), ca2)!=mapGO2parentID.at(goID[ca1]).end()){
                    v_intersection.erase(v_intersection.begin()+j);
                    j--;
                    n = n-1;
                  }else if(std::find(mapGO2parentID.at(goID[ca2]).begin(), mapGO2parentID.at(goID[ca2]).end(), ca1)!=mapGO2parentID.at(goID[ca2]).end()){
                    v_intersection.erase(v_intersection.begin()+i);
                    i--;
                    n = n-1;
                    break;
                  }
                 }
                 }
                for(int x : v_intersection){
                  ch1.push_back(as<std::string>(query[i]));
                  ch2.push_back(as<std::string>(query[j]));
                  ch3.push_back(as<std::string>(goID[x]));
                }
            
              }
            }

            DataFrame df = DataFrame::create(Named("vec1") = ch1,
                                             Named("vec2") = ch2,
                                             Named("DisjCommon") = ch3);
            return df;
            }')
dca.bp <- tibble(DisjunctiveCA(go.bp$GO_ID, go.bp, go.obo$ancestors))

a <-dca.bp %>% filter(vec1== DisjCommon) %>% select(vec2, DisjCommon) %>% rename(id = vec2)
b<-dca.bp  %>% select(vec1, DisjCommon) %>% unique()%>% rename(id = vec1)

c<-dca.bp  %>% select(vec2, DisjCommon) %>% unique()%>% rename(id = vec2)


abc <- bind_rows(b,c) %>% unique()



library(igraph)

go.t2par <-  dplyr::tibble("GO_ID"= go.obo$id, "ANC" = go.obo$parents) %>% unnest(ANC) %>%
  inner_join(go2ICMaz, by=c("ANC"="GO_ID"))
graph <- graph_from_edgelist(as.matrix(go.t2par %>% filter(GO_ID%in%go.bp$GO_ID)%>% select(GO_ID, ANC)))  

graph <- graph %>% set_vertex_attr("Ki", value = (tibble(GO_ID = V(graph)$name) %>% inner_join(go2ICRes, by=c("GO_ID"="GO_ID")) %>% mutate(kis=if_else(GO_ID=="GO:0008150",0,1/ICRes)))$kis)

library(foreach)
library(doParallel)
# To extract the longest path (LP) from term to root (GO:0008150)
no_cores <- detectCores() - 2
cl <- makeCluster(no_cores)
n = length(V(graph))
LProot <- foreach(x =1:n, .combine = "rbind") %do% {
  require(igraph)
  require(dplyr)
  (tibble("id"=V(graph)$name[x], "P"=all_simple_paths(graph, V(graph)$name[x],"GO:0008150")))
}  %>%  mutate(len = map(P, length)) %>% 
  mutate(kis = map(P, function(x){V(graph)$Ki[x]}))  %>%
  unnest(len) %>% mutate(kis=map(kis,sum)) %>%
  unnest(kis) %>% select(id, len, kis)%>% group_by(id) %>%
  arrange(id) %>% filter(kis == max(kis)) %>% unique() %>%
  mutate(wLP=kis*len) %>% select(id, wLP) %>% bind_rows(tibble("id" = "GO:0008150", "wLP"=0))
stopCluster(cl)
# To extract the shortest path (SP) from term x to y
no_cores <- detectCores() - 2
cl <- makeCluster(no_cores)
n = length(abc$id)
SPs <- foreach(x = iter(abc, by='row'), .combine = "rbind") %dopar% {
  require(igraph)
  require(dplyr)
  (tibble("id"=x$id, "DisjCommon"=x$DisjCommon, "P"=all_simple_paths(graph, x$id,x$DisjCommon)))
} %>%  mutate(len = map(P, length)) %>% mutate(SP = map(P, function(x){V(graph)$Ki[x]})) %>%
  mutate(SP=map(SP,sum)) %>% unnest(SP)%>%  unnest(len)%>% 
  group_by(id, DisjCommon) %>%  filter(SP == min(SP)) %>%  
  unique() %>% ungroup()%>% mutate(wSP = SP*len) %>%select(id, DisjCommon, wSP)  %>%
  unique()  %>% bind_rows("id" = "GO:0008150", "wLP"=0)


stopCluster(cl)




