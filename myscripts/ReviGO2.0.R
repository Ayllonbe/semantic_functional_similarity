library(ontologyIndex)
library(tidyverse)
library(Rcpp)
library(fastcluster)


# Function to extract the depth of the Gene Ontology terms
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

# Function to extract the IC mazandu of the Gene Ontology terms
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
            for(int i=0;i<query.size();i++){
              String go1 = query.at(i);
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
              std::vector<int> v1 = mapGO2parentID.at(goID[i]);
              for(int j=i+1 ;j<query.size(); j++){
                std::vector<int> v2 = mapGO2parentID.at(goID[j]);
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


# Function to get the optimal number of cluster
SilClus <- function(hc.obj,dist.obj,nc){
  require(cluster)
  asw <-c() 
  for( k in 2 : nc){ 
    sil <- silhouette(cutree(hc.obj,k = k), dist.obj)
    asw <- c(asw,summary(sil)$avg.width)
  }
  return(which.max(asw)+1)
}


# Function to get the semantic similarity of GO terms
SemSim <- function(
  micaTable,
  go2info,
  semSim="Lin" # Lin or Mazandu 
){
  
  if(semSim=="Mazandu"){
    return( micaTable %>% inner_join(go2info, by=c("vec1"="GO_ID")) %>% select(-depth) %>% rename(IC1 = ICMaz) %>%
              inner_join(go2info, by=c("vec2"="GO_ID")) %>% select(-depth) %>%rename(IC2 = ICMaz) %>%
              inner_join(go2info, by=c("Common"="GO_ID")) %>% select(-depth) %>% rename(ICC = ICMaz) %>%  mutate("SemSim"=ICC/max(IC1,IC2)) %>%
              select(vec1,vec2,SemSim))
  }
  
  if(semSim == "Lin"){
    return( micaTable %>% inner_join(go2info, by=c("vec1"="GO_ID")) %>% select(-depth) %>% rename(IC1 = ICMaz) %>%
              inner_join(go2info, by=c("vec2"="GO_ID")) %>% select(-depth) %>%rename(IC2 = ICMaz) %>%
              inner_join(go2info, by=c("Common"="GO_ID")) %>% select(-depth) %>% rename(ICC = ICMaz) %>% mutate("SemSim"=(2*ICC/(IC1+IC2))) %>%
              select(vec1,vec2,SemSim))
  }
  
  
}


# Function to get the representative term
getRep <- function(dfQuery,
                   go.obo,
                   go2info,
                   semsim="Lin", # Lin or Mazandu
                   typeRep="AlgoRep" # AlgoRep or pval or padj
){ 
  query <- dfQuery$id
  t2ic <- as.list((go2info %>% select(GO_ID,ICMaz))$ICMaz)
  names(t2ic) = go2info$GO_ID
  subOntologies=c("GO:0008150","GO:0003674","GO:0005575")
  res <- tibble()
  for(subOnt in subOntologies){
    TermsQuery <- query[query %in% get_descendants(go.obo, subOnt)]
    if(length(TermsQuery)>0){
      t2par <- go.obo$ancestors[(tibble("id"=TermsQuery,"Par"=go.obo$ancestors[TermsQuery]) %>% unnest(Par))$Par]
      # To remove the obsolete terms
      go.toanal <- (go2info %>% filter(GO_ID %in% TermsQuery))$GO_ID 
      res.mica <- tibble(CA(go.toanal, go2info%>%filter(GO_ID %in% names(t2par)), "ICMaz", t2par))
      
      lin <- SemSim(res.mica, go2info, semsim)
      
      lin.inv <- lin %>% rename(vec1="vec2",vec2="vec1") %>% select(vec1,vec2,SemSim)
      lin <- bind_rows(lin,lin.inv)
      lm <- pivot_wider(lin, names_from = vec2, values_from = SemSim,values_fill = 0) 
      M = lm %>%select(lm$vec1) %>% as.matrix()
      rownames(M) = colnames(M)
      diag(M) =1 
      d <- as.dist(1-M)
      hc <- hclust(d,method = "average")
      ncl <- SilClus(hc,d, nrow(M)-1)
      clusters <- cutree(hc, k=ncl)
      clusters <- tibble("GO_ID" = names(clusters), "cl" = clusters) %>% arrange(cl)
      for(x in unique(clusters$cl)){
        cl.t<-(clusters %>% filter(cl==x))$GO_ID
        if(length(cl.t)>1){
          t2par <- go.obo$parents[(tibble("id"=cl.t,"Par"=go.obo$ancestors[cl.t]) %>% unnest(Par))$Par]
          t2child <- go.obo$children[(tibble("id"=cl.t,"Par"=go.obo$ancestors[cl.t]) %>% unnest(Par))$Par]
          ics <-t2ic[(tibble("id"=cl.t,"Par"=go.obo$ancestors[cl.t]) %>% unnest(Par))$Par]
          if(typeRep=="AlgoRep"){
            repT <- getRepresentative(t2par,t2child,ics, cl.t,subOnt)
          }else if(typeRep=="pval"){
            df    <- dfQuery %>% filter(id %in% cl.t)
            repT  <- cl.t[which.min(df$pval)]
          }else if(typeRep=="padj"){
            df    <- dfQuery %>% filter(id %in% cl.t)
            repT  <- cl.t[which.min(df$padj)]
          }
          if(length(repT)>1){
            res <- bind_rows(res, tibble("GO_ID" = cl.t,"SubOnt"=rep(subOnt, length(cl.t)) ,"Cluster"=rep(x, length(cl.t)), "GORep" = rep(paste(repT, collapse = "; "), length(cl.t)),"GORepName" = rep(paste(go.obo$name[repT], collapse="; "), length(cl.t))))
          }else{
            res <- bind_rows(res, tibble("GO_ID" = cl.t,"SubOnt"=rep(subOnt, length(cl.t)) ,"Cluster"=rep(x, length(cl.t)), "GORep" = rep(repT, length(cl.t)),"GORepName" = rep(go.obo$name[repT], length(cl.t))))
          }
          
          
        }
        else{
          res <-  bind_rows(res, tibble("GO_ID" = cl.t,"SubOnt"=subOnt,"Cluster"=x, "GORep" = cl.t,"GORepName" = go.obo$name[cl.t]))
        }
      }
    }
  }
  obso <- query[!query%in%res$GO_ID]
  res <- bind_rows(res, tibble("GO_ID" = obso,"SubOnt"=rep("obsolete", length(obso)) ,"Cluster"=rep(NA, length(obso)), "GORep" = rep(NA, length(obso)),"GORepName" = rep(NA, length(obso))))
  return(dfQuery %>% inner_join(res,by=c("id"="GO_ID")))
}


#Preparing the Gene Ontology data (I do not like GO.db - this process is more reproducible)
ontologyFile <- "http://release.geneontology.org/2021-02-01/ontology/go-basic.obo"
go.obo <-get_OBO(ontologyFile, propagate_relationships = "is_a", extract_tags = "minimal")
go2depth <- getGOTermsDepths(go.obo)
go2ICMaz <- getGOTermsICMaz(go.obo, go2depth)
go2info <- go2ICMaz %>% inner_join(go2depth, by=c("GO_ID"="id"))  %>% 
  arrange(depth) %>%
  select(GO_ID, depth, ICMaz)

sourceCpp("Representativo.cpp")

cl.t <- c("GO:0090304", "GO:0006259", "GO:0009112", "GO:0009116",
          "GO:0090305", "GO:0006396", "GO:0010501", "GO:0016071", "GO:0000726",
          "GO:0006284", "GO:0016073", "GO:0030490", "GO:0006431")
t2par <- go.obo$parents[(tibble("id"=cl.t,"Par"=go.obo$ancestors[cl.t]) %>% unnest(Par))$Par]
t2child <- go.obo$children[(tibble("id"=cl.t,"Par"=go.obo$ancestors[cl.t]) %>% unnest(Par))$Par]
t2ic <- as.list((go2info %>% select(GO_ID,ICMaz))$ICMaz)
names(t2ic) = go2info$GO_ID
ics <-t2ic[(tibble("id"=cl.t,"Par"=go.obo$ancestors[cl.t]) %>% unnest(Par))$Par]
subOnt = "GO:0008150"
paste(getRepresentative(t2par,t2child,ics, cl.t,subOnt),collapse = ";")


load("~/Downloads/enr_go.RData")

rep.results<-lapply(enr_res, function(query) {return(getRep(query$go, go.obo, go2info, "Lin", "AlgoRep"))})

save(rep.results,file="Downloads/enr_go_rep.RData")
