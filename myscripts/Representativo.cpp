#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm> 
#include <set>
#include <iterator>
#include <cmath>
#include <map>
#include <deque> 

using namespace Rcpp;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

/*
 * REPRESENTATIVE CLASS
 */

class Representative
{
public:
  vector<string> getCombinedTerm();
  vector<string> getRepresentedTerms();
  Representative(vector<string> t, vector<string> rt);
  Representative();
  Representative(string t);
private:
  vector<string> combinedTerms;
  vector<string> representedTerms;
};

Representative::Representative(vector<string> t, vector<string> rt){
  combinedTerms = t;
  representedTerms=rt;
}
Representative::Representative(){
  
}
Representative::Representative(string t){
  combinedTerms.push_back(t);
  representedTerms.push_back(t);
  
}
vector<string> Representative::getCombinedTerm(){
  return combinedTerms;
}
vector<string> Representative::getRepresentedTerms(){
  return representedTerms;
}

/*
 * ALGOREP CLASS
 */
class AlgoRep
{
public:
  string getOntologystring();
  string getMethodHCL();
  long double getTailMin();
  long double getSemSimBT();
  long double getTermCoverage();
  vector<Representative> run(map<string, vector<string> >& t2par, map<string, vector<string> >& t2child, map<string, double >& t2ic, vector<string>& terms);
  AlgoRep();
  AlgoRep(string subOnt,long double tailmin,long double simRepFilter,long double coverage);
private:
  string subOnt;
  string mhcl;
  long double tailmin;
  long double simRepFilter;
  long double termcoverage;
  int geneSupport;
  void ANDfunction(vector<bool> boolSource, vector<bool>& boolTarget);
  void ANDNOTfunction(vector<bool> boolSource, vector<bool>& boolTarget);
  void ORfunction(vector<bool> boolSource, vector<bool>& boolTarget);
  bool SUMfunction(vector<bool>& boolS);
  Representative getRepresentative(map<string, vector<string> >& t2par, map<string, vector<string> >& t2child, map<string, double >& t2ic, vector<string>& terms);
  set<string> getManyRep(vector<string>& terms,  map<string, vector<string> >& t2par, map<string, vector<string> >& t2child, map<string, double>& t2ic, int ncombi);
  set<string> getOneRep(map<string, vector<string> >& t2child, set<string>& termSubGraph, string subOnt, int termsize, map<string,vector<bool> >& strBs);
  set<set<string> > doCombine(string& can, vector<string>& consideredC, int termsize, int limit, map<string,vector<bool> >& strBs);
  set<vector<string> > combination(vector<string> el, int k);
  string setCombinedTermInMap(set<string>& cand, map<string, vector<string> >& t2child, map<string, double>& t2ic, map<string,vector<bool> > &strBs);
};

AlgoRep::AlgoRep(){
  subOnt="GO:0008150";
  tailmin=0.1;
  simRepFilter=0.5;
  termcoverage=1.;
}

AlgoRep::AlgoRep(string ont, long double tm, long double sRF, long double cov){
  subOnt=ont;
  tailmin=tm;
  simRepFilter=sRF;
  termcoverage=cov;
}

string AlgoRep::getOntologystring(){
  return subOnt;
}
string AlgoRep::getMethodHCL(){
  return mhcl;
}
long double AlgoRep::getTailMin(){
  return tailmin;
}
long double AlgoRep::getSemSimBT(){
  return simRepFilter;
}
long double AlgoRep::getTermCoverage(){
  return termcoverage;
}

/*
 * ALGORITHM
 */
//No Mirar
vector<Representative> AlgoRep::run(map<string, vector<string> >& t2par, map<string, vector<string> >& t2child, map<string, double>& t2ic, vector<string>& terms){
  vector<Representative> res;
  Representative r;
  if(terms.size()>1){
    r= this->getRepresentative(t2par, t2child, t2ic, terms);
  } else{
    r = Representative(terms.at(0));
  }
  res.push_back(r);
  
  return res;
}



Representative AlgoRep::getRepresentative(map<string, vector<string> >& t2par, map<string, vector<string> >& t2child, map<string, double>& t2ic, vector<string>& terms){
  //long double tolerance = 0.7;
  set<string> genesInCluster;
  vector<long double> vecICs;
  vector<string> candidate;
  
  int ncombi = terms.size()>10?((int) floor(sqrt(fabs((long double)terms.size()/10.-1.)))+2):1 ; // numero que indica el numero de combinaciones deseadas
  ncombi = ncombi >3? 3:ncombi; // This is to delimit the number of combined parents we can reach
  set<string> repSet = getManyRep(terms, t2par, t2child,t2ic ,ncombi);
  for(string rep:repSet){
    candidate.push_back(rep);
    vecICs.push_back(t2ic[rep]);
  }
  if(candidate.size()>0){
    vector<long double>::iterator it = max_element(vecICs.begin(),vecICs.end());
    string termC = candidate.at(it-vecICs.begin());
    if(termC.find("_") != std::string::npos){
      vector<string> ct = t2child[termC];
      Representative r = Representative(ct,terms);
      return r;
    }else{
      vector<string> ct;
      ct.push_back(termC);
      Representative r = Representative(ct,terms);
      return r;
    }
  }else{
    return Representative();
  }
}

// No Mirar
set<string> AlgoRep::getManyRep(vector<string>& terms,  map<string, vector<string> >& t2par, map<string, vector<string> >& t2child, map<string, double>& t2ic, int ncombi){
  
  const size_t sT = terms.size();
  vector<bool> bs = vector<bool>(sT);
  deque<string> stack1;
  std::map<std::string,std::vector<bool> > strBs;
  for(int i = 0; i<terms.size();i++){
    string t = terms.at(i);
    stack1.push_back(t);
    strBs.insert(pair<string,vector<bool> >(t,bs));
    strBs[t][i] = true;
  }
  
  set<string> termSubGraph;
  while(stack1.size()>0){
    string t = stack1.at(0);
    termSubGraph.insert(t);
    vector<string> parents = t2par[t];
    for(int i=0; i < parents.size(); i++){
      string p = parents.at(i);
      if(strBs.find(p)!=strBs.end()){
        vector<bool> bsPivote = vector<bool>(sT);
        ORfunction(strBs[t],bsPivote);
        ANDNOTfunction(strBs[p],bsPivote);
        if(SUMfunction(bsPivote)){
          ORfunction(strBs[t],strBs[p]);
          stack1.push_back(p);
        }
      }else{
        strBs.insert(pair<string,vector<bool> >(p,bs));
        ORfunction(strBs[t],strBs[p]);
        stack1.push_back(p);
      }
    }
    stack1.erase(stack1.begin());   
  }
  
  set<string> candidates = getOneRep(t2child,termSubGraph, subOnt,terms.size(),strBs);
  int limit = 2;
  while(limit<=ncombi){
    set<set<string> > candidatesManyRep;
    for(string can:candidates){
      
      vector<string> consideredChildren;
      vector<string> childrens = t2child[can];
      for(int i=0; i<childrens.size();i++){
        string d = childrens.at(i);
        if(termSubGraph.find(d)!=termSubGraph.end()){
          
          vector<bool> b =  strBs[d];
          int mycount = std::count(b.begin(),b.end(),true);
          long double comparison = (long double)mycount/(long double)terms.size();
          if(comparison>=tailmin){
            consideredChildren.push_back(d);
          }
        }
      }
      set<set<string> > dC = doCombine(can,consideredChildren,terms.size(), limit,strBs);
      for(set<string> sc : dC){
        set<string> nsc;
        for(string s: sc){
          vector<bool> b = strBs[s];
          set<string> x = getOneRep(t2child,termSubGraph, s,std::count(strBs[s].begin(),strBs[s].end(),true),strBs);
          nsc.insert(x.begin(),x.end());
        }
        candidatesManyRep.insert(nsc);
      }
    }

    if(candidatesManyRep.size()>0){
      candidates.clear();
      vector<string> firstFilteredCandidates;
      vector<long double> ics;
      for(set<string> cand:candidatesManyRep){
        string t = setCombinedTermInMap(cand, t2child, t2ic, strBs); 
      
        firstFilteredCandidates.push_back(t);
        ics.push_back(t2ic[t]);
      }
      long double maxIC = ics.at(std::distance(ics.begin(),max_element(ics.begin(),ics.end())));

      for(int i=0;i<firstFilteredCandidates.size();i++){
        string ffc = firstFilteredCandidates.at(i);
        if(ics.at(i)==maxIC){
          candidates.insert(ffc);
        }
      }

    }
    limit++;
  }
  return candidates;
}


void AlgoRep::ORfunction(vector<bool> boolSource, vector<bool>& boolTarget){
  
  for(int i=0;i<boolSource.size();i++){
    boolTarget[i]= boolTarget[i] || boolSource[i];
  }
}

void AlgoRep::ANDfunction(vector<bool> boolSource, vector<bool>& boolTarget){
  
  for(int i=0;i<boolSource.size();i++){
    boolTarget[i]= boolTarget[i] && boolSource[i];
  }
}

void AlgoRep::ANDNOTfunction(vector<bool> boolSource, vector<bool>& boolTarget){
  
  for(int i=0;i<boolSource.size();i++){
    boolTarget[i]= boolTarget[i] && !boolSource[i];
  }
  
}

bool AlgoRep::SUMfunction(vector<bool>& vec){
  bool res;
  for (bool n : vec)
    res = res|| n;
  
  return res;
}


// No Mirar
string AlgoRep::setCombinedTermInMap(set<string>& cand, map<string, vector<string> >& t2child, map<string, double>& t2ic, map<string,vector<bool> > &strBs){
  string id;
  vector<bool> bs;
  long double ics = 0;
  
  vector<string> sv;
  for(string c:cand){
    sv.push_back(c);
    double ic = t2ic[c];
    ics+= ic;
    id+=c+"_";
    if(bs.size()>0){
      ORfunction(strBs[c],bs);
    }else{
      bs = vector<bool>(strBs[c].size());
    }
  }
  id = id.substr(0,id.size()-1);
  t2child.insert(pair<string, vector<string>>(id,sv));
  
  t2ic.insert(pair<string, double>(id,ics)); 
  strBs.insert(pair<string,vector<bool> >(id,bs));
  return id; 
}

//No Mirar
set<vector<string> > AlgoRep::combination(vector<string> el, int k){
  set<vector<string> > lint;
  int N = el.size();
  
  if(k>N){
    return lint;
  }
  vector<int> combination = vector<int>(k);
  int r=0;
  int index = 0;
  
  while(r>=0){
    if(index<=(N+(r-k))){
      combination[r] = index;
      if(r==k-1){
        vector<string> li;
        for(int c:combination){
          li.push_back(el.at(c));
        }
        lint.insert(li);
        index++;
      }else{
        index = combination[r]+1;
        r++;
      }
    }else{
      r--;
      if(r>0){
        index = combination[r]+1;
      }else{
        index = combination[0]+1;
      }
    }
  }
  
  return lint;
}

//No Mirar
set<set<string> > AlgoRep::doCombine(string& can, vector<string>& consideredC, int termsize, int limit, map<string,vector<bool> >& strBs){
  set<vector<string> > combiLS = combination(consideredC,limit);
  set<set<string> > combiSS;
  
  for(vector<string> ns: combiLS){
    int control =0;
    for(int i = 0;i<ns.size();i++){
      string t1 = ns.at(i);
      for(int j = i+1;j<ns.size();j++){
        string t2 = ns.at(j);
        vector<bool> inter = vector<bool>(termsize);
        ORfunction(strBs[t1],inter);
        ANDfunction(strBs[t2],inter);
        vector<bool> unionS = vector<bool>(termsize);
        ORfunction(strBs[t1],unionS);
        ORfunction(strBs[t2],unionS);
        long double jaccard = ((long double) std::count(unionS.begin(),unionS.end(),true)/(long double) std::count(inter.begin(),inter.end(),true));
        
        if(jaccard<simRepFilter){
          control++;
          break;
        }
        
      }
    }
    
    vector<bool> boolSet = vector<bool>(termsize);;
    if(control==0){
      for(int i = 0;i<ns.size();i++){
        string t1 = ns.at(i);
        ORfunction(strBs[t1],boolSet);
      }
    }
    
    if(control==0&&((long double)std::count(boolSet.begin(),boolSet.end(),true)/(long double)termsize>=termcoverage)){
      set<string> fs;
      fs.insert(ns.begin(),ns.end());
      combiSS.insert(fs);
    }
  }
  return combiSS;
}

//No Mirar
set<string> AlgoRep::getOneRep(map<string, vector<string> >& t2child, set<string>& termSubGraph, string subOnt, int termsize, map<string,vector<bool> >& strBs){
  deque<string> stack1;
  stack1.push_back(subOnt);
  set<string> candidates;
  while(stack1.size()>0){
    string t = stack1.at(0);
    int counter= 0;
    vector<string> childrens = t2child.at(t);
    for(int i=0;i<childrens.size();i++){
      string c = childrens.at(i);
      vector<bool> b =  strBs[c];
      int mycount = std::count(b.begin(),b.end(),true);
      long double comparison = (long double)mycount/(long double)termsize;
      //cout<<comparison<<" "<<termsize<<" "<<mycount<<endl;
      if(comparison>=termcoverage){
        stack1.push_back(c);
        counter++;
      }
    }
    if(counter==0){
      candidates.insert(t);
    }
    stack1.erase(stack1.begin());
  }
  return candidates;
}


// [[Rcpp::export]]
vector<string> getRepresentative( List& t2par, List& t2child, List& t2ic, CharacterVector& terms, string subOnt){
  
  double tailmin = 0.;
  double coverage = 1.;
  double simRepFilter= 0.5;
  
  vector<string> cluster_terms = as<vector<string>>(terms);
  map<string, vector<string> > term2parent;
  map<string, vector<string> > term2children;
  map<string, double > term2ic;
  vector<string> names = as<vector<string>>(t2ic.names());
  for(string n : names){
    term2parent.insert(pair<string,vector<string>>(n, t2par[n]));
    term2children.insert(pair<string,vector<string>>(n, t2child[n]));
    term2ic.insert(pair<string,double>(n, t2ic[n]));
  }
  
  AlgoRep ar = AlgoRep(subOnt,tailmin,simRepFilter,coverage);
  vector<Representative>  rep = ar.run(term2parent,term2children, term2ic, cluster_terms);
  
  vector<string> results;
  for(Representative r:rep){
    for(string t: r.getCombinedTerm()){
      results.push_back(t);
    }
  }
  return results;
  
}
