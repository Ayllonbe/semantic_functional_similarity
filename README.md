# semantic_functional_similarity

The semantic functional similarity project aimed to observe the impact of semantic similarity in multiple gene predictors that uses it. To do that we chose six gene predictors:

* ClusDCA [code](https://github.com/wangshenguiuc/clusDCA)
* HPHash [code](http://mlda.swu.edu.cn/upload/code/HPHash.zip)
* MASHUP [code](https://groups.csail.mit.edu/cb/mashup/mashup.tar.gz)
* NewGOA [code](http://mlda.swu.edu.cn/upload/code/NewGOA_code.zip)
* NMFGO [code](http://mlda.swu.edu.cn/upload/code/NMFGO_code.zip)
* PILL [code](http://mlda.swu.edu.cn/upload/code/PILL_Web.zip)

Most of these tools are developed in matlab, so we decided to prepare the data to run all of them in their matlab format. The only difference is the data was prepared using R. Since the data are matrix, we developed txt file that provide the gene network, GO network and GO diccionaries, gene hash, and semantic similarity matrix. The scripts to build these matrices are in `myscripts` folder. 

The project was organised in the following steps, splitted in multiple weeks: 

* Week 1: Learn Matlab
* Week 2: To adapt Gene Ontology and annotation with arabidopsis for 2018 and 2020
* Week 3: To include the metanetwork as input analysis
* Week 4: To modify the semantic similarity analysis to see the effect in the result.
* Week 5: To run the analysis with including the modified data
* Week 6-7: To develop a validation analysis (with existing approaches)
* Discuss the usefulness to publish it. 
* Apply this approach with Populus and Picea

This project was done until the **Week 4**. The **Week 5** must be the day that I modify the original Matlab code to includde the data that was created with R. 

I
