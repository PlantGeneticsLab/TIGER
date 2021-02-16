<img src="https://contattafiles.s3.us-west-1.amazonaws.com/tnt22006/DhRSlDOsdlFY6WL/tiger.png" height=100 align="center"> 

# Toolkits Integrated for Genetic and Evolutionary Research
*Code packages and apps to conduct efficient genetic and evolutionary analysis*

***
## Overview
TIGER integrates a collection of toolkits facilitating research analysis on genetics and evolution. It is a Java based program, combining efficiency and flexibility to realize high-performance computation of various types of analysis. 

***
## Usage
TIGER is designed to simpliy its usage for the convinience of end users. It has a list of applications (Apps). Each app has a parameter file, where app parameters can be specified in the parameter file. In parameter files, lines staring with "#" can be changed, lines starting with "@" cannot be changed, lines starting with neither "#" nor "@" are user defined papameters. Users can run apps in TIGER by selecting an app and specifify its parameters. TIGER has only two options,  
* -a, the name of the selected app  
* -p, the path of the parameter file 


Take running app FastCall as an example, the commnand line is,  
* java -Xmx100g -jar TIGER.jar -a FastCall -p ./parameter_fastcall.txt > log.txt &

TIGER will be developed in a long-term basis. Detailed usage of apps can be seen below.

***
## User Manual of Apps

* [FastCall](https://github.com/PlantGeneticsLab/TIGER/wiki/FastCall) -Superfast variant calling and genotyping for whole-genome shotgun (WGS) sequencing data.
* [HapScanner](https://github.com/PlantGeneticsLab/TIGER/wiki/HapScanner) -Superfast genotyper for whole-genome shotgun (WGS) sequencing data, based on an existing genetic variation library.
* [PopDep](https://github.com/PlantGeneticsLab/TIGER/wiki/PopDep) -Read depth profiler for a population, based on whole-genome shotgun (WGS) sequencing data.

***
## Developer Manual
* [TIGER Developer Manual](https://docs.google.com/document/d/1BU99b3joz0yItsJi2VabWbl6EyYfJmbo2oGUybl4PoM/edit?usp=sharing)

***
## Releases
* [v1.0.0](https://github.com/PlantGeneticsLab/TIGER/releases/tag/V1.0.0) 03/31/2020
* [v1.0.1](https://github.com/PlantGeneticsLab/TIGER/releases/tag/v1.0.1) 05/15/2020
* [v1.0.2](https://github.com/PlantGeneticsLab/TIGER/releases/tag/v1.0.2) 11/16/2020
* [v1.0.3](https://github.com/PlantGeneticsLab/TIGER/releases/tag/v1.0.3) 12/24/2020

***
## Citations

* FastCall  
[Punna Ramu, Williams Esuma, Robert Kawuki, Ismail Y Rabbi, Chiedozie Egesi, Jessen V Bredeson, Rebecca S Bart, Janu Verma, Edward S Buckler, Fei Lu. (2017). Cassava haplotype map highlights fixation of deleterious mutations during clonal propagation. Nature Genetics. doi: 10.1038/ng.3845!](https://www.nature.com/articles/ng.3845)
* HapScanner   
[Y. Zhou, X. Zhao, Y. Li, J. Xu, A. Bi, L. Kang, D. Xu, H. Chen, Y. Wang, Y. Wang, S. Liu, C. Jiao, H. Lu, J. Wang, C. Yin, Y. Jiao, F. Lu. (2020). Triticum population sequencing provides insights into wheat adaptation. Nature Genetics. doi:10.1038/s41588-020-00722-w](https://www.nature.com/articles/s41588-020-00722-w)
