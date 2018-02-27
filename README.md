
			 ____  ____  ____  _  _____ _____
			/ ___\/   _\/  _ \/ \/    //    /
			|    \|  /  | | \|| ||  __\|  __\
			\___ ||  \__| |_/|| || |   | |   
			\____/\____/\____/\_/\_/   \_/    

[![Build Status](https://travis-ci.org/phoenixding/scdiff.svg?branch=master)](https://travis-ci.org/phoenixding/scdiff)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


# INTRODUCTION 


SCDIFF is a python package designed to analyze the cell differentiation process 
using time-series single cell RNA-seq data. It can be used to predict the
transcription factors, which regulate the cell differentiation process and the differential genes along the differentiation process. It
also visualizes the differentiation process using a graph, in which nodes
represent different sub-population cells and edges denote the 
differentiation paths. 

![flowchart](./scdiff/img/FlowChart1.jpg)

# PREREQUISITES

* python (python 2 and python 3 are both supported)  
It was installed by default for most Linux distribution and MAC.  
If not, please check [https://www.python.org/downloads/](https://www.python.org/downloads/) for installation 
instructions. 

* Python packages dependencies:  
-- scikit-learn   
-- scipy  
-- numpy  
The python setup.py script will try to install these packages automatically.
However, please install them manually if, by any reason, the automatic 
installation fails. 

# INSTALLATION
 
cd to the package root directory 

```shell
$cd scdiff
```
run python setup to install   

```shell
$python setup.py install
```
    
Linux, Mac:

```shell  
$sudo python setup.py install 
```

or install via pip (should work in most operating systems):  

```shell
$sudo pip install --upgrade https://github.com/phoenixding/scdiff/zipball/master
```
# USAGE

```shell
usage: scdiff [-h] -i INPUT -t TF_DNA -k CLUSTERS -o OUTPUT [-s SPEEDUP] [-l LargeType] 
                   [-d DSYNC][-a VIRTUALANCESTOR]

  -h, --help            show this help message and exit
  
  -i INPUT, --input INPUT, required 
                        input single cell RNA-seq expression data
                        
  -t TF_DNA, --tf_dna TF_DNA, required
                        TF-DNA interactions used in the analysis
                        
  -k CLUSTERS, --clusters CLUSTERS, required
                        how to learn the number of clusters for each time
                        point? user-defined or auto? if user-defined, please
                        specify the configuration file path. If set as "auto"
                        scdiff will learn the parameters automatically.
                        
  -o OUTPUT, --output OUTPUT, required
                        output folder to store all results
                        
  -s SPEEDUP, --speedup SPEEDUP(1/None), optional
                        If set as 'True' or '1', SCIDFF will speedup the running
                        by reducing the iteration times.
                        (use less stringent convergence criteria).
                        
  -l LARGETYPE,  --largetype LARGETYPE (1/None), optional
						if specified as 'True' or '1', scdiff will use LargeType mode to 
						improve the running efficiency (both memory and time). 
						As spectral clustering is not scalable to large data,
						PCA+K-Means clustering was used instead. The running speed is improved 
						significantly but the performance is relatively worse. If there are
						more than 2k cells at each time point on average, it is highly 
						recommended to use this parameter to improve time and memory efficiency.
						
						
  -d DSYNC,  --dsync DSYNC (1/None), optional
						If specified as 'True' or '1', the cell synchronization will be disabled. 
						If the users believe that cells at the same time point are similar in terms of 
						differentiation/development. The synchronization can be disabled.
  
  -a VIRTUALANCESTOR, --virtualAncestor VIRTUALANCESTOR (1/None), optional
						scdiff requires a 'Ancestor' node (the starting node, 
						all other nodes are descendants).  By default, 
						the 'Ancestor' node is set as the first time point. The hypothesis behind is :  
						The cells at first time points are not differentiated yet
						( or at the very early stage of differentiation and thus no clear sub-groups, 
						all Cells at the first time point belong to the same cluster).  
						  
						If it is not the case, users can set -a as 'True' or '1' to enable
						a virtual ancestor before the first time point.  The expression of the 
						virtual ancestor is the median expression of all cells at first time point.   
						
						
```
# INPUTS AND PRE-PROCESSING

scdiff takes the two required input files (-i/--input and -t/--tf_dna) and one optional config file (-k/--clusters).

* __-i/--input__  
This specifies the single cell RNA-Seq expression data.  
If the RNA-Seq data is not processed, the instruction about how to calculate expression based on RNA-Seq raw reads can be found in many other studies, e.g (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4728800/).
For example, users can use Tophat + Cufflink to calculate the gene expression in terms of FPKM.  Please refer to corresponding tools for instructions. 
Once we get the RNA-Seq gene expression, the expression data should be transformed to log space for example by log2(x+1) where x could represent the gene expression in terms of RPKM, FPKM or TPM depending
on what tools are used to precoess the RNA-Seq expression data.  
The input file has the following formatting requirements:
	* __Header Row__  
	First 3 columns are "Cells","Time","Label" and the remaining columns are gene names.   
	* __Data Rows__  
		* __1st column__: Cell ID, represents the ID for the cell.
		* __2nd column__: Cell time, Integer, represents the measurement time of the cell. 
		* __3rd column__: Cell label, represents the label of the cell (if known). In most cases, we don't have any prior knowledge of the cell label. In this case, use "NA" instead.
		Please note that the cell label is only used for analyzing the result, we did not use this information in scdiff method. 
		* __4th- columns__: Gene expression values.  
	
	Example input:     
	[example.E](example/example.E)

* __-t/--tf_dna__  
This specifies the TF-gene interaction data.  In other words, it specifies the TF targets. 
Under the tf_dna directory, we provided a human TF-gene interaction file inferred using the strategy in our previous study (https://www.ncbi.nlm.nih.gov/pubmed/20219943). 
Although this TF-gene interaction file is collected in human, it should be also able to apply to other close species such as mouse.   
Besides, in our previous work DREM (http://sb.cs.cmu.edu/drem/), we did collected the TF-gene interactions for common species including human, mouse, fry, E.coli, yeast, Arabidopsis. 
Please refer to  http://sb.cs.cmu.edu/drem/DREMmanual.pdf appendix B for complete details. 
Those TF-gene interaction files can be downloaded from our DREM software (https://github.com/phoenixding/idrem/tree/master/TFInput).
You might need to unzip and re-format the file to satisfy the requirements. The TF-gene interaction file has the following formatting requirements:  
 
	* __Header Row__  
	```
	TF	Gene	Input
	```
	* __Data Rows__  
		* __1st column__: TF ID (gene symbol)
		* __2rd column__: gene ID (gene symbol)
		* __3rd column__: Input, optional, the interaction strength between TF and target gene. If missing, by default it's 1.  
		This column is not used in scdiff. 
		 	
	Example file:   
	[example TF gene interaction file](tf_dna/human_predicted_100.txt.update)

* __-k/--cluster__  
  This specifies the clustering parameter (String).   
  It's need to be either 'auto' or path to the 'config' file. Here, 'auto' denotes the clustering parameters will be learned automatically.  
  The path to the 'config' file specifies the file with customized initial clustering parameters. This should be used when users have some prior knowledge.  
  For example, if we know there are how many sub-populations within each time, 
  we can just directly specify the clustering parameters (K, # of clusters) using the 'config' file.
  Please note that this config file is optional, users can set -k as "auto" and the scdiff will learn the clustering parameters automatically. 
  Config file format:   
 
  * __1st column__: Time point
  * __2nd column__: # of clusters(sub cell-types) at this time point.  
  ```
  14  1  
  16  2  
  18  5
  ``` 
  please note if using the customized config file for clustering, the first line (represents the first time point) second column (# of sub cell-types) must be 1.
  The first line denotes the ancestor time point and all cells at later time points are descendants. However, if the first time point is already differentiated and there
  are already multiple sub-types, users need to use VirtualAncestor option (-a 1) to generate a virtual ancestor node (virtual root node). In the config file,
  a virtual ancestor line also needs to be added before the first time point. The virtual ancestor time point should be FirstTimePoint(Integer)-1.  For example,
  if the first time point is 14, the virtual ancestor time point would be 13.   
  Example  config file: 
  [example.config](example/example.config)
  
For the scdiff options, please refer to the [usage](#usage) section. 

# RESULTS AND VISUALIZATION

The results are given under the specified directory. 
The predicted model was provided as a json file, which is visualized by the
provided JavaScript. Please use *Chrome/FireFox/Safari* browser for best experience.   

![example_out_fig](./images/example_out.png)

The following is the manual for the visualization page. 

**Visualization Config (Left panel)**:
  
* **RESET**: It restores all configs to default.   
The slider resizes the visualization. 

* **CELL PLOTS** :
	* **Plot Cells**: show PCA/TSNE plots for all cells (use the radio button to select the dimension reduction method : PCA/T-SNE for visualization).

* **TF CONFIG** :  
	* **Show/Hide TF** : display/hide the regulating TFs for each path.  
	* **Explore  TF** : highlight the regulating paths for given TF; show the expression profile of the input TF and the average fold change of the TF targets.  

* **GENE CONFIG**:  
	* **Show/Hide DE **: display/hide the differentially expressed genes for each edge.  
	* **Explore DE gene** : highlight the paths with those DE genes and also show the expression profile of the input gene.     
	* **Find DE genes between** : find the differentially expressed genes between two specified nodes. Use the dropdown menu to specify two nodes for comparison. 

* **Download**: 
	* **Generate Figure**: download the visualization figure.
	* **Generate Json download file**: download the json file, which contains all information of the model.   
		
		```
		Json format:   
		{
			"GeneList": [List of gene symbols],
			"CellList": [                       // list of all cells
				{	
					// a cell instance
					"TE": [x0,y0], // cell coordinates in 2D space projected by T-SNE 
					"PE": [x1,y1], // cell coordinates in 2D space projected by PCA
					"typeLabel": cell Label, // label of the cell if known
					"ID": cell ID,
					"T": cell measure time
					
				},
				...
			],
			
			"NodeList":[                      // list of all nodes in Graph
				{	
					// a node instance
					"E": [list of gene expression in node, please refer to GeneList for gene symbols],
					"parent": parent node index,
					"children": list of children node indexes,
					"CELL": list of cell indexes, // all the cells assigned to the node
					"T": node time //estimated node level
					"ID": node ID 
				},
				...
			],
			"EdgeList":[                   // list of all edges in Graph
				{
					// an edge instance
					"to": edge target node index,
					"from": edge source node index,
					"de": List of differential genes along the edge,
					"etf": List of regulating TFs for the edge
					
				},
				...
			]
			
		}
		```  
	* **Generate TF download File** : download regulating TFs for all paths.
	* **Generate DE download file**: download DE genes for all paths. 

**Visualization Canvas (Right Panel)**:

* mouse over each node to show the pie chart of cell types within each node. 
* left click each node to show:
	* Cell IDs (with labels) within the node. 
	* Regulating TFs for the edge ending at the node.
	* DE genes for the edge ending at the node.
	* Regulating TFs for the whole path ending at the node. 
	
# EXAMPLES
 
Run scdiff on given time-series single cell RNA-seq data.  
An example script [exampleRun.py](example/exampleRun.py) is provided under the example directory.
  
**1) Run with automatic config**

```shell
$ scdiff -i example.E -t example.tf_dna -k auto -o example_out
```
* **-i/--input**:   
**example.E** is the single cell RNA-seq dataset with following format (tab delimited)

	```
	cell_ID	time	cell_label	ex_gene1	ex_gene2	...	ex_geneN
	```

	* cell_ID: ID of the cell.
	* time: measure time of the cell RNA-seq.
	* cell_lable: known label for the cell (if available) ,if not , denote as NA.
	* ex_genei: expression of gene i (log2 gene expression value). Gene expression can be FPKM or RPM or any other acceptable gene expression measurements.   
	
	Please read **example.E** for an example of acceptable time-series single cell dataset format.  

* **-t/--tfdna**:  
**example.tf_dna** provides the TF-DNA interaction information used in this study (TF target inforamtion) with tab delimited format.
	
	```
	TF	TF_target	Score
	```
	For example: 
	
	```
	ZNF238	TIMM8B	0.9
	SOX9	TIMM8B	0.9
	ZEB1	TIMM8B	0.9
	GATA4	TIMM8B	0.9
	CEBPA	RAB30	0.9
	NKX2-1	RAB30	0.9
	SRF	RAB30	0.9
	SOX5	RAB30	0.9
	SRY	RAB30	0.9
	POU1F1	RAB30	0.9
	POU2F1	RAB30	0.9
	NFKB1	KRI1	0.9
	E2F1	C11ORF35	0.9
	DSP	C11ORF35	0.9
	ELSPBP1	C11ORF35	0.9
	EGR2	C11ORF35	0.9
	EGR1	C11ORF35	0.9
	NR2F2	C11ORF35	0.9
	LMO2	C11ORF35	0.9
	ESR2	C11ORF35	0.9
	HNF1A	C11ORF35	0.9
	EGR3	C11ORF35	0.9
	```
The TF-DNA directory provides the TF-DNA interaction file used in this study. 
* **-k/--clusters**:  
This specifies the clustering parameter (String). It's need to be either 'auto' or  path to the 'config' file.
Here, 'auto' denotes the clustering parameters will be learned automatically. 

* **-o/--output**:  
**example_out** is the specified output directory. 


**2) Run with user-defined config** 

```shell
$scdiff -i example.E  -t example.tf_dna -k example.config -o example_out
```

The format of example.E and example.tf_dna are the same as described above. 

**example.config** specifies the custom initial clustering parameters. This was used when we have some prior knowledge.
For example, if we know they are how many sub-populations within each time, we can just directly specify the clustering parameters using
the example.config file, which provides better performance. 

example.config format(tab delimited)

```
time	#_of_clusters
```
For example:  

```
14  1  
16  2  
18  5  
```
However, if we don't have any prior knowledge about the sub-populations within each time point. We will just use the automatic initial clustering. 
:-k auto.

**3) Run scdiff on large single cell dataset** 


```shell
$scdiff -i example.E -t example.tf_dna -k auto -o example_out -l True -s True
```

-i, -t, -k, -o parameters were discussed above. 

* **-l/--large (optional)**  
String, if specified as 'True' or '1', scdiff will use LargeType mode to improve the running efficiency (both memory and time).
The performance will be sacrificed to some extent. K-means will be used for Clustering instead of Spectral clustering. 

* **-s/--speedup (optional)**  
Speedup the convergence, it will reduce the running time significantly. This is highly recommended for large dataset. 
Based on testing on lung single cell dataset (Treutlein 2014), the speedup performance is just slightly worse (2 more cells were miss-assigned )



**(4) Run scdiff on large single cell dataset with synchronization disabled and virtual ancestor**
```shell 
$scdiff -i example.E -t example.tf_dna -k auto -o example_out -l True -s True -d True -a True
```  
-i, -t , -k, -o, -l ,-s parameters were defined above. 

* **-d/--dsync (optional)**  
If set as 'True' or '1', the cell synchronization will be disabled.  By default, the cell synchronization is enabled. 
For large dataset, the users can disable the synchronization to speedup. If the authors have prior knowledge, the synchronization of cells are 
relatively good, users can also disable the synchronization.  

* **-a/--virtualAncestor (optional)**  
If set as 'True' or '1', the virtual ancestor node will be built. By default, the ancestor node is the first time point (all cells at the first time point).

**5) example running result**  

The following link present the results for an example running.   
[example_out](http://www.cs.cmu.edu/~jund/scdiff/result/treutlein2014_lung/treutlein2014_lung.html)


# MODULES  & FUNCTIONS

## scdiff module
This python module is used to perform the single cell differentiation analysis and it builds a graph (differentiation). Users can use the modules by
importing scdiff package in their program.  Besides the description below, we also provided a module testing example inside the example directory under the name [moduleTestExample.py](example/moduleTestExample.py). 

**[scdiff.Cell(Cell_ID, TimePoint, Expression,typeLabel,GeneList)](#cell)<a id="cell"></a>**    
This class defines the cell. 

**Parameters**:  

* **Cell_ID**: String  
The ID  of the cell.
* **TimePoint**: Integer  
Measurement TimePoint of the cell, Integer.  
* **Expression**: List of float  
Expression of all genes.
* **Cell_Label**: String  
The true label for the cell if available, 'NA' if not available. (Note, we don't need  this information for the  model, but it's  useful when
analyzing the result).
* **GeneList** : List of String  
List of gene symbols expressed in the cell.  

**Output**:  
A Cell class instance  (with all information regarding to  a cell)

**Attributes**:
* **ID **: String  
Cell ID  
* **T**: Integer  
Cell Time  
* **GL**: List of String  
List of gene names  
* **E** : List of float
List of gene expression
 
**Example**:

```python
import scdiff
from scdiff.scdiff import *

# reading example cells ...
AllCells=[]
print("reading cells...")
with open("example.E","r") as f:
	line_ct=0
	for line in f:
		if line_ct==0:
			GL=line.strip().split("\t")[3:]
		else:
			line=line.strip().split("\t")
			iid=line[0]
			ti=float(line[1])
			li=line[2]
			ei=[round(float(item),2) for item in line[3:]]
			ci=scdiff.Cell(iid,ti,ei,li,GL)
			AllCells.append(ci)
		line_ct+=1
		print('cell:'+str(line_ct))
```


**[scdiff.Graph(Cells, tfdna, kc, largeType=None, dsync=None, virtualAncestor=None)](#graph) <a id="graph"></a>**  
This class defines the differentiation graph. 

**Parameters**:  

* **Cells**:  Cell instances     
Please read [cell](#cell) Class definition for details.   
* **tfdna**: String  
It specifies the path to the TF-gene interaction file. 
* **kc**: String   
clustering config. It's a string with either 'auto' or clustering configure file path (-k parameter).  
* **largeType**: None(default) or String     
whether the single cell data is a 'largeType' (largeType denotes the number of cells in datasets is very large  (typically >2k cells). 
In such case, the performance will be scarified to improve the running efficiency (e.g. using K-means instead of spectral clustering). 
If not set (**None**), the dataset will be regarded as normal, if set as 'True', the dataset will be treated as largeType. 
* **dsync**: None(default) or String ('True' or '1')  
whether disable the cell time synchronization. By default, the cell time synchronization is enabled. If dsync set as "1" or "True",
this function will be disabled. No cell time synchronization will be made.  

* **virtualAncestor**: None (default) or String ('True' or '1')   
By default, all cells at the first time will be regarded as the starting ancestor for all cells at later time points. 
If users believe that the cells are already differentiated significantly and there are already more than 1 group at the first time point.  
Then, a virtual ancestor needs to be used by setting virtualAncestor as "True" or "1" . 

**Output**:  
A graph instance with all nodes and edges, which represents the differentiation structure for given inputs. 

**Attributes**:
* **Cells**: List of Cell instances  

* **Nodes**: List of Cluster instances (each cluster represents a node), all nodes in the graph. 

* **Edges**: List of Path instances (each represents an edge), all edges in the graph.   

* **dTD,dTG,dMb**:     
They are all dictionaries about TF-gene interactions.   
dTD-> key: TF, value: List of target genes  
dTG-> key: gene, value: List of regulating TFs    
dMb-> key: TF, value: List of target genes, which are expressing (non-zero) in the given single cell expression dataset.   

**Example**:

```python 
import scdiff
from scdiff.scdiff import *

print("testing scdiff.Graph module ...")
# creating graph using scdiff.Graph module and examples cells build above
g1=scdiff.Graph(AllCells,"example.tf_dna",'auto')
```

**[scdiff.Clustering(Cells, kc,largeType=None)](#graph)**   
This class represents the clustering.  

**Parameters**:

* **Cells**:  List of Cell  
Please read [Cell](#cell) Class definition for details. 
* **kc**: String   
clustering config. It's a string with either 'auto' or clustering configure file path (-k parameter).
* **largeType**: None(default) or String   
whether the single cell is a 'largeType' (largeType denotes the number of cells in datasets is very large  (typically >2k cells). 
In such case, the performance will be scarified to improve the running efficiency (e.g. using PCA + K-means instead of spectral clustering). 
If not set (**None**), the dataset will be regarded as normal, if set as 'True' or '1', the dataset will be treated as largeType. 

**Method**: **[getClusteringPars()](#clustering_getClusteringPars)**  
 
* **Output**:
Parameters needed for clustering-dCK and dBS. This function can be used to learn the
clustering parameters.  
* **dCK**: dictionary   
key：timePoint, value: K (number of clusters, Integer) , e.g {14:1, 16:2, 18:5}  
number of clusters for clustering at each time point.    

* **dBS**: dictionary  
key: timePoint, value: seed (Integer), e.g. {14:0, 16:0, 18:1}  
clustering seed for each time point

*  **Example**:

```python
import scdiff
from scdiff import *
Clustering_example=scdiff.Clustering(AllCells,'auto',None)
[dCK,dBS]=Clustering_example.getClusteringPars()
```
**Method**: **[performClustering()](#clustering_performClustering)**  

* **Output**: Clusters instances (Clustering results), please check [Cluster](#cluster) 
for details. This function is used to cluster all the nodes into clusters(Graph nodes).

* **Example**:   

```python
import scdiff 
from scdiff import *
Clustering_example=scdiff.Clustering(AllCells,'auto',None)
Clusters=Clustering_example.performClustering()

```

**[scdiff.Cluster(Cells,TimePoint,Cluster_ID)](#cluster)<a id="cluster"></a>**  
This class defines the node in the differentiation graph. 

**Parameters**:  

* **Cells**: List of Cell  
[Cell](#cell) instances.   
* **TimePoint**: Integer  
Initial Time Point for Cluster, it's the dominant measurement time for 
all cells within the cluster.   
* **Cluster_ID**: String   
Cluster ID.  

**Output**:  List of float, this function calculates the average 
gene expression of all cells in cluster. 

**Attributes**: 

* **cells**: List of Cell instances  
* **T**: Cluster time (Integer/float)   
* **ID**: Cluster ID (String)  
* **P**: Parent Cluster (Cluster instance)
* **C**: Children Clusters (List of Cluster instances)
* **E**: Mean gene expression of the Cluster  (List of float)  
* **R**: Gene expression variance of the Cluster (List of float)
* **GL**: List of gene names (List of String  )

**Example**:

```python
import scdiff 
from scdiff import *
cluster1=scdiff.Cluster([item for item in AllCells if item.T==14],14,'C1')
```

**[scdiff.Path(fromNode,toNode,Nodes,dTD,dTG,dMb)](#path)**

This class defines the edge in the differentiation graph. 

**Parameters**:

* **fromNode**: Cluster  
The beginning end  of an edge, Cluster instance

* **toNode**: Cluster
The ending end of an edge, Cluster instance 

* **Nodes**: List of Cluster  
All Nodes in Graph. 

* **dTD,dTG,dMb**:
The same as described in the scdiff.Graph class.  

**Output**:
Graph edge instance. 

**Attributes**: 
* **fromNode**: Cluster instance (source node of the edge).
* **toNode**: Cluster instance (target node of the edge). 
* **diffG**: differentially expressed genes on the edge. 
* **atf**: regulating TFs on the edge. 
                               
**Example**:

```python
import scdiff 
from scdiff import *
g1=scdiff.Graph(AllCells,"example.tf_dna",'auto')
p1=scdiff.Path(g1.Nodes[0],g1.Nodes[1],g1.Nodes,g1.dTD,g1.dTG,g1.dMb)
```


##  viz module 
This module is designed to visualize the differentiation graph structure using JavaScript.   

**[scdiff.viz(exName,Graph,output)](#viz)**

**Parameters**:

* **exName**: String
The name of the output visualization result. 

* **Graph**: Graph 
Graph instance, please refer [Graph](#graph).
* **output**: output directory 

**Output**:
a visualization folder with HTML page, JavaScript Code and Graph Structure in JSON format. 


**Example**:

```python 
import os
import scdiff
from scdiff import *
print ("testing scdiff.viz module ...")
# visualizing graph using scdiff.viz module 
os.mkdir("e1_out")
scdiff.viz("example",g1,"e1_out")

```
Then, you will find the visualized result page in HTML under 'e1_out' directory.




# CREDITS
 
This software was developed by ZIV-system biology group @ Carnegie Mellon University.  
Implemented by Jun Ding


# LICENSE 
 
This software is under MIT license.  
see the LICENSE.txt file for details. 


# CONTACT

zivbj at cs.cmu.edu  
jund  at cs.cmu.edu




                                 
                                 
                                 
                                 
                                 

                                                     
