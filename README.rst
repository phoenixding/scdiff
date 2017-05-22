::

             ____  ____  ____  _  _____ _____
            / ___\/   _\/  _ \/ \/    //    /
            |    \|  /  | | \|| ||  __\|  __\
            \___ ||  \__| |_/|| || |   | |   
            \____/\____/\____/\_/\_/   \_/    

INTRODUCTION
============

SCDIFF was designed to analyze the cell differentiation process using
time-series single cell RNA-seq data. It can be used to predict the
transcription factors, which regulate the cell differentiation process.
It also visualizes the differentiation process using a graph, in which
nodes represent different sub-population cells and edges denote the
differentiation paths.

.. figure:: http://www.cs.cmu.edu/~jund/scdiff/images/FlowChart.jpg
   :scale: 50
   :alt: flowchart


PREREQUISITES
=============
-  | python (python 2 and python 3 are both supported)
   | It was installed by default for most Linux distribution and MAC If
   | not, please check https://www.python.org/downloads/ for installation
   | instructions.

-  | Python packages dependencies:
   | -- scikit-learn
   | -- scipy
   | -- numpy
   | The python setup.py script will try to install these packages
   | automatically. However, please install them manually if, by any
   | reason, the automatic installation fails.

-  If using Graphic User Interface(GUI), please instal the Tkinter
   package. Please check
   http://tkinter.unpythonic.net/wiki/How_to_install_Tkinter to install
   Tkinter package.

INSTALLATION
============

cd to the package root directory

.. code:: shell

    $cd scdiff

run python setup to install

.. code:: shell

    $python setup.py install

Linux, Mac:

.. code:: shell

    $sudo python setup.py install 

USAGE
=====

.. code:: shell

    usage: scdiff [-h] -i INPUT -t TF_DNA -k CLUSTERS -o OUTPUT [-s SPEEDUP] [-l LargeType]

      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            input single cell RNA-seq expression data
      -t TF_DNA, --tf_dna TF_DNA
                            TF-DNA interactions used in the analysis
      -k CLUSTERS, --clusters CLUSTERS
                            how to learn the number of clusters for each time
                            point? user-defined or auto? if user-defined, please
                            specify the configuration file path
      -o OUTPUT, --output OUTPUT
                            output folder to store all results
                            
      -s SPEEDUP, --speedup SPEEDUP(1/0)
                            speedup version with less stringent convergence
                            criteria
      -l LARGETYPE  --largetype LARGETYPE (1/0)

MODULES & FUNCTIONS
===================

scdiff module
-------------

This module is used to perform the single cell differentiation analysis
and it builds a graph (differentiation) based on the analysis result.

`scdiff.Graph(Cells, kc,largeType=None) <#graph>`__

| This class defines the differentiation graph.

**Parameters**:

-  **Cells**: Cell instances
   Please read `cell <#cell>`__ Class definition for details.
-  **kc**: String
   clustering config. It's a string with either 'auto' or clustering
   configure file path (-k parameter).
-  **largeType**: None(default) or String
   whether the single cell is a 'largeType' (largeType denotes the
   number of cells in datasets is very large (typically >2k cells). In
   such case, the performance will be scarified to improve the running
   efficiency (e.g. using K-means instead of spectral clustering). If
   not set (**None**), the dataset will be regarded as normal, if set as
   'True', the dataset will be treated as largeType.

| **Output**:
| A graph instance with all nodes and edges, which represents the
| differentiation structure for given inputs.

**Example**:

.. code:: python

    import scdiff
    from scdiff import *
    graph1=scdiff.Graph(Cells,'auto',None)  #Cells: List of Cell instances 
     

| `scdiff.Cell(Cell\_ID, TimePoint,Expression,typeLabel) <#cell>`__\ 
| This class defines the cell.

**Parameters**:

-  **Cell\_ID**: String
   The ID of the cell.
-  **TimePoint**: Integer
   Measurement TimePoint of the cel, Integer.
-  **Expression**: List of float
   Expression of all genes.
-  **Cell\_Label**: String
   The true label for the cell if available, 'NA' if not available.
   (Note, we don't need this information for the model, but it's useful
   when analyzing the result).

| **Output**:
| A Cell class instance (with all information regarding to a cell)

**Example**:

.. code:: python

    import scdiff
    from scdiff.scdiff import *
    c1=Cell('C1',1,[0.1,4.2,....,3.6],'AT1')
    c2=Cell('C2',1,[0.1,4.2,....,3.6],'AT1')
    AllCells=[c1,c2]

`scdiff.Clustering(Cells, kc,largeType=None) <#graph>`__

**Parameters**:

-  **Cells**: List of Cell Please read `Cell <#cell>`__ Class definition
   for details.
-  **kc**: String
   clustering config. It's a string with either 'auto' or clustering
   configure file path (-k parameter).
-  **largeType**: None(default) or String
   whether the single cell is a 'largeType' (largeType denotes the
   number of cells in datasets is very large (typically >2k cells). In
   such case, the performance will be scarified to improve the running
   efficiency (e.g. using K-means instead of spectral clustering). If
   not set (**None**), the dataset will be regarded as normal, if set as
   'True', the dataset will be treated as largeType.

| **Method**:
| `getClusteringPars() <#clustering_getClusteringPars>`__  

| This class represents the clustering.

-  **Output**: Parameters needed for clustering-dCK and dBS. This
   function can be used to learn the clustering parameters.
-  | **dCK**: dictionary
   | key：timePoint, value: K (number of clusters, Integer) , e.g {14:1,
   | 16:2, 18:5}
   | number of clusters for clustering at each time point.

-  | **dBS**: dictionary
   | key: timePoint, value: seed (Integer), e.g. {14:0, 16:0, 18:1}
   | clustering seed for each time point

-  **Example**:

.. code:: python

    import scdiff
    from scdiff import *
    Clustering_example=scdiff.Clustering(Cells,'auto',None)
    [dCK,dBS]=Clustering_example.getClusteringPars()

**Method**: `performClustering() <#clustering_performClustering>`__

-  **Output**: Clusters instances (Clustering results), please check
   `Cluster <#cluster>`__ for details. This function is used to cluster
   all the nodes into clusters(Graph nodes).

-  **Example**:

.. code:: python

    import scdiff 
    from scdiff import *
    Clustering_example=scdiff.Clustering(Cells,'auto',None)
    Clusters=Clustering_example.performClustering()

| `scdiff.Cluster((Cells,TimePoint,Cluster\_ID)) <#cluster>`__\ 
| This class defines the node in the differentiation graph.

**Parameters**:

-  **Cells**: List of Cell `Cell <#cell>`__ instances.
-  **TimePoint**: Integer
   Initial Time Point for Cluster, it's the dominant measurement time
   for all cells within the cluster.
-  **Cluster\_ID**: String
   Cluster ID.

**Method**: `getAvgEx() <#clusters_getAvgEx>`__:

-  **Output**: List of float, this function calculates the average gene
   expression of all cells in cluster.
-  **Example**:

.. code:: python

    import scdiff 
    from scdiff import *
    cluster1=scdiff.Cluster(Cells,14,'C1')
    AvgEx=cluster1.getAvgEx()

`scdiff.Path(fromNode,toNode,Nodes) <#path>`__

This class defines the edge in the differentiation graph.

**Parameters**:

-  | **fromNode**: Cluster
   | The beginning end of an edge, Cluster instance

-  **toNode**: Cluster The ending end of an edge, Cluster instance

-  | **Nodes**: List of Cluster
   | All Nodes in Graph.

**Method**: `getFC() <#>`__:

-  **Output**: Get the log fold change (2D List):
   [[foldchange,geneIndex,fromEx,toEx],...] along the edge(Path), from
   fromNode to toNode.

**Method**: `getetf() <#>`__


-  **Output**: Get the potential regulating TFs **[etf,dtf]** for given
   edge (fromNode->toNode).
-  **etf**: 2D List: [[p-value,TF\_name],...], which represents the
   regulating TFs and p-values for each edge
-  **dtf**: 2D List: [[p-value,TF\_name],...], sub-list of etf, which
   represents the regulating TFs (and p-value) with different expression
   for different differentiation paths.

-  **example**:

.. code:: python

    import scdiff 
    from scdiff import *
    C1=scdiff.Cluster(Cells_1,14,'C1')
    C2=scdiff.Cluster(Cells_2,16,'C2')
    AC=[C1,C2,...,Cn]
    p1=scdiff.Path(C1,C2,AC)
    fc=p1.getFC()
    [etf,dtf]=p1.getetf()

viz module
----------

This module is designed to visualize the differentiation graph structure
using JavaScript.

`viz.viz(exName,Graph,GeneNames,dTD,output) <#viz>`__

**Parameters**:

-  **exName**: String The name of the output visualization result.

-  **Graph**: Graph Graph instance, please refer `Graph <#graph>`__.

-  **GeneNames**: List of String

-  **dTD**: Dictionary
   Key: String, TF name.
   Value: List of String, Gene name. e.g. {'NKX2-1':
   ['SNX13','RAB30',...]}
-  | **output**: String
   | output directory name

**Output**: a visualization folder with HTML page, JavaScript Code and
Graph Structure in JSON format.

**Example**:

.. code:: python

    import scdiff 
    from scdiff import *
    exName='example'
    output_directory='example_out'
    GeneNames=['SNX13','RAB30',...]
    dTD={'NKX2-1': ['SNX13','RAB30',...]}
    g1=Graph(Cells,'auto',None)
    viz.viz(exName,g1,GeneNames,dTD,output_directory)

Then, you will find the visualized result page in HTML under
'example\_out' directory.

EXAMPLES
========

Run scdiff on given time-series single cell RNA-seq data

**1) Run with automatic config**

.. code:: shell

    $ scdiff -i <example.E> -t <example.tf_dna> -k auto -o <example_out>

-  | **-i/--input**:
   | **example.E** is the single cell RNA-seq dataset with following
   | format (tab delimited)

   ::

       cell_ID time    cell_label  ex_gene1    ex_gene2    ... ex_geneN

   -  cell\_ID: ID of the cell.
   -  time: measure time of the cell RNA-seq.
   -  cell\_lable: known label for the cell (if available) ,if not ,
      denote as NA.
   -  ex\_genei: expression of gene i (log2 gene expression value). Gene
      expression can be FPKM or RPM or any other acceptable gene
      expression measurements.

   Please read **example.E** for an example of acceptable time-series
   single cell dataset format.

-  | **-t/--tfdna**:
   | **example.tf\_dna** provides the TF-DNA interaction information
   | used in this study (TF target inforamtion) with tab delimited format.

   ::

       TF  TF_target   Score

   For example:

   ::

       ZNF238  TIMM8B  0.9
       SOX9    TIMM8B  0.9
       ZEB1    TIMM8B  0.9
       GATA4   TIMM8B  0.9
       CEBPA   RAB30   0.9
       NKX2-1  RAB30   0.9
       SRF RAB30   0.9
       SOX5    RAB30   0.9
       SRY RAB30   0.9
       POU1F1  RAB30   0.9
       POU2F1  RAB30   0.9
       NFKB1   KRI1    0.9
       E2F1    C11ORF35    0.9
       DSP C11ORF35    0.9
       ELSPBP1 C11ORF35    0.9
       EGR2    C11ORF35    0.9
       EGR1    C11ORF35    0.9
       NR2F2   C11ORF35    0.9
       LMO2    C11ORF35    0.9
       ESR2    C11ORF35    0.9
       HNF1A   C11ORF35    0.9
       EGR3    C11ORF35    0.9

   The TF-DNA directory provides the TF-DNA interaction file used in
   this study.
-  | **-k/--clusters**:
   | This specifies the clustering parameter (String). It's need to be
   | either 'auto' or path to the 'config' file. Here, 'auto' denotes the
   | clustering parameters will be learned automatically.

-  | **-o/--output**:
   | **example\_out** is the specified output directory.

**2) Run with user-defined config**

.. code:: shell

    $scdiff -i <example.E>  -t <example.tf_dna> -k <example.config> -o <example_out>

The format of example.E and example.tf\_dna are the same as described
above.

**example.config** specifies the custom initial clustering parameters.
This was used when we have some prior knowledge. For example, if we know
they are how many sub-populations within each time, we can just directly
specify the clustering parameters using the example.config file, which
provides better performance.

example.config format(tab delimited)

::

    time    #_of_clusters

For example:

::

    14  1  
    16  2  
    18  5  

However, if we don't have any prior knowledge about the sub-populations
within each time point. We will just use the automatic initial
clustering. :-k auto.

**3) Run scdiff on large single cell dataset**

.. code:: shell

    $ scdiff -i <example.E> -t <example.tf_dna> -k auto -o <example_out> -l True -s True

-i, -t, -k parameters were discussed above.

-  | **-l/--large (optional)**
   | String, if specified as 'True' or '1', scdiff will use LargeType
   | mode to improve the running efficiency (both memory and time). The
   | performance will be sacrificed to some extent. K-means will be used
   | for Clustering instead of Spectral clustering.

-  | **-s/--speedup (optional)**
   | Speedup the convergence, it will reduce the running time
   | significantly. This is highly recommended for large dataset. Based on
   | testing on lung single cell dataset (Treutlein 2014), the speedup
   | performance is just slightly worse (2 more cells were miss-assigned )

**4) example running result**

`example\_out <http://www.cs.cmu.edu/~jund/scdiff/result/treutlein2014_lung/treutlein2014_lung.html>`__

.. figure:: ./scdiff/img/example_out.jpg
   :alt: example\_out\_fig


CREDITS
=======

| This software was developed by ZIV-system biology group @ Carnegie
| Mellon University.
| Implemented by Jun Ding

LICENSE
=======

| This software is under MIT license.
| see the LICENSE.txt file for details.

CONTACT
=======

| zivbj at cs.cmu.edu
| jund at cs.cmu.edu
