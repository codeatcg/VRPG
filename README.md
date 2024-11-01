# VRPG
an interactive web viewer for reference-projected pangenome graph

|  |  |
| --- | --- |
| <img src="https://github.com/codeatcg/VRPG/blob/main/static/images/window.svg" width=600px/> | <img src="https://github.com/codeatcg/VRPG/blob/main/static/images/description.svg" width=400px/> |


# Description  

VRPG is an interactive web viewer for reference-projected pangenome graph. It naturally supports graphs in reference Graph Fragment Assembly (<a href="https://github.com/lh3/gfatools/blob/master/doc/rGFA.md">rGFA</a>) format and for graphs in Graph Fragment Assembly (<a href="https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md">GFAv1</a>) format VRPG provides a command-line tool named gfa2view to transform the GFA files to a rGFA-like format. VRPG implements a block index system to support navigating the large and complex pangenome upon hundreds of whole genome assemblies in real time. The information about coordinate and copy number of each segment among the graph was stored in an efficient way and can be queried with almost no delay. VRPG aligns the reference nodes along the center line of the viewport, which make the reference genome easy to be recognized. VRPG also provides an intuitive way for genome comparison by highlighting the path of a particular assembly and its orientation on the rendered graph. A website shipping four pangenome graphs (one for yeast and three for human) is available at https://www.evomicslab.org/app/vrpg/. The Saccharomyces cerevisiae pangenome graph was generated using 163 assemblies and The three Homo sapiens pangenome graphs were constructed by <a href="https://humanpangenome.org/">HPRC</a>  by three different pipelines (Minigraph, Minigraph-CACTUS and PGGB) upon the same dataset with 90 whole genome assemblies. Users can also deploy the web application and view their own data.  

**Note:** 
The released version 0.1.3 is not the latest. The latest version of VRPG added a new tool named 'GraphAnno', which can be used to create indexed annotation files for reference gene track plot and interactive view of genes which a node overlaps with.

For graph in GFA format the overlap field in link line (overlap between segments)  should be specified (in graphs created by Minigraph-CACTUS and PGGB the overlap is generally specified as 0M), or the value will be set to 0 by VRPG-gfa2view. Although whether the overlap is specified doesn't affect the visualization of the graph it may affect the determination of the coordinate of the segment. 

The graphs created by Minigraph-CACTUS and PGGB include large amounts of SNPs and INDELs. The structure variations may be covered up by these small variants. VRPG (version >=0.1.3)  supplied functions to simplify the graph, i.e. remove nodes related to the small variants (with size < 50 bp). The simplification related option 'non-ref' in combobox means simplifying non-reference nodes. 'all node' means simplifying all nodes including reference and non-reference. 'none' means not simplifying the graph. Now VRPG can be used to visualize variants at different scales and find their coordinates relative to the reference  conveniently. Furthermore, the function that serves several pangenomes were added back in the latest version of VRPG.

For cola layout the node size is more proportional to the segment sequence size. But it may take a little longer time to stabilize. When the number of nodes in a window is small cola layout can be tested.  

# Installation  
Python3 (>=3.6) and pip environment are required.  

```
# For installing a historical version please access https://github.com/codeatcg/VRPG/releases and download the source code.

# install the latest version
# zlib
# gcc >= 4.9
pip install Django==3.2.4  pybind11
git clone https://github.com/codeatcg/VRPG --recursive  
cd VRPG/module
make

# By default the javascript packages that VRPG depends on are loaded from CDN. Users can also host the packages locally. 
python create.local.py
sh host.jslib.local.sh local
# switch to load packages from CDN
sh host.jslib.local.sh cdn

```

# Prepare your data  

The naming scheme of assembly should follow <a href="https://github.com/pangenome/PanSN-spec">PanSN prefix naming pattern</a>. Briefly, the assembly's name consists of sample name, delimiter, and haplotype name, e.g., sampleA#0. But it's a little looser in VRPG. It's not required that the haplotype name must be numeric, characters are also allowed. When indexing the graph users can define the search depth (VRPG version > 0.1.2) by option ‘--xDep’. In general, the default value can work well. A small value for this option may cause some big bubbles on the rendered graph uncompleted. Owing to the linearity of the reference genome on graph rendered by VRPG the uncompleted bubble and its approximate location relative to the reference genome can still be recognized generally. 

## rGFA format graph  

### Pangenome graph already exists  

The assemblies to graph mapping files are required. If these files do not exist the assembly can't be highlighted in the drawing. These files can be generated by minigraph by using command '-cxasm --vc'. Then run the following command to get files required by VRPG.  

```
Python script/vrpg_preprocess.py --rGFA all.gfa --gafList gaf_file.list --outDir out_folder --index
```

#### gaf_file.list file is formatted as follows: 

sample1#H1	sample1.H1.gaf  
sample2#0	sample2.0.gaf  
sample3#1	sample3.1.gaf  
sample3#2	sample3.2.gaf  

### Pangenome graph not exists  

Run the following command to create pangenome graph and generate files required by VRPG.  

```
Python script/vrpg_preprocess.py –-minigraph '/software/minigraph' --assList ass_file.list –-outDir out_folder --index
```

#### ass_file.list file is formatted as follows:  
sample1#H1	sample1.H1.fa  
sample2#0	sample2.0.fa  
sample3#1	sample3.1.fa  
sample3#2	sample3.2.fa  

**Note**, '/software/minigraph' represents the absolute path of minigraph executable file. Assembly in first line in file ass_file.list will be taken as reference.  

## GFA format graph

For graphs in GFA format that can be processed by VRPG segment names should be numeric. Fortunately, graphs generated by Minigraph-CACUTUS and PGGB have this feature. If the segment names are not numeric users need to modify the graph first. Also notice that all path names in the graph should follow <a href="https://github.com/pangenome/PanSN-spec">PanSN prefix naming pattern</a>. If the path names don’t obey the rule the graph needs to be modified. This can be avoided by using proper assembly names before constructing the graph. If the graph satisfied the conditions described above run the following command to get files required by VRPG.  

```
module/gfa2view --GFA in.gfa --ref refName --outDir output_dir --index --range 2000 --thread 10

# gfa2view is flexible. Users can also split the process into two steps.
# step 1: transform and calculate coverage
# This step can’t be paralleled.
module/gfa2view --GFA in.gfa --ref refName --outDir output_dir

# step 2: index
# This step can be paralleled.
module/gfa2view --outDir output_dir --index --range 2000 --thread 10 

```

By two steps users can test different options and parameters to index the graph, while avoiding to transform the graph repetitively. But note that the previous indexing results will be covered.

**Note**, For the current version of 'gfa2view' memory consumption is proportional to number of threads. A trade-off between speed and and memory consumption needs to be considered.  

If only a particular set of reference chromosomes or contigs are considered for view the option ‘--refChr’ can be used to save running time. The option only affects the process of indexing. If this option is specified, a file containing the expected chromosomes/contigs with one chromosome/contig per line is required.   

## Annotation

```
# Create files for reference gene track plot
# run 'GraphAnno addRef --help' for help 
GraphAnno addRef --inGFF gffFile --chrTrans chrTransFile --upDir upload

# Create files for interactive view of genes with which a node overlaps
# run 'GraphAnno nodeGene --help' for help
GraphAnno nodeGene --gffList gffListFile --upDir upload

```

# Execution  
## Local server or personal computer with Linux/Ubuntu operating system  

1. Move all output files in directory 'upload' generated during data preparation to the empty folder 'upload' of VRPG. If more than one pangenomes are available users can rename the results directory 'upload' generated during data preparation and then move the renamed directory into the VRPG 'upload' folder.

2. Start the development server of Django  

```
python3 manage.py runserver  

If all is well you will see the output:  
Django version 3.2.4, using settings 'primers_project.settings'  
Starting development server at http://127.0.0.1:8000/  
```

3. Then open http://127.0.0.1:8000/app/vrpg/ in your browser and you will see the drawing of the pangenome graph.  

**Note**, For large pangenome graph it's better to prepare data in a computing server and then transfer the data to local 'upload' folder of VRPG.

## Remote server  

If server is running on a different machine start the server by running

```
python3 manage.py runserver 0.0.0.0:8000
```
Please make sure the firewall is closed. Then open <a>http:\<IP of server\>:8000/app/vrpg/</a>.

If you are familiar with nginx or apache you can also deploy VRPG by using any of them.

# Tips  
1. The edge can be highlighted by clicking on it. Remove the highlight by just click on it again.
2. For graphs generated by PGGB the duplicate sequence in assigned reference assembly may be collapsed. For the collapsed region new nodes and edges will be inserted into the graph by vrpg to make the reference keeping coordinate system same to the traditional single linear reference. But the new nodes and edges will be not inserted into the path. If there are only reference nodes in a region and the highlighted reference path don’t pass these nodes the region may indicate a collapse. Another feature of these regions is that the numeric identifier of the segment that the node represents is generally much larger than the identifiers of the reference segments in the surrounding region.  
<img src="https://github.com/codeatcg/VRPG/blob/main/static/images/ref_collapse3.png"/>
3. VRPG uses 1-based coordinate system, i.e. coordinate is denoted as [start,end]. If you are interest in a particular segment (node) and the genome annotation file (GFF3 format) is available you can use <a href="https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html">bedtools</a> and the coordinate displayed by VRPG (within Segment Source:) to find the annotation of the segment conveniently. For example:
  
  
```  
bedtools intersect -wa -wb -a node.pos -b genome.gff3
``` 
  
# More about this work 

https://www.biorxiv.org/content/10.1101/2023.01.20.524991v3


 
