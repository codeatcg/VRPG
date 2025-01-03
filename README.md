# VRPG
an interactive web viewer for reference-projected pangenome graph

|  |  |
| --- | --- |
| <img src="https://github.com/codeatcg/VRPG/blob/main/static/images/window.svg" width=600px/> | <img src="https://github.com/codeatcg/VRPG/blob/main/static/images/description.svg" width=400px/> |


# Description  

VRPG is an open-source web application for fast access and interactive analysis of pangenome regions. Both <a href="https://github.com/lh3/gfatools/blob/master/doc/rGFA.md">rGFA</a> and <a href="https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md">GFAv1</a> pangenome graphs are supported. VRPG implements a block index system to support navigating the large and complex pangenomes created based on hundreds of whole genome assemblies, fluently. VRPG layouts the graph with maintaining the linearity of the primary reference. This makes VRPG capable of connecting the abundant genome annotation data based on the coordinate system of the primary linear reference to the graph. VRPG provides plenty of functions to help users to explore and understand the pangenome graph interactively, such as highlighting one or multiple assembly paths in the graph, finding the node depth, searching a sequence, simplifying the pangenome graph. We also provide an example website at https://www.evomicslab.org/app/vrpg/. Through the example website users can quickly access most of regions of the human pangenomes created by <a href="https://humanpangenome.org/">HPRC</a> based on 90 assemblies and the saccharomyces cerevisiae pangenome created by <a href="https://www.evomicslab.org/">Evomics Lab</a> based on 163 assemblies.

For detailed information about VRPG functions, please read the document at https://evomicslab.org/app/vrpg/manual/.

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

# install the dependence for sequence-to-graph mapping
cd VRPG
mkdir bin
cd bin
wget -c https://github.com/lh3/minigraph/releases/download/v0.20/minigraph-0.20_x64-linux.tar.bz2
tar -jxf minigraph-0.20_x64-linux.tar.bz2
mv minigraph-0.20_x64-linux/minigraph ./

```

# Prepare your data  

The naming scheme of assembly should follow <a href="https://github.com/pangenome/PanSN-spec">PanSN prefix naming pattern</a>. Briefly, the assembly's name consists of sample name, delimiter, and haplotype name, e.g., sampleA#0. But it's a little looser in VRPG. It's not required that the haplotype name must be numeric, characters are also allowed. When indexing the graph users can define the search depth (VRPG version > 0.1.2) by option ‘--xDep’. A small value for this option may cause some big bubbles on the rendered graph uncompleted. Owing to the linearity of the primary reference genome on graph rendered by VRPG the uncompleted bubble and its approximate location relative to the primary reference genome can still be recognized generally. In practice, users can set the value to 100.

## rGFA format graph  

### Pangenome graph already exists  

The assemblies to graph mapping files are required. If these files do not exist the assembly can't be highlighted in the graph. These files can be generated by minigraph by using command '-cxasm --vc'. Then run the following command to get files required by VRPG.  

```
Python script/vrpg_preprocess.py --rGFA all.gfa --gafList gaf_file.list --outDir out_folder --index --xDep 100
```

#### gaf_file.list file is formatted as follows: 

sample1#H1	sample1.H1.gaf  
sample2#0	sample2.0.gaf  
sample3#1	sample3.1.gaf  
sample3#2	sample3.2.gaf  

### Pangenome graph not exists  

Run the following command to create pangenome graph and generate files required by VRPG.  

```
Python script/vrpg_preprocess.py –-minigraph '/software/minigraph' --assList ass_file.list –-outDir out_folder --index --xDep 100
```

#### ass_file.list file is formatted as follows:  
sample1#H1	sample1.H1.fa  
sample2#0	sample2.0.fa  
sample3#1	sample3.1.fa  
sample3#2	sample3.2.fa  

**Note**, '/software/minigraph' represents the absolute path of minigraph executable file. Assembly in first line in file ass_file.list will be taken as reference.  

## GFA format graph

For graphs in GFA format that can be processed by VRPG segment names should be numeric. Fortunately, graphs generated by Minigraph-Cactus and PGGB have this feature. If the segment names are not numeric users need to modify the graph first. Also notice that all path names in the graph should follow <a href="https://github.com/pangenome/PanSN-spec">PanSN prefix naming pattern</a>. If the path names don’t obey the rule the graph needs to be modified. This can be avoided by using proper assembly names before constructing the graph. If the graph satisfied the conditions described above run the following command to get files required by VRPG.  

```
module/gfa2view --GFA in.gfa --ref refName --outDir output_dir --index --xDep 100 --range 2000 --thread 10

# gfa2view is flexible. Users can also split the process into two steps.
# step 1: transform and calculate coverage
# This step can’t be paralleled.
module/gfa2view --GFA in.gfa --ref refName --outDir output_dir

# step 2: index
# This step can be paralleled.
module/gfa2view --outDir output_dir --index --xDep 100 --range 2000 --thread 10

```

By two steps users can test different options and parameters to index the graph, while avoiding to transform the graph repetitively. But note that the previous indexing results will be covered.

**Note**, For the current version of 'gfa2view' memory consumption is proportional to number of threads. A trade-off between speed and and memory consumption needs to be considered.  

If only a particular set of reference chromosomes or contigs are considered to visualize the option ‘--refChr’ can be used to save running time. The option only affects the process of indexing. If this option is specified, a file containing the expected chromosomes/contigs with one chromosome/contig per line is required. For primary linear reference with too many small contigs, it's suggested to specify the option '--refChr'.

## Annotation

Overlapping annotation items are placed on separate layers. For gene annotation track the maximum number of layers allowed by VRPG is 255. For each track defined in a [BED](https://asia.ensembl.org/info/website/upload/bed.html) file the maximum number of layers allowed by VRPG is 50 by default, which can be overridden by using option ‘--layer’. Note: for track with 'score' column (total number of columns > 3) defined in a BED file, all the annotation items (including overlapping ones) will be placed on the same layer. To separate the overlapping items, please define the track with 3 columns. Currently, VRPG can utilize at most the first 6 columns of a track defined in a BED file and other columns will be ignored.

```
# Add primary reference gene track from GFF file
# The GFF file must correspond to the primary reference
# run 'GraphAnno addRef --help' for help 
GraphAnno addRef --inGFF gffFile --chrTrans chrTransFile --upDir upload

# Add annotation track from BED file
# About the BED format, please refer to https://asia.ensembl.org/info/website/upload/bed.html
# Track line is needed for adding annotation track
GraphAnno addBed --inBed track.bed --upDir upload

# Add annotation for all nodes (reference and non-reference). By clicking on the node the genes that the node overlaps with will be showed in the 'Node information' viewport. 
# run 'GraphAnno nodeGene --help' for help
GraphAnno nodeGene --gffList gffListFile --upDir upload

```

## Node sequence

```
# Extract node sequence and build index
nodeSeq --gfaFile graph.gfa --upDir upload

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

3. Open http://127.0.0.1:8000/app/vrpg/ in your browser and visit any pangenome regions you are interested in.  

**Note**, For large pangenome graph it's better to prepare data in a computing server and then transfer the data to local 'upload' folder of VRPG.

## Remote server  

If server is running on a different machine start the server by running

```
python3 manage.py runserver 0.0.0.0:8000
```
Please make sure the firewall is closed. Then open <a>http:\<IP of server\>:8000/app/vrpg/</a>.

If you are familiar with nginx or apache you can also deploy VRPG by using any of them.


# Test  
Enter the directory [test](https://github.com/codeatcg/VRPG/tree/main/test) and see [TEST_README.md](https://github.com/codeatcg/VRPG/blob/main/test/TEST_README.md)


# Tips  
1. The edge can be highlighted by clicking on it. Remove the highlight by just clicking on it again.
2. For graphs generated by PGGB the duplicate sequence in the primary reference assembly may be collapsed. For the collapsed region new nodes and edges will be inserted (into the graph, add new L line and S line) by vrpg to make the primary reference maintain the linearity. But no nodes and edges will be inserted into the path (P line in the graph). If there are only reference nodes in a visited region and the highlighted reference path doesn’t pass these nodes the region may indicate a collapse. Another feature of these regions is that the numeric identifier of the node is generally much larger than the identifiers of the node in the flanking regions.  
<img src="https://github.com/codeatcg/VRPG/blob/main/static/images/ref_collapse3.png"/>
3. VRPG uses 1-based coordinate system, i.e. coordinate is denoted as [start,end]. If you are interest in a particular segment (node) and the genome annotation file (GFF3 format) is available you can use <a href="https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html">bedtools</a> and the coordinate displayed by VRPG (within Segment Source:) to find the annotation of the segment conveniently. For example:
  
```  
bedtools intersect -wa -wb -a node.pos -b genome.gff3
``` 
  
# More about this work 

https://www.biorxiv.org/content/10.1101/2023.01.20.524991v3


 
