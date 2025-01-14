# VRPG
An interactive visualization and interpretation framework of reference-projected pangenome graphs  

|  |  |
| --- | --- |
| <img src="https://github.com/codeatcg/VRPG/blob/main/static/images/window.svg" width=600px/> | <img src="https://github.com/codeatcg/VRPG/blob/main/static/images/description.svg" width=400px/> |


# Description  

VRPG is a web-based interactive Visualization and interpretation framework for linear-Reference-projected Pangenome Graphs. VRPG provides efficient and intuitive supports for exploring and annotating pangenome graphs along a linear-genome-based coordinate system (e.g., that of a primary linear reference genome). Moreover, VRPG offers many unique features such as in-graph path highlighting for graph-constituent input assemblies, copy number characterization for graph-embedding nodes, graph-based mapping for query sequences, all of which are highly valuable for researchers working with pangenome graphs. Additionally, VRPG enables side-by-side visualization between the graph-based pangenome representation and the conventional primary-linear-reference-genome-based feature annotations, therefore seamlessly bridging the graph and linear genomic contexts.

Pangenome graphs encoded in both <a href="https://github.com/lh3/gfatools/blob/master/doc/rGFA.md">rGFA</a> and <a href="https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md">GFAv1</a> formats are supported. 

To further demonstrate its functionality and scalability, we provide a demonstration website for VRPG at https://evomicslab.org/app/vrpg/ to showcase its capacity and functionality with the cutting-edge yeast and human reference pangenome graphs derived from hundreds of high-quality genome assemblies.

For more details regarding VRPG's functionalities and usage, please check out the documentation at https://evomicslab.org/app/vrpg/manual/.


# Citation 

Zepu Miao, Jia-Xing Yue*. (2025) Interactive visualization and interpretation of pangenome graphs by linear-reference-based coordinate projection and annotation integration. Genome Research, (in press; doi: 10.1101/gr.279461.124; demonstration available at https://www.evomicslab.org/app/vrpg/; software available at https://github.com/codeatcg/VRPG) [[LINK](https://genome.cshlp.org/content/early/2025/01/13/gr.279461.124)]

# Installation  
Dependencies:  
`git`, `wget`, `bash`, `gcc/g++ >= 4.9`, `make`, `python3 (>=3.6)`, `pip`, `zlib`, `zlib-devel`  

```

# To install a historical version of VRPG, please access https://github.com/codeatcg/VRPG/releases and download the source code.
# To install the latest version of VRPG, please type
pip install Django==3.2.4  pybind11
git clone https://github.com/codeatcg/VRPG --recursive  
cd VRPG/module
make

# By default, the javascript packages that VRPG depends on are loaded from CDN. Users can optionally host these packages locally by typing: 
python create.local.py
sh host.jslib.local.sh local
# In case that users want to switch back to loading packages from CDN, type:
sh host.jslib.local.sh cdn

# To install the dependence (minigraph) for sequence-to-graph mapping, type:
cd VRPG
mkdir bin
cd bin
wget -c https://github.com/lh3/minigraph/releases/download/v0.20/minigraph-0.20_x64-linux.tar.bz2
tar -jxf minigraph-0.20_x64-linux.tar.bz2
mv minigraph-0.20_x64-linux/minigraph ./

```

# Input data preparation 

The naming scheme of assembly should follow <a href="https://github.com/pangenome/PanSN-spec">PanSN prefix naming specification</a>. Briefly, the assembly's name consists of sample name, delimiter, and haplotype name, e.g., sampleA#0. The haplotype name can be numeric (e.g. '0', '1', or '2') or characters (e.g., 'collapsed' or 'phased') or both (e.g., 'h1' or 'h2').  

When indexing the pangenome graphs, users can set the hard limit for maximal search depth  via the option ‘--xDep’ (for VRPG versions >=0.1.3). The default value for this option is 100. Setting this value too small may cause incomplete graph traversing for large bubbles.  

## For rGFA-formatted pangenome graphs  

### When already having the rGFA-formatted graph file  

Just type the following command to feed the graph file to VRPG: 

```
Python ./script/vrpg_preprocess.py --rGFA all.gfa --outDir out_folder --index --xDep 100
```

That said, if users want to use VRPG's assembly-to-path highlighting function, additional assembly-to-graph mapping files (in the GAF format) are needed, which can be generated by minigraph using the command '-cxasm --vc' (see [TEST_README.md](https://github.com/codeatcg/VRPG/blob/main/test/TEST_README.md) for more details). Once this is done, prepare a tab-delimited gaf_file.list file formatted as follows:

```
sample1#H1	sample1.H1.gaf  
sample2#0	sample2.0.gaf  
sample3#1	sample3.1.gaf  
sample3#2	sample3.2.gaf  
```

Then run the following command to prepare all the input information for VRPG.  

```
Python ./script/vrpg_preprocess.py --rGFA all.gfa --gafList gaf_file.list --outDir out_folder --index --xDep 100
```


### When not having the rGFA-formatted graph file

First prepare a tab-delimited asm_file.list file formatted as follows:

```
sample1#H1	sample1.H1.fa  
sample2#0	sample2.0.fa  
sample3#1	sample3.1.fa  
sample3#2	sample3.2.fa  
```
**Note**: The assembly defined on the first line in this file will be used as the primary linear reference genome for rGFA graph construction.

Then run the following command to create the rGFA-formatted pangenome graph and to generate all needed files for VRPG.  

```
Python ./script/vrpg_preprocess.py –-minigraph '/software/minigraph' --asmList asm_file.list –-outDir out_folder --index --xDep 100
``` 

**Note**: Here '/software/minigraph' represents the absolute path of the minigraph executable file.   

## For GFAv1-formatted pangenome graphs

For graphs encoded in GFAv1 format, VRPG requires its node (i.e., segment) name to be numeric. Graphs generated by Minigraph-Cactus and PGGB both satisfy this requirement. When dealing with a pangenome graph file of which the node/segment names are not numeric, users need to modify the graph first. Also notice that all path names defined in the graph should follow <a href="https://github.com/pangenome/PanSN-spec">PanSN prefix naming specification</a>. This requirement can be met by using proper assembly names before constructing the graph. Once these checks are all passed, run the following command to prepare the graph for VRPG:


```
./module/gfa2view --GFA in.gfa --ref refName --outDir output_dir --index --xDep 100 --range 2000 --thread 10
```

**Note**: The gfa2view step can be further divided into two steps when needed. This can be useful for testing different indexing parameters, but note that the index files generated in the new run will automatically overwrite the old ones. To run the gfa2view step in two steps, type:
```
# Step 1: format conversion and assembly-to-graph mapping depth calculation. This step cannot be parallelized.

./module/gfa2view --GFA in.gfa --ref refName --outDir output_dir

# Step 2: graph indexing. This step can be parallelized.

./module/gfa2view --outDir output_dir --index --xDep 100 --range 2000 --thread 10

```

**Note**: For now, the memory consumption of 'gfa2view' is proportional to number of threads. A trade-off between speed and and memory consumption needs to be considered.  

**Note**: For a primary linear reference genome with many small contigs (that will not be used in graph visualization), it is recommended to specify the '--refChr' option for 'gfa2view' to let VRPG only consider major chromosomes/contigs for graph indexing. This will help to substantially reduce the graph indexing time. For this option, a text file containing the names of the specified chromosomes/contigs (one line per chromosome/contig name) is needed. 


## Annotation tracks

A unique feature of VRPG lies in its native support for a side-by-side visualization between the pangenome graph and multiple annotation tracks based on the same primary linear reference coordinate system. For now, VRPG accepts annotation tracks defined in GFF3 and BED formats.

For the gene annotation track defined in GFF3 format, the maximal number of tiled layers (i.e., overlapped genes) to display is 255 and the maximal number of RNA isoforms per gene to display is 20 by default, which can be further adjusted using the '--rnaMax' option with the 'GraphAnno' command in VRPG. 

For other annotation tracks defined in the BED format, different tracks are specified via the track line notation defined in the [BED](https://asia.ensembl.org/info/website/upload/bed.html) file. The maximal number of tiled layers to display is 50 by default, which can be further adjusted using the ‘--layer’ option with the 'GraphAnno' command in VRPG. When the 5th column (i.e., score) of the BED file is defined, VRPG will assume this annotation track is score-based and render all records on the same layer. Otherwise, VRPG will try to separate overlapping records from the same track into different layers.   


```
# To add the gene annotation track based on the GFF3 file (based on the primary linear reference)
# run 'GraphAnno addRef --help' for more help information 
./module/GraphAnno addRef --inGFF gffFile --chrTrans chrTransFile --upDir upload

# To add other annotation tracks from BED files 
./module/GraphAnno addBed --inBed track.bed --upDir upload

# To add annotation for all nodes (reference and non-reference).
# run 'GraphAnno nodeGene --help' for more help information
./module/GraphAnno nodeGene --gffList gffListFile --upDir upload

```

## Node sequence extraction
Optionally, VRPG can be configured to report the exact genomic sequence for each node. 
```
# To extract the node sequences and build their index files for VRPG
./module/nodeSeq --gfaFile graph.gfa --upDir upload

```


# Deployment 
## For local server or personal computer running a Linux-based operating system (e.g., CentOS, Ubuntu)

1. Moving all files generated in directory 'upload' during data preparation to the empty directory named 'upload'. If users want to visualize more than one pangenome graphs, please see [TEST_README.md](https://github.com/codeatcg/VRPG/blob/main/test/TEST_README.md) for more details.

2. Starting the Django development server

```
python3 manage.py runserver  

If all is well you will see the output:  
Django version 3.2.4, using settings 'primers_project.settings'  
Starting development server at http://127.0.0.1:8000/  
```

3. Open http://127.0.0.1:8000/app/vrpg/ in your web browser to access VRPG's web portal.  

**Note**: For large pangenome graphs, it is recommended to use a computing server for the data preparation of VRPG and then transfer the preprocessed data to the 'upload' directory of the local server or personal computer. 

## For a remote server running a Linux-based operating system (e.g., CentOS, Ubuntu)

1. Logging in to the remote server and start the Django development server by running

```
python3 manage.py runserver 0.0.0.0:8000
```
2. On a local computer, using <a>http:\<IP of the remote server\>:8000/app/vrpg/</a> to use VRPG.

**Note**: The port 8000 can be any port number that has not occupied by any other process and is allowed by the firewall for outside visiting. If the user is familiar with nginx or apache, VRPG can be further configured using any of them. 


# Run VRPG with the testing example 
Enter the directory [test](https://github.com/codeatcg/VRPG/tree/main/test) and follow the instructions from [TEST_README.md](https://github.com/codeatcg/VRPG/blob/main/test/TEST_README.md)


# Additional Tips  
1. For graphs generated by PGGB, the duplicate sequence in the primary linear reference assembly might have been collapsed. VRPG will re-insert shadow nodes and edges into the graph (by adding new L and S lines) to restore the linearity during graph indexing (See the associated manuscript for more technical details). When primary-linear-reference path highlighting is enabled, these inserted shadow nodes can be recognized by the fact that no highlighted path goes through them. Also, the node ID of these shadow nodes should be much larger compared with the real primary-linear-reference nodes aside to them. 
<img src="https://github.com/codeatcg/VRPG/blob/main/static/images/ref_collapse3.png"/>

2. For now, VRPG does not automatically clean the directories created by running 'Sequence-to-graph mapping' jobs (i.e., '/upload/*/mapping/task_*') . To save the disk space of the hosting server,  we suggest users (or the system manager) to clean these files periodically.




  



 
