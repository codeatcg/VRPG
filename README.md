# VRPG
an interactive web viewer for reference pangenome graphs


# Description  
VRPG is an interactive web viewer for pangenome graph. It supports reference Graph Fragment Assembly format, which is taken by minigraph. It can be used to scan the large and complex pangenome graph with many assemblies. The reference nodes stand in one line in the middle of view windows and are encompassed by non-reference nodes. It’s beneficial to check the structure variations. A website shipping Saccharomyces cerevisiae pangenome graph with 163 assemblies and Homo sapiens pangenome graph with 90 assemblies can be visited at https://www.evomicslab.org/app/vrpg/. Users can also deploy the web application and view their own data by themselves.  

# Installation  
Python3 (>=3.6) and pip environment are required.  

```
git clone https://github.com/codeatcg/VRPG --recursive  
pip install Django==3.2.4  
```

# Prepare your data  

## Pangenome graph already exists  
The naming scheme of assembly should follow <a href="https://github.com/pangenome/PanSN-spec">PanSN prefix naming pattern</a>. Briefly, the assembly’s name consists of sample name, delimiter, and haplotype name, e.g., sampleA#0. But it’s a little looser in VRPG. It’s not required that the haplotype name must be numeric, characters are also allowed. The assemblies to graph mapping files are also required. If these files do not exist the assembly can’t be highlighted in the drawing. These files can be generated by minigraph by using command ‘-cxasm --vc’. Then run the following command to get files required by VRPG.  

```
Python vrpg_preprocess.py --rGFA all.gfa --gafList gaf_file.list --outDir out_folder
```

### gaf_file.list file is formatted as follows: 

sample1#H1	sample1.H1.gaf  
sample2#0	sample2.0.gaf  
sample3#1	sample3.1.gaf  
sample3#2	sample3.2.gaf  

## Pangenome graph not exists  

Run the following command to create pangenome graph and generate files required by VRPG.  

```
Python vrpg_preprocess.py –-minigraph ‘/software/minigraph’ --assList ass_file.list –-outDir out_folder  
```

### ass_file.list file is formatted as follows:  
sample1#H1	sample1.H1.fa  
sample2#0	sample2.0.fa  
sample3#1	sample3.1.fa  
sample3#2	sample3.2.fa  

*Note*, ‘/software/minigraph’ represents the absolute path of minigraph executable file. Assembly in first line in file ass_file.list will be taken as reference.  


# Execution  
1. Move all output files in directory ‘upload’ generated during data preparation to the empty folder ‘upload’ of VRPG.
2. Start the development server of Django  

```
python3 manage.py runserver  

If all is well you will see the output:  
Django version 3.2.4, using settings 'primers_project.settings'  
Starting development server at http://127.0.0.1:8000/  
```

3. Then open http://127.0.0.1:8000/app/vrpg/ in your browser and you will see the drawing of the pangenome graph.  

If server is running on different machine start the server by running

```
python3 manage.py runserver 0.0.0.0:8000
```
Please make sure the firewall is closed. Then open <a>http:\<IP of server\>:8000/app/vrpg/</a>.

If you are familiar with nginx or apache you can also deploy VRPG by using any of them.



