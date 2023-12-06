
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <regex>
#include <unistd.h>
//#include <algorithm>
#include "vgraph.h"
#include "gz.h"

using namespace std;

/*
#define FIELDSIZE 64
typedef struct NdGene{
    int node;
    int reStart;
    int len;
    char name[FIELDSIZE];
    uint8_t layer;
    char strand;
} NodeGene;
*/

typedef struct{
    int start;
    int end; 
    string id;
    string name;
    char strand;
} GeneInfo;
// id name strand chr start end
bool byPos(GeneInfo &a,GeneInfo &b){
    if(a.start < b.start){
        return true;
    }
    if(a.start == b.start){
        if(a.end < b.end){
            return true;
        }
    }
    return false;
}

//
// chr start end id name strand layer
// node gene_name gene_start gene_end gene_layer
// -----------
// | a |---
//     | b |
// node gene_name relative_start gene_len gene_layer


// nchr chr1 {offset nline} chr2 {offset nline} GeneLayer
// nchr chr1 {offset nline} chr2 {offset nline} {start end offset}    // dxFile
// NodeGene                // node gene
// ref_draw_pos   vector<int> [start,end,start,end]
// ref_draw_layer vector<int>
// ref_draw_gene_name vector<string>

typedef struct{
    int start;
    int end;
    string name;
    uint8_t layer;
    char strand;
} SimGene;

void chrNameTrans(char *chrMapFile,map<string,string> &transMap){
    ifstream in(chrMapFile);
    if(! in){
        cerr<<"Error: file open failed. "<<chrMapFile<<endl;
        exit(1);
    }
    string line;
    stringstream strStream;
    while(getline(in,line)){
        if(line[0] == '#'){
            continue;
        }
        strStream << line;
        string gffChr,graphChr;
        strStream >> gffChr;
        strStream >> graphChr;
        strStream.clear();
        strStream.str("");
        if(transMap.find(gffChr) == transMap.end()){
            transMap.emplace(gffChr,graphChr);
        }else{
            cout<<"Warning: duplicated chromosome. "<<line<<endl;
        }
    }
    
    in.close();
}

void reduceGFF(char *inGFF,char *chrMapFile,string &outFile){
    igzstream in(inGFF);
    if(! in){
        cerr<<"Error: file open failed. "<<inGFF<<endl;
        exit(1);
    }
    ofstream out(outFile.c_str());
    if(! out){
        cerr<<"Error: file open failed. "<<outFile<<endl;
        exit(1);
    }
    //
    map<string,string> transMap;
    bool trans = false;
    if(chrMapFile != nullptr){
        trans = true;
        chrNameTrans(chrMapFile,transMap);
    }
    //
    regex pat0("\\t");
    regex pat1("ID=(.+?)(;|$)");
    regex pat2("Name=(.+?)(;|$)");
    regex pat3("gene$");
    string seqID,type,strand,attr,geneID,geneName;
    string preSeqID = "";
    
    string line;
    map<string,vector<GeneInfo> > allChrGene;
    vector<string> chrVec;
    while(getline(in,line)){
        if(line[0] != '#' && line[0] != '\0'){
            sregex_token_iterator pos(line.begin(),line.end(),pat0,-1);
            sregex_token_iterator pend;
            int i = 0;
            bool flag = false;
            int start = 0,end = 0;
            while(pos != pend){
                switch(i){
                    case 0:
                        seqID = *pos;
                        break;
                    case 2:
                        type = *pos;
                        break;
                    case 3:
                        start = atoi((*pos).str().c_str());
                        break;
                    case 4:
                        end = atoi((*pos).str().c_str());
                        break;
                    case 6:
                        strand = *pos;
                        break;
                    case 8:
                        attr = *pos;
                        break;
                    default:
                        break;
                }
                //
                if(i == 2){
                    if(regex_search(type,pat3)){
                        flag = true;
                    }else{
                        break;
                    }
                }
                ++pos;
                ++i;
            }
            //
            if(flag){
                if(start < 1 || end < 1){
                    cerr<<"Error: failed to get feature position. "<<line<<endl;
                    exit(1);
                }

                smatch s;
                if(regex_search(attr,s,pat1)){
                    geneID = s[1];
                }else{
                    geneID = "unknown";
                }
                if(regex_search(attr,s,pat2)){
                    geneName = s[1];
                }else{
                    geneName = "unknown";
                }
                GeneInfo tgene = {start,end,geneID,geneName,strand[0]};
                if(trans){
                    if(transMap.find(seqID) != transMap.end()){
                        seqID = transMap[seqID];
                    }
                }
                if(seqID != preSeqID){
                    if(allChrGene.find(seqID) != allChrGene.end()){
                        allChrGene[seqID].push_back(tgene);
                    }else{
                        vector<GeneInfo> tvec;
                        tvec.push_back(tgene);
                        allChrGene.emplace(seqID,tvec);
                        //
                        chrVec.push_back(seqID);
                    }
                }else{
                    allChrGene[seqID].push_back(tgene);
                }
                preSeqID = seqID;
            }
            
        }
    }
    //
    //map<int,int> geneLayer;
    for(string &tchr : chrVec ){
        vector<GeneInfo> &chrGene = allChrGene[tchr];
        sort(chrGene.begin(),chrGene.end(),byPos);
        //
        int maxLayer = 1;
        size_t sPos = 0;
        bool overlap = false;
        vector<int> layVec;
        for(size_t k = 0; k < chrGene.size(); ++k){
            vector<int> layerSta(maxLayer,0);
            for(size_t m = sPos; m < k; ++m){
                //chrGene[k].start
                //chrGene[k].end
                if(chrGene[k].start <= chrGene[m].end && chrGene[k].end >= chrGene[m].start){
                    layerSta[layVec[m - sPos] - 1] = 1;
                    overlap = true; 
                }
            }
            if(overlap){
                bool fir = true;
                for(int n = 0; n < maxLayer; ++n){
                    //if(fir){
                        if(layerSta[n] == 0){
                            layVec.push_back(n+1);
                            fir = false;
                            break;
                        }
                    //}
                }
                if(fir){
                    ++maxLayer;
                    layVec.push_back(maxLayer);
                }
            }else{
                for(size_t j = sPos; j < k; ++j){
                    out<<tchr<<"\t"<<chrGene[j].start<<"\t"<<chrGene[j].end<<"\t"<<chrGene[j].id<<"\t"<<chrGene[j].name<<"\t"<<chrGene[j].strand<<"\t"<<layVec[j - sPos]<<endl;
                }
                //
                sPos = k;
                maxLayer = 1;
                layVec.clear();
                layVec.push_back(1);
                overlap = false;
            }
            
        }
        //
        for(size_t j = sPos; j < chrGene.size(); ++j){
            out<<tchr<<"\t"<<chrGene[j].start<<"\t"<<chrGene[j].end<<"\t"<<chrGene[j].id<<"\t"<<chrGene[j].name<<"\t"<<chrGene[j].strand<<"\t"<<layVec[j - sPos]<<endl;
        }
        
    }
    in.close();
    out.close();
}

// ldChrFile
// if(ldChrFile){ ldSet } if(ldSet)
void getDxRef(string &chrListFile,map<string,int> &refChrMap){
    ifstream in(chrListFile.c_str());
    int i = 0;
    string chrLine;
    while(getline(in,chrLine)){
        int tpos = chrLine.find("\t");
        string tchr = chrLine.substr(0,tpos);
        refChrMap.emplace(tchr,i);
        ++i;
    }
    in.close();
    if(refChrMap.empty()){
        cerr<<"Error: file is empty. "<<chrListFile<<endl;
        exit(1);
    }
}

void getGeneMap(string &geneFile,map<string,int> &refChrMap,map<int,vector<SimGene> > &allGeneLay){
    ifstream gf(geneFile.c_str());
    if(! gf){
        cerr<<"Error: file open failed. "<<geneFile<<endl;
        exit(1);
    }
    string line;
    stringstream strStream;
    //[start end layer name]
    string tchr,id,name;
    int start,end;
    char strand;
    int layer;
    string preChr = "";
    string sName;
    bool flag = false;
    int chrpos = 0;
    while(getline(gf,line)){
        strStream << line;
        strStream >> tchr;
        strStream >> start;
        strStream >> end;
        strStream >> id;
        strStream >> name;
        strStream >> strand;
        strStream >> layer;
        strStream.clear();
        strStream.str("");
        if(name != "unknown"){
            sName = name;
        }else{
            sName = id;
        }
        SimGene tsim = {start,end,sName,(uint8_t)layer,strand};
        if(tchr != preChr){
            if(refChrMap.find(tchr) != refChrMap.end()){
                flag = true;
                chrpos = refChrMap[tchr];
                //
                vector<SimGene> tvec;
                tvec.push_back(tsim);
                allGeneLay.emplace(chrpos,tvec);
            }else{
                flag = false;
            } 
        }else{
            if(flag){
                allGeneLay[chrpos].push_back(tsim);
            }
        }
        preChr = tchr;
    }
    gf.close();
}

void indexNodeGene(string &rndFile,string &rndDxFile,string &geneFile,string &chrListFile,string &ovFile,string &gDxFile){
    ofstream ov(ovFile.c_str());
    if(! ov){
        cerr<<"Error: file open failed. "<<ovFile<<endl;
        exit(1);
    }
    ofstream odx(gDxFile.c_str());
    if(! odx){
        cerr<<"Error: file open failed. "<<gDxFile<<endl;
        exit(1);
    }
    //
    map<string,int> refChrMap;
    getDxRef(chrListFile,refChrMap);
    //
    map<int,vector<SimGene> > allGeneLay;
    getGeneMap(geneFile,refChrMap,allGeneLay);
    //---------------------------------------------
    ifstream rxfh(rndDxFile.c_str());
    if(! rxfh){
        cerr<<"Error: file open failed. "<<rndDxFile<<endl;
        exit(1);
    }
    ifstream rnfh(rndFile.c_str());
    if(! rnfh){
        cerr<<"Error: file open failed. "<<rndFile<<endl;
        exit(1);
    }
    int intSize = sizeof(int);
    int crSize = sizeof(ChrRange);
    int nchr = 0;
    rxfh.read((char *)&nchr,intSize);
    odx.write((char *)&nchr,intSize);
    
    vector<int> allchr;
    allchr.reserve(nchr);
    map<int,ChrRange> chrRanMap;
    for(int t = 0; t < nchr; ++t){
        int tchr;
        ChrRange cRange;
        rxfh.read((char *)&tchr,intSize);
        rxfh.read((char *)&cRange,crSize);
        
        odx.write((char *)&tchr,intSize);
        odx.write((char *)&cRange,crSize);
        
        allchr.push_back(tchr);
        chrRanMap.emplace(tchr,cRange);
    }
    //
    int oneSize = sizeof(OneRange);
    //int llSize = sizeof(long long);
    //int dxByte = intSize + (intSize + crSize) * refChrMap.size();
    long long ndByte = 0,ndUnit = sizeof(NodeGene);
    for(int xchr : allchr){
        ChrRange cRange = chrRanMap[xchr];
        //
        //vector<OneRange> acrVec;
        //acrVec.reserve(cRange.ranNum);
        bool fdChr = false;
        map<int,vector<SimGene> >::iterator it;
        it = allGeneLay.find(xchr);
        if(it != allGeneLay.end()){
            fdChr = true;
        }else{
            cout<<"Warning: reference chromosome in the 'chr.list' can't be found in the annotation file. "<<xchr<<endl;
        }
        size_t sPos = 0;
        for(int k = 0; k < cRange.ranNum; ++k){
            OneRange aRange;
            rxfh.read((char *)&aRange,oneSize);
            // one region
            int num = 0;
            for(int j = 0; j < aRange.ranNum; ++j){
                int node,ndStart,ndEnd;
                rnfh.read((char *)&node,intSize);
                rnfh.read((char *)&ndStart,intSize);
                rnfh.read((char *)&ndEnd,intSize);
                // search gene
                if(fdChr){
                    bool flag = true;
                    size_t fir = sPos;
                    for(size_t i = sPos; i < (it->second).size(); ++i){
                        if((it->second)[i].start < ndStart){
                            if((it->second)[i].end >= ndStart){
                                if(flag){
                                    fir = i;
                                    flag = false;
                                }
                                //
                                NodeGene ndGene;
                                ndGene.node = node;
                                ndGene.reStart = 0;
                                ndGene.layer = (it->second)[i].layer;
                                ndGene.strand = (it->second)[i].strand;
                                
                                if((it->second)[i].end <= ndEnd){
                                    ndGene.len = (it->second)[i].end - ndStart + 1;
                                }else{
                                    ndGene.len = ndEnd - ndStart + 1;
                                }
                                
                                ndGene.name[FIELDSIZE-1] = '\0';
                                strncpy(ndGene.name,(it->second)[i].name.c_str(),FIELDSIZE-1);
                                
                                ov.write((char *)&ndGene,ndUnit);
                                ++num;
                            }
                            //else continue
                        }else{
                            if((it->second)[i].start <= ndEnd){
                                if(flag){
                                    fir = i;
                                    flag = false;
                                }
                                //
                                NodeGene ndGene;
                                ndGene.node = node;
                                ndGene.reStart = (it->second)[i].start - ndStart;
                                ndGene.layer = (it->second)[i].layer;
                                ndGene.strand = (it->second)[i].strand;
                                
                                if((it->second)[i].end <= ndEnd){
                                    ndGene.len = (it->second)[i].end - (it->second)[i].start + 1;
                                }else{
                                    ndGene.len = ndEnd - (it->second)[i].start + 1;
                                }
                                
                                ndGene.name[FIELDSIZE-1] = '\0';
                                strncpy(ndGene.name,(it->second)[i].name.c_str(),FIELDSIZE-1);
                                
                                ov.write((char *)&ndGene,ndUnit);
                                ++num;
                            }else{
                                if(flag){
                                    fir = i;
                                    flag = false;
                                }
                                //
                                break;
                            }
                        }
                    }
                    //
                    sPos = fir;
                }
            }
            //
            aRange.offByte = ndByte;
            aRange.ranNum = num;
            odx.write((char *)&aRange,oneSize);
            //
            ndByte += ndUnit * num;
        }
    }
    
    rxfh.close();
    rnfh.close();
    ov.close();
    odx.close();
}

//outDir
void dxRefNodeGene(char *inGFF,char *chrMapFile,string &upDir){
    
    string geneFile = upDir + "/simplify.gene";
    string rndFile = upDir + "/node.ref.bw";
    string rndDxFile = upDir + "/node.ref.bdx";
    
    string ovFile = upDir + "/gene.ref.bw";
    string gDxFile = upDir + "/gene.ref.bdx";
    
    string chrListFile = upDir + "/chr.list";
    //string chrListFile = upDir + "/load.chr.list";
    //if(access(chrListFile.c_str(),F_OK) != 0){
    //    chrListFile = upDir + "/chr.list";
    //}
    
    //if(access(geneFile.c_str(),F_OK) != 0){
        reduceGFF(inGFF,chrMapFile,geneFile);
    //}
    indexNodeGene(rndFile,rndDxFile,geneFile,chrListFile,ovFile,gDxFile);
}

// Chr start end feature strand 
void addRef_usage(){
    cout<<"Usage:   GraphAnno addRef --inGFF <gff_file> --chrTrans <name_trans_file> --upDir <upload_dir>"<<endl;
    cout<<"--inGFF       <File>    Input GFF file."<<endl;
    //  To Do 
/*    
    cout<<"--refAnno     <File>    Create files for gene track plot from this file instead of GFF file."
                                   "Five columns separated by whitespaces for each line (chromosome/contig name, start, end, feature, strand)."<<endl;
*/                                   
    
    cout<<"--chrTrans    <File>    Input file containing two columns separated by whitespaces for each line (chromosome/contig name in GFF file and chromosome/contig name in Graph file)."
                                   "This file is used to make the chromosome/contig names in GFF file matching that in graph file."
                                   "If chromosome name transformation is not needed this file can be ignored."<<endl;                                 
    
    cout<<"--upDir       <Dir>     'upload' directory which including files generated by 'gfa2view' or 'vrpg_preprocess.py'."<<endl;
}

int addRef_main(int argc,char **argv){
    string upDir = "";
    char *inGFF = nullptr,*chrMapFile = nullptr;
    if(argc == 1){
        addRef_usage();
        return 1;
    }
    for(int i = 1; i < argc; ++i){
        if(strcmp(argv[i],"--inGFF") == 0 || strcmp(argv[i],"-inGFF") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --inGFF is used but there is no value."<<endl;
                return 1;
            }
            inGFF = argv[i];
        }else if(strcmp(argv[i],"--upDir") == 0 || strcmp(argv[i],"-upDir") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --upDir is used but there is no value."<<endl;
                return 1;
            }
            upDir = argv[i];
        }else if(strcmp(argv[i],"--chrTrans") == 0 || strcmp(argv[i],"-chrTrans") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --chrTrans is used but there is no value."<<endl;
                return 1;
            }
            chrMapFile = argv[i];
        }else if(strcmp(argv[i],"--help") == 0 || strcmp(argv[i],"-help") == 0){
            addRef_usage();
            return 1;
        }else{
            cerr<<"Error: undefined parameter."<<endl;
            addRef_usage();
            return 1;
        }
    }
    //
    dxRefNodeGene(inGFF,chrMapFile,upDir);
    return 0;
}





