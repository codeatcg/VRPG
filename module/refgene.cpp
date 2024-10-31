
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <regex>
#include <unistd.h>
#include <algorithm>
#include "vgraph.h"
#include "gz.h"

using namespace std;

// replace unknown with U

typedef struct{
    int start;
    int end; 
    string id;
    string name;
    char strand;
    uint8_t num;
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

typedef struct{
    int start;
    int end;
    string name;
    char type;
    char strand;
    uint8_t layer;
    uint8_t num;
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

/*
    Output: The value in last column is 0 for gene and non zero for exon or CDS. 
*/
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
    regex pat4("Parent=(.+?)(;|$)");
    string seqID,type,strand,attr,geneID,geneName,parent;
    string preSeqID = "";
    string preGeneID = "";
    
    string line;
    map<string,vector<GeneInfo> > allChrGene;
    vector<string> chrVec;
    
    map<string,vector<GeneInfo> > emap,cmap,gmap;
    vector<string> rParent;
    int gStart = 0,gEnd = 0;
    map<string,int> rnaCtmap;
    // no exon or no CDS
    while(getline(in,line)){
        if(line[0] != '#' && line[0] != '\0'){
            sregex_token_iterator pos(line.begin(),line.end(),pat0,-1);
            sregex_token_iterator pend;
            int i = 0;
            bool flag = false, eflag = false, cflag = false;
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
                    }else if(type == "exon"){
                        eflag = true;
                    }else if(type == "CDS"){
                        cflag = true;
                    }else{
                        flag = false;
                        break;
                    }
                }
                ++pos;
                ++i;
            }
            //
            if(flag){
                if(preGeneID != ""){
                    rnaCtmap.emplace(preGeneID,rParent.size());
                    //
                    vector<GeneInfo> transcript;
                    for(string &par : rParent){
                        int x = 1;
                        if(emap.find(par) != emap.end()){
                            sort(emap[par].begin(),emap[par].end(),byPos);
                            for(GeneInfo &bk : emap[par]){
                                bk.num = x;
                                transcript.push_back(bk);
                                //
                                ++x;
                            }
                        }
                        //
                        x = 1;
                        if(cmap.find(par) != cmap.end()){
                            sort(cmap[par].begin(),cmap[par].end(),byPos);
                            for(GeneInfo &bk : cmap[par]){
                                bk.num = x;
                                transcript.push_back(bk);
                                //
                                ++x;
                            }
                        }
                    }
                    //
                    gmap.emplace(preGeneID,transcript);
                }
                //
                emap.clear();
                cmap.clear();
                rParent.clear();
                //
                if(start < 1 || end < 1){
                    cerr<<"Error: failed to get feature position. "<<line<<endl;
                    exit(1);
                }

                smatch s;
                if(regex_search(attr,s,pat1)){
                    geneID = s[1];
                }else{
                    cerr<<"Error: attribute 'ID' was not defined. "<<line<<endl;
                    exit(1);
                }
                if(regex_search(attr,s,pat2)){
                    geneName = s[1];
                }else{
                    geneName = "U";
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
                //
                preGeneID = geneID;
                gStart = start;
                gEnd = end;
            }else if(eflag){
                if(start < 1 || end < 1){
                    cerr<<"Error: failed to get feature position. "<<line<<endl;
                    exit(1);
                }
                //
                if(start >= gStart && end <= gEnd){
                    smatch s;
                    if(regex_search(attr,s,pat4)){
                        parent = s[1];
                        if(emap.find(parent) == emap.end()){
                            vector<GeneInfo> tg;
                            tg.push_back({start,end,parent,"E",strand[0]});
                            emap.emplace(parent,tg);
                            rParent.push_back(parent);
                        }else{
                            emap[parent].push_back({start,end,parent,"E",strand[0]});
                        }
                    }else{
                        cerr<<"Error: attribute 'Parent' was not defined. "<<line<<endl;
                        exit(1);
                    }
                }else{
                    cerr<<"Warning: exon is out of the range of the gene or the gene is not defined. "<<line<<endl;
                }
            }else if(cflag){
                if(start < 1 || end < 1){
                    cerr<<"Error: failed to get feature position. "<<line<<endl;
                    exit(1);
                }
                //
                if(start >= gStart && end <= gEnd){
                    smatch s;
                    if(regex_search(attr,s,pat4)){
                        parent = s[1];
                        if(cmap.find(parent) == cmap.end()){
                            vector<GeneInfo> tg;
                            tg.push_back({start,end,parent,"C",strand[0]});
                            cmap.emplace(parent,tg);
                            if(emap.find(parent) == emap.end()){
                                rParent.push_back(parent);
                            }
                        }else{
                            cmap[parent].push_back({start,end,parent,"C",strand[0]});
                        }
                    }else{
                        cerr<<"Error: attribute 'Parent' was not defined. "<<line<<endl;
                        exit(1);
                    }
                }else{
                    cerr<<"Warning: CDS is out of the range of the gene or the gene is not defined. "<<line<<endl;
                }
            }
            
        }
    }
    //
    if(preGeneID != ""){
        rnaCtmap.emplace(preGeneID,rParent.size());
        //
        vector<GeneInfo> transcript;
        for(string &par : rParent){
            int x = 1;
            if(emap.find(par) != emap.end()){
                sort(emap[par].begin(),emap[par].end(),byPos);
                for(GeneInfo &bk : emap[par]){
                    bk.num = x;
                    transcript.push_back(bk);
                    //
                    ++x;
                }
            }
            //
            x = 1;
            if(cmap.find(par) != cmap.end()){
                sort(cmap[par].begin(),cmap[par].end(),byPos);
                for(GeneInfo &bk : cmap[par]){
                    bk.num = x;
                    transcript.push_back(bk);
                    //
                    ++x;
                }
            }
        }
        //
        gmap.emplace(preGeneID,transcript);
    }
    //
    emap.clear();
    cmap.clear();
    rParent.clear();
    for(string &tchr : chrVec ){
        vector<GeneInfo> &chrGene = allChrGene[tchr];
        sort(chrGene.begin(),chrGene.end(),byPos);
        //
        int maxLayer = 1;
        size_t sPos = 0;
        bool overlap = false;
        vector<int> layVec;        
        for(size_t k = 0; k < chrGene.size(); ++k){
            // Is current gene has a overlap with previous genes ?
            // no overlap; new start
            vector<int> layerSta(maxLayer,0);
            overlap = false;
            for(size_t m = sPos; m < k; ++m){
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
                //   -------       ---------
                //       ---------------------
                //        max rna count
                // glayer --------------> layCum
                if(k > 0){
                    map<int,int> layMax, layCum;
                    for(size_t p = sPos; p < k; ++p){                    
                        int rnaNum = rnaCtmap[chrGene[p].id];
                        int glayer = layVec[p - sPos];
                        
                        if(layMax.find(glayer) == layMax.end()){
                            layMax.emplace(glayer,rnaNum);
                        }else{
                            if(rnaNum > layMax[glayer]){
                                layMax[glayer] = rnaNum;
                            }
                        }
                    }
                    //
                    if(! layMax.empty()){
                        int tnum = 0;
                        for(auto &itm : layMax){
                            layCum.emplace(itm.first,tnum+1);
                            tnum += itm.second + 1;
                        }
                        //
                        for(size_t j = sPos; j < k; ++j){
                            int tlayer = layCum[layVec[j - sPos]];
                            out<<tchr<<"\t"<<chrGene[j].start<<"\t"<<chrGene[j].end<<"\t"<<chrGene[j].id<<"\t"<<chrGene[j].name<<"\t"<<chrGene[j].strand<<"\t"<<tlayer<<"\t0"<<endl;
                            //
                            if(gmap.find(chrGene[j].id) != gmap.end()){
                                string preID = "";
                                for(GeneInfo &tinfo : gmap[chrGene[j].id]){
                                    if(tinfo.id != preID){
                                        tlayer++;
                                    }
                                    out<<tchr<<"\t"<<tinfo.start<<"\t"<<tinfo.end<<"\t"<<tinfo.id<<"\t"<<tinfo.name<<"\t"<<tinfo.strand<<"\t"<<tlayer<<"\t"<<(int)tinfo.num<<endl;
                                    preID = tinfo.id;
                                } 
                            }
                        }
                    }
                }
                //
                sPos = k;
                maxLayer = 1;
                layVec.clear();
                layVec.push_back(1);
            }
            
        }
        // last gene has not been output, whether it has a overlap with previous genes
        if(! chrGene.empty()){
            map<int,int> layMax, layCum;
            for(size_t p = sPos; p < chrGene.size(); ++p){                    
                int rnaNum = rnaCtmap[chrGene[p].id];
                int glayer = layVec[p - sPos];
                
                if(layMax.find(glayer) == layMax.end()){
                    layMax[glayer] = rnaNum;
                }else{
                    if(rnaNum > layMax[glayer]){
                        layMax[glayer] = rnaNum;
                    }
                }
            }
            //
            if(! layMax.empty()){
                int tnum = 0;
                for(auto &itm : layMax){
                    layCum.emplace(itm.first,tnum+1);
                    tnum += itm.second + 1;
                }
                //
                for(size_t j = sPos; j < chrGene.size(); ++j){
                    int tlayer = layCum[layVec[j - sPos]];
                    out<<tchr<<"\t"<<chrGene[j].start<<"\t"<<chrGene[j].end<<"\t"<<chrGene[j].id<<"\t"<<chrGene[j].name<<"\t"<<chrGene[j].strand<<"\t"<<tlayer<<"\t0"<<endl;
                    //
                    if(gmap.find(chrGene[j].id) != gmap.end()){
                        string preID = "";
                        for(GeneInfo &tinfo : gmap[chrGene[j].id]){
                            if(tinfo.id != preID){
                                tlayer++;
                            }
                            out<<tchr<<"\t"<<tinfo.start<<"\t"<<tinfo.end<<"\t"<<tinfo.id<<"\t"<<tinfo.name<<"\t"<<tinfo.strand<<"\t"<<tlayer<<"\t"<<(int)tinfo.num<<endl;
                            preID = tinfo.id;
                        }
                    }
                }
            }
        }
        
    }
    in.close();
    out.close();
}

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

/*
    Gene name is used as priority. When gene name does not exist gene ID is used.
*/
void getGeneMap(string &geneFile,map<string,int> &refChrMap,map<int,vector<SimGene> > &allGeneLay){
    ifstream gf(geneFile.c_str());
    if(! gf){
        cerr<<"Error: file open failed. "<<geneFile<<endl;
        exit(1);
    }
    string line;
    stringstream strStream;
    string tchr,id,name;
    int start,end;
    char strand;
    string preChr = "";
    bool flag = false;
    int chrpos = 0;
    int layer;
    int num;
    while(getline(gf,line)){
        strStream << line;
        strStream >> tchr;
        strStream >> start;
        strStream >> end;
        strStream >> id;
        strStream >> name;
        strStream >> strand;
        strStream >> layer;
        strStream >> num;
        strStream.clear();
        strStream.str("");
        
        if(num > 255){
            cerr<<"Warning: exon count > 255. The item ("<<tchr<<" "<<start<<" "<<end<<" "<<id<<") will be skipped."<<endl;
            continue;
        }
        if(layer > 255){
            cerr<<"Warning: layer > 255. Too many RNA or too many overlapping genes. The item ("<<tchr<<" "<<start<<" "<<end<<" "<<id<<") will be skipped."<<endl;
            continue;
        }
        SimGene tsim;
        tsim.start = start;
        tsim.end = end;
        tsim.strand = strand;
        tsim.layer = (uint8_t)layer;
        tsim.num = (uint8_t)num;
        
        if(name == "E" || name == "C"){
            tsim.type = name[0];
            tsim.name = id;
        }else{
            tsim.type = 'G';
            if(name == "U"){
                tsim.name = id;
            }else{
                tsim.name = name;   
            }
        }

        if(tchr != preChr){
            if(refChrMap.find(tchr) != refChrMap.end()){
                flag = true;
                chrpos = refChrMap[tchr];
                //
                if(allGeneLay.find(chrpos)  != allGeneLay.end()){
                    allGeneLay[chrpos].push_back(tsim);
                }else{
                    vector<SimGene> tvec;
                    tvec.push_back(tsim);
                    allGeneLay.emplace(chrpos,tvec);
                }
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
    long long ndByte = 0,ndUnit = sizeof(GeneNode);
    for(int xchr : allchr){
        ChrRange cRange = chrRanMap[xchr];
        //
        vector<RanPos> acrVec;
        acrVec.reserve(cRange.ranNum);
        int chrNdNum = 0;
        for(int k = 0; k < cRange.ranNum; ++k){
            OneRange aRange;
            rxfh.read((char *)&aRange,oneSize);
            RanPos tpos = {aRange.ranStart,aRange.ranEnd};
            acrVec.push_back(tpos);
            chrNdNum += aRange.ranNum;
        }
        //
        bool fdChr = false;
        map<int,vector<SimGene> >::iterator it;
        it = allGeneLay.find(xchr);
        if(it != allGeneLay.end()){
            fdChr = true;
        }else{
            cout<<"Warning: reference chromosome in the 'chr.list' can't be found in the annotation file. "<<xchr<<endl;
        }
        //unordered_set<int> ntNode;
        vector<RNode> chrRnode;
        chrRnode.reserve(chrNdNum);
        for(int j = 0; j < chrNdNum; ++j){
            int node,ndStart,ndEnd;
            rnfh.read((char *)&node,intSize);
            rnfh.read((char *)&ndStart,intSize);
            rnfh.read((char *)&ndEnd,intSize);
            //
            RNode trnode = {node,ndStart,ndEnd};
            chrRnode.push_back(trnode);
            //ntNode.insert(node);
        }
        //--------------------------
        int sPos = 0;
        size_t aPos = 0;
        if(fdChr){
            map<size_t,vector<GeneNode> > gdCutMap;
            //
            for(size_t i = 0; i < (it->second).size(); ++i){
                GeneNode gnode;
                bool fnd1 = false, fnd2 = false;
                gnode.node1 = -1;
                gnode.node2 = -1;
                // feature boundary
                for(int x = sPos; x < chrNdNum; ++x){
                    if((it->second)[i].start >= chrRnode[x].start){
                        if((it->second)[i].start <= chrRnode[x].pend){
                            gnode.node1 = chrRnode[x].node;
                            gnode.reStart1 = (it->second)[i].start - chrRnode[x].start;
                            fnd1 = true;
                            //
                            if((it->second)[i].type == 'G'){
                                sPos = x;
                            }
                            //
                            if((it->second)[i].end <= chrRnode[x].pend){
                                gnode.node2 = chrRnode[x].node;
                                gnode.reStart2 = (it->second)[i].end - chrRnode[x].start;
                                fnd2 = true;
                                break;
                            }
                        }else{
                            // check
                            if((it->second)[i].type == 'G'){
                                sPos = x + 1;
                            }
                        }
                    }else{
                        if((it->second)[i].end >= chrRnode[x].start){
                            if((it->second)[i].end <= chrRnode[x].pend){
                                gnode.node2 = chrRnode[x].node;
                                gnode.reStart2 = (it->second)[i].end - chrRnode[x].start;
                                fnd2 = true;
                                break;
                            }
                        }else{
                            break;
                        }
                    
                    }
                }
                //
                if(fnd1 || fnd2){
                    gnode.name[FIELDSIZE-1] = '\0';
                    strncpy(gnode.name,(it->second)[i].name.c_str(),FIELDSIZE-1);
                    gnode.layer = (it->second)[i].layer;
                    gnode.strand = (it->second)[i].strand;
                    gnode.type = (it->second)[i].type;
                    gnode.num = (it->second)[i].num;
                    // assign feature to segmentation
                    for(size_t k = aPos; k < acrVec.size(); ++k){
                        if((it->second)[i].start < acrVec[k].start){
                            if((it->second)[i].end >= acrVec[k].start){
                                if(gdCutMap.find(k) == gdCutMap.end()){
                                    vector<GeneNode> gvec;
                                    gvec.push_back(gnode);
                                    gdCutMap.emplace(k,gvec);
                                }else{
                                    gdCutMap[k].push_back(gnode);
                                }
                            }else{
                                break;
                            }
                        }else{
                            if((it->second)[i].start <= acrVec[k].pend){
                                if(gdCutMap.find(k) == gdCutMap.end()){
                                    vector<GeneNode> gvec;
                                    gvec.push_back(gnode);
                                    gdCutMap.emplace(k,gvec);
                                }else{
                                    gdCutMap[k].push_back(gnode);
                                }
                            }else{
                                if((it->second)[i].type == 'G'){
                                    aPos = k + 1;
                                }
                            }
                        }
                    }
                }
            }
            //
            
            for(size_t k = 0; k < acrVec.size(); ++k){
                int num = 0;
                if(gdCutMap.find(k) != gdCutMap.end()){
                    for(GeneNode &gd : gdCutMap[k]){
                        ov.write((char *)&gd,ndUnit);
                    }
                    num = gdCutMap[k].size();
                }
                OneRange aRange;
                aRange.ranStart = acrVec[k].start;
                aRange.ranEnd = acrVec[k].pend;
                aRange.offByte = ndByte;
                aRange.ranNum = num;
                odx.write((char *)&aRange,oneSize);
                //
                ndByte += ndUnit * num;
            }
        }else{
            for(size_t k = 0; k < acrVec.size(); ++k){
                int num = 0;
                OneRange aRange;
                aRange.ranStart = acrVec[k].start;
                aRange.ranEnd = acrVec[k].pend;
                aRange.offByte = ndByte;
                aRange.ranNum = num;
                odx.write((char *)&aRange,oneSize);
            }
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
    
    //if(access(geneFile.c_str(),F_OK) != 0){
        reduceGFF(inGFF,chrMapFile,geneFile);
    //}
    indexNodeGene(rndFile,rndDxFile,geneFile,chrListFile,ovFile,gDxFile);
}

// Chr start end feature strand 
void addRef_usage(){
    cout<<"Usage:   GraphAnno addRef --inGFF <gff_file> --chrTrans <name_trans_file> --upDir <upload_dir>"<<endl;
    cout<<"--inGFF       <File>    Input GFF file."<<endl;                                  
    
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





