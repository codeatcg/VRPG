#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <regex>
#include <unistd.h>
#include <algorithm>
#include <cstdint>
#include "vgraph.h"
//#include "gz.h"

using namespace std;

typedef struct{
    int start;
    int end;

    string name;
    float score;
    char strand;
} BedInfo;

// BedInfo,SimBed
template <typename T>
bool bedByPos(T &a, T &b){
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

typedef struct{
    int start;
    int end;
    string name;
    float score;
    char strand;
    uint16_t layer;
} SimBed;

//By default, features with layer > 50 for each track will not be output. Change the default behaviour by '--layer'
int getTrackLayer(int kLayer,int layStart,vector<string> &chrVec,map<string,vector<BedInfo> > &allChrGene,ofstream &out){
    int gLayMax = 1;
    for(string &tchr : chrVec ){
        vector<BedInfo> &chrGene = allChrGene[tchr];
        sort(chrGene.begin(),chrGene.end(),bedByPos<BedInfo>);
        //
        int maxLayer = 1;
        size_t sPos = 0;
        bool overlap = false;
        vector<int> layVec;        
        for(size_t k = 0; k < chrGene.size(); ++k){
            // Does current gene has an overlap with previous genes ?
            // no overp; new start
            vector<int> layerSta(maxLayer,0);
            overlap = false;
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
                    if(maxLayer > gLayMax){
                        gLayMax = maxLayer;
                    }
                }
            }else{
                //   -------       ---------
                //       ---------------------
                for(size_t j = sPos; j < k; ++j){
                    int tlayer = layVec[j - sPos] + layStart;
                    if(tlayer <= kLayer){
                        out<<tchr<<"\t"<<chrGene[j].start<<"\t"<<chrGene[j].end<<"\t"<<chrGene[j].name<<"\t"<<chrGene[j].score<<"\t"<<chrGene[j].strand<<"\t"<<tlayer<<endl;
                    }
                }
                //
                sPos = k;
                layVec.clear();
                layVec.push_back(1);
            }
            
        }
        // last gene has not been output, whether it has a overlap with previous genes
        if(! chrGene.empty()){
            for(size_t j = sPos; j < chrGene.size(); ++j){
                int tlayer = layVec[j - sPos] + layStart;
                if(tlayer <= kLayer){
                    out<<tchr<<"\t"<<chrGene[j].start<<"\t"<<chrGene[j].end<<"\t"<<chrGene[j].name<<"\t"<<chrGene[j].score<<"\t"<<chrGene[j].strand<<"\t"<<tlayer<<endl;
                }
            }
        }
        
    }
    //
    if(gLayMax > kLayer){
        gLayMax = kLayer;
    }
    return gLayMax;
}

void setTrackLayer(int layStart,vector<string> &chrVec,map<string,vector<BedInfo> > &allChrGene,ofstream &out){
    int tlayer = layStart + 1;
    for(string &tchr : chrVec ){
        vector<BedInfo> &chrGene = allChrGene[tchr];
        for(size_t k = 0; k < chrGene.size(); ++k){
            out<<tchr<<"\t"<<chrGene[k].start<<"\t"<<chrGene[k].end<<"\t"<<chrGene[k].name<<"\t"<<chrGene[k].score<<"\t"<<chrGene[k].strand<<"\t"<<tlayer<<endl;
        }
    }
}

void reduceBed(int kLayer,char *inBed,string &tkDesFile,string &outFile){
    ifstream in(inBed);
    if(! in){
        cerr<<"Error: file open failed. "<<inBed<<endl;
        exit(1);
    }
    ofstream kfh(tkDesFile.c_str());
    if(! kfh){
        cerr<<"Error: file open failed. "<<tkDesFile<<endl;
        exit(1);
    }
    ofstream out(outFile.c_str());
    if(! out){
        cerr<<"Error: file open failed. "<<outFile<<endl;
        exit(1);
    }
    //
    regex pat0("\\t");
    regex pat1("^track");
    regex pat2("name=\"(.+?)\"");
    regex pat3("description=\"(.+?)\"");
    int trackNum = 0;

    //char stat = '0';
    string seqID;
    string preSeqID = "";
    int start,end;
    string name = "*";
    float score = 0;
    char strand = '+';
    string line;
    map<string,vector<BedInfo> > allChrGene;
    vector<string> chrVec;
    string trackName = "*";
    string trackDes = "*";
    int tnum = 0;
    int cumLayMax = 0;
    while(getline(in,line)){
        if(line[0] == '#'){
            continue;
        }
        if(regex_search(line,pat1)){
            if(trackNum > 0){
                int layStart = cumLayMax;
                if(tnum > 4){
                    setTrackLayer(layStart,chrVec,allChrGene,out);
                    cumLayMax += 1;
                }else{
                    cumLayMax += getTrackLayer(kLayer,layStart,chrVec,allChrGene,out);
                }
                
                allChrGene.clear();
                chrVec.clear();
                //#nColumns layer   track_name  description
                kfh<<tnum<<"\t"<<cumLayMax<<"\t"<<trackName<<"\t"<<trackDes<<endl;
            }
            //
            ++trackNum;
            smatch tkname,tkdes;
            if(regex_search(line,tkname,pat2)){
                trackName = tkname.str(1);
            }else{
                trackName = "track" + to_string(trackNum);
            }
            if(regex_search(line,tkdes,pat3)){
                trackDes = tkdes.str(1);
            }else{
                trackDes = "*";
            }
            //
            preSeqID = "";
        }else{
            sregex_token_iterator pos(line.begin(),line.end(),pat0,-1);
            sregex_token_iterator pend;
            vector<string> vec;
            while(pos != pend){
                vec.push_back(pos->str());
                ++pos;
            }
            //
            tnum = vec.size();
            if(tnum < 3){
                cerr<<"Error: bed format error. The total number of columns should be >= 3. "<<line<<endl;
                exit(1);
            }
            //
            start = atoi(vec[1].c_str());
            end = atoi(vec[2].c_str());
            if(end <= start){
                cerr<<"Error: bed format error. chromEnd should be > chromStart. "<<line<<endl;
                exit(1);                
            }
            //
            start += 1;
            if(tnum == 3){
                //stat = '0';
                name = "*";
                score = 0;
                strand = '+';
            }else if(tnum == 4){
                //stat = '1';
                name = vec[3];
                score = 0;
                strand = '+';
            }else if(tnum == 5){
                //stat = '2';
                name = vec[3];
                score = atof(vec[4].c_str());
                strand = '+';
            }else{
                name = vec[3];
                score = atof(vec[4].c_str());
                if(vec[5] != "+" && vec[5] != "-"){
                    cerr<<"Error: bed format error. Strand should be '+' or '-'. "<<line<<endl;
                    exit(1); 
                }
                strand = vec[5][0];
            }
            //
            BedInfo tgene = {start,end,name,score,strand};
            if(vec[0] != preSeqID){
                if(allChrGene.find(vec[0]) != allChrGene.end()){
                    allChrGene[vec[0]].push_back(tgene);
                }else{
                    vector<BedInfo> tvec;
                    tvec.push_back(tgene);
                    allChrGene.emplace(vec[0],tvec);
                    //
                    chrVec.push_back(vec[0]);
                }
            }else{
                allChrGene[vec[0]].push_back(tgene);
            }
            //
            preSeqID = vec[0];
        }
    }
    //
    if(trackNum > 0){
        int layStart = cumLayMax;
        if(tnum > 4){
            setTrackLayer(layStart,chrVec,allChrGene,out);
            cumLayMax += 1;
        }else{
            cumLayMax += getTrackLayer(kLayer,layStart,chrVec,allChrGene,out);
        }
        //#nColumns layer   track_name  description
        kfh<<tnum<<"\t"<<cumLayMax<<"\t"<<trackName<<"\t"<<trackDes<<endl;
    }else{
        cerr<<"Error: format error. Track was not defined."<<endl;
        exit(1);
    }
    
    in.close();
    out.close();
    kfh.close();    
}

void getGeneMap(string &tkDesFile,string &geneFile,map<string,int> &refChrMap,map<int,vector<SimBed> > &allGeneLay){
    ifstream gf(geneFile.c_str());
    if(! gf){
        cerr<<"Error: file open failed. "<<geneFile<<endl;
        exit(1);
    }
    string line;
    //stringstream strStream;
    //[start end layer name]
    string tchr,name;
    int start,end;
    float score;
    char strand;
    string preChr = "";
    bool flag = false;
    int chrpos = 0;
    int layer;
    regex pat0("\\t");
    
    while(getline(gf,line)){
        sregex_token_iterator pos(line.begin(),line.end(),pat0,-1);
        sregex_token_iterator pend;
        vector<string> vec;
        vec.reserve(7);
        while(pos != pend){
            vec.push_back(pos->str());
            ++pos;
        }
        tchr = vec[0];
        start = atoi(vec[1].c_str());
        end = atoi(vec[2].c_str());
        name = vec[3];
        score = atof(vec[4].c_str());
        strand = vec[5][0];
        layer = atoi(vec[6].c_str());

        if(layer > 65535){
            cerr<<"Warning: layer > 65535. Too many RNA or too many overlapping genes. The item ("<<tchr<<" "<<start<<" "<<end<<" "<<name<<") will be skipped."<<endl;
            continue;
        }
        SimBed tsim;
        tsim.start = start;
        tsim.end = end;
        tsim.name = name;
        tsim.score = score;
        tsim.strand = strand;
        tsim.layer = (uint16_t)layer;

        if(tchr != preChr){
            if(refChrMap.find(tchr) != refChrMap.end()){
                flag = true;
                chrpos = refChrMap[tchr];
                //
                if(allGeneLay.find(chrpos) != allGeneLay.end()){
                    allGeneLay[chrpos].push_back(tsim);
                }else{
                    vector<SimBed> tvec;
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
    //
    int tkNum = 0;
    ifstream kfh(tkDesFile.c_str());
    while(getline(kfh,line)){
        tkNum++;
    }
    kfh.close();
    //
    if(tkNum > 1){
        for(auto &tm : allGeneLay){
            sort(tm.second.begin(),tm.second.end(),bedByPos<SimBed>);
        }
    }
}

void getBdDxRef(string &chrListFile,map<string,int> &refChrMap){
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

void indexNodeBed(string &tkDesFile,string &rndFile,string &rndDxFile,string &geneFile,string &chrListFile,string &ovFile,string &gDxFile){
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
    getBdDxRef(chrListFile,refChrMap);
    //
    map<int,vector<SimBed> > allGeneLay;
    getGeneMap(tkDesFile,geneFile,refChrMap,allGeneLay);
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
    long long ndByte = 0,ndUnit = sizeof(BedNode);
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
        map<int,vector<SimBed> >::iterator it;
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
            map<size_t,vector<BedNode> > gdCutMap;
            //
            for(size_t i = 0; i < (it->second).size(); ++i){
                BedNode gnode;
                bool fnd1 = false, fnd2 = false;
                gnode.node1 = 0;
                gnode.node2 = 0;
                // feature boundary
                for(int x = sPos; x < chrNdNum; ++x){
                    if((it->second)[i].start >= chrRnode[x].start){
                        if((it->second)[i].start <= chrRnode[x].pend){
                            gnode.node1 = chrRnode[x].node;
                            gnode.reStart1 = (it->second)[i].start - chrRnode[x].start;
                            fnd1 = true;
                            //
                            sPos = x;
                            //
                            if((it->second)[i].end <= chrRnode[x].pend){
                                gnode.node2 = chrRnode[x].node;
                                gnode.reStart2 = (it->second)[i].end - chrRnode[x].start;
                                fnd2 = true;
                                break;
                            }
                        }else{
                            // check
                            sPos = x + 1;
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
                    gnode.score = (it->second)[i].score;
                    // assign feature to segmentation
                    for(size_t k = aPos; k < acrVec.size(); ++k){
                        if((it->second)[i].start < acrVec[k].start){
                            if((it->second)[i].end >= acrVec[k].start){
                                if(gdCutMap.find(k) == gdCutMap.end()){
                                    vector<BedNode> gvec;
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
                                    vector<BedNode> gvec;
                                    gvec.push_back(gnode);
                                    gdCutMap.emplace(k,gvec);
                                }else{
                                    gdCutMap[k].push_back(gnode);
                                }
                            }else{
                                aPos = k + 1;
                            }
                        }
                    }
                }
            }
            //
            
            for(size_t k = 0; k < acrVec.size(); ++k){
                int num = 0;
                if(gdCutMap.find(k) != gdCutMap.end()){
                    for(BedNode &gd : gdCutMap[k]){
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
void dxRefNodeBed(int kLayer,char *inBed,char *chrMapFile,string &upDir){
    
    string bedFile = upDir + "/simplify.bed";
    string tkDesFile = upDir + "/track.info";
    string rndFile = upDir + "/node.ref.bw";
    string rndDxFile = upDir + "/node.ref.bdx";
    
    string ovFile = upDir + "/bed.ref.bw";
    string gDxFile = upDir + "/bed.ref.bdx";
    
    string chrListFile = upDir + "/chr.list";    
    
    //if(access(geneFile.c_str(),F_OK) != 0){
        reduceBed(kLayer,inBed,tkDesFile,bedFile);
    //}
    indexNodeBed(tkDesFile,rndFile,rndDxFile,bedFile,chrListFile,ovFile,gDxFile);
}

// Chr start end feature strand 
void addBed_usage(){
    cout<<"Usage:   GraphAnno addBed --inBed <bed_file> --upDir <upload_dir>"<<endl;
    cout<<"--inBed       <File>    Input bed file."<<endl;                                
    
    cout<<"--upDir       <Dir>     'upload' directory which including files generated by 'gfa2view' or 'vrpg_preprocess.py'"<<endl;
    cout<<"--layer       <Int>     maximum layers for each track, by default: 50"<<endl;
}

int addBed_main(int argc,char **argv){
    string upDir = "";
    char *inBed = nullptr,*chrMapFile = nullptr;
    int kLayer = 50;
    if(argc == 1){
        addBed_usage();
        return 1;
    }
    for(int i = 1; i < argc; ++i){
        if(strcmp(argv[i],"--inBed") == 0 || strcmp(argv[i],"-inBed") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --inBed is used but there is no value."<<endl;
                return 1;
            }
            inBed = argv[i];
        }else if(strcmp(argv[i],"--upDir") == 0 || strcmp(argv[i],"-upDir") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --upDir is used but there is no value."<<endl;
                return 1;
            }
            upDir = argv[i];
        }else if(strcmp(argv[i],"--layer") == 0 || strcmp(argv[i],"-layer") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --layer is used but there is no value."<<endl;
                return 1;
            }
            kLayer = atoi(argv[i]);
        }else if(strcmp(argv[i],"--help") == 0 || strcmp(argv[i],"-help") == 0){
            addBed_usage();
            return 1;
        }else{
            cerr<<"Error: undefined parameter."<<endl;
            addBed_usage();
            return 1;
        }
    }
    if(kLayer < 1){
        kLayer = 1;
    }
    //
    dxRefNodeBed(kLayer,inBed,chrMapFile,upDir);
    return 0;
}












