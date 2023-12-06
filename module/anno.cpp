


#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <regex>
//#include <map>
//#include <set>
#include <sys/stat.h>
#include <dirent.h>
#include <algorithm>

#include "gz.h"
#include "vgraph.h"

using namespace std;
/*
typedef struct{
    char seqid[FIELDSIZE];
    int start;
    int end;
    char geneID[FIELDSIZE];
    char geneName[FIELDSIZE];
    char strand;
} AnnoLine;

typedef struct{
    long long offset;
    int num;
} AnnoDx;
*/
//----------------------
typedef struct{
    string seqid;
    int start;
    int end;
    string geneID;
    string geneName;
    string strand;
} SimAnno;
typedef struct{
    int start;
    int end;
} NdPos;

// NdPos ChrAnnoPos

bool byPos(SimAnno &a,SimAnno &b){
    if(a.seqid < b.seqid){
        return true;
    }
    if(a.seqid > b.seqid){
        return false;
    }
    
    if(a.start < b.start){
        return true;
    }
    if(a.start > b.start){
        return false;
    }
    
    if(a.end < b.end){
        return true;
    }
    return false;
}

string getSep(string &sepFile){
    ifstream in(sepFile.c_str());
    if(! in){
        cerr<<"Error: file open failed. "<<sepFile<<endl;
        exit(1);
    }
    string sep;
    in >> sep;
    in.close();
    return sep;
}

// assChrPos
void getAssNode(string &nodeFile,string &sep,map<string,map<string,vector<NdPos> > > &nodeChrPos,string &refStr){
    ifstream in(nodeFile.c_str());
    string line;
    string node,r_chr;
    int start,end,len,refOr;
    stringstream strStream;
    string preAss = "";
    string preChr = "";
    bool isRef = false;
    getline(in,line);
    while(getline(in,line)){
        strStream << line;
        strStream >> node;
        strStream >> r_chr;
        strStream >> start;
        strStream >> end;
        strStream >> len;
        strStream >> refOr;
        
        strStream.clear();
        strStream.str("");
        //
        string tName,t_hap,tchr;
        assSplit(r_chr,sep,tName,t_hap,tchr);
        string assStr = tName + sep + t_hap;
        
        if(assStr != preAss){
            if(refOr == 0){
                refStr = assStr;
                isRef = true;
            }else{
                isRef = false;
            }
            //
            if(! isRef){
                map<string,vector<NdPos> > tChrPos;
                vector<NdPos> tvec;
                NdPos tndpos{start,end};
                tvec.push_back(tndpos);
                tChrPos.emplace(tchr,tvec);
                nodeChrPos.emplace(assStr,tChrPos);
            }
        }else{
            if(! isRef){
                NdPos tndpos{start,end};
                if(tchr != preChr){
                    vector<NdPos> tvec;
                    tvec.push_back(tndpos);
                    nodeChrPos[assStr].emplace(tchr,tvec);
                }else{
                    nodeChrPos[assStr][tchr].push_back(tndpos);
                }
            }
        }
        //
        preChr = tchr;
        preAss = assStr;
    }
    in.close();
}

void rgSearch(vector<NdPos> &chrNode,vector<SimAnno> &assAnnoVec,int usize,ofstream &out,int &total){
    size_t sPos = 0;
    //bool flag = true;
    for(size_t i = 0; i < assAnnoVec.size(); ++i){
        size_t tPos = sPos;
        for(size_t s = sPos; s < chrNode.size(); ++s){
            if(assAnnoVec[i].end < chrNode[s].start){
                tPos = s;
                break;
            }else{
                if(assAnnoVec[i].start <= chrNode[s].end){
                    AnnoLine anno;
                    anno.seqid[FIELDSIZE-1] = '\0';
                    anno.geneID[FIELDSIZE-1] = '\0';
                    anno.geneName[FIELDSIZE-1] = '\0';
                    
                    strncpy(anno.seqid,assAnnoVec[i].seqid.c_str(),FIELDSIZE-1);
                    strncpy(anno.geneID,assAnnoVec[i].geneID.c_str(),FIELDSIZE-1);
                    strncpy(anno.geneName,assAnnoVec[i].geneName.c_str(),FIELDSIZE-1);
                    anno.strand = assAnnoVec[i].strand[0];
                    anno.start = assAnnoVec[i].start;
                    anno.end = assAnnoVec[i].end;
                    
                    out.write((char *)&anno,usize);
                    
                    ++total;
                    //
                    tPos = s;
                    break;
                }
            }
        }
        sPos = tPos;
    }
}

void simpGFF(string &gffFile,string &chrNameFile,map<string,vector<NdPos> > &assChrPos,bool isRef,string &outFile){
    igzstream in(gffFile.c_str());
    ofstream out(outFile.c_str());
    if(! in){
        cerr<<"Error: file open failed. "<<gffFile;
        exit(1);
    }
    
    if(! out){
        cerr<<"Error: file open failed. "<<outFile;
        exit(1);
    }
    //
    stringstream strStream;
    string line;
    string pName,gName;
    map<string,string> chrMap;
    
    bool flag = false;
    if(chrNameFile != "NA"){
        ifstream rn(chrNameFile.c_str());
        if(! rn){
            cerr<<"Error: file open failed. "<<chrNameFile;
            exit(1);
        }
        while(getline(rn,line)){
            strStream << line;
            strStream >> pName;
            strStream >> gName;
            chrMap.emplace(pName,gName);
            strStream.clear();
            strStream.str("");
        }
        rn.close();
        //
        flag = true;
    }
    //
    regex pat0("\\t");
    regex pat1("ID=(.+?)(;|$)");
    regex pat2("Name=(.+?)(;|$)");
    regex pat3("gene$");
    
    
    string preSeqID = "";
    string gSeqID = "";
    string preGid = "";
    int usize = sizeof(AnnoLine);
    string seqID,type,strand,attr,geneID,geneName;

    bool chrFind = false;
    
    vector<SimAnno> assAnnoVec;
    //map<string,vector<SimAnno> > allGene;
    
    int total = 0;
    out.write((char *)&total,sizeof(int));
    while(getline(in,line)){
        if(line[0] != '#'){
            sregex_token_iterator pos(line.begin(),line.end(),pat0,-1);
            sregex_token_iterator pend;
            int i = 0;
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
                }
                //
                ++pos;
                ++i;
            }
            //
            if(start < 1 || end < 1){
                cerr<<"Warning: failed to get feature position. "<<line<<endl;
                continue;
            }
            if(! regex_search(type,pat3)){
                continue;
            }
            //------------------------------------
            if(seqID != preSeqID){
                gSeqID = seqID;
                if(flag){
                    if(chrMap.find(seqID) != chrMap.end()){
                        gSeqID = chrMap[seqID];
                    }
                }
            }
            // ref  {gSeqID,start,end,geneID,geneName,strand};
            if(isRef){
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
                //
                AnnoLine anno;
                anno.seqid[FIELDSIZE-1] = '\0';
                anno.geneID[FIELDSIZE-1] = '\0';
                anno.geneName[FIELDSIZE-1] = '\0';
                
                strncpy(anno.seqid,gSeqID.c_str(),FIELDSIZE-1);
                strncpy(anno.geneID,geneID.c_str(),FIELDSIZE-1);
                strncpy(anno.geneName,geneName.c_str(),FIELDSIZE-1);
                anno.strand = strand[0];
                anno.start = start;
                anno.end = end;
                
                out.write((char *)&anno,usize);
                
                ++total;
                
            }else{
                if(seqID != preSeqID){
                    if(chrFind){
                        if(! assAnnoVec.empty()){
                            sort(assAnnoVec.begin(),assAnnoVec.end(),byPos);
                            vector<NdPos> &chrNode = assChrPos[preGid];
                            rgSearch(chrNode,assAnnoVec,usize,out,total);
                        }
                    }
                    assAnnoVec.clear();
                    //
                    if(assChrPos.find(gSeqID) != assChrPos.end()){
                        chrFind = true;
                    }else{
                        chrFind = false;
                    }
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
                
                SimAnno tAnno = {gSeqID,start,end,geneID,geneName,strand};
                assAnnoVec.push_back(tAnno);
            }
            preSeqID = seqID;
            preGid = gSeqID;
        }
    }
    //
    if(! isRef){
        if(chrFind){
            if(! assAnnoVec.empty()){
                sort(assAnnoVec.begin(),assAnnoVec.end(),byPos);
                vector<NdPos> &chrNode = assChrPos[preGid];
                rgSearch(chrNode,assAnnoVec,usize,out,total);
            }
        }
    }
    //
    out.seekp(0,ios::beg);
    out.write((char *)&total,sizeof(int));
    //
    in.close();
    out.close();
}

void simAllGFF(const char *assListFile,const char *gffListFile,string &nodeFile,string &sep,string &annoDir){
    ifstream af(assListFile);
    ifstream gf(gffListFile);
    
    string line;
    stringstream strStream;
    int i = 0;
    map<string,int> assMap;
    while(getline(af,line)){
        assMap.emplace(line,i);
        ++i;
    }
    
    map<string,map<string,vector<NdPos> > > nodeChrPos;
    string refStr = "";
    getAssNode(nodeFile,sep,nodeChrPos,refStr);
    
    string ass,gffFile,chrNameFile;
    while(getline(gf,line)){
        strStream << line;
    
        strStream >> ass;
        strStream >> gffFile;
        strStream >> chrNameFile;
        
        string outFile = annoDir + "/" + to_string(assMap[ass]) + ".anno.bw";
        cout<<"Process -- "<<gffFile<<endl;
        
        if(ass == refStr){
            map<string,vector<NdPos> > assChrPos;
            bool isRef = true;
            simpGFF(gffFile,chrNameFile,assChrPos,isRef,outFile);
        }else{
            if(nodeChrPos.find(ass) == nodeChrPos.end()){
                cout<<"Warning: "<<ass<<" can't be found in the graph. Please check the assembly name in gffList file."<<endl;
            }else{
                bool isRef = false;
                map<string,vector<NdPos> > &assChrPos = nodeChrPos[ass];
                simpGFF(gffFile,chrNameFile,assChrPos,isRef,outFile);
            }
        }
        strStream.clear();
        strStream.str("");
    }
    //
    af.close();
    gf.close();
    
}

// anno.bdx, anno.bw

int getMax(const char *nodeFile){
    ifstream nfh(nodeFile);
    stringstream strStream;
    string line;
    getline(nfh,line);
    int node;
    int i = 0;
    int maxNode = 0;
    while(getline(nfh,line)){
        strStream << line;
        strStream >> node;
        strStream.clear();
        strStream.str("");
        //
        if(i == 0){
            maxNode = node;
        }else{
            if(node > maxNode){
                maxNode = node;
            }
        }
        ++i;
    }
    nfh.close();
    return maxNode;
}

void findGene(const char *nodeFile,const char *assListFile,string &annoDir,string &sep,const char *bdxFile,const char *numFile){
    ifstream af(assListFile);
    ifstream nfh(nodeFile);
    ofstream bdx(bdxFile);
    ofstream mfh(numFile);
    
    int i = 0;
    map<string,int> assMap;
    string line;
    while(getline(af,line)){
        assMap.emplace(line,i);
        ++i;
    }
    af.close();    
    //
    stringstream strStream;
    getline(nfh,line);
    string r_chr;
    int node,start,end;
    string preSeq = "";
    AnnoLine *assAnno = nullptr;
    bool flag = false;
    int maxNode = getMax(nodeFile);
    AnnoDx *dxArr = new AnnoDx[maxNode];
    //AnnoDx zero = {0LL,0};
    long long offset = 0LL;
    int intSize = sizeof(int);
    int uSize = sizeof(AnnoDx);
    string tName = "",t_hap = "",tchr = "";
    int nline = 0;
    //int pSearchStart = 0;
    int searchStart = 0,searchEnd = 0,moveStart = 0;
        
    map<string,NdPos> sChrMap;
    string p_tchr = "";
    bool assFind = false;
    while(getline(nfh,line)){
        strStream << line;
        strStream >> node;
        strStream >> r_chr;
        strStream >> start;
        strStream >> end;
        strStream.clear();
        strStream.str("");
        //
        tName = "",t_hap = "",tchr = "";
        assSplit(r_chr,sep,tName,t_hap,tchr);
        string assStr = tName + sep + t_hap;
        if(assStr != preSeq){
            if(flag){
                delete []assAnno;
                flag = false;
                sChrMap.clear();
            }
            //
            string annoFile = annoDir + "/" + to_string(assMap[assStr]) + ".anno.bw";
            ifstream in(annoFile.c_str());
            if(! in){
                cerr<<"Warning: Annotation of "<<assStr<<" can't be found."<<endl;
                assFind = false;
            }else{
                assFind = true;
            }
            if(assFind){
                in.read((char *)&nline,sizeof(int));
                
                assAnno = new AnnoLine[nline];
                in.read((char *)assAnno,sizeof(AnnoLine)*nline);
                flag = true;
                in.close();
                //
                string preChr = "";
                int tstart = 0;
                for(int k = 0; k < nline; ++k){
                    if(strcmp(assAnno[k].seqid,preChr.c_str()) != 0){
                        if(k > 0){
                            if(sChrMap.find(preChr) == sChrMap.end()){
                                NdPos txpos = {tstart,k-1};
                                sChrMap.emplace(preChr,txpos);
                            }
                            //else cerr<<"Error: items in GFF was not ordered by chromosome/scaffold."<<endl;
                        }
                        //
                        tstart = k;
                        preChr = assAnno[k].seqid;
                    }
                }
                if(sChrMap.find(preChr) == sChrMap.end()){
                    NdPos txpos = {tstart,nline-1};
                    sChrMap.emplace(preChr,txpos);
                }
            }
            //
            //pSearchStart = 0;
            searchStart = 0;
            searchEnd = 0;
            p_tchr = "";
        }
        vector<int> posVec;
        if(assFind){
            bool fir = true;
            bool chrFind = true;
            if(tchr != p_tchr){
                if(sChrMap.find(tchr) != sChrMap.end()){
                    searchStart = sChrMap[tchr].start;
                    searchEnd = sChrMap[tchr].end;
                    chrFind = true;
                    moveStart = searchStart;
                }else{
                    chrFind = false;
                }
            }else{
                if(chrFind){
                    searchStart = moveStart;
                }
            }
            //for(int i = pSearchStart; i < nline; ++i){
            if(chrFind){
                for(int i = searchStart; i <= searchEnd; ++i){    
                    //if(strcmp(assAnno[i].seqid,tchr.c_str()) == 0){
                        if(assAnno[i].end >= start){
                            if(assAnno[i].start <= end){
                                posVec.push_back(i);
                                if(fir){
                                    moveStart = i;
                                    fir = false;
                                }
                            }else{
                                break;
                            }
                        }
                    //}else{
                    //    break;
                    //}
                }
            }
        }
       // pSearchStart = searchStart;
        for(auto tpos : posVec){
            mfh.write((char *)&tpos,intSize);
        }
        
        int tnum = posVec.size();
        //AnnoDx tdx = {offset,tnum};
        dxArr[node-1].offset = offset;
        dxArr[node-1].num = tnum;
        offset += posVec.size() * intSize;
        //
        preSeq = assStr;
        p_tchr = tchr;
    }
    bdx.write((char *)dxArr,uSize*maxNode);
    delete []dxArr;
    if(flag){
        delete []assAnno;
    }
    //
    nfh.close();
    mfh.close();
    bdx.close();
}


void ndg_usage(){
    cout<<"Usage: GraphAnno nodeGene --gffList <gff_list_file> --upDir <upload_dir>"<<endl;
    cout<<"--gffList    <File>     Input file containing three columns separated by whitespace(assemlby name, absolute file path of GFF file, absolute file path of translation table file)."
                                   "Translation table is used to make the chromosome/contig names in GFF file matching that in graph file. This file contain two columns "
                                   "separated by whitespace (chromosome/contig name in GFF file, chromosome/contig name in graph file). If chromosome name transformation "
                                   "is not needed the file path of the translation table can be set to NA in the gffList file."<<endl;
                                   
    cout<<"--upDir      <Dir>      'upload' directory which including files generated by 'gfa2view' or 'vrpg_preprocess.py'."<<endl;
}

int ndg_main(int argc,char **argv){
    string outDir;
    char *gffListFile = nullptr;
    if(argc == 1){
        ndg_usage();
        return 1;
    }
    for(int i = 1; i < argc; ++i){
        if(strcmp(argv[i],"--upDir") == 0 || strcmp(argv[i],"-upDir") == 0){
            i++;
            if(i == argc){
                cerr<<"Error: --upDir is used but there is no value"<<endl;
                return 1;
            }
            outDir = argv[i];
        }else if(strcmp(argv[i],"--gffList") == 0 || strcmp(argv[i],"-gffList") == 0){
            i++;
            if(i == argc){
                cerr<<"Error: --gffList is used but there is no value"<<endl;
                return 1;
            }
            gffListFile = argv[i];
        }else if(strcmp(argv[i],"--help") == 0 || strcmp(argv[i],"-help") == 0){
            ndg_usage();
            return 1;
        }else{
            cerr<<"Error: undefined parameter"<<endl;
            ndg_usage();
            return 1;
        }
    }
    
    if(gffListFile == nullptr){
        cerr<<"Error: --gffList is required."<<endl;
        return 1;
    }
    string sepFile = outDir + "/sep.info";
    string sep = getSep(sepFile);
    string annoDir = outDir + "/anno";
    if(! opendir(annoDir.c_str())){
        if(mkdir(annoDir.c_str(),S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)){
            cerr<<"Error: failed to create directory. "<<annoDir<<endl;
            return 1;
        }
    }
    string nodeFile = outDir + "/node.info";
    string assFile = outDir + "/node.ass.list";
    string assFile2 = outDir + "/ass.list";
    
    if(access(assFile.c_str(),F_OK) != 0){
        assFile =  assFile2;
    }
    string annoDxFile = outDir + "/anno.bdx";
    string annoNumFile = outDir + "/anno.num";
    //
    simAllGFF(assFile.c_str(),gffListFile,nodeFile,sep,annoDir);
    findGene(nodeFile.c_str(),assFile.c_str(),annoDir,sep,annoDxFile.c_str(),annoNumFile.c_str());
    return 0;
}









