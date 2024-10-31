
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <stdlib.h>
#include <string.h>
#include "vgraph.h"

using namespace std;

//
//
void readEdge(char *edgeFile,unordered_map<NodeType,vector<NodeType> > &oedge,unordered_map<NodeType,vector<NodeType> > &iedge){
    ifstream in(edgeFile);
    if(! in){
        cerr<<"Error: file open failed. "<<edgeFile<<endl;
        exit(1);
    }
    string strLine;
    getline(in,strLine);
    stringstream strStream;
    unordered_map<NodeType,vector<NodeType> >::iterator it;
    while(getline(in,strLine)){
        strStream << strLine;
        NodeType node1,node2;
        char sign1,sign2;
        strStream >> node1;
        strStream >> node2;
        strStream >> sign1;
        strStream >> sign2;
        
        it = oedge.find(node1);
        if(it != oedge.end()){
            (it->second).push_back(node2);
        }else{
            vector<NodeType> tvec{node2};
            oedge.emplace(node1,tvec);
        }
        
        it = iedge.find(node2);
        if(it != iedge.end()){
            (it->second).push_back(node1);
        }else{
            vector<NodeType> tvec{node1};
            iedge.emplace(node2,tvec);
        }
        
        strStream.clear();
        strStream.str("");
    }
    in.close();
}

//
void readRnode(char *rndFile,char *rndDxFile,unordered_map<NodeType,ANode> &refNdMap){    
    ifstream dfh(rndDxFile);
    if(! dfh){
        cerr<<"Error: file open failed. "<<rndDxFile<<endl;
        exit(1);
    }
    
    ifstream nfh(rndFile);
    if(! nfh){
        cerr<<"Error: file open failed. "<<rndFile<<endl;
        exit(1);
    }
    
    int nchr = 0;
    int intSize = sizeof(int);
    dfh.read((char *)&nchr,intSize);

    ChrRange crRange;
    int crSize = sizeof(ChrRange);
    int refChr = -1;
    //
    vector<int> chrCnVec;
    chrCnVec.reserve(nchr);
    
    vector<int> chrVec;
    chrVec.reserve(nchr);
    
    vector<int> chrNdVec;
    chrNdVec.reserve(nchr);
    for(int x = 0; x < nchr; ++x){
        dfh.read((char *)&refChr,intSize);      
        dfh.read((char *)&crRange,crSize);
        chrVec.push_back(refChr);
        chrCnVec.push_back(crRange.ranNum);
    }
    
    int oneSize = sizeof(OneRange);
    for(int nc : chrCnVec){
        int n = 0;
        for(int k = 0; k < nc; ++k){
            OneRange tRange;
            dfh.read((char *)&tRange,oneSize);
            n += tRange.ranNum;
        }
        chrNdVec.push_back(n);
    }
    //
    int r_node,r_start,r_end;
    int j = 0;
    for(int nd : chrNdVec){
        int tchr = chrVec[j];
        for(int k = 0; k < nd; ++k){
            nfh.read((char *)&r_node,intSize);
            nfh.read((char *)&r_start,intSize);
            nfh.read((char *)&r_end,intSize);
            ANode tn = {r_start,r_end,tchr};
            refNdMap.emplace(r_node,tn);
        }
        ++j;
    }
    //
    dfh.close();
    nfh.close();
}

void readChr(char *rChrFile,vector<string> &chrVec){
    ifstream cfh(rChrFile);
    if(! cfh){
        cerr<<"Error: file open failed. "<<rChrFile<<endl;
        exit(1);
    }
    string line;
    //int pos = 0;
    while(getline(cfh,line)){
        int tpos = line.find("\t");
        string tchr = line.substr(0,tpos);
        chrVec.push_back(tchr);
    }
    cfh.close();
}

void searchPos(unordered_map<NodeType,ANode> &refNdMap,unordered_map<NodeType,vector<NodeType> > &oedge,unordered_map<NodeType,vector<NodeType> > &iedge,vector<NodeType> &selNode,map<int,RanPos> &chrRanMap){    
    unordered_map<NodeType,vector<NodeType> >::iterator it;
    int queryDep = 500;
    
    for(NodeType tnode : selNode){
        it = oedge.find(tnode);
        vector<NodeType> tNref;
        vector<int> deep;
        bool findRef = false;
        
        set<NodeType> travNode;
        travNode.insert(tnode);
        if(it != oedge.end()){
            for(NodeType &o_node : it->second){
                if(travNode.find(o_node) == travNode.end()){
                    if(refNdMap.find(o_node) != refNdMap.end()){
                        ANode &tRef = refNdMap[o_node];
                        if(chrRanMap.find(tRef.achr) != chrRanMap.end()){
                            RanPos &tRan = chrRanMap[tRef.achr];
                            if(tRan.start < tRef.start){
                                tRan.pend = tRef.pend;
                            }else{
                                tRan.start = tRef.start;
                            }
                        }else{
                            RanPos rp = {tRef.start,tRef.pend};
                            chrRanMap.emplace(tRef.achr,rp);
                        }
                        findRef = true;
                        break;
                    }
                    tNref.push_back(o_node);
                    deep.push_back(0);
                    //
                    travNode.insert(tnode);
                }
            }
        }
        
        if(findRef){
            break;
        }
        it = iedge.find(tnode);
        if(it != iedge.end()){
            for(NodeType &o_node : it->second){
                if(travNode.find(o_node) == travNode.end()){
                    if(refNdMap.find(o_node) != refNdMap.end()){
                        ANode &tRef = refNdMap[o_node];
                        if(chrRanMap.find(tRef.achr) != chrRanMap.end()){
                            RanPos &tRan = chrRanMap[tRef.achr];
                            if(tRan.start < tRef.start){
                                tRan.pend = tRef.pend;
                            }else{
                                tRan.start = tRef.start;
                            }
                        }else{
                            RanPos rp = {tRef.start,tRef.pend};
                            chrRanMap.emplace(tRef.achr,rp);
                        }
                        findRef = true;
                        break;
                    }
                    tNref.push_back(o_node);
                    deep.push_back(0);
                    //
                    travNode.insert(tnode);
                }
            }
        }
        if(findRef){
            break;
        }    
        if(! tNref.empty()){
            size_t i = 0;
            while(i < tNref.size()){
                if(deep[i] > queryDep){
                    break;   
                }
                
                it = oedge.find(tNref[i]);
                for(NodeType &o_node : it->second){
                    if(travNode.find(o_node) == travNode.end()){
                        if(refNdMap.find(o_node) != refNdMap.end()){
                            ANode &tRef = refNdMap[o_node];
                            if(chrRanMap.find(tRef.achr) != chrRanMap.end()){
                                RanPos &tRan = chrRanMap[tRef.achr];
                                if(tRan.start < tRef.start){
                                    tRan.pend = tRef.pend;
                                }else{
                                    tRan.start = tRef.start;
                                }
                            }else{
                                RanPos rp = {tRef.start,tRef.pend};
                                chrRanMap.emplace(tRef.achr,rp);
                            }
                            findRef = true;
                            break;
                        }
                        tNref.push_back(o_node);
                        deep.push_back(deep[i]+1);
                        //
                        travNode.insert(tnode);
                    }
                }
                
                if(findRef){
                    break;
                }
                
                it = iedge.find(tNref[i]);
                for(NodeType &o_node : it->second){
                    if(travNode.find(o_node) == travNode.end()){
                        if(refNdMap.find(o_node) != refNdMap.end()){
                            ANode &tRef = refNdMap[o_node];
                            if(chrRanMap.find(tRef.achr) != chrRanMap.end()){
                                RanPos &tRan = chrRanMap[tRef.achr];
                                if(tRan.start < tRef.start){
                                    tRan.pend = tRef.pend;
                                }else{
                                    tRan.start = tRef.start;
                                }
                            }else{
                                RanPos rp = {tRef.start,tRef.pend};
                                chrRanMap.emplace(tRef.achr,rp);
                            }
                            findRef = true;
                            break;
                        }
                        tNref.push_back(o_node);
                        deep.push_back(deep[i]+1);
                        //
                        travNode.insert(tnode);
                    }
                }
                
                if(findRef){
                    break;
                }
                
                i += 1;
            }
        }
    }
    //
}

void writeBed(vector<string> &chrVec,map<int,RanPos> &chrRanMap,char *rChrFile,ofstream &ofh){
    int n = chrVec.size();
    for(auto &aran : chrRanMap){
        if(aran.first < n){
            ofh<<chrVec[aran.first]<<"\t"<<aran.second.start<<"\t"<<aran.second.pend<<endl;
        }
    }
}

bool pathPos(unordered_map<NodeType,ANode> &refNdMap,string &path,vector<NodeType> &selNode,map<int,RanPos> &chrRanMap){
    int snode = 0, firNode = 0, lastNode = 0;
    //int pmin = 0, pmax = 0;
    unordered_map<NodeType,ANode>::iterator it;
    bool findRef = false;
    bool fir = true;
    
    for(size_t i = 0; i < path.length(); ++i){
        if(path[i] == '>' || path[i] == '<'){
            if(i > 0){
                if(fir){
                    firNode = snode;
                    fir = false;
                }
                it = refNdMap.find(snode);
                if(it != refNdMap.end()){
                    findRef = true;
                    //
                    if(chrRanMap.find((it->second).achr) != chrRanMap.end()){
                        RanPos &tRan = chrRanMap[(it->second).achr];
                        if(tRan.start < (it->second).start){
                            tRan.pend = (it->second).pend;
                        }else{
                            tRan.start = (it->second).start;
                        }
                    }else{
                        RanPos rp = {(it->second).start,(it->second).pend};
                        chrRanMap.emplace((it->second).achr,rp);
                    }
                }
                snode = 0;
            }
        }else{
            snode = snode * 10 + (path[i] - '0');
        }
    }
    //
    if(fir){
        firNode = snode;
    }
    lastNode = snode;
    it = refNdMap.find(snode);
    if(it != refNdMap.end()){
        findRef = true;
        //
        if(chrRanMap.find((it->second).achr) != chrRanMap.end()){
            RanPos &tRan = chrRanMap[(it->second).achr];
            if(tRan.start < (it->second).start){
                tRan.pend = (it->second).pend;
            }else{
                tRan.start = (it->second).start;
            }
        }else{
            RanPos rp = {(it->second).start,(it->second).pend};
            chrRanMap.emplace((it->second).achr,rp);
        }
    }
    //
    if(! findRef){
        selNode.push_back(firNode);
        if(lastNode != firNode){
            selNode.push_back(lastNode);
        }
    }
    return findRef;
    
}

void gaf2rbed(char *rChrFile,char *rndFile,char *rndDxFile,char *edgeFile,char *pathFile,char *outFile){
    ifstream pfh(pathFile);
    if(! pfh){
        cerr<<"Error: file open failed. "<<pathFile<<endl;
        exit(1);
    }
    
    ofstream ofh(outFile);
    if(! ofh){
        cerr<<"Error: file open failed. "<<outFile<<endl;
        exit(1);
    }
    
    bool useEdge = true;
    //
    unordered_map<NodeType,ANode> refNdMap;
    readRnode(rndFile,rndDxFile,refNdMap);
    vector<string> chrVec;
    readChr(rChrFile,chrVec);
    unordered_map<NodeType,vector<NodeType> > oedge,iedge;
    //
    string pLine,pName,path;
    stringstream strStream;
    while(getline(pfh,pLine)){
        strStream << pLine;
        strStream >> pName;
        strStream >> path;
        
        strStream.clear();
        strStream.str("");
        //
        vector<NodeType> selNode;
        map<int,RanPos> chrRanMap;
        bool findRef = pathPos(refNdMap,path,selNode,chrRanMap);
        if(! findRef){
            if(useEdge){
                readEdge(edgeFile,oedge,iedge);
                useEdge = false;
            }
            searchPos(refNdMap,oedge,iedge,selNode,chrRanMap);
        }
        //
        writeBed(chrVec,chrRanMap,rChrFile,ofh);
    }
    //
    pfh.close();
    ofh.close();
}

void ga2bd_usage(){
    cout<<"Usage: gaf2bed --chr <chr.list> --rnode <node.> --dxnode <node_index_file> --edge <edge.info> --path <query.path> --out <out_file>"<<endl;
    cout<<"--chr       <FILE>      chromosome list of reference genome."<<endl;
    cout<<"--rnode     <FILE>      graph nodes from reference genome."<<endl;
    cout<<"--dxnode    <FILE>      index of '--rnode' file."<<endl;
    cout<<"--edge      <FILE>      edge file."<<endl;
    cout<<"--path      <FILE>      query path."<<endl;
    cout<<"--out       <FILE>      output file."<<endl;
}

int main(int argc,char **argv){
    char *rChrFile=nullptr,*rndFile=nullptr,*rndDxFile=nullptr,*edgeFile=nullptr,*pathFile=nullptr,*outFile=nullptr;
    for(int i = 1; i < argc; ++i){
        if(strcmp(argv[i],"--chr") == 0 || strcmp(argv[i],"-chr") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --chr is used but there is no value."<<endl;
                return 1;
            }
            rChrFile = argv[i];
        }else if(strcmp(argv[i],"--rnode") == 0 || strcmp(argv[i],"-rnode") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --rnode is used but there is no value."<<endl;
                return 1;
            }
            rndFile = argv[i]; 
        }else if(strcmp(argv[i],"--dxnode") == 0 || strcmp(argv[i],"-dxnode") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --dxnode is used but there is no value."<<endl;
                return 1;
            }
            rndDxFile = argv[i];
        }else if(strcmp(argv[i],"--edge") == 0 || strcmp(argv[i],"-edge") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --edge is used but there is no value."<<endl;
                return 1;
            }
            edgeFile = argv[i];
        }else if(strcmp(argv[i],"--path") == 0 || strcmp(argv[i],"-path") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --path is used but there is no value."<<endl;
                return 1;
            }
            pathFile = argv[i];
        }else if(strcmp(argv[i],"--out") == 0 || strcmp(argv[i],"-out") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --out is used but there is no value."<<endl;
                return 1;
            }
            outFile = argv[i];
        }else if(strcmp(argv[i],"--help") == 0 || strcmp(argv[i],"-help") == 0){
            ga2bd_usage();
            return 1;
        }else{
            cerr<<"Error: undefined parameter"<<endl;
            ga2bd_usage();
            return 1;
        }
    }
    //
    gaf2rbed(rChrFile,rndFile,rndDxFile,edgeFile,pathFile,outFile);
}










