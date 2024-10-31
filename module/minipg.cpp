

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>
#include <algorithm>
#include <cmath>

#ifdef PYMODULE
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#endif

#include "vgraph.h"

using namespace std;

#ifdef PYMODULE
namespace py = pybind11;
#endif


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

GraphRange::GraphRange(string &t_upDir,int index):upDir(t_upDir),indexFlag(index){
    
    asmFile = upDir + "/asm.list";
    chrFile = upDir + "/chr.list";
    comChrFile = upDir + "/complete.chr.list";
    sepFile = upDir + "/sep.info";
    pathDir = upDir + "/path";
    nodeFile = upDir + "/node.info";
    edgeFile = upDir + "/edge.info";
    
    sep = getSep(sepFile);
    draw_node.reserve(64);
    draw_pos.reserve(64);
    draw_edge.reserve(64);
    genome.reserve(64);
    nnames.reserve(64);
    hnGroup.reserve(64);
    hLinks.reserve(64);
}

void GraphRange::conformEdge(NodeType &node1,NodeType &node2,char mark,unordered_map<NodeType,vector<ENode> > &iedge,unordered_map<NodeType,vector<ENode> > &oedge){
    ENode x = {node2,mark};
    unordered_map<NodeType,vector<ENode> >::iterator it;
    it = oedge.find(node1);
    if(it != oedge.end()){
        (it->second).push_back(x);
    }else{
        
        vector<ENode> eV;
        eV.push_back(x);
        oedge.emplace(node1,eV);
    }
    
    ENode y = {node1,mark};
    it = iedge.find(node2);
    if(it != iedge.end()){
        (it->second).push_back(y);
    }else{
        vector<ENode> eV;
        eV.push_back(y);
        iedge.emplace(node2,eV);
    }
}

void GraphRange::parseEdge(unordered_map<NodeType,vector<ENode> > &iedge,unordered_map<NodeType,vector<ENode> > &oedge){
    ifstream in(edgeFile.c_str());
    if(! in){
        cerr<<"Error: file open failed. "<<edgeFile<<endl;
        return;
        //exit(1);
    }
    string strLine;
    getline(in,strLine);
    stringstream strStream;
    while(getline(in,strLine)){
        strStream << strLine;
        NodeType node1,node2;
        char sign1,sign2;
        strStream >> node1;
        strStream >> node2;
        strStream >> sign1;
        strStream >> sign2;
        
        char mark = getMark(sign1,sign2);
        conformEdge(node1,node2,mark,iedge,oedge);
        
        strStream.clear();
        strStream.str("");
    }
    in.close();
}
//
void GraphRange::parseNode(string &sChr,int sStart,int sEnd,int ex,vector<NodeType> &rangeNode,unordered_map<NodeType,char> &exNode,unordered_map<NodeType,LenAsm> &info,int &realLen){
    int eStart = sStart - ex > 0 ? sStart - ex : 1;
    int eEnd = sEnd + ex;
    
    string nodeLine,asmLine;
    stringstream strStream;
    unordered_map<string,int> asmMap;
    
    ifstream afh(asmFile.c_str());
    if(! afh){
        cerr<<"Error: file open failed. "<<asmFile<<endl;
        return;
        //exit(1);
    }
    int pos = 0;
    while(getline(afh,asmLine)){
        asmMap.emplace(asmLine,pos);
        ++pos;
    }
    afh.close();
    //
    string jAsm = "Jump" + sep + "H";
    string uAsm = "Un" + sep + "H";
    asmMap.emplace(jAsm,pos);
    ++pos;
    asmMap.emplace(uAsm,pos);
    //
    
    ifstream in(nodeFile.c_str());
    if(! in){
        cerr<<"Error: file open failed. "<<nodeFile<<endl;
        return;
        //exit(1);
    }
    
    getline(in,nodeLine);
    while(getline(in,nodeLine)){
        strStream << nodeLine;
        NodeType r_node;
        string r_chr;
        int r_start,r_end,r_len,r_ref;
        strStream >> r_node;
        strStream >> r_chr;
        strStream >> r_start;
        strStream >> r_end;
        strStream >> r_len;
        strStream >> r_ref;
        
        string tName = "",t_hap = "",tchr = "";
        asmSplit(r_chr,sep,tName,t_hap,tchr);
        string asmStr = tName + sep + t_hap;
        int asmb = asmMap[asmStr];
        if(r_ref == 0){
            if(tchr == sChr){
                if(r_start <= eStart){
                    if(r_end >= eStart){
                        //exNode.insert(r_node);
                        LenAsm tx = {r_len,asmb};
                        info.emplace(r_node,tx);
                        //asmDict.emplace(r_node,tName);
                        if(r_end >= sStart){
                            rangeNode.push_back(r_node);
                            realLen += r_len;
                            //
                            exNode.emplace(r_node,'0');
                        }else{
                            exNode.emplace(r_node,'2');
                        }
                    }
                }else{
                    if(r_start <= eEnd){
                        //exNode.insert(r_node);
                        //info.emplace(r_node,r_len);
                        //asmDict.emplace(r_node,tName);
                        LenAsm tx = {r_len,asmb};
                        info.emplace(r_node,tx);
                        if(r_start <= sStart){
                            if(r_end >= sStart){
                                rangeNode.push_back(r_node);
                                realLen += r_len;
                                //
                                exNode.emplace(r_node,'0');
                            }else{
                                exNode.emplace(r_node,'2');
                            }
                        }else{
                            if(r_start <= sEnd){
                                rangeNode.push_back(r_node);
                                realLen += r_len;
                                //
                                exNode.emplace(r_node,'0');
                            }else{
                                exNode.emplace(r_node,'1');
                            }
                        }
                        
                    }
                }
            }
        }else{
            LenAsm tx = {r_len,asmb};
            info.emplace(r_node,tx);
        }
        
        strStream.clear();
        strStream.str("");
    }
    
    in.close();
}

// [posStart,posEnd)
string GraphRange::subCigar(string &rgCigar,int ndStart,int ndEnd){
    int n = rgCigar.size();
    int posStart = 0, posEnd = 0;
    int count = 0;
    string outStr = "";
    for(int i = 0; i < n; ++i){
        if(rgCigar[i] == '>' || rgCigar[i] == '<'){
            if(ndStart == count){
                posStart = i;
            }else if(ndEnd == count){
                posEnd = i;
                break;
            }            
            ++count;
        }
    }
    //
    if(posEnd == 0){
        outStr = rgCigar.substr(posStart);
    }else{
        outStr = rgCigar.substr(posStart,posEnd-posStart);
    }
    return outStr;
}    

/*
    When option 'Simplify all nodes' is taken new edges that link two reference nodes may be created.
    Variable 'tIn' is used to record the status of all nodes.
*/
void GraphRange::eAsmFind(bool visCigar,bool refSim,vector<char> &orient,vector<NodeType> &nodes,vector<int> &qPosVec,string &rgCigar,string &rgName,map<NEdge,int> &r_edge_dict,unordered_map<NodeType,Nid> &nid_dict,unordered_set<int> &ndGroup){
    char preOri = '\0';
    int preIn = 0;
    NodeType preNode = 0;
    int tIn = 0;
    //
    int ndStart = -1, ndEnd = 0;
    string tpath = "";
    
    for(size_t j = 0; j < orient.size(); ++j){
        //if(refSkip.find(nodes[j]) != refSkip.end()){
        //    continue;
        //}
        if(nid_dict.find(nodes[j]) != nid_dict.end()){
            if(ndGroup.find(nodes[j]) == ndGroup.end()){
                hnGroup.push_back(nid_dict[nodes[j]].gNum);
                ndGroup.insert(nodes[j]);
            }
            if(preIn > 0){
                char mark = getMark(preOri,orient[j],'>');
                NEdge sym = {preNode,nodes[j],mark};               
                if(r_edge_dict.find(sym) != r_edge_dict.end()){                   
                    //if(edGroup.find(sym) == edGroup.end()){
                        hDir.push_back('1');
                        //
                        if(mark == '2'){
                            if(r_edge_dict[sym] == 0){
                                mark = '1';
                            }
                        }                       
                        switch(mark){
                            case '1':
                                hLinks.push_back(to_string(nid_dict[preNode].e_nid) + "_" + to_string(nid_dict[nodes[j]].s_nid) + "_" + string(1,mark));                          
                                break;
                            case '2':
                                hLinks.push_back(to_string(nid_dict[preNode].e_nid) + "_" + to_string(nid_dict[nodes[j]].s_nid) + "_" + string(1,mark));
                                break;
                            case '3':
                                hLinks.push_back(to_string(nid_dict[preNode].e_nid) + "_" + to_string(nid_dict[nodes[j]].e_nid) + "_" + string(1,mark));
                                break;
                            case '4':
                                hLinks.push_back(to_string(nid_dict[preNode].s_nid) + "_" + to_string(nid_dict[nodes[j]].s_nid) + "_" + string(1,mark));
                                break;
                            default:
                                hLinks.push_back(to_string(nid_dict[preNode].s_nid) + "_" + to_string(nid_dict[nodes[j]].e_nid) + "_" + string(1,mark));
                        }
                        //
                        //edGroup.insert(sym);
                    //}
                }else{
                    char rmark = revMark(mark);
                    NEdge sym = {nodes[j],preNode,rmark};
                    
                    if(r_edge_dict.find(sym) != r_edge_dict.end()){
                        //if(edGroup.find(sym) == edGroup.end()){
                            hDir.push_back('0');
                            //
                            if(rmark == '2'){
                                if(r_edge_dict[sym] == 0){
                                    rmark = '1';   
                                }
                            }
                            switch(rmark){
                                case '1':
                                    hLinks.push_back(to_string(nid_dict[nodes[j]].e_nid) + "_" + to_string(nid_dict[preNode].s_nid) + "_" + string(1,rmark));
                                    break;
                                case '2':
                                    hLinks.push_back(to_string(nid_dict[nodes[j]].e_nid) + "_" + to_string(nid_dict[preNode].s_nid) + "_" + string(1,rmark));
                                    break;
                                case '3':
                                    hLinks.push_back(to_string(nid_dict[nodes[j]].e_nid) + "_" + to_string(nid_dict[preNode].e_nid) + "_" + string(1,rmark));
                                    break;
                                case '4':
                                    hLinks.push_back(to_string(nid_dict[nodes[j]].s_nid) + "_" + to_string(nid_dict[preNode].s_nid) + "_" + string(1,rmark));
                                    break;
                                default:
                                    hLinks.push_back(to_string(nid_dict[nodes[j]].s_nid) + "_" + to_string(nid_dict[preNode].e_nid) + "_" + string(1,rmark));
                            }
                            //
                            //edGroup.insert(sym);
                        //}
                    }
                    
                }
            }
            //
            if(visCigar){
                if(tIn){
                    tpath += string(1,orient[j]) + to_string(nodes[j]);
                }else{
                    tpath = string(1,orient[j]) + to_string(nodes[j]);
                    ndStart = j;                   
                }
            }
            //
            preOri = orient[j];
            preNode = nodes[j];
            preIn = 1;
            tIn = 1;
        }else{
            if(! refSim){
                preOri = orient[j];
                preNode = nodes[j];
                preIn = 0;
            }
            
            if(visCigar){
                if(tIn > 0){
                    ndEnd = j;
                    qChr.push_back(rgName);
                    qStart.push_back(qPosVec[ndStart]);
                    qEnd.push_back(qPosVec[ndEnd-1]);
                    qPath.push_back(tpath);
                    // * -> not contain CIGAR
                    if(rgCigar != "*"){
                        //[ndStart,ndEnd)
                        string tCigar = subCigar(rgCigar,ndStart,ndEnd);                       
                        qCigar.push_back(tCigar);
                    }
                }
            }
            tIn = 0;
        }
        
    }
    //
    if(visCigar){
        if(tIn > 0){
            ndEnd = nodes.size();
            qChr.push_back(rgName);
            qStart.push_back(qPosVec[ndStart]);
            qEnd.push_back(qPosVec[ndEnd-1]);
            qPath.push_back(tpath);              
            if(rgCigar != "*"){
                string tCigar = subCigar(rgCigar,ndStart,ndEnd);
                qCigar.push_back(tCigar);
            }
        }
    }
}

void GraphRange::eAsmFind2(bool refSim,int asmCode,vector<char> &orient,vector<NodeType> &nodes,unordered_map<NodeType,Nid> &nid_dict,unordered_set<int> &ndGroup,map<NEdge,int> &h_edge_dict){
    char preOri = '\0';
    int preIn = 0;
    NodeType preNode = 0;
    int tIn = 0;
    for(size_t j = 0; j < orient.size(); ++j){
        //if(refSkip.find(nodes[j]) != refSkip.end()){
        //    continue;
        //}
        if(nid_dict.find(nodes[j]) != nid_dict.end()){
            tIn = 1;
            if(ndGroup.find(nodes[j]) == ndGroup.end()){
                hnGroup.push_back(nid_dict[nodes[j]].gNum);
                ndGroup.insert(nodes[j]);
            }
            if(preIn > 0){
                char mark = getMark(preOri,orient[j],'>');
                NEdge sym = {preNode,nodes[j],mark};
                if(h_edge_dict.find(sym) == h_edge_dict.end()){
                    h_edge_dict.emplace(sym,asmCode);
                }
            }
            //
            preOri = orient[j];
            preIn = tIn;
            preNode = nodes[j]; 
        }else{
            if(! refSim){
                tIn = 0;
                preOri = orient[j];
                preIn = tIn;
                preNode = nodes[j];
            }
        }
        
    }
}

bool GraphRange::f_cigar2pos(int qStart,int rStart,unordered_map<NodeType,LenAsm> &info,string &cigar,vector<int> &nodes,vector<string> &qCigarVec,vector<int> &qPosVec){    
    int tvalue = 0;
    int rMove = 0, qMove = 0;
    int qPreMove = 0;
    int nextPos = 0;
    // firLen,rStart
    int firLen;  
    if(info.find(nodes[0]) != info.end()){
        firLen = info[nodes[0]].len;
    }else{
        return false;
    }
    int cumLen = firLen - rStart;
    //
    qPosVec.push_back(qStart);
    string tCigar = "";
    int i = 0;
    int tlen = 0;
    int preCum = cumLen;
    int n = nodes.size();
    for(char x : cigar){
        if(x >= '0' && x <= '9'){
            tvalue = tvalue * 10 + (x - '0');
        }else{
            switch(x){
                case 'M':
                    qMove += tvalue;
                    rMove += tvalue;
                    break;
                case 'X':
                    qMove += tvalue;
                    rMove += tvalue;
                    break;
                case '=':
                    qMove += tvalue;
                    rMove += tvalue;
                    break;
                case 'I':
                    qMove += tvalue;
                    break;
                case 'S':
                    qMove += tvalue;
                    break;
                case 'D':
                    rMove += tvalue;
                    break;
                case 'N':
                    rMove += tvalue;
                    break;
            }
            if(rMove < cumLen){
                tCigar += to_string(tvalue) + string(1,x);
            }else if(rMove == cumLen){
                tCigar += to_string(tvalue) + string(1,x);
                qCigarVec.push_back(tCigar);
                tCigar = "";
                //
                ++i;
                if(i < n){
                    if(info.find(nodes[i]) != info.end()){
                        tlen = info[nodes[i]].len;
                    }else{
                        return false;
                    }
                    cumLen += tlen;
                    //
                    nextPos = qStart + qMove;
                    qPosVec.push_back(nextPos);                    
                }
            }else{
                // |-------|-------|---------|
                //      |----------------------|
                int nd = tvalue - (rMove - cumLen);
                //if(nd > 0){
                    tCigar += to_string(nd) + string(1,x);
                //}
                qCigarVec.push_back(tCigar);
                tCigar = "";
                // for last node cumLen >= rMove
                preCum = cumLen;
                ++i;
                //
                int qd = nd;
                if(i < n){
                    if(qMove > qPreMove){
                        nextPos = qStart + qMove + nd;
                    }else{
                        nextPos = qStart + qMove;
                    }
                    qPosVec.push_back(nextPos);
                    //
                    if(info.find(nodes[i]) != info.end()){
                        tlen = info[nodes[i]].len;
                    }else{
                        return false;
                    }
                    cumLen += tlen;
                    qd += tlen;
                }
                //if(i < nodes.size(){
                    while(1){
                        if(rMove < cumLen){
                            //int nd = tlen - (cumLen - rMove);
                            int nd = rMove - preCum;
                            tCigar = to_string(nd) + string(1,x);
                            break;
                        }else if(rMove == cumLen){
                            tCigar = to_string(tlen) + string(1,x);
                            qCigarVec.push_back(tCigar);
                            tCigar = "";
                            //
                            preCum = cumLen;
                            ++i;
                            if(i < n){
                                nextPos = qStart + qMove;
                                qPosVec.push_back(nextPos);
                                //
                                if(info.find(nodes[i]) != info.end()){
                                    tlen = info[nodes[i]].len;
                                }else{
                                    return false;
                                }
                                cumLen += tlen;
                            }
                            break;
                        }else{
                            tCigar = to_string(tlen) + string(1,x);
                            qCigarVec.push_back(tCigar);
                            tCigar = "";
                            //
                            preCum = cumLen;
                            ++i;
                            if(i < n){
                                if(qMove > qPreMove){
                                    nextPos = qStart + qPreMove + qd;
                                }else{
                                    nextPos = qStart + qMove;
                                }
                                qPosVec.push_back(nextPos);
                                //
                                if(info.find(nodes[i]) != info.end()){
                                    tlen = info[nodes[i]].len;
                                }else{
                                    return false;
                                }
                                cumLen += tlen;
                                qd += tlen;
                            }
                        }
                    }
                //}
            }
            //
            tvalue = 0;
            qPreMove = qMove;
        }
    }
    //
    if(tCigar != ""){
        qCigarVec.push_back(tCigar);
    }
    //
    return true;
}

bool GraphRange::f_path2pos(int rStart,unordered_map<NodeType,LenAsm> &info,vector<int> &nodes,vector<int> &qPosVec){
    int tpos = rStart;
    for(int x : nodes){        
        qPosVec.push_back(tpos);
        if(info.find(x) != info.end()){
            tpos += info[x].len;
        }else{
            return false;
        }
    }
    return true;
}

void GraphRange::hAsmNode(string &taskDir,bool refSim,unordered_map<NodeType,LenAsm> &info,int asmNum,map<NEdge,int> &r_edge_dict,unordered_map<NodeType,Nid> &nid_dict){
    string pathFile = "";
    //
    string rgCigar = "*";
    bool formR = true;
    //
    if(taskDir != ""){
        pathFile = taskDir + "/query.path";
    }else{
        pathFile = pathDir + "/" + to_string(asmNum) + ".path";
        //
        string formFile = upDir + "/form.info";
        if(access(formFile.c_str(),F_OK) == 0){
            formR = false;
        }
    }
    ifstream in(pathFile.c_str());
    if(! in){
        cerr<<"Error: file open failed. "<<pathFile<<endl;
        return;
        //exit(1);
    }
    //   
    string pathLine;
    stringstream strStream;
    unordered_set<int> ndGroup;
    //
    while(getline(in,pathLine)){
        strStream << pathLine;
        string r_chr,path;
        strStream >> r_chr;
        strStream >> path;
        
        //string tName = "",t_hap = "",tchr = "";
        //asmSplit(r_chr,sep,tName,t_hap,tchr);
        //string tAsm = tName + sep + t_hap;
        
        vector<char> orient;
        vector<NodeType> nodes;
        //
        vector<int> qPosVec;
        //
        NodeType snode = 0;
        for(size_t i = 0; i < path.length(); ++i){
            if(path[i] == '>' || path[i] == '<'){
                orient.push_back(path[i]);
                if(i > 0){
                    nodes.push_back(snode);
                    snode = 0;
                }
            }else{
                snode = snode * 10 + (path[i] - '0');
            }
        }
        nodes.push_back(snode);
        // 
        if(formR){
            int qStart,rStart;
            string cigar;
            vector<string> qCigarVec;
            vector<int> qPosVec;
            strStream >> qStart >> rStart >> cigar;          
            bool visCigar = f_cigar2pos(qStart,rStart,info,cigar,nodes,qCigarVec,qPosVec);
            if(visCigar){
                rgCigar = "";
                for(size_t i = 0; i < orient.size(); ++i){
                    rgCigar += string(1,orient[i]) + qCigarVec[i];
                } 
            }
            eAsmFind(visCigar,refSim,orient,nodes,qPosVec,rgCigar,r_chr,r_edge_dict,nid_dict,ndGroup);
        }else{
            char pType;
            int rStart;
            strStream >> pType;
            if(pType == 'W'){
                strStream >>  rStart;
            }else{
                rStart = 0;
            }
            bool visCigar = f_path2pos(rStart,info,nodes,qPosVec);
            //
            eAsmFind(visCigar,refSim,orient,nodes,qPosVec,rgCigar,r_chr,r_edge_dict,nid_dict,ndGroup);
        }
        strStream.clear();
        strStream.str("");
    }   
    //
    in.close();
}

// hEdgeAsm
void GraphRange::hEdge2fig(unordered_map<NodeType,Nid> &nid_dict,map<NEdge,int> &h_edge_dict,map<NEdge,int> &r_edge_dict){
    for(auto &tedge : h_edge_dict){
        if(r_edge_dict.find(tedge.first) != r_edge_dict.end()){                   
            hDir.push_back('1');
            char mark = tedge.first.mark;
            hEdgeAsm.push_back(tedge.second);
            //
            if(mark == '2'){
                if(r_edge_dict[tedge.first] == 0){
                    mark = '1';
                }
            }                       
            switch(mark){
                case '1':
                    hLinks.push_back(to_string(nid_dict[tedge.first.node1].e_nid) + "_" + to_string(nid_dict[tedge.first.node2].s_nid) + "_" + string(1,mark));                          
                    break;
                case '2':
                    hLinks.push_back(to_string(nid_dict[tedge.first.node1].e_nid) + "_" + to_string(nid_dict[tedge.first.node2].s_nid) + "_" + string(1,mark));
                    break;
                case '3':
                    hLinks.push_back(to_string(nid_dict[tedge.first.node1].e_nid) + "_" + to_string(nid_dict[tedge.first.node2].e_nid) + "_" + string(1,mark));
                    break;
                case '4':
                    hLinks.push_back(to_string(nid_dict[tedge.first.node1].s_nid) + "_" + to_string(nid_dict[tedge.first.node2].s_nid) + "_" + string(1,mark));
                    break;
                default:
                    hLinks.push_back(to_string(nid_dict[tedge.first.node1].s_nid) + "_" + to_string(nid_dict[tedge.first.node2].e_nid) + "_" + string(1,mark));
            }
        }else{
            char rmark = revMark(tedge.first.mark);
            NEdge sym = {tedge.first.node2,tedge.first.node1,rmark};
            
            if(r_edge_dict.find(sym) != r_edge_dict.end()){
                    hDir.push_back('0');
                    hEdgeAsm.push_back(tedge.second);
                    //
                    if(rmark == '2'){
                        if(r_edge_dict[sym] == 0){
                            rmark = '1';   
                        }
                    }
                    switch(rmark){
                        case '1':
                            hLinks.push_back(to_string(nid_dict[tedge.first.node2].e_nid) + "_" + to_string(nid_dict[tedge.first.node1].s_nid) + "_" + string(1,rmark));
                            break;
                        case '2':
                            hLinks.push_back(to_string(nid_dict[tedge.first.node2].e_nid) + "_" + to_string(nid_dict[tedge.first.node1].s_nid) + "_" + string(1,rmark));
                            break;
                        case '3':
                            hLinks.push_back(to_string(nid_dict[tedge.first.node2].e_nid) + "_" + to_string(nid_dict[tedge.first.node1].e_nid) + "_" + string(1,rmark));
                            break;
                        case '4':
                            hLinks.push_back(to_string(nid_dict[tedge.first.node2].s_nid) + "_" + to_string(nid_dict[tedge.first.node1].s_nid) + "_" + string(1,rmark));
                            break;
                        default:
                            hLinks.push_back(to_string(nid_dict[tedge.first.node2].s_nid) + "_" + to_string(nid_dict[tedge.first.node1].e_nid) + "_" + string(1,rmark));
                    }
            }
            
        }
    }
}

void GraphRange::hAsmNode2(bool refSim,vector<int> &asmNumVec,map<NEdge,int> &r_edge_dict,unordered_map<NodeType,Nid> &nid_dict){
    int asmCode = 0;
    map<NEdge,int> h_edge_dict;
    //
    for(int asmNum : asmNumVec){
        string pathFile = pathDir + "/" + to_string(asmNum) + ".path";
        ifstream in(pathFile.c_str());
        if(! in){
            cerr<<"Error: file open failed. "<<pathFile<<endl;
            //exit(1);
        }        
        string pathLine;
        stringstream strStream;
        unordered_set<int> ndGroup;
        while(getline(in,pathLine)){
            strStream << pathLine;
            string r_chr,path;
            strStream >> r_chr;
            strStream >> path;
            
            string tName = "",t_hap = "",tchr = "";
            asmSplit(r_chr,sep,tName,t_hap,tchr);
            
            string tAsm = tName + sep + t_hap;
            vector<char> orient;
            vector<NodeType> nodes;
            
            NodeType snode = 0;
            for(size_t i = 0; i < path.length(); ++i){
                if(path[i] == '>' || path[i] == '<'){
                    orient.push_back(path[i]);
                    if(i > 0){
                        nodes.push_back(snode);
                        snode = 0;
                    }
                }else{
                    snode = snode * 10 + (path[i] - '0');
                }
            }
            nodes.push_back(snode);
            //
            eAsmFind2(refSim,asmCode,orient,nodes,nid_dict,ndGroup,h_edge_dict);
            strStream.clear();
            strStream.str("");
        }
        in.close();
        //
        ++asmCode;
    }
    //
    hEdge2fig(nid_dict,h_edge_dict,r_edge_dict);
}

//------------------
void GraphRange::parseIndex(int chrNum,int sStart,int sEnd,unordered_map<NodeType,vector<ENode> > &iedge,unordered_map<NodeType,vector<ENode> > &oedge){
    string eIndexFile = upDir + "/edge.bdx";
    ifstream efh(eIndexFile.c_str());
    if(! efh){
        cerr<<"Error: file open failed. "<<eIndexFile<<endl;
        return;
        //exit(1);
    }
    
    string bEdgeFile = upDir + "/edge.bw";
    ifstream xfh(bEdgeFile.c_str());
    if(! xfh){
        cerr<<"Error: file open failed. "<<bEdgeFile<<endl;
        return;
        //exit(1);
    }
    
    string strLine;
    stringstream strStream;

    int esize = sizeof(CEdge);
    set<NEdge> redgeSet;
    
    int nchr = 0;
    int intSize = sizeof(int);
    //int llSize = sizeof(long long);
    efh.read((char *)&nchr,intSize);
    
    //int refChr = -1,offByte,nLine;
    int refChr = -1;
    bool findChr = false;
    
    ChrRange crRange;
    int crSize = sizeof(ChrRange);
    for(int x = 0; x < nchr; ++x){
        efh.read((char *)&refChr,intSize);
        //efh.read((char *)&offByte,intSize);
        //efh.read((char *)&nLine,intSize);
        efh.read((char *)&crRange,crSize);
        if(refChr == chrNum){
            findChr = true;
            break;
        }
    }
    
    vector<OneRange> acrVec;
    int oneSize = sizeof(OneRange);
    int total = 0;
    bool flag = true;
    long long edOff = 0LL;
    if(findChr){
        efh.clear();
        efh.seekg(crRange.rByte,ios::beg);
        for(int k = 0; k < crRange.ranNum; ++k){
            OneRange aRange;
            efh.read((char *)&aRange,oneSize);
            if(aRange.ranStart <= sStart){
                if(aRange.ranEnd >= sStart){
                    if(flag){
                        edOff = aRange.offByte;
                        flag = false;
                    }
                    total += aRange.ranNum;
                }
            }else{
                if(aRange.ranStart <= sEnd){
                    total += aRange.ranNum;
                }else{
                    break;
                }
            }
        }
    }
    efh.close();
    //
    if(flag){
        cerr<<"Error: edges in the interval can't be found. "<<chrNum<<" : "<<sStart<<" - "<<sEnd<<endl;
        return;
    }
    xfh.seekg(edOff,ios::beg);
    for(int i = 0; i < total; ++i){
        CEdge xedge;
        xfh.read((char *)&xedge,esize);
        NEdge tedge = {xedge.node1,xedge.node2,xedge.mark};
        if(redgeSet.find(tedge) == redgeSet.end()){
            conformEdge(xedge.node1,xedge.node2,xedge.mark,iedge,oedge);
            redgeSet.insert(tedge); 
        } 
    }
    //
    xfh.close();
}

int GraphRange::getChrNum(string &sChr){
    ifstream cfh(chrFile.c_str());
    if(! cfh){
        cerr<<"Error: file open failed. "<<chrFile<<endl;
        //exit(1);
        return -1;
    }
    string chrLine;
    int chrNum = -1,pos = 0;
    while(getline(cfh,chrLine)){
        int tpos = chrLine.find("\t");
        string tchr = chrLine.substr(0,tpos);
//cout<<tchr<<" "<<sChr<<" @@@"<<endl;        
        if(tchr == sChr){
            chrNum = pos;
            break;
        }
        ++pos;
    }
    cfh.close();
    return chrNum;
}

int GraphRange::getAsmNum(string &asmb){
    ifstream in(asmFile.c_str());
    if(! in){
        cerr<<"Error: file open failed. "<<asmFile<<endl;
        exit(1);
    }
    int pos = 0,tpos = -1;
    string asmLine;
    while(getline(in,asmLine)){
        if(asmLine == asmb){
            tpos = pos;
            break;
        }
        ++pos;
    }
    in.close();
    return tpos;
}

void GraphRange::getAsmNum2(set<string> &nameSet,vector<int> &asmNumVec){
    ifstream in(asmFile.c_str());
    if(! in){
        cerr<<"Error: file open failed. "<<asmFile<<endl;
        exit(1);
    }
    int pos = 0;
    string asmLine;
    while(getline(in,asmLine)){
        if(nameSet.find(asmLine) != nameSet.end()){
            asmNumVec.push_back(pos);
        }
        ++pos;
    }
    in.close();
}


void GraphRange::getExNode(int chrNum,int sStart,int sEnd,int ex,vector<NodeType> &rangeNode,unordered_map<NodeType,char> &exNode,unordered_map<NodeType,LenAsm> &info,int &realLen,int &realStart){
    
    //int pos = eIndexFile.find_last_of("/");
    //string outDir = eIndexFile.substr(0,pos + 1);
    
    string rndFile = upDir + "/node.ref.bw";
    string nrdFile = upDir + "/node.nonref.bw";
    string mgDxFile = upDir + "/node.merge.bdx";
    
    ifstream mgfh(mgDxFile.c_str());
    ifstream rnfh(rndFile.c_str());
    ifstream ndfh(nrdFile.c_str());
 
    //
    int nchr = 0;
    int intSize = sizeof(int);

    mgfh.read((char *)&nchr,intSize);

    int refChr = -1;
    bool findChr = false;

    ChrRange crRange;
    int crSize = sizeof(ChrRange);
    for(int x = 0; x < nchr; ++x){
        mgfh.read((char *)&refChr,intSize);
        
        mgfh.read((char *)&crRange,crSize);      
        if(refChr == chrNum){
            findChr = true;
            break;
        }
    }
    
    //
    long long refOff = 0LL,nRefOff = 0LL;
    int rt = 0, nrt = 0;
    bool flag = true, nflag = true;
    int oneSize = sizeof(OneRange);
    int secSize = sizeof(EdRang);
    int eStart = sStart - ex > 0 ? sStart - ex : 1;
    int eEnd = sEnd + ex;
    
    
    if(findChr){
        mgfh.clear();
        mgfh.seekg(crRange.rByte,ios::beg);
        for(int k = 0; k < crRange.ranNum; ++k){
            OneRange tRange;
            EdRang sRange;
            mgfh.read((char *)&tRange,oneSize);
            mgfh.read((char *)&sRange,secSize);
            //
            if(tRange.ranStart <= sStart){
                if(tRange.ranEnd >= sStart){
                    if(nflag){
                        nRefOff = sRange.rOffsize;
                        nflag = false;
                    }
                    nrt += sRange.rCount;
                }
            }else{
                if(tRange.ranStart <= sEnd){
                    if(nflag){
                        nRefOff = sRange.rOffsize;
                        nflag = false;
                    }
                    nrt += sRange.rCount;
                }
            }
            //
            if(tRange.ranStart <= eStart){
                if(tRange.ranEnd >= eStart){
                    if(flag){
                        refOff = tRange.offByte;
                        flag = false;
                    }
                    rt += tRange.ranNum;
                }
                
            }else{
                if(tRange.ranStart <= eEnd){
                    if(flag){
                        refOff = tRange.offByte;
                        flag = false;
                    }
                    rt += tRange.ranNum;
                }else{
                    break;
                }
                
            }
        }
    }
    mgfh.close();
    if(flag){
        cerr<<"Error: nodes can't be found in the interval. "<<chrNum<<" : "<<sStart<<" - "<<sEnd<<endl;
        return;
        //exit(1);
    }
    // ref node consecutive
    rnfh.seekg(refOff,ios::beg);
    for(int m = 0; m < rt; ++m){
        int r_node,r_start,r_end,r_len;
        rnfh.read((char *)&r_node,intSize);
        rnfh.read((char *)&r_start,intSize);
        rnfh.read((char *)&r_end,intSize);
        
        r_len = r_end - r_start + 1;
        LenAsm tx = {r_len,0};
        if(info.find(r_node) == info.end()){
            info.emplace(r_node,tx);
        }
        //
        if(r_start <= eStart){
            if(r_end >= eStart){
                //
                if(r_end >= sStart){
                    rangeNode.push_back(r_node);
                    realLen += r_len;
                    //
                    realStart = r_start;
                    //
                    exNode.emplace(r_node,'0');
                }else{
                    exNode.emplace(r_node,'2');
                }
            }
        }else{
            if(r_start <= eEnd){
                //
                if(r_start <= sStart){
                    if(r_end >= sStart){
                        rangeNode.push_back(r_node);
                        realLen += r_len;
                        //
                        realStart = r_start;
                        //
                        exNode.emplace(r_node,'0');
                    }else{
                        exNode.emplace(r_node,'2');
                    }
                }else{
                    if(r_start <= sEnd){
                        rangeNode.push_back(r_node);
                        realLen += r_len;
                        //
                        exNode.emplace(r_node,'0');
                    }else{
                        exNode.emplace(r_node,'1');
                    }
                }
            }else{
                break;
            }
            
        }
    }
    rnfh.close();
    //
    ndfh.seekg(nRefOff,ios::beg);
    for(int n = 0; n < nrt; ++n){
        int t_node,t_len,t_asm;
        ndfh.read((char *)&t_node,intSize);
        ndfh.read((char *)&t_len,intSize);
        ndfh.read((char *)&t_asm,intSize);
        //
        LenAsm tx = {t_len,t_asm};
        if(info.find(t_node) == info.end()){
            info.emplace(t_node,tx);
        }
    }

    ndfh.close();

}

void GraphRange::queryDbPath(bool formR,int asmNum,int chrNum,int sStart,int sEnd,unordered_map<NodeType,LenAsm> &info,vector<vector<char> > &oriMulti,vector<vector<int> > &nodeMulti,vector<vector<int> > &qPosMulti, vector<string> &cigarMulti,vector<string> &nameMulti){
    string bPathFile = pathDir + "/" + to_string(asmNum) + ".path.bw";
    string dxPathFile = pathDir + "/" + to_string(asmNum) + ".path.bdx";    
    
    ifstream pfh(bPathFile);
    if(! pfh){
        cerr<<"Error: file open failed. "<<bPathFile<<endl;
        exit(1);
    }
    ifstream xfh(dxPathFile);
    if(! xfh){
        cerr<<"Error: file open failed. "<<dxPathFile<<endl;
        exit(1);
    }
    
    int nchr = 0;
    int intSize = sizeof(int);
    xfh.read((char *)&nchr,intSize);

    int refChr;
    ChrRange crRange;
    int crSize = sizeof(ChrRange);
    //
    bool flag = false;
    for(int x = 0; x < nchr; ++x){
        xfh.read((char *)&refChr,intSize);
        xfh.read((char *)&crRange,crSize);
        if(refChr == chrNum){
            flag = true;
            break;
        }
    }
    
    int oneSize = sizeof(OneRange);
    long long tOff;
    int total = 0;
    bool tflag = true;
    if(flag){
        xfh.clear();
        xfh.seekg(crRange.rByte,ios::beg);
        for(int i = 0; i < crRange.ranNum; ++i){
            OneRange tRange;
            xfh.read((char *)&tRange,oneSize);
            //
            if(tRange.ranStart <= sStart){
                if(tRange.ranEnd >= sStart){
                    if(tflag){
                        tOff = tRange.offByte;
                        tflag = false;
                    }
                    total += tRange.ranNum;
                }
            }else{
                if(tRange.ranStart <= sEnd){
                    if(tflag){
                        tOff = tRange.offByte;
                        tflag = false;
                    }
                    total += tRange.ranNum;
                }else{
                    break;
                }
            }
        }
    }else{
        cerr<<"Error: chromosome can't be found (path search). "<<endl;
        //exit(1);
    }
    
    int lndSize = sizeof(LagNode);
    if(total > 0){
        pfh.seekg(tOff,ios::beg);
        int tnum = 0;
        
        vector<char> oriVec;
        vector<int> ndVec;
        //
        string nameFile = pathDir + "/" + to_string(asmNum) + ".name";
        ifstream mfh(nameFile);
        if(! mfh){
            cerr<<"Error: file open failed. "<<nameFile<<endl;
            exit(1);
        }
        vector<int> qPosVec;
        string cigarStr = "";
        string nameLine;
        vector<string> allName;
        while(getline(mfh,nameLine)){
            allName.push_back(nameLine);  
        }
        string firName = "";
        mfh.close();
        //
       
        if(formR){
            while(1){
                int pathNum;
                int firNode;
                char firOri;
                int lagNum;   
                
                pfh.read((char *)&pathNum,intSize);
                pfh.read((char *)&firNode,intSize);
                pfh.read(&firOri,1);
                pfh.read((char *)&lagNum,intSize);

                if(firOri == '>' || firOri == '<'){
                    if(! oriVec.empty()){
                        oriMulti.push_back(oriVec);
                        nodeMulti.push_back(ndVec);
                        qPosMulti.push_back(qPosVec);
                        cigarMulti.push_back(cigarStr);
                        //
                        nameMulti.push_back(firName);

                        oriVec.clear();
                        ndVec.clear();
                        qPosVec.clear();
                        cigarStr = "";
                    }
                    firName = allName[pathNum];
                }else{
                    if(firOri == '1'){
                        firOri = '>';
                    }else{
                        firOri = '<';
                    }
                }
                
                oriVec.push_back(firOri);
                ndVec.push_back(firNode);
                for(int i = 0; i < lagNum; ++i){
                    LagNode tLag;
                    pfh.read((char *)&tLag,lndSize);
                    ndVec.push_back(firNode + tLag.diff);
                    oriVec.push_back(tLag.ori);
                  
                }
             
                //
                int tpos;
                for(int j = 0; j <= lagNum; ++j){
                    pfh.read((char *)&tpos,intSize);
                    qPosVec.push_back(tpos);
                   
                }

                int cigarSize;
                pfh.read((char *)&cigarSize,intSize);
                char *tCigar = new char[cigarSize + 1];
                pfh.read(tCigar,cigarSize);
                tCigar[cigarSize] = '\0';
                if(cigarStr == ""){
                    cigarStr = tCigar;
                }else{
                    string tStr = tCigar;
                    cigarStr += tStr;
                }
                delete []tCigar;
                //
                tnum += (lagNum + 1);
           
                if(tnum == total){
                    break;
                }
            }
            if(! oriVec.empty()){
                oriMulti.push_back(oriVec);
                nodeMulti.push_back(ndVec);
                qPosMulti.push_back(qPosVec);
                cigarMulti.push_back(cigarStr);
                //
                nameMulti.push_back(firName);

                oriVec.clear();
                ndVec.clear();
                qPosVec.clear();
                cigarStr = "";
            }
        }else{
            while(1){
                int pathNum;
                int firNode;
                char firOri;
                int lagNum;   
                
                pfh.read((char *)&pathNum,intSize);
                pfh.read((char *)&firNode,intSize);
                pfh.read(&firOri,1);
                pfh.read((char *)&lagNum,intSize);
                if(firOri == '>' || firOri == '<'){
                    if(! oriVec.empty()){
                        oriMulti.push_back(oriVec);
                        nodeMulti.push_back(ndVec);
                        qPosMulti.push_back(qPosVec);
                        //
                        nameMulti.push_back(firName);

                        oriVec.clear();
                        ndVec.clear();
                        qPosVec.clear();
                    }
                    //
                    firName = allName[pathNum];
                }else{
                    if(firOri == '1'){
                        firOri = '>';
                    }else{
                        firOri = '<';
                    }
                }
                // index of first node
                int curPos = ndVec.size();
                oriVec.push_back(firOri);
                ndVec.push_back(firNode);
                for(int i = 0; i < lagNum; ++i){
                    LagNode tLag;
                    pfh.read((char *)&tLag,lndSize);
                    ndVec.push_back(firNode + tLag.diff);
                    oriVec.push_back(tLag.ori);
                }
                //
                int tpos;
                pfh.read((char *)&tpos,intSize);
                qPosVec.push_back(tpos);            
                bool nf = false;
                // -1 --> Some graph nodes were not stored during graph indexing. Coordinate is uncertain.
                for(int j = 0; j < lagNum; ++j){
                    if(nf){
                        qPosVec.push_back(-1);
                    }else{
                        if(info.find(ndVec[curPos + j]) != info.end()){
                            tpos += info[ndVec[curPos + j]].len;
                            qPosVec.push_back(tpos);                          
                        }else{
                            nf = true;
                            qPosVec.push_back(-1);
                        }
                    }
                }
                //
                tnum += (lagNum + 1);
                if(tnum == total){
                    break;
                }
            }
            if(! oriVec.empty()){
                oriMulti.push_back(oriVec);
                nodeMulti.push_back(ndVec);
                qPosMulti.push_back(qPosVec);
                //
                nameMulti.push_back(firName);

                oriVec.clear();
                ndVec.clear();
                qPosVec.clear();
            }
        }
        //

    }else{
        cerr<<"Warning: path can't be found. "<<endl;
    }
    pfh.close();
    xfh.close();
}

void GraphRange::dxAsmNode(bool refSim,int asmNum,int chrNum,int sStart,int sEnd,unordered_map<NodeType,LenAsm> &info,map<NEdge,int> &r_edge_dict,unordered_map<NodeType,Nid> &nid_dict){
    vector<vector<char> > oriMulti;
    vector<vector<int> > nodeMulti;
    //
    vector<vector<int> > qPosMulti;
    vector<string> cigarMulti;
    vector<string> nameMulti;
    
    bool formR = true;
    string formFile = upDir + "/form.info";
	if(access(formFile.c_str(),F_OK) == 0){
        formR = false;
    }
    
    queryDbPath(formR,asmNum,chrNum,sStart,sEnd,info,oriMulti,nodeMulti,qPosMulti,cigarMulti,nameMulti);
    unordered_set<int> ndGroup;
    //
    bool visCigar = true;
    if(formR){
        for(size_t i =0; i < oriMulti.size(); ++i){
            eAsmFind(visCigar,refSim,oriMulti[i],nodeMulti[i],qPosMulti[i],cigarMulti[i],nameMulti[i],r_edge_dict,nid_dict,ndGroup);
        }
    }else{
        string rgCigar = "*";
        for(size_t i =0; i < oriMulti.size(); ++i){
            eAsmFind(visCigar,refSim,oriMulti[i],nodeMulti[i],qPosMulti[i],rgCigar,nameMulti[i],r_edge_dict,nid_dict,ndGroup);
        }
    }
}
//----------------------------------
/*
    When multiple assemblies are selected to highlight the detailed information on the highlighted paths will not be displayed 
*/
void GraphRange::dxAsmNode2(bool refSim,vector<int> &asmVec,int chrNum,int sStart,int sEnd,unordered_map<NodeType,LenAsm> &info,map<NEdge,int> &r_edge_dict,unordered_map<NodeType,Nid> &nid_dict){
    int asmCode = 0;
    map<NEdge,int> h_edge_dict;
    for(int asmNum : asmVec){
        vector<vector<char> > oriMulti;
        vector<vector<int> > nodeMulti;
        //
        vector<vector<int> > qPosMulti;
        vector<string> cigarMulti;
        vector<string> nameMulti;
        
        bool formR = true;
        string formFile = upDir + "/form.info";
        if(access(formFile.c_str(),F_OK) == 0){
            formR = false;
        }
        queryDbPath(formR,asmNum,chrNum,sStart,sEnd,info,oriMulti,nodeMulti,qPosMulti,cigarMulti,nameMulti);
        unordered_set<int> ndGroup;
        for(size_t i =0; i < oriMulti.size(); ++i){
            eAsmFind2(refSim,asmCode,oriMulti[i],nodeMulti[i],nid_dict,ndGroup,h_edge_dict);
        }
        //
        ++asmCode;
    }
    //
    hEdge2fig(nid_dict,h_edge_dict,r_edge_dict);
}

//--------------------------
void GraphRange::readRefGene(string &ovFile,string &gDxFile,int chrNum,int sStart,int sEnd,unordered_set<NodeType> &retainID,unordered_map<NodeType,char> &exNode,vector<GeneNode> &refGeneNode,unordered_set<string> midExon){
    ifstream gf(ovFile.c_str());
    if(! gf){
        cerr<<"Warning: file open failed. "<<ovFile<<endl;
        return;
        //exit(1);
    }
    ifstream gdx(gDxFile.c_str());
    if(! gdx){
        cerr<<"Warning: file open failed. "<<gDxFile<<endl;
        return;
        //exit(1);
    }
    //
    int nchr = 0;
    int intSize = sizeof(int);
    gdx.read((char *)&nchr,intSize);
    int refChr = -1;
    bool findChr = false;
    //int total = 0;
    ChrRange crRange;
    int crSize = sizeof(ChrRange);
    for(int x = 0; x < nchr; ++x){
        gdx.read((char *)&refChr,intSize);
        gdx.read((char *)&crRange,crSize);
        if(refChr == chrNum){
            findChr = true;
            break;
        }
    }
    
    long long refOff = 0LL;
    int rt = 0;
    bool flag = true;
    int oneSize = sizeof(OneRange);

    if(findChr){
        gdx.clear();
        gdx.seekg(crRange.rByte,ios::beg);
        for(int k = 0; k < crRange.ranNum; ++k){
            OneRange tRange;
            gdx.read((char *)&tRange,oneSize);
            if(tRange.ranStart <= sStart){
                if(tRange.ranEnd >= sStart){
                    if(flag){
                        refOff = tRange.offByte;
                        flag = false;
                    }
                    rt += tRange.ranNum;
                }
                
            }else{
                if(tRange.ranStart <= sEnd){
                    if(flag){
                        refOff = tRange.offByte;
                        flag = false;
                    }
                    rt += tRange.ranNum;
                }else{
                    break;
                }
                
            }
        }
    }
    gdx.close();
    if(flag){
        cerr<<"Warning: annotation (GFF) can't be found in the interval. "<<chrNum<<" : "<<sStart<<" - "<<sEnd<<endl;
        return;
        //exit(1);
    }
    //
    gf.seekg(refOff,ios::beg);
    int unitSize = sizeof(GeneNode);
    unordered_map<NodeType,char>::iterator itl,itr;
    string rnaName = "";
    for(int m = 0; m < rt; ++m){
        GeneNode ndg;
        gf.read((char *)&ndg,unitSize);
        if(retainID.find(ndg.node1) != retainID.end() || retainID.find(ndg.node2) != retainID.end()){
            refGeneNode.push_back(ndg);
            //
            if(ndg.type == 'E'){
                rnaName = ndg.name;
            }
        }else{
            if(ndg.type == 'E'){
                string tRName = ndg.name;
                if(tRName == rnaName){
                    midExon.insert(rnaName);
                }
            }
            //
            if(ndg.type == 'G'){
                itl = exNode.find(ndg.node1);
                if(itl != exNode.end()){
                    if(itl->second == '2'){
                        itr = exNode.find(ndg.node2);
                        if(itr != exNode.end()){
                            if(itr->second == '1'){
                                refGeneNode.push_back(ndg);
                            }
                        }
                    }
                }
            }
        }
    }
    gf.close();
}

/*
    The bound of a searching region may be not the boundary of a gene
*/
void GraphRange::getFigGene(string &bwGeneFile,string &gDxFile,int chrNum,int sStart,int sEnd,float wPerK,unordered_map<NodeType,char> &exNode){

    unordered_set<NodeType> retainID;
    map<NodeType,float> rNodeStart;
    //

    for(size_t i = 0; i < draw_pos.size(); i+=2){
        NodeType id = nnames[draw_node[i]["group"]];
        retainID.insert(id);
        rNodeStart.emplace(id,draw_pos[i]);
    }
    //
    // read range node gene -> refGeneNode
    // sort refGeneNode
    //
    vector<GeneNode> refGeneNode;
    unordered_set<string> midExon;
    readRefGene(bwGeneFile,gDxFile,chrNum,sStart,sEnd,retainID,exNode,refGeneNode,midExon);    
    if(refGeneNode.empty()){
        return;
    }    
    
    // merge gene
    //
    map<string,FigGene> geneMap;
    map<string,FigGene>::iterator it;
    
    map<string,vector<FigEle> > exonMap;
    map<string,vector<FigEle> >::iterator eit;
    
    map<string,vector<FigEle> > cdsMap;
    //
    uint8_t eFir = 1;
    for(size_t i = 0; i < refGeneNode.size(); ++i){
        string geneName = refGeneNode[i].name;
        if(refGeneNode[i].type == 'G'){
            FigGene tf;
            tf.margin = '0';
            tf.layer = refGeneNode[i].layer;
            tf.strand = refGeneNode[i].strand;
            if(rNodeStart.find(refGeneNode[i].node1) != rNodeStart.end()){
                float ndFigStart = rNodeStart[refGeneNode[i].node1];
                tf.start = ndFigStart + refGeneNode[i].reStart1 * wPerK;
            }else{
                float ndFigStart = draw_pos.front();
                tf.start = ndFigStart;
                if(refGeneNode[i].strand == '-'){
                    tf.margin = '2';
                }
            }
            
            if(rNodeStart.find(refGeneNode[i].node2) != rNodeStart.end()){
                float ndFigStart = rNodeStart[refGeneNode[i].node2];
                tf.end = ndFigStart + refGeneNode[i].reStart2 * wPerK;
            }else{
                tf.end = draw_pos.back();
                if(refGeneNode[i].strand == '+'){
                    tf.margin = '1';
                }
            }  
            geneMap.emplace(geneName,tf);
        }else if(refGeneNode[i].type == 'E'){
            FigEle te;
            te.num = refGeneNode[i].num;
            te.layer = refGeneNode[i].layer;
            te.strand = refGeneNode[i].strand;
            //
            if(rNodeStart.find(refGeneNode[i].node1) != rNodeStart.end()){
                float ndFigStart = rNodeStart[refGeneNode[i].node1];
                te.start = ndFigStart + refGeneNode[i].reStart1 * wPerK;
            }else{
                float ndFigStart = draw_pos.front();
                te.start = ndFigStart;
            }
            
            if(rNodeStart.find(refGeneNode[i].node2) != rNodeStart.end()){
                float ndFigStart = rNodeStart[refGeneNode[i].node2];
                te.end = ndFigStart + refGeneNode[i].reStart2 * wPerK;
            }else{
                te.end = draw_pos.back();
            }
            //
            eit = exonMap.find(geneName);
            if(eit != exonMap.end()){
                (eit->second).push_back(te);
            }else{
                vector<FigEle> tv;
                tv.push_back(te);
                //
                exonMap.emplace(geneName,tv);
            }
            
        }else{
            FigEle te;
            //te.num = refGeneNode[i].num;
            te.layer = refGeneNode[i].layer;
            te.strand = refGeneNode[i].strand;
            //
            if(rNodeStart.find(refGeneNode[i].node1) != rNodeStart.end()){
                float ndFigStart = rNodeStart[refGeneNode[i].node1];
                te.start = ndFigStart + refGeneNode[i].reStart1 * wPerK;
            }else{
                float ndFigStart = draw_pos.front();
                te.start = ndFigStart;
            }
            
            if(rNodeStart.find(refGeneNode[i].node2) != rNodeStart.end()){
                float ndFigStart = rNodeStart[refGeneNode[i].node2];
                te.end = ndFigStart + refGeneNode[i].reStart2 * wPerK;
            }else{
                te.end = draw_pos.back();
            }
            //
            eit = cdsMap.find(geneName);
            if(eit != cdsMap.end()){
                (eit->second).push_back(te);
            }else{
                vector<FigEle> tv;
                tv.push_back(te);
                //
                cdsMap.emplace(geneName,tv);
            }
        }
    }
    for(auto &tg : geneMap){
        geneVec.push_back(tg.first);
        
        ndGenePos.push_back(tg.second.start);
        ndGenePos.push_back(tg.second.end);
        layerVec.push_back(tg.second.layer);
        strandVec.push_back(tg.second.strand);
        mgFlagVec.push_back(tg.second.margin);
    }
    
    for(auto &te : exonMap){
        rnaVec.push_back(te.first);
        for(auto &nd : te.second){
            ndExonPos.push_back(nd.start);
            ndExonPos.push_back(nd.end);
            eLayerVec.push_back(nd.layer);
            eStrandVec.push_back(nd.strand);
        }
        eNumVec.push_back(te.second.size());
        //
        if(te.second.front().num != eFir){
            if(midExon.find(te.first) != midExon.end()){
                eFlagVec.push_back('3');
            }else{
                eFlagVec.push_back('2');
            }
        }else{
            if(midExon.find(te.first) != midExon.end()){
                eFlagVec.push_back('1');
            }else{
                eFlagVec.push_back('0');
            }
        }
    }
    
    for(auto &tc : cdsMap){
        cdsVec.push_back(tc.first);
        for(auto &nd : tc.second){
            ndCDSPos.push_back(nd.start);
            ndCDSPos.push_back(nd.end);
            cLayerVec.push_back(nd.layer);
        }
        cNumVec.push_back(tc.second.size());
    }
}
//------------------

void GraphRange::readRefBed(string &ovFile,string &gDxFile,int chrNum,int sStart,int sEnd,unordered_set<NodeType> &retainID,unordered_map<NodeType,char> &exNode,vector<BedNode> &refBedNode){
    ifstream gf(ovFile.c_str());
    if(! gf){
        cerr<<"Warning: file open failed. "<<ovFile<<endl;
        return;
        //exit(1);
    }
    ifstream gdx(gDxFile.c_str());
    if(! gdx){
        cerr<<"Warning: file open failed. "<<gDxFile<<endl;
        return;
        //exit(1);
    }
    //
    int nchr = 0;
    int intSize = sizeof(int);
    gdx.read((char *)&nchr,intSize);
    int refChr = -1;
    bool findChr = false;
    //
    ChrRange crRange;
    int crSize = sizeof(ChrRange);
    for(int x = 0; x < nchr; ++x){
        gdx.read((char *)&refChr,intSize);
        gdx.read((char *)&crRange,crSize);
        if(refChr == chrNum){
            findChr = true;
            break;
        }
    }
    
    long long refOff = 0LL;
    int rt = 0;
    bool flag = true;
    int oneSize = sizeof(OneRange);

    if(findChr){
        gdx.clear();
        gdx.seekg(crRange.rByte,ios::beg);
        for(int k = 0; k < crRange.ranNum; ++k){
            OneRange tRange;
            gdx.read((char *)&tRange,oneSize);
            if(tRange.ranStart <= sStart){
                if(tRange.ranEnd >= sStart){
                    if(flag){
                        refOff = tRange.offByte;
                        flag = false;
                    }
                    rt += tRange.ranNum;                   
                }
                
            }else{
                if(tRange.ranStart <= sEnd){
                    if(flag){
                        refOff = tRange.offByte;
                        flag = false;
                    }
                    rt += tRange.ranNum;        
                }else{
                    break;
                }
                
            }
        }
    }
    gdx.close();
    if(flag){
        cerr<<"Warning: annotation (BED) can't be found in the interval. "<<chrNum<<" : "<<sStart<<" - "<<sEnd<<endl;
        return;
        //exit(1);
    }
    //
    gf.seekg(refOff,ios::beg);
    int unitSize = sizeof(BedNode);
    unordered_map<NodeType,char>::iterator itl,itr;
    for(int m = 0; m < rt; ++m){
        BedNode ndg;
        gf.read((char *)&ndg,unitSize);      
        if(retainID.find(ndg.node1) != retainID.end() || retainID.find(ndg.node2) != retainID.end()){
            refBedNode.push_back(ndg);
        }else{
            itl = exNode.find(ndg.node1);
            if(itl != exNode.end()){
                if(itl->second == '2'){
                    itr = exNode.find(ndg.node2);
                    if(itr != exNode.end()){
                        if(itr->second == '1'){
                            refBedNode.push_back(ndg);
                        }
                    }
                }
            }
        }
    }
    gf.close();
}

// #nColumns layer   track_name  description
/*
   ncol: 3 --> (chr,start,end); 4 --> (chr,start,end,name); 5 --> (chr,start,end,name,score); 6 --> (chr,start,end,name,score,strand)
*/
void GraphRange::getFigBed(string &tkDesFile,string &bwGeneFile,string &gDxFile,int chrNum,int sStart,int sEnd,float wPerK,unordered_map<NodeType,char> &exNode){
    ifstream kfh(tkDesFile.c_str());
    if(! kfh){
        cerr<<"Warning: file open failed. "<<tkDesFile<<endl;
        return;
    }
    string line;
    while(getline(kfh,line)){
        int fir = line.find("\t");
        int sec = line.find("\t",fir+1);
        int thi = line.find("\t",sec+1);
        int ncol = atoi(line.substr(0,fir).c_str());
        int cum = atoi(line.substr(fir+1,sec-fir-1).c_str());
        string name = line.substr(sec+1,thi-sec-1);
        string des = line.substr(thi+1);
        tkNameVec.push_back(name);
        tkDesVec.push_back(des);
        tkColVec.push_back(ncol);
        tkCumVec.push_back(cum);
    }
    kfh.close();
    //
    unordered_set<NodeType> retainID;
    map<NodeType,float> rNodeStart;
    //
    for(size_t i = 0; i < draw_pos.size(); i+=2){
        NodeType id = nnames[draw_node[i]["group"]];
        retainID.insert(id);
        rNodeStart.emplace(id,draw_pos[i]);
    }
    //
    vector<BedNode> refBedNode;   
    readRefBed(bwGeneFile,gDxFile,chrNum,sStart,sEnd,retainID,exNode,refBedNode);    
    if(refBedNode.empty()){
        return;
    }
    //
    map<int,vector<FigBed> > geneMap;
    for(size_t j = 0; j < refBedNode.size(); ++j){
        FigBed fb;
        int layer = refBedNode[j].layer;
        fb.name = refBedNode[j].name;
        fb.strand = refBedNode[j].strand;
        fb.score = refBedNode[j].score;
        if(rNodeStart.find(refBedNode[j].node1) != rNodeStart.end()){
            float ndFigStart = rNodeStart[refBedNode[j].node1];
            fb.start = ndFigStart + refBedNode[j].reStart1 * wPerK;
        }else{
            float ndFigStart = draw_pos.front();
            fb.start = ndFigStart;
        }
        
        if(rNodeStart.find(refBedNode[j].node2) != rNodeStart.end()){
            float ndFigStart = rNodeStart[refBedNode[j].node2];
            fb.end = ndFigStart + refBedNode[j].reStart2 * wPerK;
        }else{
            fb.end = draw_pos.back();
        }
        //
        if(geneMap.find(layer) != geneMap.end()){
            geneMap[layer].push_back(fb);
        }else{
            vector<FigBed> vec{fb};
            geneMap.emplace(layer,vec);
        }
    }
    //
    int layStart = 1;
    // tnum --> number of items for each track
    for(size_t k = 0; k < tkColVec.size(); ++k){
        int tnum = 0;
        if(tkColVec[k] == 3){
            // start end
            for(int x = layStart; x <= tkCumVec[k]; ++x){
                for(auto &gm : geneMap[x]){
                    rBedPos.push_back(gm.start);
                    rBedPos.push_back(gm.end);
                    rBedLayer.push_back(x);
                }
                tnum += geneMap[x].size();
            }
        }else if(tkColVec[k] == 4){
            // start end name
            for(int x = layStart; x <= tkCumVec[k]; ++x){
                for(auto &gm : geneMap[x]){
                    rBedPos.push_back(gm.start);
                    rBedPos.push_back(gm.end);
                    rBedLayer.push_back(x);
                    //
                    rBedName.push_back(gm.name);
                }
                tnum += geneMap[x].size();
            }
        }else if(tkColVec[k] == 5){
            // start end (name; not store for the present) score
            for(int x = layStart; x <= tkCumVec[k]; ++x){
                for(auto &gm : geneMap[x]){
                    rBedPos.push_back(gm.start);
                    rBedPos.push_back(gm.end);
                    rBedLayer.push_back(x);
                    //
                    rBedName.push_back(gm.name);
                    rBedScore.push_back(gm.score);
                }
                tnum += geneMap[x].size();
            }
        }else{
            // start end (name; not store for the present) score  strand
            for(int x = layStart; x <= tkCumVec[k]; ++x){
                for(auto &gm : geneMap[x]){
                    rBedPos.push_back(gm.start);
                    rBedPos.push_back(gm.end);
                    rBedLayer.push_back(x);
                    //
                    rBedName.push_back(gm.name);
                    rBedScore.push_back(gm.score);
                    rBedStrand.push_back(gm.strand);
                }
                tnum += geneMap[x].size();
            }
        }
        //
        tkItem.push_back(tnum);
        layStart = tkCumVec[k] + 1;
    }
}

//-----------------
void GraphRange::formatGraph(string &asmb,string &sChr,int sStart,int sEnd,int ex,int wStart,int wWidth,int wCut,int wY,int queryDep,int varLen,bool sim,bool refSim){
    unordered_map<NodeType,vector<ENode> > iedge;
    unordered_map<NodeType,vector<ENode> > oedge;
    vector<NodeType> rangeNode;
    unordered_map<NodeType,char> exNode;
    unordered_map<NodeType,LenAsm> info;
    int realLen = 0;
    int realStart = 0;
    //
    int chrNum = -1;
    if(indexFlag == 1){
        chrNum = getChrNum(sChr);
        if(chrNum < 0){
            cerr<<"Error: chromosome can't be found. "<<sChr<<endl;
            return;
            //exit(1);
        }
        getExNode(chrNum,sStart,sEnd,ex,rangeNode,exNode,info,realLen,realStart);
        parseIndex(chrNum,sStart,sEnd,iedge,oedge);
    }else{
        parseNode(sChr,sStart,sEnd,ex,rangeNode,exNode,info,realLen);
        parseEdge(iedge,oedge);
    }
    //
    int rcount = rangeNode.size();
    int baseWidth = wWidth - 80 * (rcount-1);
    if(baseWidth > 0){
        int spLimit = 80 * realLen / baseWidth + 1;
        if(wCut > spLimit){
            wCut = spLimit;
        }
    }
    
    float wPerK = (float)wWidth / (realLen + (rcount-1) * wCut);
    //
    float x_wCut = wPerK * wCut;   
    figScale = x_wCut;
    //
    int gNum = 0;
    int iNum = 0;    
    //
    Ndic firNode;
    firNode.insert(Ndic::value_type("id",0));
    firNode.insert(Ndic::value_type("group",0));
    draw_pos.push_back(wStart);
    draw_node.push_back(firNode);

    float node_pre = wStart;
    unordered_set<NodeType> nRefNode;
    unordered_map<NodeType,Nid> nid_dict;
    int s_nid = 0;
    int e_nid = 0;
    NodeType pre_node = 0;
    unordered_map<NodeType,int> range_set;
    int nx = 0;
    for(NodeType &rnode: rangeNode){
        range_set.emplace(rnode,nx);
        ++nx;
    }
    
    map<NEdge,int> r_edge_dict;
    unordered_map<NodeType,std::vector<ENode> >::iterator it,itx;
    unordered_set<NodeType> subMap;
    int mnx = nx - 1;
    if(sim){
        int point = 0;
        int lpoint = 0,rpoint = 0;
        while(point < nx){
            lpoint = point;
            rpoint = point;

            int tnode = rangeNode[point];
            it = oedge.find(tnode);
            if(it != oedge.end()){
                vector<NodeType> tNref;
                vector<int> deep;
                vector<char> sign;
                //
                char lnrMark = '2',rnrMark = '2';
                char lstMark = '2',rstMark = '2';
                char lOri = 'o',rOri = 'o';

                
                int maxLen = 0;
                //bool flag = true;
                for(ENode &o_node : it->second){
                    if(exNode.find(o_node.node) == exNode.end()){
                        if(subMap.find(o_node.node) == subMap.end()){
                            tNref.push_back(o_node.node);
                            deep.push_back(0);
                            int tlen = info[o_node.node].len;
                            if(tlen > maxLen){
                                maxLen = tlen;
                            }
                            //
                            sign.push_back(o_node.mark);
                        }
                        //else{
                        //    flag = false;
                        //    break;
                        //}
                    }
                }
                //
                //if(flag){
                    if(maxLen < varLen){
                        if(! tNref.empty()){
                            size_t i = 0;
                            int aMaxLen = 0;
                            int preDeep = 1;
                            
                            unordered_set<NodeType> tset;
                            while(i < tNref.size()){
                                if(deep[i] != preDeep){
                                    maxLen += aMaxLen;
                                    if(maxLen >= varLen){
                                        break;
                                    }
                                    //
                                    aMaxLen = 0;
                                }
                                if(deep[i] > queryDep){
                                    break;   
                                }
                                
                                it = oedge.find(tNref[i]);
                                itx = iedge.find(tNref[i]); 
                                if(it != oedge.end()){
                                    for(ENode &to_node : it->second){
                                        if(exNode.find(to_node.node) == exNode.end()){
                                            //if(subMap.find(to_node.node) == subMap.end()){
                                            if(tset.find(to_node.node) == tset.end()){
                                                tNref.push_back(to_node.node);
                                                deep.push_back(deep[i]+1);
                                                tset.insert(to_node.node);
                                                //
                                                if(info[to_node.node].len > aMaxLen){
                                                    aMaxLen = info[to_node.node].len;
                                                }
                                                //
                                                sign.push_back(sign[i]);
                                            }
                                            //}
                                        }else{
                                            if(range_set.find(to_node.node) != range_set.end()){
                                                int tpoint = range_set[to_node.node];
                                                if(tpoint > rpoint){
                                                    rpoint = tpoint;
                                                    //
                                                    rnrMark = to_node.mark;
                                                    rstMark = sign[i];
                                                    rOri = 'o';
                                                }else{
                                                    if(tpoint < lpoint){
                                                        lpoint = tpoint;
                                                        //
                                                        lnrMark = to_node.mark;
                                                        lstMark = sign[i];
                                                        lOri = 'o';
                                                    }
                                                }
                                                
                                            }
                                        }
                                    }
                                }
                                //
                                if(itx != iedge.end()){
                                    for(ENode &to_node : itx->second){
                                        if(exNode.find(to_node.node) == exNode.end()){
                                            //if(subMap.find(to_node.node) == subMap.end()){
                                            if(tset.find(to_node.node) == tset.end()){
                                                tNref.push_back(to_node.node);
                                                deep.push_back(deep[i]+1);
                                                tset.insert(to_node.node);
                                                //
                                                if(info[to_node.node].len > aMaxLen){
                                                    aMaxLen = info[to_node.node].len;
                                                }
                                                //
                                                sign.push_back(sign[i]);
                                            }
                                            //}
                                        }else{
                                            if(range_set.find(to_node.node) != range_set.end()){
                                                int tpoint = range_set[to_node.node];
                                                if(tpoint > rpoint){
                                                    rpoint = tpoint;
                                                    rnrMark = to_node.mark;
                                                    rstMark = sign[i];
                                                    rOri = 'i';
                                                }else{
                                                    if(tpoint < lpoint){
                                                        lpoint = tpoint;
                                                        //
                                                        lnrMark = to_node.mark;
                                                        lstMark = sign[i];
                                                        lOri = 'i';
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                //
                                preDeep = deep[i];
                                ++i;
                            }
                            //
                            maxLen += aMaxLen;
                        }
                    }
                
                    //
                    //if(rpoint > lpoint); remove end point
                    if(maxLen < varLen){
                        int rlen = 0,llen = 0;
                        bool inv = false;
                        int lStart = lpoint + 1,rEnd = rpoint;
                        if(point > lpoint){
                            
                            if(lstMark < '4'){
                                // + out
                                inv = true;
                            }else{
                                if(lOri == 'o'){
                                    if(lnrMark == '2' || lnrMark == '4'){
                                        lStart = lpoint;
                                    }
                                }else{
                                    if(lnrMark == '4' || lnrMark == '5'){
                                        lStart = lpoint;
                                    }
                                }
                            }

                            if(! inv){
                                for(int j = lStart; j < point; ++j){
                                    llen += info[rangeNode[j]].len;
                                }
                                if(llen < varLen){
                                    if(rpoint > point){
                                        if(rstMark < '4'){
                                            if(rOri == 'o'){
                                                // start -> end(-)
                                                if(rnrMark == '3' || rnrMark == '5'){
                                                    rEnd = rpoint + 1;
                                                }
                                            }else{
                                                // start(+) -> end
                                                if(rnrMark == '2' || rnrMark == '3'){
                                                    rEnd = rpoint + 1;
                                                }
                                            }
                                        }else{
                                            // - out
                                            inv = true;
                                        }
                                        //
                                        if(! inv){
                                            for(int j = point + 1; j < rEnd; ++j){
                                                rlen += info[rangeNode[j]].len;
                                            }
                                            if(rlen < varLen){
                                                for(auto &xnode : tNref){
                                                    subMap.insert(xnode);
                                                }
                                            }
                                        }
                                    }else{
                                        for(auto &xnode : tNref){
                                            subMap.insert(xnode);
                                        }
                                    }
                                }
                            }
                        }else{
                            if(rpoint > point){
                                if(rstMark < '4'){
                                    if(rOri == 'o'){
                                        if(rnrMark == '3' || rnrMark == '5'){
                                            rEnd = rpoint + 1;
                                        }
                                    }else{
                                        if(rnrMark == '2' || rnrMark == '3'){
                                            rEnd = rpoint + 1;
                                        }
                                    }
                                }else{
                                    inv = true;
                                }
                                //
                                if(! inv){
                                    for(int j = point + 1; j < rEnd; ++j){
                                        rlen += info[rangeNode[j]].len;
                                    }
                                    if(rlen < varLen){
                                        for(auto &xnode : tNref){
                                            subMap.insert(xnode);
                                        }
                                    }
                                }
                            }
                        }
                    }
                //}
            }
            //
            //if(rpoint == point){
                ++point;
            //}else{
            //    point = rpoint;   
            //}
        }
    }
    //---------------------------------
    float pTickPos = wStart - 100;
    float curStart = wStart;
    float curEnd = wStart;
    int preCum = 0;
    int refLen = 0;
    //----------------------------------
    NodeType simNode = rangeNode.front();
    bool preOte = false;
    for(NodeType &tnode: rangeNode){    
        int node_len = info[tnode].len;
        nnames.push_back(tnode);
        genome.push_back(info[tnode].asmb);
        //
        refLen += node_len;
        //
        it = oedge.find(tnode);
        bool emp1 = true,emp2 = true;
        bool alink1 = true,alink2 = true;
        //
        vector<NodeType> tNref;
        vector<int> deep;
        if(it != oedge.end()){   
            for(ENode &o_node : it->second){
                if(exNode.find(o_node.node) != exNode.end()){
                    if(range_set.find(o_node.node) != range_set.end()){
                        bool flag = true;
                        if(sim){
                            if(range_set[o_node.node] > range_set[tnode]){
                                int vlen = 0;
                                for(int k = range_set[tnode] + 1; k < range_set[o_node.node]; ++k){
                                    vlen += info[rangeNode[k]].len;
                                }
                                if(vlen < varLen){
                                    if(o_node.mark == '2'){
                                        flag = false;
                                    }
                                }
                            }else{
                                int vlen = 0;
                                for(int k = range_set[o_node.node] + 1; k < range_set[tnode]; ++k){
                                    vlen += info[rangeNode[k]].len;
                                }
                                if(vlen < varLen){
                                    if(o_node.mark == '5'){
                                        flag = false;
                                    }
                                }
                            }
                        }
                        if(flag){
                            NEdge sym = {tnode,o_node.node,o_node.mark};
                            if(r_edge_dict.find(sym) == r_edge_dict.end()){
                                r_edge_dict.emplace(sym,2);
                            }
                            emp1 = false;
                        }
                    }
                }else{
                    if(nRefNode.find(o_node.node) == nRefNode.end()){
                        if(sim){
                            if(subMap.find(o_node.node) != subMap.end()){
                                continue;
                            }
                        }
                        
                        tNref.push_back(o_node.node);
                        deep.push_back(0);
                        nRefNode.insert(o_node.node);                        
                    }
                    //
                    NEdge sym = {tnode,o_node.node,o_node.mark};
                    r_edge_dict.emplace(sym,2);
                    emp1 = false;
                    //
                    alink1 = false;
                }
            }
        }
        it = iedge.find(tnode);
        if(it != iedge.end()){
            for(ENode &o_node : it->second){
                if(exNode.find(o_node.node) != exNode.end()){
                    if(sim){
                        NEdge sym = {o_node.node,tnode,o_node.mark};
                        if(r_edge_dict.find(sym) != r_edge_dict.end()){
                            emp2 = false;
                        }else{
                            if(range_set.find(o_node.node) != range_set.end()){
                                bool flag = true;
                                
                                if(range_set[o_node.node] > range_set[tnode]){
                                    int vlen = 0;
                                    for(int k = range_set[tnode] + 1; k < range_set[o_node.node]; ++k){
                                        vlen += info[rangeNode[k]].len;
                                    }
                                    if(vlen < varLen){
                                        if(o_node.mark == '5'){
                                            flag = false;
                                        }
                                    }
                                }else{
                                    int vlen = 0;
                                    for(int k = range_set[o_node.node] + 1; k < range_set[tnode]; ++k){
                                        vlen += info[rangeNode[k]].len;
                                    }
                                    if(vlen < varLen){
                                        if(o_node.mark == '2'){
                                            flag = false;
                                        }
                                    }
                                }
                                
                                if(flag){
                                    //NEdge sym = {o_node.node,tnode,o_node.mark};
                                    //if(r_edge_dict.find(sym) == r_edge_dict.end()){
                                        r_edge_dict.emplace(sym,2);
                                    //}
                                    emp2 = false;
                                }
                            }
                        }
                    }
                }else{
                    if(nRefNode.find(o_node.node) == nRefNode.end()){
                        if(sim){
                            if(subMap.find(o_node.node) != subMap.end()){
                                continue;
                            }
                        }
                        //
                        tNref.push_back(o_node.node);
                        deep.push_back(0);
                        nRefNode.insert(o_node.node);
                        //vRefNode.push_back(o_node.node);
                    }
                    //
                    emp2 = false;
                    //
                    NEdge sym = {o_node.node,tnode,o_node.mark};
                    r_edge_dict.emplace(sym,2);
                    //
                    alink2 = false;
                }
            }
        }
            
        if(! tNref.empty()){
            size_t i = 0;
            while(i < tNref.size()){             
                if(deep[i] > queryDep){
                    break;   
                }
                
                it = oedge.find(tNref[i]);
                if(it != oedge.end()){
                    for(ENode &to_node : it->second){
                        if(exNode.find(to_node.node) != exNode.end()){
                            if(range_set.find(to_node.node) != range_set.end()){
                                NEdge sym = {tNref[i],to_node.node,to_node.mark};
                                r_edge_dict.emplace(sym,2);
                            }
                        }else{
                            if(nRefNode.find(to_node.node) == nRefNode.end()){
                                tNref.push_back(to_node.node);
                                deep.push_back(deep[i]+1);
                                nRefNode.insert(to_node.node);
                                //vRefNode.push_back(to_node.node);
                            }
                            NEdge sym = {tNref[i],to_node.node,to_node.mark};
                            r_edge_dict.emplace(sym,2);
                        }
                    }
                }

                it = iedge.find(tNref[i]);
                if(it != iedge.end()){
                    for(ENode &to_node : it->second){
                        if(exNode.find(to_node.node) != exNode.end()){
                            if(range_set.find(to_node.node) != range_set.end()){
                                NEdge sym = {to_node.node,tNref[i],to_node.mark};
                                r_edge_dict.emplace(sym,2);
                            }
                        }else{
                            if(nRefNode.find(to_node.node) == nRefNode.end()){
                                tNref.push_back(to_node.node);
                                deep.push_back(deep[i]+1);
                                nRefNode.insert(to_node.node);
                                //vRefNode.push_back(to_node.node);
                            }
                            NEdge sym = {to_node.node,tNref[i],to_node.mark};
                            r_edge_dict.emplace(sym,2);
                        }
                    }
                }
                i += 1;
            }
        }
        //
        if(refSim){
            bool dRef = false;
            if(emp1 && emp2){
                dRef = true;
            }
            //
            bool curRes = false;
            bool linkOend = false,linkOte = false;
            if(! dRef || (gNum == mnx)){
                linkOend = true;
            }
            if(linkOend || (gNum == 0)){
                linkOte = true;
                curRes = true;
            }else{
                if(preOte && alink1 && alink2){
                    curRes = true;
                }
            }
            // alink1 && alink2 => has no non-ref link
            //
            if(gNum > 0){
                node_pre += x_wCut;
                //
                curStart = node_pre;
                //
                //if(! dRef || (gNum == mnx)){
                if(curRes){
                    iNum += 1;
                    //
                    Ndic tdNode;
                    tdNode.insert(Ndic::value_type("id",iNum));
                    tdNode.insert(Ndic::value_type("group",gNum));
                    //tdNode.insert(Ndic::value_type("fx",node_pre));
                    //tdNode.insert(Ndic::value_type("fy",wY));
                    draw_pos.push_back(node_pre);
                    draw_node.push_back(tdNode);
                    
                    map<string,int> tdEdge;
                    tdEdge.insert(map<string,int>::value_type("source",iNum-1));
                    tdEdge.insert(map<string,int>::value_type("target",iNum));
                    tdEdge.insert(map<string,int>::value_type("type",1));
                    draw_edge.push_back(tdEdge);
                    s_nid = iNum;
                    NEdge tsym = {simNode,tnode,'2'};
                    if(r_edge_dict.find(tsym) != r_edge_dict.end()){
                        if(simNode == pre_node){
                            r_edge_dict[tsym] = 0;
                        }
                    }else{
                        r_edge_dict.emplace(tsym,0);
                    }
                    //
                    simNode = tnode;
                }
            }
            
            float trans = node_len * wPerK;
            node_pre += trans;
            //
            curEnd = node_pre;
            //
            //if(! dRef || (gNum == 0) || (gNum == mnx)){
            if(curRes){    
                iNum += 1;
                //
                Ndic tdNode;
                tdNode.insert(Ndic::value_type("id",iNum));
                tdNode.insert(Ndic::value_type("group",gNum));
                draw_pos.push_back(node_pre);
                draw_node.push_back(tdNode);
                
                map<string,int> tdEdge;
                tdEdge.insert(map<string,int>::value_type("source",iNum-1));
                tdEdge.insert(map<string,int>::value_type("target",iNum));
                tdEdge.insert(map<string,int>::value_type("type",0));
                draw_edge.push_back(tdEdge);
                //
                e_nid = iNum;
                Nid tnid = {s_nid,e_nid,gNum};            
                nid_dict.emplace(tnode,tnid);
                //
                
            }
            //else{
            //    refSkip.insert(tnode);
            //}
            gNum += 1;
            pre_node = tnode;
            preOte = linkOte;
            //----------------------------
            if(curStart - pTickPos >= 80.0f){
                int tValue = preCum + realStart;
                tickValue.push_back(preCum + realStart);
                tickPos.push_back(curStart);
                pTickPos = curStart;
                float tPos = curStart;
                while(tPos <= curEnd){            
                    tPos += 80.0f;
                    if(tPos <= curEnd){
                        tValue += int(80.0f / (curEnd - curStart) * node_len);
                        tickValue.push_back(tValue);
                        tickPos.push_back(tPos);
                        pTickPos = tPos;
                    }
                }
            }else{
                if(curEnd - pTickPos >= 80){
                    float tPos = pTickPos + 80;
                    int tValue = int((tPos - curStart) / (curEnd - curStart) * node_len) + preCum + realStart - 1;
                    tickValue.push_back(tValue);
                    tickPos.push_back(tPos);
                    pTickPos = tPos;
                    while(tPos <= curEnd){            
                        tPos += 80;
                        if(tPos <= curEnd){
                            tValue += int(80 / (curEnd - curStart) * node_len);
                            tickValue.push_back(tValue);
                            tickPos.push_back(tPos);
                            pTickPos = tPos;
                        }
                    }
                }    
            }
            //---------------------------------
        }else{
            if(gNum > 0){
                node_pre += x_wCut;
                //
                curStart = node_pre;
                //
                iNum += 1;
                
                Ndic tdNode;
                tdNode.insert(Ndic::value_type("id",iNum));
                tdNode.insert(Ndic::value_type("group",gNum));
                //tdNode.insert(Ndic::value_type("fx",node_pre));
                //tdNode.insert(Ndic::value_type("fy",wY));
                draw_pos.push_back(node_pre);
                draw_node.push_back(tdNode);
                
                map<string,int> tdEdge;
                tdEdge.insert(map<string,int>::value_type("source",iNum-1));
                tdEdge.insert(map<string,int>::value_type("target",iNum));
                tdEdge.insert(map<string,int>::value_type("type",1));
                draw_edge.push_back(tdEdge);
                s_nid = iNum;
                NEdge tsym = {pre_node,tnode,'2'};
                if(r_edge_dict.find(tsym) != r_edge_dict.end()){
                    r_edge_dict[tsym] = 0;
                }else{
                    r_edge_dict.emplace(tsym,0);
                }
            }
            
            float trans = node_len * wPerK;
            node_pre += trans;
            //
            curEnd = node_pre;
            //
            iNum += 1;
            
            Ndic tdNode;
            tdNode.insert(Ndic::value_type("id",iNum));
            tdNode.insert(Ndic::value_type("group",gNum));
            draw_pos.push_back(node_pre);
            draw_node.push_back(tdNode);
            
            map<string,int> tdEdge;
            tdEdge.insert(map<string,int>::value_type("source",iNum-1));
            tdEdge.insert(map<string,int>::value_type("target",iNum));
            tdEdge.insert(map<string,int>::value_type("type",0));
            draw_edge.push_back(tdEdge);
            //
            e_nid = iNum;
            Nid tnid = {s_nid,e_nid,gNum};
            nid_dict.emplace(tnode,tnid);
            
            gNum += 1;
            pre_node = tnode;
            //----------------------------
            if(curStart - pTickPos >= 80.0f){
                int tValue = preCum + realStart;
                tickValue.push_back(preCum + realStart);
                tickPos.push_back(curStart);
                pTickPos = curStart;
                float tPos = curStart;
                while(tPos <= curEnd){            
                    tPos += 80.0f;
                    if(tPos <= curEnd){
                        tValue += int(80.0f / (curEnd - curStart) * node_len);
                        tickValue.push_back(tValue);
                        tickPos.push_back(tPos);
                        pTickPos = tPos;
                    }
                }
            }else{
                if(curEnd - pTickPos >= 80){
                    float tPos = pTickPos + 80;
                    int tValue = int((tPos - curStart) / (curEnd - curStart) * node_len) + preCum + realStart - 1;
                    tickValue.push_back(tValue);
                    tickPos.push_back(tPos);
                    pTickPos = tPos;
                    while(tPos <= curEnd){            
                        tPos += 80;
                        if(tPos <= curEnd){
                            tValue += int(80 / (curEnd - curStart) * node_len);
                            tickValue.push_back(tValue);
                            tickPos.push_back(tPos);
                            pTickPos = tPos;
                        }
                    }
                }    
            }
            //---------------------------------
        }
        //
        preCum = refLen;
    }   
    //-----------------------------
    for(auto &nk : nRefNode){       
        nnames.push_back(nk);
        int node_len = info[nk].len;
        genome.push_back(info[nk].asmb);
        //
        if(node_len <= wCut){
            iNum += 1;
            s_nid = iNum;
            Ndic tdNode;
            tdNode.insert(Ndic::value_type("id",iNum));
            tdNode.insert(Ndic::value_type("group",gNum));
            draw_node.push_back(tdNode);
            
            iNum += 1;
            e_nid = iNum;
            
            Ndic tdNode2;
            tdNode2.insert(Ndic::value_type("id",iNum));
            tdNode2.insert(Ndic::value_type("group",gNum));
            draw_node.push_back(tdNode2);
            
            map<string,int> tdEdge;
            tdEdge.insert(map<string,int>::value_type("source",iNum-1));
            tdEdge.insert(map<string,int>::value_type("target",iNum));
            tdEdge.insert(map<string,int>::value_type("type",0));
            draw_edge.push_back(tdEdge);
            
            float trans = node_len * wPerK;
            dnode_len.push_back(trans);
        }else{
            int preLen = wCut;
            iNum += 1;
            s_nid = iNum;
            
            Ndic tdNode;
            tdNode.insert(Ndic::value_type("id",iNum));
            tdNode.insert(Ndic::value_type("group",gNum));
            draw_node.push_back(tdNode);
            while(preLen < node_len){
                iNum += 1;
                Ndic tdNode;
                tdNode.insert(Ndic::value_type("id",iNum));
                tdNode.insert(Ndic::value_type("group",gNum));
                draw_node.push_back(tdNode);
                
                map<string,int> tdEdge;
                tdEdge.insert(map<string,int>::value_type("source",iNum-1));
                tdEdge.insert(map<string,int>::value_type("target",iNum));
                tdEdge.insert(map<string,int>::value_type("type",0));
                draw_edge.push_back(tdEdge);
                
                dnode_len.push_back(x_wCut);
                preLen += wCut;
            }
            iNum += 1;
            e_nid = iNum;
            Ndic tdNode2;
            tdNode2.insert(Ndic::value_type("id",iNum));
            tdNode2.insert(Ndic::value_type("group",gNum));
            draw_node.push_back(tdNode2);
            
            map<string,int> tdEdge;
            tdEdge.insert(map<string,int>::value_type("source",iNum-1));
            tdEdge.insert(map<string,int>::value_type("target",iNum));
            tdEdge.insert(map<string,int>::value_type("type",0));
            draw_edge.push_back(tdEdge);
            
            float trans = (wCut + node_len - preLen) * wPerK;
            dnode_len.push_back(trans);
        }        
              
        Nid tnid = {s_nid,e_nid,gNum};
        nid_dict.emplace(nk,tnid);
        gNum += 1;
    }
    
    for(auto &ed : r_edge_dict){
        if(ed.second > 1){     
            map<string,int> tdEdge;
            switch(ed.first.mark){
                case '2':
                    tdEdge.insert(map<string,int>::value_type("source",nid_dict[ed.first.node1].e_nid));
                    tdEdge.insert(map<string,int>::value_type("target",nid_dict[ed.first.node2].s_nid));
                    tdEdge.insert(map<string,int>::value_type("type",2));
                    draw_edge.push_back(tdEdge);
                    break;
                case '3':
                    tdEdge.insert(map<string,int>::value_type("source",nid_dict[ed.first.node1].e_nid));
                    tdEdge.insert(map<string,int>::value_type("target",nid_dict[ed.first.node2].e_nid));
                    tdEdge.insert(map<string,int>::value_type("type",3));
                    draw_edge.push_back(tdEdge);
                    break;
                case '4':
                    tdEdge.insert(map<string,int>::value_type("source",nid_dict[ed.first.node1].s_nid));
                    tdEdge.insert(map<string,int>::value_type("target",nid_dict[ed.first.node2].s_nid));
                    tdEdge.insert(map<string,int>::value_type("type",4));
                    draw_edge.push_back(tdEdge);
                    break;
                default:
                    tdEdge.insert(map<string,int>::value_type("source",nid_dict[ed.first.node1].s_nid));
                    tdEdge.insert(map<string,int>::value_type("target",nid_dict[ed.first.node2].e_nid));
                    tdEdge.insert(map<string,int>::value_type("type",5));
                    draw_edge.push_back(tdEdge);
                
            }
        }
    }
    if(asmb != ""){
        if(asmb[0] == '!'){
            string taskID = asmb.substr(1);
            string taskDir = upDir + "/mapping/" + taskID;
            int asmNum = -1;          
            hAsmNode(taskDir,refSim,info,asmNum,r_edge_dict,nid_dict);
        }else{
            if(asmb.find(",") != string::npos){
                set<string> nameSet;
                vector<int> asmNumVec;
                int t = 0;
                for(size_t i = 0; i < asmb.size(); ++i){
                    if(asmb[i] == ','){
                        int tlen = i - t;
                        nameSet.insert(asmb.substr(t,tlen));
                        t = i + 1;
                    }
                }
                nameSet.insert(asmb.substr(t));
                //
                getAsmNum2(nameSet,asmNumVec);
                //
                if(asmNumVec.empty()){
                    cerr<<"Error: assemblies can't be found. "<<asmb<<endl;
                }
                //
                if(indexFlag){
                    dxAsmNode2(refSim,asmNumVec,chrNum,sStart,sEnd,info,r_edge_dict,nid_dict);
                }else{
                    hAsmNode2(refSim,asmNumVec,r_edge_dict,nid_dict);
                }
            }else{
                int asmNum = getAsmNum(asmb);
                if(asmNum < 0){
                    cerr<<"Error: assembly can't be found. "<<asmb<<endl;
                    //exit(1);
                }else{
                    if(indexFlag){
                        dxAsmNode(refSim,asmNum,chrNum,sStart,sEnd,info,r_edge_dict,nid_dict);
                    }else{
                        string taskDir = "";
                        hAsmNode(taskDir,refSim,info,asmNum,r_edge_dict,nid_dict);
                    }
                }
            }
        }
    }
    //
    if(! refSim){
        string bwGeneFile = upDir + "/gene.ref.bw";
        string gDxFile = upDir + "/gene.ref.bdx";
        getFigGene(bwGeneFile,gDxFile,chrNum,sStart,sEnd,wPerK,exNode);
        
        string tkDesFile = upDir + "/track.info";
        string bwBedFile = upDir + "/bed.ref.bw";
        string bedDxFile = upDir + "/bed.ref.bdx";
        getFigBed(tkDesFile,bwBedFile,bedDxFile,chrNum,sStart,sEnd,wPerK,exNode);
    }
    
}

//---------------
void GraphRange::splitRange(int rangeNum,unordered_map<string,int> &chrMap,unordered_map<string,int> &refChrMap,string &rndDxFile,string &rndFile,string &nspecFile,string &snFile){
    ofstream dfh(rndDxFile.c_str());
    if(! dfh){
        cerr<<"Error: file open failed. "<<rndDxFile<<endl;
        exit(1);
    }
    
    ofstream nfh(rndFile.c_str());
    if(! nfh){
        cerr<<"Error: file open failed. "<<rndFile<<endl;
        exit(1);
    }
    
    ofstream nspfh(nspecFile.c_str());
    if(! nspfh){
        cerr<<"Error: file open failed. "<<nspecFile<<endl;
        exit(1);
    }
    
    bool sn = true;
    if(access(snFile.c_str(),F_OK) == 0){
        sn = false;
    }
    
    ifstream in(nodeFile.c_str());
    if(! in){
        cerr<<"Error: file open failed. "<<nodeFile<<endl;
        exit(1);
    }
    //
    string nodeLine;
    stringstream strStream;
    getline(in,nodeLine);
    int num = 0;
    string preChr = "";
    int rangeStart = 1,rangeEnd = 1;
    int intSize = sizeof(int);

    int oneSize = sizeof(OneRange);
    int crSize = sizeof(ChrRange);
    long long ndByte = 0,ndUnit = intSize * 3;

    int dxByte = intSize + (intSize + crSize) * refChrMap.size();
    int dxUnit = sizeof(OneRange);
    int nchr = refChrMap.size();
    dfh.write((char *)&nchr,intSize);
    //
    char *dummy = new char[dxByte-intSize];
    memset(dummy,0,dxByte-intSize);
    dfh.write(dummy,dxByte-intSize);
    delete []dummy;
    
    map<string,ChrRange> dxChrMap;
    int nLine = 0;
    //
    list<SANode> allNode;
    vector<string> chrVec;
    bool flag = false;
    //int rCount = 0;
    int nsp = 0;
    nspfh.write((char *)&nsp,intSize);
    while(getline(in,nodeLine)){
        strStream << nodeLine;
        NodeType r_node;
        string r_chr;
        int r_start,r_end,r_len,r_ref;
        strStream >> r_node;
        strStream >> r_chr;
        strStream >> r_start;
        strStream >> r_end;
        strStream >> r_len;
        strStream >> r_ref;
        
        if(sn){
            SANode xNode = {r_node,r_start,r_end,chrMap[r_chr]};
            allNode.push_back(xNode);
        }
        if(r_ref == 0){
            string tName = "",t_hap = "",tchr = "";
            asmSplit(r_chr,sep,tName,t_hap,tchr);
            //
            if(tchr != preChr){
                if(refChrMap.find(tchr) != refChrMap.end()){
                    chrVec.push_back(tchr);
                    if(preChr != ""){
                        //int tnum = num - 1;
                        int tnum = num;
                        /*
                        dfh.write((char *)&rangeStart,intSize);
                        dfh.write((char *)&rangeEnd,intSize);
                        dfh.write((char *)&ndByte,llSize);
                        dfh.write((char *)&tnum,intSize);
                        */
                        OneRange aRange = {rangeStart,rangeEnd,ndByte,tnum};
                        dfh.write((char *)&aRange,oneSize);
                        ++nLine;
                        ChrRange chrinfo = {dxByte,nLine};
                        dxChrMap.emplace(preChr,chrinfo);
                        // next start
                        dxByte += nLine * dxUnit;
                        ndByte += tnum * ndUnit;
                        //ndByte = ndUnit * rCount;
                    }
                    rangeStart = r_start;
                    rangeEnd = r_end;
                    num = 1;
                    //
                    nLine = 0;
                    //
                    flag = true;
                    preChr = tchr;
                }else{
                    flag = false;
                }
            }else{
                if(flag){
                    ++num;
                    if(num == rangeNum){
                        //ofh<<tchr<<"\t"<<rangeStart<<"\t"<<r_end<<endl;
                        /*
                        dfh.write((char *)&rangeStart,intSize);
                        dfh.write((char *)&rangeEnd,intSize);
                        dfh.write((char *)&ndByte,llSize);
                        dfh.write((char *)&num,intSize);
                        */
                        rangeEnd = r_end;
                        OneRange aRange = {rangeStart,rangeEnd,ndByte,num};
                        dfh.write((char *)&aRange,oneSize);
                        //
                        ++nLine;
                        //
                        ndByte += num * ndUnit;
                        //ndByte = ndUnit * (rCount + 1);
                        num = 0;
                    }else{
                        if(num == 1){
                            rangeStart = r_start;   
                        }
                        rangeEnd = r_end;
                    }
                }
            }
            //
            if(flag){
                nfh.write((char *)&r_node,intSize);
                nfh.write((char *)&r_start,intSize);
                nfh.write((char *)&r_end,intSize);
                //
                //++rCount;
            }else{
                nspfh.write((char *)&r_node,intSize);
                ++nsp;
            }
            //preChr = tchr;
        }
        strStream.clear();
        strStream.str("");
    }
    
    if(nsp > 0){
        nspfh.seekp(0,ios::beg);
        nspfh.write((char *)&nsp,intSize);
    }
    
    if(num > 0){
        /*
        dfh.write((char *)&rangeStart,intSize);
        dfh.write((char *)&rangeEnd,intSize);
        dfh.write((char *)&ndByte,llSize);
        dfh.write((char *)&num,intSize);
        */
        OneRange aRange = {rangeStart,rangeEnd,ndByte,num};
        dfh.write((char *)&aRange,oneSize);
        //
        ++nLine;
        ChrRange chrinfo = {dxByte,nLine};
        dxChrMap.emplace(preChr,chrinfo);
    }       
    dfh.seekp(intSize,ios::beg);
    for(auto &echr : chrVec){
        dfh.write((char *)&refChrMap[echr],intSize);
        dfh.write((char *)&dxChrMap[echr],crSize);
    }

    if(sn){
        ofstream sfh(snFile.c_str());
        if(! sfh){
            cerr<<"Error: file open failed. "<<snFile<<endl;
            exit(1);
        }
        allNode.sort();
        //sfh.seekp(0,ios::beg);
        int total = allNode.back().node;
        sfh.write((char *)&total,intSize);
        int pre = 0;
        int ndSize = sizeof(ANode);
        for(auto &gnode : allNode){
            /*
            sfh.write((char *)&gnode.start,intSize);
            sfh.write((char *)&gnode.pend,intSize);
            sfh.write((char *)&gnode.achr,intSize);
            */
            for(int i = pre + 1; i < gnode.node; ++i){
                ANode axNode = {0,0,0};
                sfh.write((char *)&axNode,ndSize);    
            }
            ANode axNode = {gnode.start,gnode.pend,gnode.achr};
            sfh.write((char *)&axNode,ndSize);
            pre = gnode.node;
        }
        sfh.close();
    }
    
    in.close();
    dfh.close();
    nfh.close();
    nspfh.close();
}

void GraphRange::getNrefEdge(string &rndFile,string &nspecFile,vector<NEdge> &resEdge){
    //string rndFile = pathDir + "/node.ref.bw";
    ifstream nfh(rndFile.c_str());
    if(! nfh){
        cerr<<"Error: file open failed. "<<nodeFile<<endl;
        exit(1);
    }

    unordered_set<int> ntNode;
    int r_node,r_start,r_end;
    int intSize = sizeof(int);
    while(nfh.read((char *)&r_node,intSize)){
        nfh.read((char *)&r_start,intSize);
        nfh.read((char *)&r_end,intSize);
        ntNode.insert(r_node);
    }
    nfh.close();
    
    ifstream nspfh(nspecFile.c_str());
    int nsp = 0;
    if(! nspfh){
        cerr<<"Warning: files in folder upload may be generated by old version gfa2view."<<endl;
    }else{
        nspfh.read((char *)&nsp,intSize);
        for(int k = 0; k < nsp; ++k){
            nspfh.read((char *)&r_node,intSize);
            ntNode.insert(r_node);
        }        
    }
    nspfh.close();
    //
    ifstream efh(edgeFile.c_str());
    if(! efh){
        cerr<<"Error: file open failed. "<< edgeFile<<endl;
        exit(1);
    }
    string strLine;
    stringstream strStream;
    getline(efh,strLine);
    while(getline(efh,strLine)){
        strStream << strLine;
        NodeType node1,node2;
        char sign1,sign2;
        strStream >> node1;
        strStream >> node2;
        strStream >> sign1;
        strStream >> sign2;
        if(ntNode.find(node1) == ntNode.end() && ntNode.find(node2) == ntNode.end()){
            char mark = getMark(sign1,sign2);
            NEdge tedge = {node1,node2,mark};
            resEdge.push_back(tedge);
        }
        
        strStream.clear();
        strStream.str("");
    }
    efh.close();
}

void GraphRange::getChrRmEdge(unordered_set<int> &ntNode,vector<NEdge> &chrRmEdge){
    ifstream efh(edgeFile.c_str());
    if(! efh){
        cerr<<"Error: file open failed. "<< edgeFile<<endl;
        exit(1);
    }
    string strLine;
    stringstream strStream;
    getline(efh,strLine);
    while(getline(efh,strLine)){
        strStream << strLine;
        NodeType node1,node2;
        char sign1,sign2;
        strStream >> node1;
        strStream >> node2;
        strStream >> sign1;
        strStream >> sign2;
        if(ntNode.find(node1) != ntNode.end() || ntNode.find(node2) != ntNode.end()){
            char mark = getMark(sign1,sign2);
            NEdge tedge = {node1,node2,mark};
            chrRmEdge.push_back(tedge);
        }
        
        strStream.clear();
        strStream.str("");
    }
    efh.close();
}

void GraphRange::parseRange(vector<RNode> &chrRnode,vector<OneRange> &arcVec,int sStart,int sEnd,int ex,vector<NodeType> &rangeNode,unordered_set<NodeType> &exNode){
    
    int eStart = sStart - ex > 0 ? sStart - ex : 1;
    int eEnd = sEnd + ex;
    //int rtotal = 0;
    bool flag = true;
    
    int pos = 0,posStart = 0;
    for(auto &tRange : arcVec){
        if(tRange.ranStart <= eStart){
            if(tRange.ranEnd >= eStart){
                if(flag){
                    posStart = pos;
                    flag = false;
                }
            }
        }else{
            if(tRange.ranStart <= eEnd){
                if(flag){
                    posStart = pos;
                    flag = false;
                }
            }else{
                break;
            }
            
        }
        pos += tRange.ranNum;
    }
    //
    if(flag){
        cerr<<"Error: nodes can't be found in the interval."<<sStart<<" - "<<sEnd<<endl;
        exit(1);
    }
    int r_node;
    int r_start,r_end;
    for(int s = posStart; s < pos; ++s){
        r_node = chrRnode[s].node;
        r_start = chrRnode[s].start;
        r_end = chrRnode[s].pend;
        
        if(r_start <= eStart){
            if(r_end >= eStart){
                exNode.insert(r_node);
                if(r_end >= sStart){
                    rangeNode.push_back(r_node);
                }
            }
        }else{
            if(r_start <= eEnd){
                exNode.insert(r_node);
                if(r_start <= sStart){
                    if(r_end >= sStart){
                        rangeNode.push_back(r_node);
                    }
                }else{
                    if(r_start <= sEnd){
                        rangeNode.push_back(r_node);
                    }
                }
            }else{
                break;
            }
            
        }
    }
}

void GraphRange::edgeRange(vector<RNode> &chrRnode,vector<OneRange> &arcVec,int sStart,int sEnd,int ex,int nocross,int storeDep,vector<NEdge> &chrRmEdge,unordered_map<NodeType,vector<ENode> > &iedge,unordered_map<NodeType,vector<ENode> > &oedge,set<NEdge> &r_edge_dict,unordered_set<NodeType> &nRefNode){
    vector<NodeType> rangeNode;
    unordered_set<NodeType> exNode;
    parseRange(chrRnode,arcVec,sStart,sEnd,ex,rangeNode,exNode);
    
    if(nocross){
        for(auto &tedge : chrRmEdge){
            if(exNode.find(tedge.node1) != exNode.end() || exNode.find(tedge.node2) != exNode.end()){
                conformEdge(tedge.node1,tedge.node2,tedge.mark,iedge,oedge);
            }
        }
    }
    
    unordered_set<NodeType> range_set;
    for(NodeType &rnode : rangeNode){
        range_set.insert(rnode);
    }
    unordered_map<NodeType,std::vector<ENode> >::iterator it;
    for(NodeType &tnode: rangeNode){
        vector<NodeType> tNref;
        vector<int> deep;
        it = oedge.find(tnode);
        if(it != oedge.end()){
            for(ENode &o_node : it->second){
                if(exNode.find(o_node.node) != exNode.end()){
                    //if(range_set.find(o_node.node) != range_set.end()){
                        NEdge sym = {tnode,o_node.node,o_node.mark};
                        
                        if(r_edge_dict.find(sym) == r_edge_dict.end()){
                            r_edge_dict.insert(sym);
                        }
                    //}
                }else{
                    if(nRefNode.find(o_node.node) == nRefNode.end()){
                        tNref.push_back(o_node.node);
                        deep.push_back(0);
                        nRefNode.insert(o_node.node);
                    }
                    NEdge sym = {tnode,o_node.node,o_node.mark};
                    r_edge_dict.insert(sym);
                }
            }
        }    
        //
        it = iedge.find(tnode);
        if(it != iedge.end()){
            for(ENode &o_node : it->second){
                if(exNode.find(o_node.node) == exNode.end()){
                    if(nRefNode.find(o_node.node) == nRefNode.end()){
                        tNref.push_back(o_node.node);
                        deep.push_back(0);
                        nRefNode.insert(o_node.node);
                    }
                    NEdge sym = {o_node.node,tnode,o_node.mark};
                    r_edge_dict.insert(sym);
                }else{
                    NEdge sym = {o_node.node,tnode,o_node.mark};
                    if(r_edge_dict.find(sym) == r_edge_dict.end()){
                        r_edge_dict.insert(sym);
                    }
                }
            }
        }
        //
        if(! tNref.empty()){
            size_t i = 0;
            while(i < tNref.size()){
                if(deep[i] > storeDep){
                    break;   
                }
                it = oedge.find(tNref[i]);
                if(it != oedge.end()){
                    for(ENode &to_node : it->second){
                        if(exNode.find(to_node.node) != exNode.end()){
                            if(range_set.find(to_node.node) != range_set.end()){
                                NEdge sym = {tNref[i],to_node.node,to_node.mark};
                                r_edge_dict.insert(sym);
                            }
                        }else{
                            if(nRefNode.find(to_node.node) == nRefNode.end()){
                                tNref.push_back(to_node.node);
                                deep.push_back(deep[i]+1);
                                nRefNode.insert(to_node.node);
                            }
                            NEdge sym = {tNref[i],to_node.node,to_node.mark};
                            r_edge_dict.insert(sym);
                        }
                    }
                }
                    
                it = iedge.find(tNref[i]);
                if(it != iedge.end()){
                    for(ENode &to_node : it->second){
                        if(exNode.find(to_node.node) != exNode.end()){
                            if(range_set.find(to_node.node) != range_set.end()){
                                NEdge sym = {to_node.node,tNref[i],to_node.mark};
                                r_edge_dict.insert(sym);
                            }
                        }else{
                            if(nRefNode.find(to_node.node) == nRefNode.end()){
                                tNref.push_back(to_node.node);
                                deep.push_back(deep[i]+1);
                                nRefNode.insert(to_node.node);
                            }
                            NEdge sym = {to_node.node,tNref[i],to_node.mark};
                            r_edge_dict.insert(sym);
                        }
                    }
                }

                i += 1;
            }
        }
    }
}

void GraphRange::fillNode(string &comChrFile,string &ndAsmFile,string &nrNodeFile,string &nrNumFile,string &snFile,string &nrdFile){
    ofstream nfh(nrdFile.c_str());
    if(! nfh){
        cerr<<"Error: file open failed. "<<nrdFile<<endl;
        exit(1);
    }
    
    ifstream acfh(comChrFile.c_str());
    ifstream mfh(ndAsmFile.c_str());
    ifstream rdfh(nrNodeFile.c_str());
    ifstream rnfh(nrNumFile.c_str());
    ifstream sfh(snFile.c_str());
    
    int total = 0;
    int intSize = sizeof(int);
    //int unit = intSize * 3;
    int unit = sizeof(ANode);
    
    string chrLine;
    unordered_map<int,string> chrMap;
    unordered_map<string,string> chrAsmMap;
    unordered_map<string,int> asmMap;
    
    int pos = 0;
    while(getline(acfh,chrLine)){
        chrMap.emplace(pos,chrLine);
        string tName = "",t_hap = "",tchr = "";
        asmSplit(chrLine,sep,tName,t_hap,tchr);
        chrAsmMap.emplace(chrLine,tName + sep + t_hap);
        ++pos;
    }
    acfh.close();
    //
    string jAsm = "Jump" + sep + "H";
    string uAsm = "Un" + sep + "H";
    string jComChr = jAsm + sep + "1";
    string uComChr = uAsm + sep + "1";
    chrMap.emplace(pos,jComChr);
    chrAsmMap.emplace(jComChr,jAsm);
    ++pos;
    chrMap.emplace(pos,uComChr);
    chrAsmMap.emplace(uComChr,uAsm);
    //
    pos = 0;
    while(getline(mfh,chrLine)){
        asmMap.emplace(chrLine,pos);
        ++pos;
    }
    mfh.close();
    //
    asmMap.emplace(jAsm,pos);
    ++pos;
    asmMap.emplace(uAsm,pos);
    //
    
    sfh.read((char *)&total,intSize);
    ANode *allNode = new ANode[total];
    sfh.read((char *)allNode,total * unit);
    
    total = 0;
    rnfh.read((char *)&total,intSize);
    for(int i = 0; i < total; ++i){
        int num = 0;
        rnfh.read((char *)&num,intSize);
        for(int k = 0; k < num; ++k){
            int tnode = 0,tstart,tpend,tasm;
            rdfh.read((char *)&tnode,intSize);
            tstart = allNode[tnode-1].start;
            tpend = allNode[tnode-1].pend;
            
            tasm = asmMap[chrAsmMap[chrMap[allNode[tnode-1].achr]]];
            
            int tlen = tpend - tstart + 1;
            nfh.write((char *)&tnode,intSize);
            nfh.write((char *)&tlen,intSize);
            nfh.write((char *)&tasm,intSize);
        }
    }
    
    delete []allNode;
    nfh.close();
    rdfh.close();
    rnfh.close();
    sfh.close();
}

void GraphRange::mergeDx(string &rndDxFile,string &nrNumFile,string &mgDxFile){
    ofstream mfh(mgDxFile.c_str());
    if(! mfh){
        cerr<<"Error: file open failed. "<<mgDxFile<<endl;
        exit(1);
    }
    //
    ifstream rdfh(rndDxFile.c_str());
    if(! rdfh){
        cerr<<"Error: file open failed. "<<rndDxFile<<endl;
        exit(1);
    }
    
    ifstream rnfh(nrNumFile.c_str());
    if(! rnfh){
        cerr<<"Error: file open failed. "<<nrNumFile<<endl;
        exit(1);
    }
    
    int intSize = sizeof(int);
    int nchr = 0;
    rdfh.read((char *)&nchr,intSize);
    
    mfh.write((char *)&nchr,intSize);
    //int refChr,offByte,nLine;
    int refChr;
    ChrRange crRange;
    int crSize = sizeof(ChrRange);
    if(nchr > 0){
        rdfh.read((char *)&refChr,intSize);
        rdfh.read((char *)&crRange,crSize);
        
        mfh.write((char *)&refChr,intSize);
        mfh.write((char *)&crRange,crSize);
    }else{
        cerr<<"Error: number of reference chromosome is 0."<<endl;
        exit(1);
    }
    int oneSize = sizeof(OneRange);
    int secSize =  sizeof(EdRang);
    int lineSize = oneSize + secSize;
    int preNum = crRange.ranNum;
    int preSize = crRange.rByte;
    for(int i = 1; i < nchr; ++i){
        rdfh.read((char *)&refChr,intSize);
        rdfh.read((char *)&crRange,crSize);
        
        crRange.rByte = preSize + preNum * lineSize;
        mfh.write((char *)&refChr,intSize);
        mfh.write((char *)&crRange,crSize);
        preNum += crRange.ranNum;
    }

    int total = 0;
    rnfh.read((char *)&total,intSize);
    long long nrByte = 0LL;
    
    OneRange aRange;
    EdRang sRange;
    for(int j = 0; j < total; ++j){
        int num = 0;
        rnfh.read((char *)&num,intSize);
        rdfh.read((char *)&aRange,oneSize);
        
        if(j == 0){
            mfh.write((char *)&aRange,oneSize);
            long long tOff = 0LL;
            sRange = {tOff,num};
            mfh.write((char *)&sRange,secSize);
        }else{
            mfh.write((char *)&aRange,oneSize);
            sRange = {nrByte,num};
            mfh.write((char *)&sRange,secSize);
        }
        nrByte += (long long)num * intSize * 3;
    }
    mfh.close();
    rdfh.close();
    rnfh.close();
}

// nodeA, nodeB not in the same range, add nodeB to the range the nodeA belong to
void GraphRange::rangePath(bool formR,int num,vector<string> &qCigarVec,vector<int> &qPosVec,vector<char> &orient,vector<NodeType> &nodes,unordered_map<NodeType,vector<int> > &ndCutMap,list<PathRang> &allPaRa){
    unordered_map<int,vector<int> > fragPos;
    unordered_map<NodeType,vector<int> >::iterator it;
    int preMin,preMax,curMin,curMax;
    bool flag = false;
    for(size_t j = 0; j < orient.size(); ++j){
        it = ndCutMap.find(nodes[j]);
        if(it != ndCutMap.end()){
            curMin = (it->second).front();
            curMax = (it->second).back();
            if(flag){
                if(curMin > preMax || curMax < preMin){
                    if(fragPos.find(curMin) != fragPos.end()){
                        fragPos[curMin].push_back(j - 1);
                    }else{
                        int frag = j - 1;
                        vector<int> fr{frag};
                        fragPos.emplace(curMin,fr);
                    } 
                }
            }
            for(auto x : it->second){
                if(fragPos.find(x) != fragPos.end()){
                    fragPos[x].push_back(j);
                }else{
                    vector<int> fr{(int)j};
                    fragPos.emplace(x,fr);
                }
            }
            //
            preMin = curMin;
            preMax = curMax;
            flag = true;
        }else{
            flag = false;   
        }
    }
    
    int firNode = 0;
    char firOri = '\0';
    list<LagNode> nodeCons;
    list<int> posList;
    list<string> cigarList;
    for(auto &nd : fragPos){
        int pre = -2;
        for(auto pos : nd.second){
            if(pos - pre > 1){
                //if(! nodeCons.empty()){
                if(pre >= 0){
                    PathRang tpRan = {nd.first,num,firNode,firOri,nodeCons,posList,cigarList};
                    allPaRa.push_back(tpRan);
                    
                    nodeCons.clear();
                    //
                    posList.clear();
                    if(formR){
                        cigarList.clear();
                        cigarList.push_back(qCigarVec[pos]);
                    }
                }else{
                    if(formR){
                        cigarList.push_back(qCigarVec[pos]);
                    }
                }
                firNode = nodes[pos];
                firOri = orient[pos];
                //
                posList.push_back(qPosVec[pos]);
                
            }else{
                int diff = nodes[pos] - firNode;
                if(diff > 127 || diff < -128){
                    PathRang tpRan = {nd.first,num,firNode,firOri,nodeCons,posList,cigarList};
                    allPaRa.push_back(tpRan);
                    nodeCons.clear();
                    
                    //
                    posList.clear();
                    if(formR){
                        cigarList.clear();
                        cigarList.push_back(qCigarVec[pos]);
                    }
                    //
                    firNode = nodes[pos];
                    firOri = orient[pos];
                    if(firOri == '>'){
                        firOri = '1';
                    }else{
                        firOri = '2';
                    }
                    //
                    posList.push_back(qPosVec[pos]);
                    //cigarList.push_back(qCigarVec[pos]);
                }else{
                    LagNode tLag = {(char)diff,orient[pos]};
                    nodeCons.push_back(tLag);
                    //
                    if(formR){
                        posList.push_back(qPosVec[pos]);
                        cigarList.push_back(qCigarVec[pos]);
                    }
                }
            }
            pre = pos;
        }
        //
        //if(! nodeCons.empty()){
        if(pre >= 0){
            PathRang tpRan = {nd.first,num,firNode,firOri,nodeCons,posList,cigarList};
            allPaRa.push_back(tpRan);
            
            nodeCons.clear();
            //
            posList.clear();
            if(formR){
                cigarList.clear();
            }
        }
    }
}

// allLen
void GraphRange::pthTask(bool formR,vector<int> &allLen,unordered_map<NodeType,vector<int> > &ndCutMap,vector<RanPos> &allpos,char *header,int dxByte,int frStart,int frEnd,vector<ifstream> &pthVec,vector<ofstream> &xpthVec,vector<ofstream> &wpthVec){
    
    int oneSize = sizeof(OneRange);
    int intSize = sizeof(int);
    int lndSize = sizeof(LagNode);
    string fullName,path;
    for(int k = frStart; k < frEnd; ++k){
        list<PathRang> allPaRa;
        int num = 0;
        while(pthVec[k] >> fullName){
            pthVec[k] >> path;
            //
            vector<char> orient;
            vector<NodeType> nodes;
            int snode = 0;
            for(size_t i = 0; i < path.length(); ++i){
                if(path[i] == '>' || path[i] == '<'){
                    orient.push_back(path[i]);
                    if(i > 0){
                        nodes.push_back(snode);
                        snode = 0;
                    }
                }else{
                    snode = snode * 10 + (path[i] - '0');
                }
            }
            nodes.push_back(snode);
            vector<string> qCigarVec;
            vector<int> qPosVec;
            // formR
            if(formR){
                int qStart,rStart;
                string cigar;
                pthVec[k] >> qStart >> rStart >> cigar;
                cigar2pos(qStart,rStart,allLen,cigar,nodes,qCigarVec,qPosVec);
            }else{
                char pType;
                int rStart;
                pthVec[k] >> pType;
                if(pType == 'W'){
                    pthVec[k] >>  rStart;
                }else{
                    rStart = 0;
                }
                path2pos(rStart,allLen,nodes,qPosVec);
            }
            rangePath(formR,num,qCigarVec,qPosVec,orient,nodes,ndCutMap,allPaRa);
            //
            ++num;
        }
        //
        int nr = 0;
        xpthVec[k].write(header,dxByte);
        
        if(! allPaRa.empty()){
            allPaRa.sort();
            long long offByte = 0LL, preOff = 0LL;
            int preFrag = allPaRa.front().frag;
            int num = 0;
            for(int x = 0; x < preFrag; ++x){
                OneRange abRange = {allpos[nr].start,allpos[nr].pend,preOff,num};
                xpthVec[k].write((char *)&abRange,oneSize);
                ++nr;
            }
            //
            if(formR){
                for(auto &tp : allPaRa){
                    list<string>::iterator it = tp.cigarList.begin();
                    string outStr;
                    if(tp.firOri == '1'){
                        outStr = ">";
                    }else if(tp.firOri == '2'){
                        outStr = "<";
                    }else{
                        outStr = string(1,tp.firOri);
                    }
                    outStr += (*it);
                    ++it;
                    
                    list<LagNode>::iterator itx = tp.lag.begin();
                    while(it != tp.cigarList.end()){
                        outStr += string(1,itx->ori) + (*it);
                        ++it;
                        ++itx;
                    }
                    int cigarSize = outStr.size();
                    //
                    if(tp.frag != preFrag){
                        OneRange abRange = {allpos[nr].start,allpos[nr].pend,preOff,num};
                        xpthVec[k].write((char *)&abRange,oneSize);
                        // next range start
                        num = 1 + tp.lag.size();
                        preOff = offByte;
                        //offByte += (intSize * 2 + 1 + tp.lag.size() * 2);
                        // int: pathNum, firNode, tp.lag.size(); char: firOri; 2 chars: tp.lag
                        offByte += (intSize * 3 + 1 + tp.lag.size() * 2 + intSize * (tp.posList.size() + 1) + cigarSize);
                        ++nr;
                        
                        for(int x = preFrag + 1; x < tp.frag; ++x){
                            int tnum = 0;
                            OneRange abRange = {allpos[nr].start,allpos[nr].pend,preOff,tnum};
                            xpthVec[k].write((char *)&abRange,oneSize);
                            ++nr;
                        }
                    }else{
                        num += (1 + tp.lag.size());
                        //offByte += (intSize * 2 + 1 + tp.lag.size() * 2);
                        offByte += (intSize * 3 + 1 + tp.lag.size() * 2 + intSize * (tp.posList.size() + 1) + cigarSize);
                    }
                    //
                    wpthVec[k].write((char *)&tp.pathNum,intSize);
                    wpthVec[k].write((char *)&tp.firNode,intSize);
                    wpthVec[k].write(&tp.firOri,1);
                    int tori = tp.lag.size();
                    wpthVec[k].write((char *)&tori,intSize);
                    for(auto &x : tp.lag){
                        wpthVec[k].write((char *)&x,lndSize);
                    }
                    //
                    for(int tpos : tp.posList){
                        wpthVec[k].write((char *)&tpos,intSize);
                    }
                    wpthVec[k].write((char *)&cigarSize,intSize);
                    wpthVec[k].write(outStr.c_str(),cigarSize);
                    //
                    preFrag = tp.frag;
                }
            }else{
                for(auto &tp : allPaRa){
                    if(tp.frag != preFrag){
                        OneRange abRange = {allpos[nr].start,allpos[nr].pend,preOff,num};
                        xpthVec[k].write((char *)&abRange,oneSize);
                        // next range start
                        num = 1 + tp.lag.size();
                        preOff = offByte;
                        //offByte += (intSize * 2 + 1 + tp.lag.size() * 2);
                        offByte += (intSize * 3 + 1 + tp.lag.size() * 2 + intSize * tp.posList.size());
                        ++nr;
                        
                        for(int x = preFrag + 1; x < tp.frag; ++x){
                            int tnum = 0;
                            OneRange abRange = {allpos[nr].start,allpos[nr].pend,preOff,tnum};
                            xpthVec[k].write((char *)&abRange,oneSize);
                            ++nr;
                        }
                    }else{
                        num += (1 + tp.lag.size());
                        //offByte += (intSize * 2 + 1 + tp.lag.size() * 2);
                        offByte += (intSize * 3 + 1 + tp.lag.size() * 2 + intSize * tp.posList.size());
                    }
                    //
                    wpthVec[k].write((char *)&tp.pathNum,intSize);
                    wpthVec[k].write((char *)&tp.firNode,intSize);
                    wpthVec[k].write(&tp.firOri,1);
                    int tori = tp.lag.size();
                    wpthVec[k].write((char *)&tori,intSize);
                    for(auto &x : tp.lag){
                        wpthVec[k].write((char *)&x,lndSize);
                    }
                    //
                    for(int tpos : tp.posList){
                        wpthVec[k].write((char *)&tpos,intSize);
                    }
                    //
                    preFrag = tp.frag;
                }
            }
            
            //
            OneRange abRange = {allpos[nr].start,allpos[nr].pend,preOff,num};
            xpthVec[k].write((char *)&abRange,oneSize);
            ++nr;
            //
            num = 0;
            for(size_t m = nr; m < allpos.size(); ++m){
                OneRange abRange = {allpos[nr].start,allpos[nr].pend,preOff,num};
                xpthVec[k].write((char *)&abRange,oneSize);
                ++nr;
            }
        }else{
            xpthVec[k].seekp(0,ios::beg);
            int tnchr = 0;
            xpthVec[k].write((char *)&tnchr,intSize);
        }
        
    }
    
    //efh.close();
    //bfh.close();
}

//
// ref, nonref len
// all node len    int *allLen = ; delete []allLen;
void GraphRange::getAllLen(string &snFile,vector<int> &allLen){    
    ifstream sfh(snFile.c_str());
    if(! sfh){
        cerr<<"Error: file open failed. "<<snFile<<endl;
        exit(1);
    }
    int total = 0;
    int intSize = sizeof(int);
    int unit = sizeof(ANode);
    sfh.read((char *)&total,intSize);
    ANode eachNode;
    allLen.reserve(total);
    for(int i = 0; i < total; ++i){
        sfh.read((char *)&eachNode,unit);
        int tlen = eachNode.pend - eachNode.start + 1;
        allLen.push_back(tlen);
    }
    sfh.close();
}

// [MIS=X] q move, [MDN=X] r move,   [H,P] not move
// vector<string> qCigarVec;vector<int> qPosVec;
void GraphRange::cigar2pos(int qStart,int rStart,vector<int> &allLen,string &cigar,vector<int> &nodes,vector<string> &qCigarVec,vector<int> &qPosVec){    
    int tvalue = 0;
    int rMove = 0, qMove = 0;
    int qPreMove = 0;
    int nextPos = 0;
    int firLen = allLen[nodes[0] - 1];
    int cumLen = firLen - rStart;
    //
    qPosVec.push_back(qStart);
    string tCigar = "";
    int i = 0;
    int tlen = 0;
    int preCum = cumLen;
    int n = nodes.size();
    for(char x : cigar){
        if(x >= '0' && x <= '9'){
            tvalue = tvalue * 10 + (x - '0');
        }else{
            switch(x){
                case 'M':
                    qMove += tvalue;
                    rMove += tvalue;
                    break;
                case 'X':
                    qMove += tvalue;
                    rMove += tvalue;
                    break;
                case '=':
                    qMove += tvalue;
                    rMove += tvalue;
                    break;
                case 'I':
                    qMove += tvalue;
                    break;
                case 'S':
                    qMove += tvalue;
                    break;
                case 'D':
                    rMove += tvalue;
                    break;
                case 'N':
                    rMove += tvalue;
                    break;
            }
            if(rMove < cumLen){
                tCigar += to_string(tvalue) + string(1,x);
            }else if(rMove == cumLen){
                tCigar += to_string(tvalue) + string(1,x);
                qCigarVec.push_back(tCigar);
                tCigar = "";
                //
                ++i;
                if(i < n){
                    tlen = allLen[nodes[i] - 1];
                    cumLen += tlen;
                    //
                    nextPos = qStart + qMove;
                    qPosVec.push_back(nextPos);
                }
            }else{
                // |-------|-------|---------|
                //      |----------------------|
                int nd = tvalue - (rMove - cumLen);
                //if(nd > 0){
                    tCigar += to_string(nd) + string(1,x);
                //}
                qCigarVec.push_back(tCigar);
                tCigar = "";
                // for last node cumLen >= rMove
                preCum = cumLen;
                ++i;
                // query changes relative to previous position
                int qd = nd;
                if(i < n){
                    if(qMove > qPreMove){
                        nextPos = qStart + qMove + nd;
                    }else{
                        nextPos = qStart + qMove;
                    }
                    qPosVec.push_back(nextPos);
                    //
                    tlen = allLen[nodes[i] - 1];
                    cumLen += tlen;
                    qd += tlen;
                }
                //if(i < nodes.size(){
                    while(1){
                        if(rMove < cumLen){
                            int td = rMove - preCum;
                            tCigar = to_string(td) + string(1,x);
                            break;
                        }else if(rMove == cumLen){
                            tCigar = to_string(tlen) + string(1,x);
                            qCigarVec.push_back(tCigar);
                            tCigar = "";
                            //
                            preCum = cumLen;
                            ++i;
                            if(i < n){
                                nextPos = qStart + qMove;
                                qPosVec.push_back(nextPos);
                                //
                                tlen = allLen[nodes[i] - 1];
                                cumLen += tlen;
                                //qd += tlen;
                            }
                            break;
                        }else{
                            tCigar = to_string(tlen) + string(1,x);
                            qCigarVec.push_back(tCigar);
                            tCigar = "";
                            //
                            preCum = cumLen;
                            ++i;
                            if(i < n){
                                if(qMove > qPreMove){
                                    nextPos = qStart + qPreMove + qd;
                                }else{
                                    nextPos = qStart + qMove;
                                }
                                qPosVec.push_back(nextPos);
                                //
                                tlen = allLen[nodes[i] - 1];
                                cumLen += tlen;
                                qd += tlen;
                            }
                        }
                    }
                //}
            }
            //
            tvalue = 0;
            qPreMove = qMove;
        }
    }
    //
    if(tCigar != ""){
        qCigarVec.push_back(tCigar);
    }
}


void GraphRange::path2pos(int rStart,vector<int> &allLen,vector<int> &nodes,vector<int> &qPosVec){
    int tpos = rStart;
    for(int x : nodes){        
        qPosVec.push_back(tpos);
        tpos += allLen[x-1];
    }
    
}
//
void GraphRange::indexPath(bool formR,string &asmFile,string &eIndexFile,string &bEdgeFile,string &snFile,int nthread){
    ifstream in(asmFile.c_str());
    if(! in){
        cerr<<"Error: file open failed. "<<asmFile<<endl;
        exit(1);
    }
    string asmLine;
    int nAsm = 0;
    while(getline(in,asmLine)){
        ++nAsm;
    }
    in.close();
    //
    ifstream efh(eIndexFile.c_str());
    if(! efh){
        cerr<<"Error: file open failed. "<<eIndexFile<<endl;
        exit(1);
    }
    ifstream bfh(bEdgeFile.c_str());
    if(! bfh){
        cerr<<"Error: file open failed. "<<bEdgeFile<<endl;
        exit(1);
    }
    //
    int nchr = 0;
    int intSize = sizeof(int);
    efh.read((char *)&nchr,intSize);

    //
    int refChr;
    ChrRange crRange;
    int crSize = sizeof(ChrRange);
    int total = 0;
    for(int x = 0; x < nchr; ++x){
        efh.read((char *)&refChr,intSize);
        efh.read((char *)&crRange,crSize);
        total += crRange.ranNum;
    }
    
    unordered_map<NodeType,vector<int> > ndCutMap;
    unordered_map<NodeType,vector<int> >::iterator it1,it2;
    
    vector<RanPos> allpos;
    allpos.reserve(total);
    int oneSize = sizeof(OneRange);
    int esize = sizeof(CEdge);
    for(int j = 0; j < total; ++j){
        OneRange aRange;
        efh.read((char *)&aRange,oneSize);

        RanPos trpos = {aRange.ranStart,aRange.ranEnd};
        allpos.push_back(trpos);
        
        unordered_set<NodeType> cutSet;
        for(int k = 0; k < aRange.ranNum; ++k){
            CEdge xedge;
            bfh.read((char *)&xedge,esize);
            if(cutSet.find(xedge.node1) == cutSet.end()){
                it1 = ndCutMap.find(xedge.node1);
                if(it1 != ndCutMap.end()){
                    (it1->second).push_back(j);
                }else{
                    vector<int> ndFr{j};
                    ndCutMap.emplace(xedge.node1,ndFr);
                }
                cutSet.insert(xedge.node1);
            }
            
            if(cutSet.find(xedge.node2) == cutSet.end()){
                it2 = ndCutMap.find(xedge.node2);
                if(it2 != ndCutMap.end()){
                    (it2->second).push_back(j);
                }else{
                    vector<int> ndFr{j};
                    ndCutMap.emplace(xedge.node2,ndFr);
                }
                cutSet.insert(xedge.node2);
            }
        }
    }
    bfh.close();
    //
    int dxByte = intSize + (intSize + crSize) * nchr;
    char *header = new char[dxByte];
    efh.clear();
    efh.seekg(0,ios::beg);
    efh.read(header,dxByte);
    efh.close();
    //
    int cthread = nthread;
    if(nthread > nAsm){
        cthread = nAsm;
    }
    int etask = nAsm / cthread;
    int redis = nAsm % cthread;

    vector<ifstream> pthVec;
    vector<ofstream> xpthVec;
    vector<ofstream> wpthVec;
    for(int i = 0; i < nAsm; ++i){
        string tpFile = pathDir + "/" + to_string(i) + ".path";
        string txpFile = pathDir + "/" + to_string(i) + ".path.bdx";
        string bwpFile = pathDir + "/" + to_string(i) + ".path.bw";
        
        pthVec.push_back(ifstream(tpFile.c_str()));
        xpthVec.push_back(ofstream(txpFile.c_str()));
        wpthVec.push_back(ofstream(bwpFile.c_str()));
    }
    vector<thread> pdVec;
    //
    vector<int> allLen;
    getAllLen(snFile,allLen);
    //
    for(int n = 0; n < cthread; ++n){
        int frStart,frEnd;
        if(n < redis){
            frStart = (etask + 1) * n;
            frEnd = frStart + (etask + 1);
        }else{
            if(n > 0){
                frStart = (etask + 1) * redis + etask * (n - redis);
            }else{
                frStart = 0;
            }
            frEnd = frStart + etask;
        }
        //
        pdVec.push_back(thread(&GraphRange::pthTask,this,formR,ref(allLen),ref(ndCutMap),ref(allpos),header,dxByte,frStart,frEnd,ref(pthVec),ref(xpthVec),ref(wpthVec))); 
    }
    //
    for(auto &td : pdVec){
        td.join();
    }
    delete []header;
    for(int j = 0; j < nAsm; ++j){
        pthVec[j].close();
        xpthVec[j].close();
        wpthVec[j].close();
    }
}

void GraphRange::oneTask(unordered_map<NodeType,vector<ENode> > &iedge,unordered_map<NodeType,vector<ENode> > &oedge,vector<RNode> &chrRnode,vector<OneRange> &acrVec,vector<NEdge> &chrRmEdge,
             int ex,int nocross,int frStart,int frEnd,int storeDep,ofstream &tndfh,ofstream &tbfh,int *frNrefNum,int *frEdgeNum
    ){
    unordered_map<NodeType,vector<ENode> > iRanEdge;
    unordered_map<NodeType,vector<ENode> > oRanEdge;
    int intSize = sizeof(int);
    int esize = sizeof(CEdge);
    for(int k = frStart; k < frEnd; ++k){
        set<NEdge> r_edge_dict;
        unordered_set<NodeType> nRefNode;
        if(nocross){
            iRanEdge = iedge;
            oRanEdge = oedge;
            // rStart, rEnd
            edgeRange(chrRnode,acrVec,acrVec[k].ranStart,acrVec[k].ranEnd,ex,nocross,storeDep,chrRmEdge,iRanEdge,oRanEdge,r_edge_dict,nRefNode);
        }else{
            edgeRange(chrRnode,acrVec,acrVec[k].ranStart,acrVec[k].ranEnd,ex,nocross,storeDep,chrRmEdge,iedge,oedge,r_edge_dict,nRefNode);
        }
        //
        int tnum = nRefNode.size();
        //nufh.write((char *)&tnum,intSize);
        frNrefNum[k] = tnum;
        for(auto &tnode : nRefNode){
            tndfh.write((char *)&tnode,intSize);
        }
        //
        for(auto &tedge : r_edge_dict){
            CEdge xedge;
            xedge.node1 = tedge.node1;
            xedge.node2 = tedge.node2;
            xedge.mark = tedge.mark;
            tbfh.write((char *)&xedge,esize);
        }
        
        int xnum = r_edge_dict.size();
        //
        frEdgeNum[k] = xnum;
    }
}

/*
  If a subset of chromosomes is selected to create index 'load.chr.list' will be created.
  For original graph in rGFA format 'node.asm.list' will be used.
*/
void GraphRange::edgeWrite(string &spChrFile,int rangeSize,int ex,int nocross,int nthread,int storeDep){
    
    
    string rndDxFile = upDir + "/node.ref.bdx";
    string rndFile = upDir + "/node.ref.bw";
    
    string nspecFile = upDir + "/node.nsp.ref.bw";
    
    string mgDxFile = upDir + "/node.merge.bdx";
    string nrdFile = upDir + "/node.nonref.bw";
    
    string snFile = upDir + "/node.sort.bw";
    
    string nrNodeFile = upDir + "/node.nref.id";
    string nrNumFile = upDir + "/node.nref.num";

    string bEdgeFile = upDir + "/edge.bw";
    string eIndexFile = upDir + "/edge.bdx";
    
    string depFile = upDir + "/index.dep";
    ofstream xdep(depFile.c_str());
    if(! xdep){
        cerr<<"Error: file open failed. "<<depFile<<endl;
        exit(1);
    }
    xdep<<storeDep<<endl;
    xdep.close();
    
    ofstream ndfh(nrNodeFile.c_str());
    if(! ndfh){
        cerr<<"Error: file open failed. "<<nrNodeFile<<endl;
        exit(1);
    }
    
    ofstream nufh(nrNumFile.c_str());
    if(! nufh){
        cerr<<"Error: file open failed. "<<nrNumFile<<endl;
        exit(1);
    }
    
    
    ofstream bfh(bEdgeFile.c_str());
    if(! bfh){
        cerr<<"Error: file open failed. "<<bEdgeFile<<endl;
        exit(1);
    }
    
    ofstream xfh(eIndexFile.c_str());
    if(! xfh){
        cerr<<"Error: file open failed. "<<eIndexFile<<endl;
        exit(1);
    }
    //
    string chrLine;
    unordered_map<string,int> chrMap;
    //unordered_map<string,string> chrAsmMap;
    int pos = 0;
    ifstream acfh(comChrFile.c_str());
    if(! acfh){
        cerr<<"Error: file open failed. "<<comChrFile<<endl;
        exit(1);
    }
    while(getline(acfh,chrLine)){
        chrMap.emplace(chrLine,pos);
        ++pos;
    }
    acfh.close();
    //
    string jAsm = "Jump" + sep + "H";
    string uAsm = "Un" + sep + "H";
    string jComChr = jAsm + sep + "1";
    string uComChr = uAsm + sep + "1";
    chrMap.emplace(jComChr,pos);
    ++pos;
    chrMap.emplace(uComChr,pos);
    //
    unordered_map<string,int> refChrMap;
    //vector<string> chrVec;
    pos = 0;
    ifstream cfh(chrFile.c_str());
    if(! cfh){
        cerr<<"Error: file open failed. "<<chrFile<<endl;
        exit(1);
    }
    //
    string loadChrFile = upDir + "/load.chr.list";
    if(spChrFile != "00000000"){
        set<string> speSet;
        ifstream spfh(spChrFile.c_str());
        if(! spfh){
            cerr<<"Error: file containing specific chromosomes for indexing open failed. "<<spChrFile<<endl;
            exit(1);
        }
        while(getline(spfh,chrLine)){
            speSet.insert(chrLine);
        }
        spfh.close();
        //
        ofstream ldfh(loadChrFile.c_str());
        while(getline(cfh,chrLine)){
            int tpos = chrLine.find("\t");
            string tchr = chrLine.substr(0,tpos);
            if(speSet.find(tchr) != speSet.end()){
                refChrMap.emplace(tchr,pos);
                ldfh<<chrLine<<endl;
            }
            //chrVec.push_back(tchr);
            ++pos;
        }
        ldfh.close();
    }else{
        while(getline(cfh,chrLine)){
            int tpos = chrLine.find("\t");
            string tchr = chrLine.substr(0,tpos);
            refChrMap.emplace(tchr,pos);
            //chrVec.push_back(tchr);
            ++pos;
        }
        if(access(loadChrFile.c_str(),F_OK) == 0){
            remove(loadChrFile.c_str());
        }
    }
    //
    cfh.close();
    if(refChrMap.empty()){
        cerr<<"Error: file is empty. "<<chrFile<<endl;
        exit(1);
    }
    
    int intSize = sizeof(int);
    int crSize = sizeof(ChrRange);
    splitRange(rangeSize,chrMap,refChrMap,rndDxFile,rndFile,nspecFile,snFile);
    //
    vector<NEdge> resEdge;
    //
    unordered_map<NodeType,vector<ENode> > iedge;
    unordered_map<NodeType,vector<ENode> > oedge;
    
    //
    if(nocross == 0){
        parseEdge(iedge,oedge);
    }else{
        resEdge.reserve(102400);
        getNrefEdge(rndFile,nspecFile,resEdge);
        for(auto &tedge : resEdge){
            conformEdge(tedge.node1,tedge.node2,tedge.mark,iedge,oedge);
        }
        resEdge.clear();
        vector<NEdge> ().swap(resEdge);
    }
    long long rSite = 0LL;
    
    int esize = sizeof(CEdge);
    
    //---------------------------------------------
    string preChr = "";
    vector<NEdge> chrRmEdge;
    chrRmEdge.reserve(10240);
    
    int total = 0;
    int nchr = 0;
    // split alread check
    ifstream rxfh(rndDxFile.c_str());
    ifstream rnfh(rndFile.c_str());
    rxfh.read((char *)&nchr,intSize);
    xfh.write((char *)&nchr,intSize);
    map<int,ChrRange> chrRanMap;
    
    vector<int> allchr;
    allchr.reserve(nchr);

    for(int t = 0; t < nchr; ++t){
        int tchr;
        ChrRange cRange;
        rxfh.read((char *)&tchr,intSize);
        rxfh.read((char *)&cRange,crSize);
        xfh.write((char *)&tchr,intSize);
        xfh.write((char *)&cRange,crSize);
        
        allchr.push_back(tchr);
        chrRanMap.emplace(tchr,cRange);
        total += cRange.ranNum;
    }
    nufh.write((char *)&total,intSize);
    int oneSize = sizeof(OneRange);
    for(int xchr : allchr){
        ChrRange cRange = chrRanMap[xchr];
        //
        vector<OneRange> acrVec;
        acrVec.reserve(cRange.ranNum);
        int chrNdNum = 0;
        for(int k = 0; k < cRange.ranNum; ++k){
            OneRange aRange;
            rxfh.read((char *)&aRange,oneSize);
            acrVec.push_back(aRange);
            chrNdNum += aRange.ranNum;
            //if(k == 0){
            //    chrNdOff = aRange.offByte;
            //}
        }
        //
        unordered_set<int> ntNode;
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
            ntNode.insert(node);
        }
        //
        if(nocross){
            chrRmEdge.clear();
            getChrRmEdge(ntNode,chrRmEdge);
        }
        
        cout<<"INFO -- index chromosome: "<<xchr<<endl;
        int fr = cRange.ranNum;
        int *frNrefNum = new int[fr];
        int *frEdgeNum = new int[fr];
        int cthread = nthread;
        if(nthread > fr){
            cthread = fr;
        }
        int etask = fr / cthread;
        int redis = fr % cthread;

        vector<thread> thVec;
        vector<ofstream> nrFhVec;
        vector<ofstream> edFhVec;
        for(int p = 0; p < cthread; ++p){
            string tNdFile = upDir + "/" + to_string(p) + ".nr.part";
            string tbFile = upDir + "/" + to_string(p) + ".ed.part";
            nrFhVec.push_back(ofstream(tNdFile.c_str()));
            edFhVec.push_back(ofstream(tbFile.c_str()));
        }
        for(int n = 0; n < cthread; ++n){
            int frStart,frEnd;
            if(n < redis){
                frStart = (etask + 1) * n;
                frEnd = frStart + (etask + 1);
            }else{
                if(n > 0){
                    frStart = (etask + 1) * redis + etask * (n - redis);
                }else{
                    frStart = 0;
                }
                frEnd = frStart + etask;
            }
            thVec.push_back(thread(&GraphRange::oneTask,this,ref(iedge),ref(oedge),ref(chrRnode),ref(acrVec),ref(chrRmEdge),ex,nocross,frStart,frEnd,storeDep,ref(nrFhVec[n]),ref(edFhVec[n]),frNrefNum,frEdgeNum));
        }

        for(auto &th : thVec){
            th.join();
        }

        for(int x = 0; x < cthread; ++x){
            nrFhVec[x].close();
            edFhVec[x].close();
            // File size for each part should be <= std::numeric_limits<int>::max, 
            // Empirically, human pangenome graph with more than 1000 assemblies can still be handled by this step. 
            string tNdFile = upDir + "/" + to_string(x) + ".nr.part";
            ifstream tf(tNdFile.c_str());
            tf.seekg(0,ios::end);
            int nrSize = tf.tellg();
            tf.seekg(0,ios::beg);
            char *nrNode = new char[nrSize];
            tf.read(nrNode,nrSize);
            ndfh.write(nrNode,nrSize);
            delete []nrNode;
            tf.close();
            
            string tbFile = upDir + "/" + to_string(x) + ".ed.part";
            ifstream kf(tbFile.c_str());
            kf.seekg(0,ios::end);
            int edSize = kf.tellg();
            kf.seekg(0,ios::beg);
            char *ed = new char[edSize];
            kf.read(ed,edSize);
            bfh.write(ed,edSize);
            delete []ed;
            kf.close();
            //
            remove(tNdFile.c_str());
            remove(tbFile.c_str());
        }

        for(int m = 0; m < fr; ++m){
            nufh.write((char *)&frNrefNum[m],intSize);
            
            OneRange abRange = {acrVec[m].ranStart,acrVec[m].ranEnd,rSite,frEdgeNum[m]};
            xfh.write((char *)&abRange,oneSize);
            rSite += esize * frEdgeNum[m];
        }
        //
        delete []frNrefNum;
        delete []frEdgeNum;
    }
    //
    rxfh.close();
    rnfh.close();
    //in.close();
    bfh.close();
    xfh.close();
    nufh.close();
    ndfh.close();
    //--------------------------------
    string ndAsmFile = upDir + "/node.asm.list";
    if(access(ndAsmFile.c_str(),F_OK) == 0 ){
        fillNode(comChrFile,ndAsmFile,nrNodeFile,nrNumFile,snFile,nrdFile);
    }else{
        fillNode(comChrFile,asmFile,nrNodeFile,nrNumFile,snFile,nrdFile);
    }
    mergeDx(rndDxFile,nrNumFile,mgDxFile);
    //
    bool formR = true;
    string formFile = upDir + "/form.info";
    if(access(formFile.c_str(),F_OK) == 0){
        formR = false;
    }
    indexPath(formR,asmFile,eIndexFile,bEdgeFile,snFile,nthread);
    //
    remove(nrNodeFile.c_str());
    remove(nrNumFile.c_str());
}

//-----------------------------------
QueryNode::QueryNode(string &t_dbDir):dbDir(t_dbDir){
    
    asmFile = dbDir + "/asm.list";
    comChrFile = dbDir + "/complete.chr.list";
    sepFile = dbDir + "/sep.info";
    
    sep = getSep(sepFile);
}

int QueryNode::countHeader(){
    ifstream in(asmFile.c_str());
    if(! in){
        cerr<<"Error: file open failed. "<<asmFile<<endl;
        exit(1);
    }
    int pos = 0;
    string asmLine;
    while(getline(in,asmLine)){
        //if(asmLine = ""){
        //    continue;
        //}
        ++pos;
    }
    in.close();
    return pos;
}

void QueryNode::getHeader(){
    ifstream in(asmFile.c_str());
    if(! in){
        cerr<<"Error: file open failed. "<<asmFile<<endl;
        exit(1);
    }
    string asmLine;
    while(getline(in,asmLine)){
        header.push_back(asmLine);
    }
    in.close();
}

void QueryNode::queryDbNode(int node){
    
    ifstream cfh(comChrFile);
    if(! cfh){
        cerr<<"Error: file open failed. "<<comChrFile<<endl;
        exit(1);
    }
    unordered_map<int,string> comChrMap;
    string chrLine;
    int pos = 0;
    while(getline(cfh,chrLine)){
        comChrMap.emplace(pos,chrLine);
        ++pos;
    }
    cfh.close();
    //
    string jAsm = "Jump" + sep + "H";
    string uAsm = "Un" + sep + "H";
    string jComChr = jAsm + sep + "1";
    string uComChr = uAsm + sep + "1";
    comChrMap.emplace(pos,jComChr);
    ++pos;
    comChrMap.emplace(pos,uComChr);
    //
    string bNodeFile = dbDir + "/node.sort.bw";
    ifstream in(bNodeFile.c_str());
    if(! in){
        cerr<<"Error: file open failed. "<<bNodeFile<<endl;
        exit(1);
    }
    int ndSize = sizeof(ANode);
    long long offByte = sizeof(int) + (long long)(node - 1) * ndSize;
    in.seekg(offByte,ios::beg);
    ANode tnode;
    in.read((char *)&tnode,ndSize);
    nodeStart = tnode.start;
    nodeEnd = tnode.pend;
    string fullName = comChrMap[tnode.achr];
    string tName = "",t_hap = "",tchr = "";
    asmSplit(fullName,sep,tName,t_hap,tchr);
    nodeAsm = tName + sep + t_hap;
    nodeChr = tchr;
    in.close();
}
//
void QueryNode::queryGene(int node,string &nodeAsm){
    string annoDxFile = dbDir + "/anno.bdx";
    string annoNumFile = dbDir + "/anno.num";
    
    string asmFile = dbDir + "/asm.list";
    string asmFile2 = dbDir + "/node.asm.list";
    if(access(asmFile2.c_str(),F_OK) == 0){
        asmFile = asmFile2;
    }
    int pos = 0;
    ifstream afh(asmFile.c_str());
    string line;
    while(getline(afh,line)){
        if(line == nodeAsm){
            break;
        }
        ++pos;
    }
    afh.close();
    
    string annoFile = dbDir + "/anno/" + to_string(pos) + ".anno.bw";     
    ifstream in(annoDxFile.c_str()); 
    if(! in){
        cerr<<"Warning: file open failed. "<<annoDxFile<<endl;
        return;
    }
    in.seekg(sizeof(AnnoDx) * (node - 1),ios::beg);
    AnnoDx adx;
    in.read((char *)&adx,sizeof(AnnoDx));
    //cout<<node<<" "<<adx.offset<<" "<<adx.num<<endl;    
    ifstream mfh(annoNumFile.c_str());
    ifstream nfh(annoFile.c_str());
    if(! nfh){
        cerr<<"Warning: file open failed. "<<annoFile<<endl;
        return;
    }
    long long offset = adx.offset;
    mfh.seekg(offset,ios::beg);
    vector<int> gPos;
    int intSize = sizeof(int);
    int uSize = sizeof(AnnoLine);   
    for(int i = 0; i < adx.num; ++i){       
        int tnum;
        mfh.read((char *)&tnum,intSize);       
        long long toffset = uSize * tnum + intSize;
        nfh.seekg(toffset,ios::beg);
        AnnoLine annoInfo;
        nfh.read((char *)&annoInfo,uSize);
        vector<string> geneStr;
        geneStr.push_back(annoInfo.seqid);
        geneStr.push_back(to_string(annoInfo.start));
        geneStr.push_back(to_string(annoInfo.end));
        geneStr.push_back(annoInfo.geneID);
        geneStr.push_back(annoInfo.geneName);
        geneStr.push_back(string(1,annoInfo.strand));        

        geneList.push_back(geneStr);

    }
    in.close();
    mfh.close();
    nfh.close();
}

void QueryNode::fetchNdSeq(int node){
    string seqFile = dbDir + "/graphSeq.fa";
    string dxFile = dbDir + "/graphSeq.fa.bdx";
    nodeSeq = "";
    ifstream in(seqFile.c_str());
    if(in){
        ifstream dxfh(dxFile.c_str());
        int uSize = sizeof(SeqDx);
        long long dxOffset = (long long)(node - 1) * uSize;
        dxfh.seekg(dxOffset);
        SeqDx tdx;
        dxfh.read((char *)&tdx,uSize);
        if(tdx.offset == 0){            
            nodeSeq.append(tdx.len,'N');
        }else{           
            in.seekg(tdx.offset);           
            char *seq = new char[tdx.len+1];
            in.read(seq,tdx.len);
            seq[tdx.len] = '\0';
            
            nodeSeq = seq;
            //
            delete []seq;
        }
        in.close();
        dxfh.close();
    }
}
//
void QueryNode::queryDbCov(int node){
    string dxCovFile = dbDir + "/cover.bdx";
    string covFile = dbDir + "/cover.bw";
    ifstream xfh(dxCovFile.c_str());
    if(! xfh){
        cerr<<"Error: file open failed. "<<dxCovFile<<endl;
        exit(1);
    }
    ifstream vfh(covFile.c_str());
    if(! vfh){
        cerr<<"Error: file open failed. "<<covFile<<endl;
        exit(1);
    }
    
    int llSize = sizeof(long long);
    int usintSize = sizeof(unsigned short int);
    long long unit = llSize + usintSize * 2;
    long long offByte = (long long)(node - 1) * unit;
    xfh.seekg(offByte);
    long long tByte;
    unsigned short int num1,num2;
    xfh.read((char *)&tByte,llSize);
    xfh.read((char *)&num1,usintSize);
    xfh.read((char *)&num2,usintSize);
    
    vfh.seekg(tByte,ios::beg);
    int asmNum = countHeader();
    ndCov.reserve(asmNum);
    
    if(num1 < asmNum){
        for(int x = 0; x < asmNum; ++x){
            ndCov.push_back(0);
        }
        for(unsigned short int i = 0; i < num1; ++i){
            unsigned short int tpos;
            vfh.read((char *)&tpos,usintSize);
            ndCov[tpos] = 1;
        }
        
        vector<unsigned short int> posVec;
        vector<unsigned short int> valVec;
        for(unsigned short int j = 0; j < num2; ++j){
            unsigned short int tpos;
            vfh.read((char *)&tpos,usintSize);
            posVec.push_back(tpos);
        }
        for(unsigned short int k = 0; k < num2; ++k){
            unsigned short int value;
            vfh.read((char *)&value,usintSize);
            valVec.push_back(value);
        }
        for(size_t w = 0; w < posVec.size(); ++w){
            ndCov[posVec[w]] = (int)valVec[w];
        }
    }else{
        for(int x = 0; x < asmNum; ++x){
            unsigned short int value;
            vfh.read((char *)&value,usintSize);
            ndCov.push_back((int)value);
        }
    }
    xfh.close();
    vfh.close();
}

void QueryNode::queryAsmCov(vector<int> &nodeVec,string &asmb){
    string dxCovFile = dbDir + "/cover.bdx";
    string covFile = dbDir + "/cover.bw";
    ifstream xfh(dxCovFile.c_str());
    if(! xfh){
        cerr<<"Error: file open failed. "<<dxCovFile<<endl;
        exit(1);
    }
    ifstream vfh(covFile.c_str());
    if(! vfh){
        cerr<<"Error: file open failed. "<<covFile<<endl;
        exit(1);
    }
    
    int llSize = sizeof(long long);
    int usintSize = sizeof(unsigned short int);
    long long unit = llSize + usintSize * 2;
    
    getHeader();
    int asmNum = header.size();
    int asmPos = 0;
    for(int t = 0; t < asmNum; ++t){
        if(header[t] == asmb){
            asmPos = t;
            break;
        }
    }
    //
    ndCov.reserve(nodeVec.size());
    bool flag = false;
    for(int node : nodeVec){
        long long offByte = (long long)(node - 1) * unit;
        xfh.clear();
        xfh.seekg(offByte);
        long long tByte;
        unsigned short int num1,num2;
        xfh.read((char *)&tByte,llSize);
        xfh.read((char *)&num1,usintSize);
        xfh.read((char *)&num2,usintSize);
        
        vfh.clear();
        vfh.seekg(tByte,ios::beg);
        //
        flag = false;
        if(num1 < asmNum){
            for(unsigned short int i = 0; i < num1; ++i){
                unsigned short int tpos;
                vfh.read((char *)&tpos,usintSize);
                if(tpos == asmPos){
                    ndCov.push_back(1);
                    flag = true;
                    break;
                }
            }
            
            if(flag){
                continue;
            }
            vector<unsigned short int> posVec;
            int p = 0;
            for(unsigned short int j = 0; j < num2; ++j){
                unsigned short int tpos;
                vfh.read((char *)&tpos,usintSize);
                if(tpos == asmPos){
                    flag = true;
                    p = j;
                }
                posVec.push_back(tpos);
            }
            if(flag){
                for(unsigned short int k = 0; k < num2; ++k){
                    unsigned short int value;
                    vfh.read((char *)&value,usintSize);
                    if(k == p){
                        ndCov.push_back((int)value);
                        break;
                    }
                }
            }else{
                cerr<<"Warning: node : "<<node<<" may be not in the assembly."<<endl;
                ndCov.push_back(0);
            }
            
        }else{
            for(int x = 0; x < asmNum; ++x){
                unsigned short int value;
                vfh.read((char *)&value,usintSize);
                if(x == asmPos){
                    ndCov.push_back((int)value);
                    break;
                }
                
            }
        }
    }
    
    xfh.close();
    vfh.close();
}

#ifdef PYMODULE
PYBIND11_MODULE(minipg,m){
    py::class_<GraphRange>(m,"GraphRange")
        .def(py::init<string &,int>())
        .def("formatGraph",&GraphRange::formatGraph)
        .def("edgeWrite",&GraphRange::edgeWrite)
        .def_readwrite("draw_node",&GraphRange::draw_node)
        .def_readwrite("draw_pos",&GraphRange::draw_pos)
        .def_readwrite("draw_edge",&GraphRange::draw_edge)
        .def_readwrite("dnode_len",&GraphRange::dnode_len)
        .def_readwrite("genome",&GraphRange::genome)
        .def_readwrite("nnames",&GraphRange::nnames)
        .def_readwrite("hnGroup",&GraphRange::hnGroup)
        .def_readwrite("hLinks",&GraphRange::hLinks)
        .def_readwrite("hDir",&GraphRange::hDir)
        .def_readwrite("hEdgeAsm",&GraphRange::hEdgeAsm)
        .def_readwrite("ndGenePos",&GraphRange::ndGenePos)
        .def_readwrite("geneVec",&GraphRange::geneVec)
        .def_readwrite("layerVec",&GraphRange::layerVec)
        .def_readwrite("strandVec",&GraphRange::strandVec)
        .def_readwrite("mgFlagVec",&GraphRange::mgFlagVec)
        .def_readwrite("ndExonPos",&GraphRange::ndExonPos)
        .def_readwrite("rnaVec",&GraphRange::rnaVec)
        .def_readwrite("eLayerVec",&GraphRange::eLayerVec)
        .def_readwrite("eStrandVec",&GraphRange::eStrandVec)
        .def_readwrite("eNumVec",&GraphRange::eNumVec)
        .def_readwrite("eFlagVec",&GraphRange::eFlagVec)
        .def_readwrite("ndCDSPos",&GraphRange::ndCDSPos)
        .def_readwrite("cdsVec",&GraphRange::cdsVec)
        .def_readwrite("cLayerVec",&GraphRange::cLayerVec)
        .def_readwrite("cNumVec",&GraphRange::cNumVec)
        .def_readwrite("tkNameVec",&GraphRange::tkNameVec)
        .def_readwrite("tkDesVec",&GraphRange::tkDesVec)
        .def_readwrite("tkColVec",&GraphRange::tkColVec)
        .def_readwrite("tkCumVec",&GraphRange::tkCumVec)
        .def_readwrite("tkItem",&GraphRange::tkItem)
        .def_readwrite("rBedPos",&GraphRange::rBedPos)
        .def_readwrite("rBedName",&GraphRange::rBedName)
        .def_readwrite("rBedLayer",&GraphRange::rBedLayer)
        .def_readwrite("rBedScore",&GraphRange::rBedScore)
        .def_readwrite("rBedStrand",&GraphRange::rBedStrand)
        .def_readwrite("tkItem",&GraphRange::tkItem)
        .def_readwrite("figScale",&GraphRange::figScale)
        .def_readwrite("tickValue",&GraphRange::tickValue)
        .def_readwrite("tickPos",&GraphRange::tickPos)
        .def_readwrite("qChr",&GraphRange::qChr)
        .def_readwrite("qStart",&GraphRange::qStart)
        .def_readwrite("qEnd",&GraphRange::qEnd)
        .def_readwrite("qPath",&GraphRange::qPath)
        .def_readwrite("qCigar",&GraphRange::qCigar);
                
    py::class_<QueryNode>(m,"QueryNode")
        .def(py::init<string &>())
        .def("fetchNdSeq",&QueryNode::fetchNdSeq)
        .def("queryDbNode",&QueryNode::queryDbNode)
        .def("queryGene",&QueryNode::queryGene)
        .def("queryDbCov",&QueryNode::queryDbCov)
        .def("queryAsmCov",&QueryNode::queryAsmCov)
        //.def_readwrite("header",&QueryNode::header)
        .def_readwrite("geneList",&QueryNode::geneList)
        .def_readwrite("ndCov",&QueryNode::ndCov)
        .def_readwrite("nodeAsm",&QueryNode::nodeAsm)
        .def_readwrite("nodeChr",&QueryNode::nodeChr)
        .def_readwrite("nodeStart",&QueryNode::nodeStart)
        .def_readwrite("nodeEnd",&QueryNode::nodeEnd)
        .def_readwrite("nodeSeq",&QueryNode::nodeSeq);
}

#endif

