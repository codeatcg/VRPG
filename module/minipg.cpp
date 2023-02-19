

#include <fstream>
#include <iostream>
#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "vgraph.h"

using namespace std;
namespace py = pybind11;


GraphRange::GraphRange(string &t_nodeFile,string &t_edgeFile,string &t_pathFile,
                       string &t_sepFile
                       ):nodeFile(t_nodeFile),edgeFile(t_edgeFile),pathFile(t_pathFile),sepFile(t_sepFile){
    indexFlag = 0;
    getSep();
    //
    draw_node.reserve(64);
    draw_pos.reserve(64);
    draw_edge.reserve(64);
    genome.reserve(64);
    nnames.reserve(64);
    hnGroup.reserve(64);
    hLinks.reserve(64);
}

GraphRange::GraphRange(string &t_nodeFile,string &t_edgeFile,string &t_pathFile,
                       string &t_sepFile,string &t_chrFile,string &t_bEdgeFile,string &t_eIndexFile
                      ):nodeFile(t_nodeFile),edgeFile(t_edgeFile),pathFile(t_pathFile),sepFile(t_sepFile),chrFile(t_chrFile),bEdgeFile(t_bEdgeFile),eIndexFile(t_eIndexFile){
    indexFlag = 1;
    getSep();
    //
    draw_node.reserve(64);
    draw_pos.reserve(64);
    draw_edge.reserve(64);
    genome.reserve(64);
    nnames.reserve(64);
    hnGroup.reserve(64);
    hLinks.reserve(64);
}

void GraphRange::getSep(){
    ifstream in(sepFile.c_str());
    getline(in,sep),
    in.close();
}

void GraphRange::conformEdge(string &node1,string &node2,char mark){
    ENode x = {node2,mark};
    if(edge.find(node1) != edge.end()){
        map<string,vector<ENode> > &t = edge[node1];
        if(t.find("o") != t.end()){
            t["o"].push_back(x);
        }else{
            vector<ENode> eV;
            eV.push_back(x);
            t.insert(map<string,vector<ENode> >::value_type("o",eV));
        }
    }else{
        vector<ENode> eV;
        eV.push_back(x);
        map<string,vector<ENode> > neo;
        neo.insert(map<string,vector<ENode> >::value_type("o",eV));
        edge.insert(map<string, map<string,vector<ENode> > >::value_type(node1,neo));
    }
    //
    ENode y = {node1,mark};
    if(edge.find(node2) != edge.end()){
        map<string,vector<ENode> > &t = edge[node2];
        if(t.find("i") != t.end()){
            t["i"].push_back(y);
        }else{
            vector<ENode> eV;
            eV.push_back(y);
            t.insert(map<string,vector<ENode> >::value_type("i",eV));
        }
    }else{
        vector<ENode> eV;
        eV.push_back(y);
        map<string,vector<ENode> > neo;
        neo.insert(map<string,vector<ENode> >::value_type("i",eV));
        edge.insert(map<string, map<string,vector<ENode> > >::value_type(node2,neo));
    }
}

void GraphRange::parseEdge(){
    ifstream in(edgeFile.c_str());
    string strLine;
    getline(in,strLine);
    stringstream strStream;
    while(getline(in,strLine)){
        strStream << strLine;
        string node1,node2;
        char sign1,sign2;
        strStream >> node1;
        strStream >> node2;
        strStream >> sign1;
        strStream >> sign2;
        
        char mark = getMark(sign1,sign2);
        conformEdge(node1,node2,mark);
        
        strStream.str("");
        strStream.clear();
        
        
    }
    in.close();
}

void GraphRange::parseIndex(string &sChr,int sStart,int sEnd,int ex){
    ifstream in(eIndexFile.c_str());
    ifstream xfh(bEdgeFile.c_str());
    string strLine;
    stringstream strStream;
    bool flag = false;
    int esize = NODESIZE * 2 + 1;
    set<NEdge> redgeSet;
    while(getline(in,strLine)){
        strStream << strLine;
        string chrName;
        int rStart,rEnd,rOffsize,rCount;
        strStream >> chrName;
        strStream >> rStart;
        strStream >> rEnd;
        strStream >> rOffsize;
        strStream >> rCount;
        
        if(chrName == sChr){
            if(sStart >= rStart){
                if(sStart <= rEnd){
                    xfh.seekg(rOffsize,ios::beg);
                    for(int i = 0; i < rCount; ++i){
                        CEdge xedge;
                        xfh.read((char *)&xedge,esize);
                        string tnode1 = xedge.node1;
                        string tnode2 = xedge.node2;
                        NEdge tedge = {tnode1,tnode2,xedge.mark};
                        if(redgeSet.find(tedge) == redgeSet.end()){
                            conformEdge(tnode1,tnode2,xedge.mark);
                            redgeSet.insert(tedge); 
                        }
                    }
                    flag = true;
                }
            }else{
                if(sEnd >= rStart){
                    xfh.seekg(rOffsize,ios::beg);
                    for(int i = 0; i < rCount; ++i){
                        CEdge xedge;
                        xfh.read((char *)&xedge,esize);
                        string tnode1 = xedge.node1;
                        string tnode2 = xedge.node2;
                        NEdge tedge = {tnode1,tnode2,xedge.mark};
                        if(redgeSet.find(tedge) == redgeSet.end()){
                            conformEdge(tnode1,tnode2,xedge.mark);
                            redgeSet.insert(tedge); 
                        }
                    }
                    flag = true;
                }
            }
        }else{
            if(flag){
                break;
            }
        }
        strStream.str("");
        strStream.clear();
    }
    in.close();
    xfh.close();
}

void GraphRange::parseNode(string &sChr,int sStart,int sEnd,int ex,vector<string> &rangeNode,set<string> &exNode,map<string,int> &info,map<string,string> &assDict,int &realLen){
    int eStart = sStart - ex > 0 ? sStart - ex : 1;
    int eEnd = sEnd + ex;
    
    ifstream in(nodeFile.c_str());
    string nodeLine;
    stringstream strStream;
    getline(in,nodeLine);

    while(getline(in,nodeLine)){
        strStream << nodeLine;
        string r_node,r_chr;
        int r_start,r_end,r_len,r_ref;
        strStream >> r_node;
        strStream >> r_chr;
        strStream >> r_start;
        strStream >> r_end;
        strStream >> r_len;
        strStream >> r_ref;
        
        string tName,t_hap,tchr;
        assSplit(r_chr,sep,tName,t_hap,tchr);
        if(r_ref == 0){
            if(tchr == sChr){
                if(r_start <= eStart){
                    if(r_end >= eStart){
                        exNode.insert(r_node);
                        info.insert(map<string,int>::value_type(r_node,r_len));
                        assDict.insert(map<string,string>::value_type(r_node,tName));
                        if(r_end >= sStart){
                            rangeNode.push_back(r_node);
                            realLen += r_len;
                        }
                    }
                }else{
                    if(r_start <= eEnd){
                        exNode.insert(r_node);
                        info.insert(map<string,int>::value_type(r_node,r_len));
                        assDict.insert(map<string,string>::value_type(r_node,tName));
                        if(r_start <= sStart){
                            if(r_end >= sStart){
                                rangeNode.push_back(r_node);
                                realLen += r_len;
                            }
                        }else{
                            if(r_start <= sEnd){
                                rangeNode.push_back(r_node);
                                realLen += r_len;
                            }
                        }
                        
                    }
                    
                }
            }
        }else{
            info.insert(map<string,int>::value_type(r_node,r_len));
            assDict.insert(map<string,string>::value_type(r_node,tName));
        }
        
        strStream.str("");
        strStream.clear();
    }
    
    in.close();
}

void GraphRange::hAssNode(string &ass,map<NEdge,int> &r_edge_dict,map<string,Nid> &nid_dict){
    
    ifstream in(pathFile.c_str());
    string pathLine;
    stringstream strStream;
    while(getline(in,pathLine)){
        strStream << pathLine;
        string r_chr,path;
        strStream >> r_chr;
        strStream >> path;
        
        string tName,t_hap,tchr;
        assSplit(r_chr,sep,tName,t_hap,tchr);
        
        string tAss = tName + sep + t_hap;
        vector<char> orient;
        vector<string> nodes;
        bool flag = false;
        if(tAss == ass){
            flag = true;
            char snode[NODESIZE];
            int k = 0;
            for(size_t i = 0; i < path.length(); ++i){
                if(path[i] == '>' || path[i] == '<'){
                    orient.push_back(path[i]);
                    if(i > 0){
                        snode[k] = '\0';
                        nodes.push_back(snode);
                        //memset(snode,0,NODESIZE);
                        k = 0;
                    }
                }else{
                    snode[k] = path[i];
                    ++k;
                }
            }
            snode[k] = '\0';
            nodes.push_back(snode);
            //
            char preOri = '\0';
            int preIn = 0;
            string preNode = "";
            int tIn = 0;
            for(size_t j = 0; j < orient.size(); ++j){
                if(nid_dict.find(nodes[j]) != nid_dict.end()){
                    tIn = 1;
                    hnGroup.push_back(nid_dict[nodes[j]].gNum);
                    if(preIn > 0){
                        char mark = getMark(preOri,orient[j],'>');
                        NEdge sym = {preNode,nodes[j],mark};
                        
                        if(r_edge_dict.find(sym) != r_edge_dict.end()){
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
                                    break;
                            }
                            
                        }else{
                            char rmark = revMark(mark);
                            NEdge sym = {nodes[j],preNode,rmark};
                            
                            if(r_edge_dict.find(sym) != r_edge_dict.end()){
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
                                
                            }
                            
                        }
                    }
                }else{
                    tIn = 0;
                }
                preOri = orient[j];
                preIn = tIn;
                preNode = nodes[j];
            }
        }else{
            if(flag){
                break;
            }
        }
        strStream.str("");
        strStream.clear();
    }
    
    in.close();
    
}

void GraphRange::formatGraph(string &ass,string &sChr,int sStart,int sEnd,int ex,int wStart,int wWidth,int wCut,int wY){
    
    if(indexFlag == 1){
        parseIndex(sChr,sStart,sEnd,ex);
    }else{
        parseEdge();   
    }
    //
    std::vector<string> rangeNode;
    std::set<string> exNode;
    std::map<string,int> info;
    std::map<string,string> assDict;
    int realLen = 0;
    
    parseNode(sChr,sStart,sEnd,ex,rangeNode,exNode,info,assDict,realLen);
    
    int rcount = rangeNode.size();
    float wPerK = (float)wWidth / (realLen + (rcount-1) * wCut);
    float x_wCut = wPerK * wCut;   
    int gNum = 0;
    int iNum = 0;
    
    //string swY = to_string(wY);
    Ndic firNode;
    firNode.insert(Ndic::value_type("id",0));
    firNode.insert(Ndic::value_type("group",0));
    //firNode.insert(Ndic::value_type("fx",wStart));
    //firNode.insert(Ndic::value_type("fy",wY));
    draw_pos.push_back(wStart);
    draw_node.push_back(firNode);

    float node_pre = wStart;
    set<string> nRefNode;
    //vector<string> vRefNode;
    map<string,Nid> nid_dict;
    int s_nid = 0;
    int e_nid = 0;
    string pre_node = "";
    set<string> range_set;
    for(string rnode: rangeNode){
        range_set.insert(rnode);
    }
    
    map<NEdge,int> r_edge_dict;
    for(string &tnode: rangeNode){
        int node_len = info[tnode];
        if(edge.find(tnode) != edge.end()){
            if(edge[tnode].find("o") != edge[tnode].end()){
                vector<string> tNref;
                vector<int> deep;
                for(ENode &o_node : edge[tnode]["o"]){
                    if(exNode.find(o_node.node) != exNode.end()){
                        if(range_set.find(o_node.node) != range_set.end()){
                            NEdge sym = {tnode,o_node.node,o_node.mark};
                            
                            if(r_edge_dict.find(sym) == r_edge_dict.end()){
                                r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                            }
                        }
                    }else{
                        if(nRefNode.find(o_node.node) == nRefNode.end()){
                            tNref.push_back(o_node.node);
                            deep.push_back(0);
                            nRefNode.insert(o_node.node);
                            //vRefNode.push_back(o_node.node);
                        }
                        NEdge sym = {tnode,o_node.node,o_node.mark};
                        r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                    }
                }
                
                if(! tNref.empty()){
                    size_t i = 0;
                    while(i < tNref.size()){
                        if(deep[i] > 10){
                            break;   
                        }
                        if(edge.find(tNref[i]) != edge.end()){
                            if(edge[tNref[i]].find("o") != edge[tNref[i]].end()){
                                for(ENode &to_node : edge[tNref[i]]["o"]){
                                    if(exNode.find(to_node.node) != exNode.end()){
                                        if(range_set.find(to_node.node) != range_set.end()){
                                            NEdge sym = {tNref[i],to_node.node,to_node.mark};
                                            r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                                        }
                                    }else{
                                        if(nRefNode.find(to_node.node) == nRefNode.end()){
                                            tNref.push_back(to_node.node);
                                            deep.push_back(deep[i]+1);
                                            nRefNode.insert(to_node.node);
                                            //vRefNode.push_back(to_node.node);
                                            
                                        }
                                        NEdge sym = {tNref[i],to_node.node,to_node.mark};
                                        r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                                    }
                                }
                            }

                            if(edge[tNref[i]].find("i") != edge[tNref[i]].end()){
                                for(ENode &to_node : edge[tNref[i]]["i"]){
                                    if(exNode.find(to_node.node) != exNode.end()){
                                        if(range_set.find(to_node.node) != range_set.end()){
                                            NEdge sym = {to_node.node,tNref[i],to_node.mark};
                                            r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                                        }
                                    }else{
                                        if(nRefNode.find(to_node.node) == nRefNode.end()){
                                            tNref.push_back(to_node.node);
                                            deep.push_back(deep[i]+1);
                                            nRefNode.insert(to_node.node);
                                            //vRefNode.push_back(to_node.node);
                                        }
                                        NEdge sym = {to_node.node,tNref[i],to_node.mark};
                                        r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                                    }
                                }
                            }
                        }
                        
                        i += 1;
                    }
                }
            }
            //
            if(edge[tnode].find("i") != edge[tnode].end()){
                vector<string> tNref;
                vector<int> deep;
                for(ENode &o_node : edge[tnode]["i"]){
                    if(exNode.find(o_node.node) == exNode.end()){
                        if(nRefNode.find(o_node.node) == nRefNode.end()){
                            tNref.push_back(o_node.node);
                            deep.push_back(0);
                            nRefNode.insert(o_node.node);
                            //vRefNode.push_back(o_node.node);
                        }
                        NEdge sym = {o_node.node,tnode,o_node.mark};
                        r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                    }
                }
                
                if(! tNref.empty()){
                    size_t i = 0;
                    while(i < tNref.size()){
                        if(deep[i] > 10){
                            break;   
                        }
                        if(edge.find(tNref[i]) != edge.end()){
                            if(edge[tNref[i]].find("o") != edge[tNref[i]].end()){
                                for(ENode &to_node : edge[tNref[i]]["o"]){
                                    if(exNode.find(to_node.node) != exNode.end()){
                                        if(range_set.find(to_node.node) != range_set.end()){
                                            NEdge sym = {tNref[i],to_node.node,to_node.mark};
                                            r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                                        }
                                    }else{
                                        if(nRefNode.find(to_node.node) == nRefNode.end()){
                                            tNref.push_back(to_node.node);
                                            deep.push_back(deep[i]+1);
                                            nRefNode.insert(to_node.node);
                                            //vRefNode.push_back(to_node.node);
                                        }
                                        NEdge sym = {tNref[i],to_node.node,to_node.mark};
                                        r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                                    }
                                }
                            }
                            
                            if(edge[tNref[i]].find("i") != edge[tNref[i]].end()){
                                for(ENode &to_node : edge[tNref[i]]["i"]){
                                    if(exNode.find(to_node.node) != exNode.end()){
                                        if(range_set.find(to_node.node) != range_set.end()){
                                            NEdge sym = {to_node.node,tNref[i],to_node.mark};
                                            r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                                        }
                                    }else{
                                        if(nRefNode.find(to_node.node) == nRefNode.end()){
                                            tNref.push_back(to_node.node);
                                            deep.push_back(deep[i]+1);
                                            nRefNode.insert(to_node.node);
                                            //vRefNode.push_back(to_node.node);
                                        }
                                        NEdge sym = {to_node.node,tNref[i],to_node.mark};
                                        r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                                    }
                                }
                            }
                        }
                        
                        i += 1;
                    }
                }
            }
        }
        //
        if(gNum > 0){
            node_pre += x_wCut;
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
            }
        }
        
        float trans = node_len * wPerK;
        node_pre += trans;
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
        nnames.push_back(tnode);
        genome.push_back(assDict[tnode]);
        e_nid = iNum;
        Nid tnid = {s_nid,e_nid,gNum};
        nid_dict.insert(map<string,Nid>::value_type(tnode,tnid));
        
        gNum += 1;
        pre_node = tnode;
        
    }
    
    for(string nk : nRefNode){
        int node_len = info[nk];
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
        
        nnames.push_back(nk);
        genome.push_back(assDict[nk]);
        
        Nid tnid = {s_nid,e_nid,gNum};
        nid_dict.insert(map<string,Nid>::value_type(nk,tnid));
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
    if(ass != "0"){
        hAssNode(ass,r_edge_dict,nid_dict);
    }
}

void GraphRange::edgeRange(string &sChr,int sStart,int sEnd,int ex,std::map<NEdge,int> &r_edge_dict){
    vector<string> rangeNode;
    set<string> exNode;
    map<string,int> info;
    map<string,string> assDict;
    int realLen = 0;
    
    parseNode(sChr,sStart,sEnd,ex,rangeNode,exNode,info,assDict,realLen);
    set<string> nRefNode;
    set<string> range_set;
    for(string rnode: rangeNode){
        range_set.insert(rnode);
    }
    
    for(string &tnode: rangeNode){

        if(edge.find(tnode) != edge.end()){
            if(edge[tnode].find("o") != edge[tnode].end()){
                vector<string> tNref;
                vector<int> deep;
                for(ENode &o_node : edge[tnode]["o"]){
                    if(exNode.find(o_node.node) != exNode.end()){
                        if(range_set.find(o_node.node) != range_set.end()){
                            NEdge sym = {tnode,o_node.node,o_node.mark};
                            
                            if(r_edge_dict.find(sym) == r_edge_dict.end()){
                                r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                            }
                        }
                    }else{
                        if(nRefNode.find(o_node.node) == nRefNode.end()){
                            tNref.push_back(o_node.node);
                            deep.push_back(0);
                            nRefNode.insert(o_node.node);
                        }
                        NEdge sym = {tnode,o_node.node,o_node.mark};
                        r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                    }
                }
                
                if(! tNref.empty()){
                    size_t i = 0;
                    while(i < tNref.size()){
                        if(deep[i] > 10){
                            break;   
                        }
                        if(edge.find(tNref[i]) != edge.end()){
                            if(edge[tNref[i]].find("o") != edge[tNref[i]].end()){
                                for(ENode &to_node : edge[tNref[i]]["o"]){
                                    if(exNode.find(to_node.node) != exNode.end()){
                                        if(range_set.find(to_node.node) != range_set.end()){
                                            NEdge sym = {tNref[i],to_node.node,to_node.mark};
                                            r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                                        }
                                    }else{
                                        if(nRefNode.find(to_node.node) == nRefNode.end()){
                                            tNref.push_back(to_node.node);
                                            deep.push_back(deep[i]+1);
                                            nRefNode.insert(to_node.node);
                                        }
                                        NEdge sym = {tNref[i],to_node.node,to_node.mark};
                                        r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                                    }
                                }
                            }
                            
                            if(edge[tNref[i]].find("i") != edge[tNref[i]].end()){
                                for(ENode &to_node : edge[tNref[i]]["i"]){
                                    if(exNode.find(to_node.node) != exNode.end()){
                                        if(range_set.find(to_node.node) != range_set.end()){
                                            NEdge sym = {to_node.node,tNref[i],to_node.mark};
                                            r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                                        }
                                    }else{
                                        if(nRefNode.find(to_node.node) == nRefNode.end()){
                                            tNref.push_back(to_node.node);
                                            deep.push_back(deep[i]+1);
                                            nRefNode.insert(to_node.node);
                                        }
                                        NEdge sym = {to_node.node,tNref[i],to_node.mark};
                                        r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                                    }
                                }
                            }
                        }
                        
                        i += 1;
                    }
                }
            }
            //
            if(edge[tnode].find("i") != edge[tnode].end()){
                vector<string> tNref;
                vector<int> deep;
                for(ENode &o_node : edge[tnode]["i"]){
                    if(exNode.find(o_node.node) == exNode.end()){
                        if(nRefNode.find(o_node.node) == nRefNode.end()){
                            tNref.push_back(o_node.node);
                            deep.push_back(0);
                            nRefNode.insert(o_node.node);
                        }
                        NEdge sym = {o_node.node,tnode,o_node.mark};
                        r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                    }
                }
                
                if(! tNref.empty()){
                    size_t i = 0;
                    while(i < tNref.size()){
                        if(deep[i] > 10){
                            break;   
                        }
                        if(edge.find(tNref[i]) != edge.end()){
                            if(edge[tNref[i]].find("o") != edge[tNref[i]].end()){
                                for(ENode &to_node : edge[tNref[i]]["o"]){
                                    if(exNode.find(to_node.node) != exNode.end()){
                                        if(range_set.find(to_node.node) != range_set.end()){
                                            NEdge sym = {tNref[i],to_node.node,to_node.mark};
                                            r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                                        }
                                    }else{
                                        if(nRefNode.find(to_node.node) == nRefNode.end()){
                                            tNref.push_back(to_node.node);
                                            deep.push_back(deep[i]+1);
                                            nRefNode.insert(to_node.node);
                                        }
                                        NEdge sym = {tNref[i],to_node.node,to_node.mark};
                                        r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                                    }
                                }
                            }
                            
                            if(edge[tNref[i]].find("i") != edge[tNref[i]].end()){
                                for(ENode &to_node : edge[tNref[i]]["i"]){
                                    if(exNode.find(to_node.node) != exNode.end()){
                                        if(range_set.find(to_node.node) != range_set.end()){
                                            NEdge sym = {to_node.node,tNref[i],to_node.mark};
                                            r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                                        }
                                    }else{
                                        if(nRefNode.find(to_node.node) == nRefNode.end()){
                                            tNref.push_back(to_node.node);
                                            deep.push_back(deep[i]+1);
                                            nRefNode.insert(to_node.node);
                                        }
                                        NEdge sym = {to_node.node,tNref[i],to_node.mark};
                                        r_edge_dict.insert(map<NEdge,int>::value_type(sym,2));
                                    }
                                }
                            }
                        }
                        
                        i += 1;
                    }
                }
            }
        }
    }
}

void GraphRange::edgeWrite(int rangeSize,int ex){
    ifstream in(chrFile.c_str());
    ofstream bfh(bEdgeFile.c_str());
    ofstream xfh(eIndexFile.c_str());
    stringstream strStream;
    string chrLine;
    parseEdge();
    int rSite = 0;
    while(getline(in,chrLine)){
        strStream << chrLine;
        string chrName;
        int chrStart,chrLen;
        strStream >> chrName;
        strStream >> chrStart;
        strStream >> chrLen;
        
        int tEnd = rangeSize;
        int rStart = chrStart;
        int rEnd = rangeSize;
        
        int esize = NODESIZE * 2 + 1;
        while(tEnd < chrLen){
            std::map<NEdge,int> r_edge_dict;
            edgeRange(chrName,rStart,rEnd,ex,r_edge_dict);
            for(auto &tedge : r_edge_dict){
                CEdge xedge;
                strcpy(xedge.node1,tedge.first.node1.c_str());
                strcpy(xedge.node2,tedge.first.node2.c_str());
                xedge.mark = tedge.first.mark;
                bfh.write((char *)&xedge,esize);
            }
            xfh << chrName << "\t" << rStart << "\t" << rEnd << "\t" << rSite << "\t" << r_edge_dict.size() << endl;
            
            tEnd += rangeSize;
            rStart = rEnd + 1;
            rEnd = tEnd;
            rSite += esize * r_edge_dict.size();
        }
        
        rEnd = chrLen;
        std::map<NEdge,int> r_edge_dict;
        edgeRange(chrName,rStart,rEnd,ex,r_edge_dict);
        for(auto &tedge : r_edge_dict){
            CEdge xedge;
            strcpy(xedge.node1,tedge.first.node1.c_str());
            strcpy(xedge.node2,tedge.first.node2.c_str());
            xedge.mark = tedge.first.mark;
            bfh.write((char *)&xedge,esize);
        }
        xfh << chrName << "\t" << rStart << "\t" << rEnd << "\t" << rSite << "\t" << r_edge_dict.size() << endl;
        
        rSite += esize * r_edge_dict.size();
        strStream.str("");
        strStream.clear();
    }
    
    in.close();
    bfh.close();
    xfh.close(); 
}

PYBIND11_MODULE(minipg,m){
    py::class_<GraphRange>(m,"GraphRange")
        .def(py::init<string &,string &,string &,string &>())
        .def(py::init<string &,string &,string &,string &,string &,string &,string &>())
        .def("formatGraph",&GraphRange::formatGraph)
        .def("edgeWrite",&GraphRange::edgeWrite)
        .def_readwrite("draw_node",&GraphRange::draw_node)
        .def_readwrite("draw_pos",&GraphRange::draw_pos)
        .def_readwrite("draw_edge",&GraphRange::draw_edge)
        .def_readwrite("dnode_len",&GraphRange::dnode_len)
        .def_readwrite("genome",&GraphRange::genome)
        .def_readwrite("nnames",&GraphRange::nnames)
        .def_readwrite("hnGroup",&GraphRange::hnGroup)
        .def_readwrite("hLinks",&GraphRange::hLinks);
}



