
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>

#include "gfa2v.h"

using namespace std;
// For 'Jump line' if the distance is unavailable (*) the value will be set to 100. 
int readGFA(string &tmpFolder,ifstream &in,string &refStr,string &sep,unordered_map<int,int> &mNodeLen,map<NEdge,int> &edgeMap,map<NEdge,int> &jumpMap,ofstream &afh,ofstream &acfh){
    string gfaLine;
    stringstream strStream;
    //int nodeCount = 0;
    int maxNode = 0;
    int node,node1,node2;
    string tag,nodeSeq,align,fullName,path;
    char sign1,sign2;
    map<string,int> gMap;
    vector<string> assVec;
    
    int tlen;
    char mark;
    NEdge tedge;
    vector<ofstream> assfh;
    int nf = 1;
    set<string> allComChr;
    
    string rpFile = tmpFolder + "/0.ass";
    ofstream rpfh(rpFile.c_str());
    if(! rpfh){
        cerr<<"Error: file open failed. "<<rpFile<<endl;
        exit(1);
    }
    
    string fixRef = "REF" + sep + "HAP";
    string assName = "",hapName = "",chrName = "";
    string seqStart,seqEnd;
    int numStart;
    string tStr;
    while(getline(in,gfaLine)){
        strStream << gfaLine;
        switch(gfaLine[0]){
            case 'S':
                strStream >> tag;
                strStream >> node;
                strStream >> nodeSeq;
                mNodeLen.emplace(node,nodeSeq.length());
                //++nodeCount;
                if(node > maxNode){
                    maxNode = node;    
                }
                break;
            case 'L':
                strStream >> tag;
                strStream >> node1;
                strStream >> sign1;
                strStream >> node2;
                strStream >> sign2;
                strStream >> align;
                tlen = parseAlign(align);
                mark = getMark(sign1,sign2);
                tedge = {node1,node2,mark};
                edgeMap.emplace(tedge,tlen);
                break;
            case 'J':
                strStream >> tag;
                strStream >> node1;
                strStream >> sign1;
                strStream >> node2;
                strStream >> sign2;
                strStream >> align;
                tlen = 0;
                if(align == "*"){
                    tlen = 100;
                }else{
                    tlen = atoi(align.c_str());
                }
                mark = getMark(sign1,sign2);
                tedge = {node1,node2,mark};
                jumpMap.emplace(tedge,tlen);
                break;
            case 'W':
                assName = "";
                hapName = "";
                chrName = "";
                strStream >> tag;
                strStream >> assName;
                strStream >> hapName;
                strStream >> chrName;
                strStream >> seqStart;
                strStream >> seqEnd;
                if(seqStart == "*"){
                    numStart = 0;
                }else{
                    numStart = atoi(seqStart.c_str());
                }
                strStream >> path;
                tStr = assName + sep + hapName;
                fullName = tStr + sep + chrName;
                //
                if(allComChr.find(fullName) == allComChr.end()){
                    allComChr.insert(fullName);
                }
                
                if(tStr != refStr){
                    if(gMap.find(tStr) == gMap.end()){
                        gMap.emplace(tStr,nf);
                        string tfile = tmpFolder + "/" + to_string(nf) + ".ass";
                        assfh.push_back(ofstream(tfile.c_str()));
                        assVec.push_back(tStr);
                        
                        ++nf;
                    }
                    assfh[gMap[tStr] - 1]<<"W\t"<<fullName<<"\t"<<numStart<<"\t"<<path<<endl;
                    
                }else{
                    if(tStr == fixRef){
                        fullName = tStr + sep + fullName;
                    }
                    rpfh<<"W\t"<<fullName<<"\t"<<numStart<<path<<endl;
                }
                //
                break;
            case 'P':
                strStream >> tag;
                strStream >> fullName;
                strStream >> path;
                
                assName = "";
                hapName = "";
                chrName = "";
                assSplit(fullName,sep,assName,hapName,chrName);
                tStr = assName + sep + hapName;
                
                if(allComChr.find(fullName) == allComChr.end()){
                    allComChr.insert(fullName);
                }
                
                if(tStr != refStr){
                    if(gMap.find(tStr) == gMap.end()){
                        gMap.emplace(tStr,nf);
                        string tfile = tmpFolder + "/" + to_string(nf) + ".ass";
                        assfh.push_back(ofstream(tfile.c_str()));
                        assVec.push_back(tStr);
                        
                        ++nf;
                    }
                    assfh[gMap[tStr] - 1]<<fullName<<"\t"<<path<<endl;
                    
                }else{
                    if(tStr == fixRef){
                        fullName = tStr + sep + fullName;
                    }
                    rpfh<<fullName<<"\t"<<path<<endl;
                }
        }
        strStream.clear();
        strStream.str("");
    }
    //
    rpfh.close();
    for(auto &fh : assfh){
        fh.close();
    }
    //
    for(auto &txChr : allComChr){
        acfh<<txChr<<endl;
    }
    
    //
    afh<<refStr<<endl;

    for(auto &tAss : assVec){
        afh<<tAss<<endl;
    }
    //return nodeCount;
    return maxNode;
}

int psRchrWalk(string &refPath,string &fullName,int refStart,unordered_map<int,int> &mNodeLen,int &neoID,unordered_set<int> &refNodeSet,
               unordered_set<int> &flipSet,unordered_map<int,int> &rCovMap,ofstream &nfh,ofstream &efh,ofstream &pfh){
    int ndVal = 0;
    int pos = refStart, prePos = 0;
    bool preNeo = false;
    pfh<<fullName<<"\t";
    char preOri = '\0';
    int preNode = 0,pNeoNode = 0;
    for(char x : refPath){
        if(x >= '0' && x <= '9'){
            ndVal = ndVal * 10 + (x - '0');
        }else{
            //
            if(x == '>' || x == '<'){       
                if(ndVal > 0){
                    prePos = pos + 1;
                    pos += mNodeLen[ndVal];
                    if(refNodeSet.find(ndVal) == refNodeSet.end()){
                        refNodeSet.insert(ndVal);
                        nfh<<ndVal<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<mNodeLen[ndVal]<<"\t0"<<endl;
                        if(preOri == '<'){
                            flipSet.insert(ndVal);
                        }
                        pfh<<">"<<ndVal;
                        rCovMap.emplace(ndVal,1);
                        //
                        if(preNeo){
                            efh<<pNeoNode<<"\t"<<ndVal<<"\t"<<"+"<<"\t"<<"+"<<endl;
                        }
                        //
                        preNeo = false;
                    }else{
                        rCovMap[ndVal] += 1;
                        ++neoID;
                        cout<<"Add reference node (L): "<<neoID<<endl;
                        rCovMap.emplace(neoID,1);
                        nfh<<neoID<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<mNodeLen[ndVal]<<"\t0"<<endl;
                        pfh<<">"<<neoID;
                        //
                        if(preNeo){
                            efh<<pNeoNode<<"\t"<<neoID<<"\t"<<"+"<<"\t"<<"+"<<endl;
                        }else{
                            efh<<preNode<<"\t"<<neoID<<"\t"<<"+"<<"\t"<<"+"<<endl;
                        }
                        preNeo = true;
                        pNeoNode = neoID;
                    }
                    //
                    preOri = x;
                    preNode = ndVal;
                    ndVal = 0;
                }
            }
        }
    }
    //
    prePos = pos + 1;
    pos += mNodeLen[ndVal];
    if(refNodeSet.find(ndVal) == refNodeSet.end()){
        refNodeSet.insert(ndVal);
        nfh<<ndVal<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<mNodeLen[ndVal]<<"\t0"<<endl;
        if(preOri == '<'){
            flipSet.insert(ndVal);
        }
        pfh<<">"<<ndVal;
        rCovMap.emplace(ndVal,1);
        //
        if(preNeo){
            efh<<pNeoNode<<"\t"<<ndVal<<"\t"<<"+"<<"\t"<<"+"<<endl;
        }
        //
        preNeo = false;
    }else{
        rCovMap[ndVal] += 1;
        ++neoID;
        cout<<"Add reference node (L): "<<neoID<<endl;
        rCovMap.emplace(neoID,1);
        nfh<<neoID<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<mNodeLen[ndVal]<<"\t0"<<endl;
        pfh<<">"<<neoID;
        //
        if(preNeo){
            efh<<pNeoNode<<"\t"<<neoID<<"\t"<<"+"<<"\t"<<"+"<<endl;
        }else{
            efh<<preNode<<"\t"<<neoID<<"\t"<<"+"<<"\t"<<"+"<<endl;
        }
    }
    //
    pfh<<endl;
    return pos;
}

//--------

void psNWalk(string &nrPath,string &ass,int walkStart,unordered_set<int> &flipSet,unordered_map<int,int> &tmap,unordered_map<int,int> &mNodeLen,unordered_set<int> &noutNode,ofstream &pfh,ofstream &nfh){
    int ndVal = 0;
    int pos = walkStart,prePos = 0;
    pfh<<ass<<"\t";
    char tag = '\0';
    char preOri = '\0';
    //
    unordered_map<int,int>::iterator ix;
    for(char x : nrPath){
        if(x >= '0' && x <= '9'){
            ndVal = ndVal * 10 + (x - '0');
        }else{
            if(x == '>' || x == '<'){
                if(ndVal > 0){
                    if(flipSet.find(ndVal) != flipSet.end()){
                        if(preOri == '>'){
                            tag = '<';
                        }else{
                            tag = '>';
                        }
                    }else{
                        if(preOri == '>'){
                            tag = '>';
                        }else{
                            tag = '<';
                        }
                    }
                    pfh<<tag<<ndVal;
                    ix = tmap.find(ndVal);
                    prePos = pos + 1;
                    pos += mNodeLen[ndVal];
                    if(ix != tmap.end()){
                        ix->second += 1;
                    }else{
                        tmap.emplace(ndVal,1);
                        if(noutNode.find(ndVal) != noutNode.end()){
                            nfh<<ndVal<<"\t"<<ass<<"\t"<<prePos<<"\t"<<pos<<"\t"<<mNodeLen[ndVal]<<"\t"<<1<<endl;
                            noutNode.erase(ndVal);
                        }
                    }
                    //
                    ndVal = 0;
                }
                //
                preOri = x;
            }
        }
    }
    //
    if(flipSet.find(ndVal) != flipSet.end()){
        if(preOri == '>'){
            tag = '<';
        }else{
            tag = '>';
        }
    }else{
        if(preOri == '>'){
            tag = '>';
        }else{
            tag = '<';
        }
    }
    pfh<<tag<<ndVal;
    ix = tmap.find(ndVal);
    prePos = pos + 1;
    pos += mNodeLen[ndVal];
    if(ix != tmap.end()){
        ix->second += 1;
        
    }else{
        tmap.emplace(ndVal,1);
        if(noutNode.find(ndVal) != noutNode.end()){
            nfh<<ndVal<<"\t"<<ass<<"\t"<<prePos<<"\t"<<pos<<"\t"<<mNodeLen[ndVal]<<"\t"<<1<<endl;
            noutNode.erase(ndVal);
        }
    }
    //
    pfh<<endl;
}

//
int psRchrPath(string &refPath,string &fullName,unordered_map<int,int> &mNodeLen,map<NEdge,int> &edgeMap,map<NEdge,int> &jumpMap,int &neoID,unordered_set<int> &refNodeSet,
               unordered_set<int> &flipSet,set<NEdge> &rEset,set<NEdge> &rJset,map<NEdge,Jnode> &jNeoMap,unordered_map<int,int> &rCovMap,ofstream &nfh,ofstream &efh,ofstream &pfh){
    int ndVal = 0, node1 = 0, node2 = 0;
    int pNeoNode = 0;
    char sign1 = '\0', sign2 = '\0';
    bool linkType = false;
    bool preNeo = false, tNeo = false;
    map<NEdge,int>::iterator it;
    map<NEdge,Jnode>::iterator ij;
    int pos = 0, prePos = 0, moveLen = 0;
    
    pfh<<fullName<<"\t";
    map<int,int>::iterator ix;
    for(char x : refPath){
        if(x >= '0' && x <= '9'){
            ndVal = ndVal * 10 + (x - '0');
        }else{
            //
            if(x == '+' || x == '-'){
                
                if(node1 == 0){
                    node1 = ndVal;
                    sign1 = x;
                    pos = mNodeLen[node1];
                    //
                    if(refNodeSet.find(node1) == refNodeSet.end()){
                        refNodeSet.insert(node1);
                        nfh<<node1<<"\t"<<fullName<<"\t1\t"<<pos<<"\t"<<pos<<"\t0"<<endl;
                        if(x == '-'){
                            flipSet.insert(node1);
                        }
                        pfh<<">"<<node1;
                        //
                        rCovMap.emplace(node1,1); 
                    }else{
                        rCovMap[ndVal] += 1;
                        ++neoID;
                        cout<<"Add reference node (L): "<<neoID<<endl;
                        rCovMap.emplace(neoID,1);   
                        //
                        if(flipSet.find(node1) == flipSet.end()){
                            if(x == '+'){
                                pfh<<">"<<node1;
                            }else{
                                pfh<<"<"<<node1;
                            }
                        }else{
                            if(x == '+'){
                                pfh<<"<"<<node1;
                            }else{
                                pfh<<">"<<node1;
                            }
                        }
                        preNeo = true;
                        pNeoNode = neoID;
                    }
                }else{
                    if(node2 != 0){
                        node1 = node2;
                        sign1 = sign2;
                        
                    }
                    node2 = ndVal;
                    sign2 = x;
                    char mark = getMark(sign1,sign2);
                    //
                    if(refNodeSet.find(node2) != refNodeSet.end()){
                        tNeo = true;
                    }else{
                        refNodeSet.insert(node2);
                        tNeo = false;
                        //
                        if(x == '-'){
                            flipSet.insert(node2);
                        }
                    }
                    //
                    NEdge tedge = {node1,node2,mark};
                    if(linkType){
                        if(tNeo){
                            rCovMap[node2] += 1;
                            ++neoID;
                            cout<<"Add reference node (L): "<<neoID<<endl;
                            rCovMap.emplace(neoID,1);
                        }else{
                            rCovMap.emplace(node2,1);
                        }
                        
                        it = edgeMap.find(tedge);
                        if(it != edgeMap.end()){
                            moveLen = mNodeLen[node2] - it->second;
                            prePos = pos + 1;
                            pos += moveLen;
                            if(preNeo){
                                if(tNeo){
                                    nfh<<neoID<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                                    efh<<pNeoNode<<"\t"<<neoID<<"\t"<<"+"<<"\t"<<"+"<<endl;
                                }else{
                                    nfh<<node2<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                                    efh<<pNeoNode<<"\t"<<node2<<"\t"<<"+"<<"\t"<<"+"<<endl;
                                }
                            }else{
                                if(tNeo){
                                    nfh<<neoID<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                                    efh<<node1<<"\t"<<neoID<<"\t"<<"+"<<"\t"<<"+"<<endl;
                                }else{
                                    nfh<<node2<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                                    efh<<node1<<"\t"<<node2<<"\t"<<"+"<<"\t"<<"+"<<endl;
                                    //
                                    rEset.insert(tedge);
                                }
                            }
                        }else{
                            char rmark = revMark(mark);
                            NEdge rtedge = {node2,node1,rmark};
                            it = edgeMap.find(rtedge);
                            if(it != edgeMap.end()){
                                moveLen = mNodeLen[node2] - it->second;
                                prePos = pos + 1;
                                pos += moveLen;
                                if(preNeo){
                                    if(tNeo){
                                        nfh<<neoID<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                                        efh<<pNeoNode<<"\t"<<neoID<<"\t"<<"+"<<"\t"<<"+"<<endl;
                                    }else{
                                        nfh<<node2<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                                        efh<<pNeoNode<<"\t"<<node2<<"\t"<<"+"<<"\t"<<"+"<<endl;
                                    }
                                }else{
                                    if(tNeo){
                                        nfh<<neoID<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                                        efh<<node1<<"\t"<<neoID<<"\t"<<"+"<<"\t"<<"+"<<endl;
                                    }else{
                                        nfh<<node2<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                                        efh<<node1<<"\t"<<node2<<"\t"<<"+"<<"\t"<<"+"<<endl;
                                        //
                                        rEset.insert(tedge);
                                    }
                                }
                            }else{
                                cerr<<"Warning: edge (L) ["<<node1<<"\t"<<node2<<"\t"<<sign1<<"\t"<<sign2<<"] in reference path can't be found in Link lines"<<endl;
                                moveLen = mNodeLen[node2];
                                prePos = pos + 1;
                                pos += moveLen;
                                if(preNeo){
                                    if(tNeo){
                                        nfh<<neoID<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                                        efh<<pNeoNode<<"\t"<<neoID<<"\t"<<"+"<<"\t"<<"+"<<endl;
                                    }else{
                                        nfh<<node2<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                                        efh<<pNeoNode<<"\t"<<node2<<"\t"<<"+"<<"\t"<<"+"<<endl;
                                    }
                                }else{
                                    if(tNeo){
                                        nfh<<neoID<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                                        efh<<node1<<"\t"<<neoID<<"\t"<<"+"<<"\t"<<"+"<<endl;
                                    }else{
                                        nfh<<node2<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                                        efh<<node1<<"\t"<<node2<<"\t"<<"+"<<"\t"<<"+"<<endl;
                                    }
                                }
                            }
                        }
                    }else{
                        it = jumpMap.find(tedge);
                        if(it != jumpMap.end()){
                            moveLen = it->second;
                            prePos = pos + 1;
                            pos += moveLen;
                            ++neoID;
                            cout<<"Add reference node (J): "<<neoID<<endl;
                            rCovMap.emplace(neoID,1);
                            nfh<<neoID<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                            //
                            ij = jNeoMap.find(tedge);
                            if(ij == jNeoMap.end()){
                                Jnode tj = {neoID,'+'};
                                jNeoMap.emplace(tedge,tj);
                                pfh<<">"<<neoID;
                                //
                                rJset.insert(tedge);
                                //
                                char tsign1 = '\0',tsign2 = '\0';
                                if(flipSet.find(node1) != flipSet.end()){
                                    if(sign1 == '+'){
                                        tsign1 = '-';
                                    }else{
                                        tsign1 = '+';
                                    }
                                }else{
                                    tsign1 = sign1;
                                }
                                
                                if(flipSet.find(node2) != flipSet.end()){
                                    if(sign2 == '+'){
                                        tsign2 = '-';
                                    }else{
                                        tsign2 = '+';
                                    }
                                }else{
                                    tsign2 = sign2;
                                }
                                efh<<node1<<"\t"<<neoID<<"\t"<<tsign1<<"\t"<<"+"<<endl;
                                efh<<neoID<<"\t"<<node2<<"\t"<<"+"<<"\t"<<tsign2<<endl;
                            }else{
                                if((ij->second).ori == '+'){
                                    pfh<<">"<<(ij->second).neoID;
                                }else{
                                    pfh<<"<"<<(ij->second).neoID;   
                                }
                                //
                                rCovMap[(ij->second).neoID] += 1;
                            }
                            //
                            if(preNeo){
                                efh<<pNeoNode<<"\t"<<neoID<<"\t"<<"+"<<"\t"<<"+"<<endl;
                            }
                            moveLen = mNodeLen[node2];
                            prePos = pos + 1;
                            pos += moveLen;
                            pNeoNode = neoID;
                            if(tNeo){
                                rCovMap[node2] += 1;
                                ++neoID;
                                cout<<"Add reference node (J): "<<neoID<<endl;
                                rCovMap.emplace(neoID,1);
                                //
                                nfh<<neoID<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                                efh<<pNeoNode<<"\t"<<neoID<<"\t"<<"+"<<"\t"<<"+"<<endl;
                            }else{
                                //
                                rCovMap.emplace(node2,1);
                                //
                                nfh<<node2<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                            }
                        }else{
                            char rmark = revMark(mark);
                            NEdge rtedge = {node2,node1,rmark};
                            it = jumpMap.find(rtedge);
                            if(it != jumpMap.end()){
                                moveLen = it->second;
                                prePos = pos + 1;
                                pos += moveLen;
                                ++neoID;
                                //
                                cout<<"Add reference node (J): "<<neoID<<endl;
                                rCovMap.emplace(neoID,1);
                                nfh<<neoID<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                                //
                                ij = jNeoMap.find(rtedge);
                                if(ij == jNeoMap.end()){
                                    Jnode tj = {neoID,'-'};
                                    jNeoMap.emplace(tedge,tj);
                                    pfh<<">"<<neoID;
                                    //
                                    rJset.insert(rtedge);
                                    //
                                    //
                                    char tsign1 = '\0',tsign2 = '\0';
                                    if(flipSet.find(node1) != flipSet.end()){
                                        if(sign1 == '+'){
                                            tsign1 = '-';
                                        }else{
                                            tsign1 = '+';
                                        }
                                    }else{
                                        tsign1 = sign1;
                                    }
                                    
                                    if(flipSet.find(node2) != flipSet.end()){
                                        if(sign2 == '+'){
                                            tsign2 = '-';
                                        }else{
                                            tsign2 = '+';
                                        }
                                    }else{
                                        tsign2 = sign2;
                                    }
                                    efh<<node1<<"\t"<<neoID<<"\t"<<tsign1<<"\t"<<"+"<<endl;
                                    efh<<neoID<<"\t"<<node2<<"\t"<<"+"<<"\t"<<tsign2<<endl;
                                }else{
                                    if((ij->second).ori == '+'){
                                        pfh<<"<"<<(ij->second).neoID;
                                    }else{
                                        pfh<<">"<<(ij->second).neoID;   
                                    }
                                    //
                                    rCovMap[(ij->second).neoID] += 1;
                                }
                                //
                                if(preNeo){
                                    efh<<pNeoNode<<"\t"<<neoID<<"\t"<<"+"<<"\t"<<"+"<<endl;
                                }
                                //
                                moveLen = mNodeLen[node2];
                                prePos = pos + 1;
                                pos += moveLen;
                                pNeoNode = neoID;
                                if(tNeo){
                                    rCovMap[node2] += 1;
                                    ++neoID;
                                    cout<<"Add reference node (J): "<<neoID<<endl;
                                    rCovMap.emplace(neoID,1);
                                    //
                                    nfh<<neoID<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                                    efh<<pNeoNode<<"\t"<<neoID<<"\t"<<"+"<<"\t"<<"+"<<endl;
                                }else{
                                    rCovMap.emplace(node2,1);
                                    //
                                    nfh<<node2<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                                }
                            }else{
                                cerr<<"Warning: edge (J) ["<<node1<<"\t"<<node2<<"\t"<<sign1<<"\t"<<sign2 <<"] in reference path can't be found in Jump lines"<<endl;
                                moveLen = mNodeLen[node2];
                                prePos = pos + 1;
                                pos += moveLen;
                                if(preNeo){
                                    if(tNeo){
                                        rCovMap[node2] += 1;
                                        ++neoID;
                                        cout<<"Add reference node (J): "<<neoID<<endl;
                                        rCovMap.emplace(neoID,1);
                                        //
                                        nfh<<neoID<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                                        efh<<pNeoNode<<"\t"<<neoID<<"\t"<<"+"<<"\t"<<"+"<<endl;
                                    }else{
                                        rCovMap.emplace(node2,1);
                                        //
                                        nfh<<node2<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                                        efh<<pNeoNode<<"\t"<<node2<<"\t"<<"+"<<"\t"<<"+"<<endl;
                                    }
                                }else{
                                    if(tNeo){
                                        rCovMap[node2] += 1;
                                        ++neoID;
                                        cout<<"Add reference node (J): "<<neoID<<endl;
                                        rCovMap.emplace(neoID,1);
                                        //
                                        nfh<<neoID<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                                        efh<<node1<<"\t"<<neoID<<"\t"<<"+"<<"\t"<<"+"<<endl;
                                    }else{
                                        rCovMap.emplace(node2,1);
                                        //
                                        nfh<<node2<<"\t"<<fullName<<"\t"<<prePos<<"\t"<<pos<<"\t"<<moveLen<<"\t0"<<endl;
                                        efh<<node1<<"\t"<<node2<<"\t"<<"+"<<"\t"<<"+"<<endl;
                                    }
                                }
                            }
                        }
                    }
                    //
                    if(flipSet.find(node2) != flipSet.end()){
                        if(x == '+'){
                            pfh<<"<"<<node2;
                        }else{
                            pfh<<">"<<node2;
                        }
                    }else{
                        if(x == '+'){
                            pfh<<">"<<node2;
                        }else{
                            pfh<<"<"<<node2;
                        }
                    }
                    //
                    pNeoNode = neoID;
                    preNeo = tNeo;
                }
                
            }else if(x == ','){
                linkType = true;
                ndVal = 0;
            }else if(x == ';'){
                linkType = false;
                ndVal = 0;
            }
        }
        
    }
    pfh<<endl;
    if(node2 == 0){
        cerr<<"Warning: the path only has one node"<<endl;
    }
    
    return pos;
}

//
void psNPath(string &nrPath,string &ass,unordered_set<int> &flipSet,map<NEdge,Jnode> &jNeoMap,unordered_map<int,int> &tmap,unordered_map<int,int> &mNodeLen,unordered_set<int> &noutNode,ofstream &pfh,ofstream &nfh){
    int ndVal = 0, node1 = 0, node2 = 0;
    char sign1 = '\0', sign2 = '\0';
    bool linkType = false;
    map<NEdge,Jnode>::iterator ij;
    int pos = 0;
    pfh<<ass<<"\t";
    char tag;
    //
    unordered_map<int,int>::iterator ix;
    for(char x : nrPath){
        if(x >= '0' && x <= '9'){
            ndVal = ndVal * 10 + (x - '0');
        }else{
            if(x == '+' || x == '-'){
                if(node1 == 0){
                    node1 = ndVal;
                    sign1 = x;
                    if(flipSet.find(node1) != flipSet.end()){
                        if(x == '+'){
                            tag = '<';
                        }else{
                            tag = '>';
                        }
                    }else{
                        if(x == '+'){
                            tag = '>';
                        }else{
                            tag = '<';
                        }
                    }
                    pfh<<tag<<node1;
                    //
                    ix = tmap.find(node1);
                    if(ix != tmap.end()){
                        ix->second += 1;
                        pos += mNodeLen[node1];
                    }else{
                        tmap.emplace(node1,1);
                        if(noutNode.find(node1) != noutNode.end()){
                            nfh<<node1<<"\t"<<ass<<"\t"<<1<<"\t"<<mNodeLen[node1]<<"\t"<<mNodeLen[node1]<<"\t"<<1<<endl;
                            pos += mNodeLen[node1];
                            noutNode.erase(node1);
                        }
                    }
                }else{
                    if(node2 == 0){
                        node2 = ndVal;
                    }else{
                        node1 = node2;
                        sign1 = sign2;
                        node2 = ndVal;
                    }
                    sign2 = x;
                    
                    ix = tmap.find(node2);
                    if(ix != tmap.end()){
                        ix->second += 1;
                        pos += mNodeLen[node2];
                    }else{
                        tmap.emplace(node2,1);
                        if(noutNode.find(node2) != noutNode.end()){
                            nfh<<node2<<"\t"<<ass<<"\t"<<pos+1<<"\t"<<pos+mNodeLen[node2]<<"\t"<<mNodeLen[node2]<<"\t"<<1<<endl;
                            pos += mNodeLen[node2];
                            noutNode.erase(node2);
                        }
                    }
                    //
                    if(! linkType){
                        char mark = getMark(sign1,sign2);
                        NEdge tedge = {node1,node2,mark};
                        ij = jNeoMap.find(tedge);
                        if(ij != jNeoMap.end()){
                            if((ij->second).ori == '+'){
                                tag = '>';
                            }else{
                                tag = '<';
                            }
                            pfh<<tag<<(ij->second).neoID;
                            //
                            ix = tmap.find((ij->second).neoID);
                            if(ix != tmap.end()){
                                ix->second += 1;
                            }else{
                                tmap.emplace((ij->second).neoID,1);
                            }
                        }else{
                            char rmark = revMark(mark);
                            NEdge rtedge = {node2,node1,rmark};
                            ij = jNeoMap.find(rtedge);
                            if(ij != jNeoMap.end()){
                                if((ij->second).ori == '+'){
                                    tag = '<';
                                }else{
                                    tag = '>';
                                }
                                pfh<<tag<<(ij->second).neoID;
                                //
                                ix = tmap.find((ij->second).neoID);
                                if(ix != tmap.end()){
                                    ix->second += 1;
                                }else{
                                    tmap.emplace((ij->second).neoID,1);
                                }
                            }else{
                                cerr<<"Warning: edge (J) ["<<node1<<"\t"<<node2<<"\t"<<sign1<<"\t"<<sign2 <<"] in path can't be found in Jump lines"<<endl;
                            }
                        }
                    }
                    if(flipSet.find(node2) != flipSet.end()){
                        if(x == '+'){
                            tag = '<';
                        }else{
                            tag = '>';
                        }
                    }else{
                        if(x == '+'){
                            tag = '>';
                        }else{
                            tag = '<';
                        }
                    }
                    pfh<<tag<<node2;
                }
                //
            }else if(x == ','){
                linkType = true;
                ndVal = 0;
            }else if(x == ';'){
                linkType = false;
                ndVal = 0;
            }
        }
        
    }
    pfh<<endl;
    if(node2 == 0){
        cerr<<"Warning: the path only has one node. "<<nrPath<<endl;
    }
}
 
void psAllPath(char *rfChrFile,string &tmpFolder,string &assFile,string &sepStr,int &neoID,unordered_map<int,int> &mNodeLen,map<NEdge,int> &edgeMap,map<NEdge,int> &jumpMap,ofstream &nfh,ofstream &efh,ofstream &covfh,ofstream &xcovfh,ofstream &cfh,string &pathDir){
    
    string gfaLine,chrLine;
    stringstream strStream;
    string fullName,path;
    //
    set<string> resChr;
    bool flag = false;
    if(rfChrFile != nullptr){
        flag = true;
        ifstream rfh(rfChrFile);
        while(getline(rfh,chrLine)){
            resChr.insert(chrLine);
        }   
        rfh.close();
    }
    //
    unordered_set<int> refNodeSet;
    unordered_set<int> flipSet;
    set<NEdge> rEset;
    set<NEdge> rJset;
    map<NEdge,Jnode> jNeoMap;
    unordered_map<int,int> rCovMap;
    map<string,int> chrSet;
    map<string,int>::iterator cit;
    
    string rpFile = tmpFolder + "/0.ass";
    ifstream rpfh(rpFile.c_str());
    
    string formRfile = pathDir + "/0.path";
    ofstream fpfh(formRfile.c_str());
    if(! fpfh){
        cerr<<"Error: file open failed. "<<formRfile<<endl;
        exit(1);
    }
    int numStart = 0;
    while(getline(rpfh,gfaLine)){
        strStream<<gfaLine;
        strStream >> fullName;
        bool walk = false;
        if(fullName == "W"){
            strStream >> fullName;
            strStream >> numStart;
            walk = true;
        }
        strStream >> path;
        strStream.clear();
        strStream.str("");
        string assName = "",hapName = "",chrName = "";
        assSplit(fullName,sepStr,assName,hapName,chrName);
        
        string tName = chrName;
        cit = chrSet.find(tName);
        if(cit != chrSet.end()){
            cit->second += 1;
            string tchr = "_" + to_string(cit->second);
            tName += tchr;
        }else{
            chrSet.emplace(tName,0);
        }
        int tpos = 0;
        if(walk){
            tpos = psRchrWalk(path,fullName,numStart,mNodeLen,neoID,refNodeSet,flipSet,rCovMap,nfh,efh,fpfh);
        }else{
            tpos = psRchrPath(path,fullName,mNodeLen,edgeMap,jumpMap,neoID,refNodeSet,flipSet,rEset,rJset,jNeoMap,rCovMap,nfh,efh,fpfh);
        }
        //
        if(flag){
            if(resChr.find(chrName) == resChr.end()){
                continue;
            }
        }
        cfh<<tName<<"\t1\t"<<tpos<<endl;
    }
    chrSet.clear();
    rpfh.close();
    fpfh.close();
    remove(rpFile.c_str());
    //
    unordered_set<int> noutNode;
    for(auto &node : mNodeLen){
        if(refNodeSet.find(node.first) == refNodeSet.end()){
            noutNode.insert(node.first);
        }
    }
    refNodeSet.clear();
    //
    for(auto &edge : edgeMap){
        if(rEset.find(edge.first) == rEset.end()){
            char rmark = revMark(edge.first.mark);
            NEdge tedge = {edge.first.node2,edge.first.node1,rmark};
            if(rEset.find(tedge) == rEset.end()){
                char sign1 = '\0',sign2 = '\0';
                bool flip1,flip2;
                if(flipSet.find(edge.first.node1) != flipSet.end()){
                    flip1 = true;
                }else{
                    flip1 = false;
                }
                
                if(flipSet.find(edge.first.node2) != flipSet.end()){
                    flip2 = true;
                }else{
                    flip2 = false;
                }
                markSign(edge.first.mark,flip1,flip2,sign1,sign2);
                efh<<edge.first.node1<<"\t"<<edge.first.node2<<"\t"<<sign1<<"\t"<<sign2<<endl;
            }
        }
    }
    edgeMap.clear();
    rEset.clear();
    //
    for(auto &jedge : jumpMap){
        if(rJset.find(jedge.first) == rJset.end()){
            char rmark = revMark(jedge.first.mark);
            NEdge tedge = {jedge.first.node2,jedge.first.node1,rmark};
            if(rJset.find(tedge) == rJset.end()){
                ++neoID;
                Jnode tnode = {neoID,'+'};
                nfh<<neoID<<"\tJump"<<sepStr<<"H"<<sepStr<<"1"<<"\t1\t"<<jedge.second<<"\t"<<jedge.second<<"\t1"<<endl;
                jNeoMap.emplace(jedge.first,tnode);
                //
                char sign1 = '\0',sign2 = '\0';
                bool flip1,flip2;
                if(flipSet.find(jedge.first.node1) != flipSet.end()){
                    flip1 = true;
                }else{
                    flip1 = false;
                }
                
                if(flipSet.find(jedge.first.node2) != flipSet.end()){
                    flip2 = true;
                }else{
                    flip2 = false;
                }
                markSign(jedge.first.mark,flip1,flip2,sign1,sign2);
                
                efh<<jedge.first.node1<<"\t"<<neoID<<"\t"<<sign1<<"\t+"<<endl;
                efh<<neoID<<"\t"<<jedge.first.node2<<"\t"<<"+\t"<<sign2<<endl;
            }
        }
    }
    jumpMap.clear();
    rJset.clear();
    //-----------------------------------------
    string drcFile = tmpFolder + "/0.cov.dx";
    ofstream drcfh(drcFile.c_str());
    //
    unordered_map<int,int>::iterator it;
    int dxNum;
    if(neoID % 8){
        dxNum = neoID / 8 + 1;
    }else{
        dxNum = neoID / 8;        
    }
    char *dxArr = new char[dxNum];
    memset(dxArr,0,dxNum);
    
    //
    string rcFile = tmpFolder + "/0.cov";
    ofstream rcfh(rcFile.c_str());          
    for(int k=1; k<=neoID; ++k){
        it = rCovMap.find(k);
        if(it != rCovMap.end()){
            rcfh.write((char *)&(it->second),sizeof(int));
            int j = (k - 1) / 8;
            int x = 7 - (k - 1) % 8;
            dxArr[j] |= (0x1 << x);
        }
    }
    drcfh.write(dxArr,dxNum);
    drcfh.close();
    rcfh.close();
    rCovMap.clear();
    //-----------------------------
    
    string assLine;
    int nf = 0;
    ifstream assfh(assFile.c_str());
    while(getline(assfh,assLine)){
        ++nf;
    }
    assfh.close();
    if(nf > 65535){
        cerr<<"Error: too many assemblies, not supported by current version of gfa2view. "<<endl;
    }
    unsigned short int snf = (unsigned short int)nf;
    
    for(int h = 1; h < nf; ++h){
        string tpFile = tmpFolder + "/" + to_string(h) + ".ass";
        ifstream tpfh(tpFile.c_str());
        
        string fpFile = pathDir + "/" + to_string(h) + ".path";
        ofstream fpfh(fpFile.c_str());
        if(! fpfh){
            cerr<<"Error: file open failed. "<<fpFile<<endl;
            exit(1);
        }
        
        string rcFile = tmpFolder + "/" + to_string(h) + ".cov";
        ofstream rcfh(rcFile.c_str());
        if(! rcfh){
            cerr<<"Error: file open failed. "<<rcFile<<endl;
            exit(1);
        }
        //
        while(getline(tpfh,gfaLine)){
            strStream<<gfaLine;
            strStream >> fullName;
            bool walk = false;
            if(fullName == "W"){
                strStream >> fullName;
                strStream >> numStart;
                walk = true;
            }
            strStream >> path;
            if(walk){
                psNWalk(path,fullName,numStart,flipSet,rCovMap,mNodeLen,noutNode,fpfh,nfh);
            }else{
                psNPath(path,fullName,flipSet,jNeoMap,rCovMap,mNodeLen,noutNode,fpfh,nfh);
            }
            strStream.clear();
            strStream.str("");
        }
        tpfh.close();
        fpfh.close();
        
        remove(tpFile.c_str());
        //
        memset(dxArr,0,dxNum);            
        for(int k=1; k<=neoID; ++k){
            it = rCovMap.find(k);
            if(it != rCovMap.end()){
                rcfh.write((char *)&(it->second),sizeof(int));
                int j = (k - 1) / 8;
                int x = 7 - (k - 1) % 8;
                dxArr[j] |= (0x1 << x);
            }
        }
        //
        string drcFile = tmpFolder + "/" + to_string(h) + ".cov.dx";
        ofstream drcfh(drcFile.c_str());
        
        drcfh.write(dxArr,dxNum);
        rcfh.close();
        drcfh.close();
        rCovMap.clear();
    }
    
    for(auto &mNode : noutNode){
        nfh<<mNode<<"\tUn"<<sepStr<<"H"<<sepStr<<"1"<<"\t1\t"<<mNodeLen[mNode]<<"\t"<<mNodeLen[mNode]<<"\t1"<<endl;
    }
    mNodeLen.clear();
    //------------------
    delete []dxArr;
    //------------------
    int tz = 1000;
    int tb = tz * 8;
    int ng = nf;
    int **covArr = new int*[tb];
    for(int x=0; x<tb; ++x){
        covArr[x] = new int[ng];   
    }
    char *tdxArr = new char[tz];
    int ptz = tz;
    
    //ifstream idrcfh(drcFile);
    vector<ifstream> alltmpfh;
    vector<ifstream> alldxfh;
    for(int i=0; i<nf; ++i){
        string trcFile = tmpFolder + "/" + to_string(i) + ".cov";
        alltmpfh.push_back(ifstream(trcFile.c_str()));

        string tdrcFile = tmpFolder + "/" + to_string(i) + ".cov.dx";
        alldxfh.push_back(ifstream(tdrcFile.c_str()));
    }
    
    int numLimit = nf / 2;
    int sintSize = sizeof(unsigned short int);
    long long cOff = 0LL;
    int llSize = sizeof(long long);
    //
    while(ptz < dxNum){
        for(int ia = 0; ia < nf; ++ia){
            memset(tdxArr,0,tz);
            alldxfh[ia].read(tdxArr,tz);
            vector<int> nonzero;
            for(int m=0; m<tb; ++m){
                int n = m / 8;
                int p = 7 - m % 8;
                if(tdxArr[n] & (0x1 << p)){
                    nonzero.push_back(m);                 
                }else{
                    covArr[m][ia] = 0;
                }
            }
            
            int nonzCount = nonzero.size();
            if(nonzCount > 0){
                int *nonZarr = new int[nonzCount];
                alltmpfh[ia].read((char *)nonZarr,sizeof(int)*nonzCount);
                
                for(int i=0; i<nonzCount; ++i){
                    covArr[nonzero[i]][ia] = nonZarr[i];                   
                }
                delete []nonZarr;
            }           
        }
        for(int u=0; u<tb; ++u){
            vector<int> vaIndex;
            vector<int> gvaIndex;
            vector<int> gva;
            
            for(int ib=0; ib<nf; ++ib){
                if(covArr[u][ib] == 1){
                    vaIndex.push_back(ib);
                }else if(covArr[u][ib] > 0){
                    gvaIndex.push_back(ib);
                    gva.push_back(covArr[u][ib]);
                }
            }
            unsigned short int vaNum = vaIndex.size();
            unsigned short int gvaNum = gvaIndex.size();
            int nz = vaIndex.size() + gvaIndex.size();
            if(nz < numLimit){
                for(auto val : vaIndex){
                    covfh.write((char *)&val,sintSize);
                }
                for(auto gv : gvaIndex){
                    covfh.write((char *)&gv,sintSize);
                }
                for(auto g : gva){
                    unsigned short int value = g > 65535 ? 65535 : g;
                    covfh.write((char *)&value,sintSize);
                }
                xcovfh.write((char *)&cOff,llSize);
                xcovfh.write((char *)&vaNum,sintSize);
                xcovfh.write((char *)&gvaNum,sintSize);
                //
                cOff += (vaNum + gvaNum * 2) * sintSize;
            }else{
                for(int ib=0; ib<nf; ++ib){
                    unsigned short int value = covArr[u][ib] > 65535 ? 65535 : covArr[u][ib];
                    covfh.write((char *)&value,sintSize);
                }
                //
                xcovfh.write((char *)&cOff,llSize);
                xcovfh.write((char *)&snf,sintSize);
                xcovfh.write((char *)&snf,sintSize);
                //
                cOff += nf * sintSize;
            }
        }
        //
        ptz += tz;
    }
    //
    int ktz = tz + dxNum - ptz;
    int res = neoID % 8;
    int ktb;
    if(res){
        ktb = (ktz - 1) * 8 + res;
    }else{
        ktb = ktz * 8;
    }
    
    for(int ia=0; ia<nf; ++ia){
        memset(tdxArr,0,tz);
        alldxfh[ia].read(tdxArr,ktz);
        vector<int> nonzero;
        for(int m=0; m<ktb; ++m){
            int n = m / 8;
            int p = 7 - m % 8;
            if(tdxArr[n] & (0x1 << p)){
                nonzero.push_back(m);
            }else{
                covArr[m][ia] = 0;
            }
        }
        int nonzCount = nonzero.size();
        if(nonzCount > 0){
            int *nonZarr = new int[nonzCount];
            alltmpfh[ia].read((char *)nonZarr,sizeof(int)*nonzCount);
            
            for(int i=0; i<nonzCount; ++i){
                covArr[nonzero[i]][ia] = nonZarr[i];
            }
            delete []nonZarr;
        }
    }
    //
    for(int u=0; u<ktb; ++u){
        vector<int> vaIndex;
        vector<int> gvaIndex;
        vector<int> gva;
        
        for(int ib=0; ib<nf; ++ib){
            if(covArr[u][ib] == 1){
                vaIndex.push_back(ib);
            }else if(covArr[u][ib] > 0){
                gvaIndex.push_back(ib);
                gva.push_back(covArr[u][ib]);
            }
        }
        //
        unsigned short int vaNum = vaIndex.size();
        unsigned short int gvaNum = gvaIndex.size();
        int nz = vaIndex.size() + gvaIndex.size();
        if(nz < numLimit){
            for(auto val : vaIndex){
                covfh.write((char *)&val,sintSize);
            }
            for(auto gv : gvaIndex){
                covfh.write((char *)&gv,sintSize);
            }
            for(auto g : gva){
                unsigned short int value = g > 65535 ? 65535 : g;
                covfh.write((char *)&value,sintSize);
            }
            xcovfh.write((char *)&cOff,llSize);
            xcovfh.write((char *)&vaNum,sintSize);
            xcovfh.write((char *)&gvaNum,sintSize);
        }else{
            for(int ib=0; ib<nf; ++ib){
                unsigned short int value = covArr[u][ib] > 65535 ? 65535 : covArr[u][ib];
                covfh.write((char *)&value,sintSize);
            }
            //
            xcovfh.write((char *)&cOff,llSize);
            xcovfh.write((char *)&snf,sintSize);
            xcovfh.write((char *)&snf,sintSize);
        }
        cOff += (vaNum + gvaNum * 2) * sintSize;
        
    }
    
    //----------------------------
    for(int ix=0; ix<tb; ++ix){
        delete []covArr[ix];        
    }
    delete []covArr;
    delete []tdxArr;
    
    int g = 0;
    for(auto &fh : alltmpfh){
        fh.close();
        //
        string trcFile = tmpFolder + "/" + to_string(g) + ".cov";
        remove(trcFile.c_str());
        g++;
    }
    
    g = 0;
    for(auto &dfh : alldxfh){
        dfh.close();
        string tdrcFile = tmpFolder + "/" + to_string(g) + ".cov.dx";
        remove(tdrcFile.c_str());
        g++;
    }
}

void gfa2view(char *rfChrFile,char *gfaFile,char *refName,char *sep,int range,int ex,bool index,int nocross,int nthread,char *outDir){
    
    string sOutDir = outDir;
    string upFolder = sOutDir + "/upload";
    //
    string edgeFile = upFolder + "/edge.info";
    string nodeFile = upFolder + "/node.info";
    string covFile = upFolder + "/cover.bw";
    string xcovFile = upFolder + "/cover.bdx";
    string sepFile = upFolder + "/sep.info";
    string assListFile = upFolder + "/ass.list";
    string chrListFile = upFolder + "/chr.list";
    
    string formatFile = upFolder + "/form.info";
    string comChrFile = upFolder + "/complete.chr.list";
    //
    if(gfaFile != nullptr){
        if(opendir(outDir)){
            cerr<<"Error: output directory already exists"<<endl;
            exit(1);
        }

        if(mkdir(outDir,S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)){
            cerr<<"Error: failed to create output directory, please make sure the parent directory of output exists"<<endl;
            exit(1);
        }
        mkdir(upFolder.c_str(),S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
        //
        string tmpFolder = sOutDir + "/tmp";
        mkdir(tmpFolder.c_str(),S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
        //
        string pathDir = upFolder + "/path";
        mkdir(pathDir.c_str(),S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
        //
        ofstream ofmfh(formatFile.c_str());
        if(! ofmfh){
            cerr<<"Error: file including format information open failed"<<formatFile<<endl; 
            exit(1);        
        }
        ofmfh<<"GFAv1"<<endl;
        ofmfh.close();
        //
        ifstream in(gfaFile);
        if(! in){
            cerr<<"Error: GFA file open failed"<<endl;
            exit(1);
        }
        
        ofstream sfh(sepFile.c_str());
        if(! sfh){
            cerr<<"Error: delimiter file open failed. "<<sepFile<<endl;
            exit(1);
        }
        sfh<<sep<<endl;
        sfh.close();
        
        ofstream afh(assListFile.c_str());
        if(! afh){
            cerr<<"Error: assembly list file open failed. "<<assListFile<<endl;
            exit(1);
        }
        ofstream acfh(comChrFile.c_str());
        if(! acfh){
            cerr<<"Error: file open failed. "<<comChrFile<<endl; 
            exit(1);        
        }
        //
        unordered_map<int,int> mNodeLen;
        map<NEdge,int> edgeMap;
        map<NEdge,int> jumpMap;
        string refStr = refName;
        string sepStr = sep;
        int maxNode = readGFA(tmpFolder,in,refStr,sepStr,mNodeLen,edgeMap,jumpMap,afh,acfh);
        in.close();
        afh.close();
        acfh.close();
        //
        ofstream nfh(nodeFile.c_str());
        if(! nfh){
            cerr<<"Error: node file open failed"<<endl;
            exit(1);
        }
        nfh<<"#Segment\tChr\tStart\tEnd\tLen\tRefOrNot"<<endl;
        
        ofstream efh(edgeFile.c_str());
        if(! efh){
            cerr<<"Error: edge file open failed"<<endl;
            exit(1);
        }
        efh<<"#Source\tTarget\tOrigin1\tOrigin2"<<endl;

        ofstream covfh(covFile.c_str());
        if(! covfh){
            cerr<<"Error: coverage file open failed. "<<covFile<<endl;
            exit(1);
        }
        
        ofstream xcovfh(xcovFile.c_str());
        if(! xcovfh){
            cerr<<"Error: coverage file open failed. "<<xcovFile<<endl;
            exit(1);
        }
        
        ofstream cfh(chrListFile.c_str());
        if(! cfh){
            cerr<<"Error: chromosome list file open failed. "<<chrListFile<<endl;
            exit(1);
        }
        psAllPath(rfChrFile,tmpFolder,assListFile,sepStr,maxNode,mNodeLen,edgeMap,jumpMap,nfh,efh,covfh,xcovfh,cfh,pathDir);
        nfh.close();
        efh.close();
        covfh.close();
        xcovfh.close();
        cfh.close();
    }else{
        if(! opendir(upFolder.c_str())){
            cerr<<"Error: upload directory not exists"<<endl;
            exit(1);
        }
    }
    //-----------------------
    if(index){
        int flag = 0;
        GraphRange gr(upFolder,flag);
        gr.edgeWrite(range,ex,nocross,nthread);
    }
}

void g2v_usage(){
    cout<<"Usage: gfa2view --GFA input.gfa --index --ref REF#HAP --outDir output_dir"<<endl;
    cout<<"--sep     <String>   Delimiter between sample and haplotype names, by default: #"<<endl;
    cout<<"--GFA     <File>     Input GFA file"<<endl;
    cout<<"--ref     <String>   Reference name (assembly_name + delimiter + haplotype)"<<endl;
    cout<<"--refChr  <File>     When index the graph only consider reference chromosomes (or contig) contained in this file (saving time for situation that there are too many contigs)."<<endl; 
    cout<<"--outDir  <Dir>      Output directory"<<endl;
    cout<<"--index              Index the graph to speed up data visualization"<<endl;
    cout<<"--range   <Int>      number of reference nodes in a chunk, which is used for indexing the graph, by default: 2000"<<endl;
    cout<<"--cross              There are crosses between reference chromosomes. It will take more running time."<<endl;
    cout<<"--thread  <Int>      number of threads."<<endl;
    cout<<"--help               "<<endl;
}

int main(int argc,char **argv){
    char p[2] = {'#'};
    char *csep = p; 
    char *gfaFile = nullptr, *outDir = nullptr, *refName = nullptr, *rfChrFile = nullptr;
    bool gIndex = false;
    int range = 2000;
    int ex = 1000000;
    int nocross = 1;
    int nthread = 1;
    for(int i = 1; i < argc; ++i){
        if(strcmp(argv[i],"--sep") == 0 || strcmp(argv[i],"-sep") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --sep is used but there is no value"<<endl;
                return 1;
            }
            csep = argv[i];
        }else if(strcmp(argv[i],"--GFA") == 0 || strcmp(argv[i],"-GFA") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --GFA is used but there is no value"<<endl;
                return 1;
            }
            gfaFile = argv[i];
            
        }else if(strcmp(argv[i],"--refChr") == 0 || strcmp(argv[i],"-refChr") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --refChr is used but there is no value"<<endl;
                return 1;
            }
            rfChrFile = argv[i];
            
        }else if(strcmp(argv[i],"--ref") == 0 || strcmp(argv[i],"-ref") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --ref is used but there is no value"<<endl;
                return 1;
            }
            refName = argv[i];
            
        }else if(strcmp(argv[i],"--outDir") == 0 || strcmp(argv[i],"-outDir") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --outDir is used but there is no value"<<endl;
                return 1;
            }
            outDir = argv[i];
        }else if(strcmp(argv[i],"--index") == 0 || strcmp(argv[i],"-index") == 0){
            gIndex = true;
        }else if(strcmp(argv[i],"--cross") == 0 || strcmp(argv[i],"-cross") == 0){
            nocross = 0;
        }else if(strcmp(argv[i],"--range") == 0 || strcmp(argv[i],"-range") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --range is used but there is no value"<<endl;
                return 1;
            }
            range = atoi(argv[i]);
        }else if(strcmp(argv[i],"--thread") == 0 || strcmp(argv[i],"-thread") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --thread is used but there is no value"<<endl;
                return 1;
            }
            nthread = atoi(argv[i]);
        }else if(strcmp(argv[i],"--help") == 0 || strcmp(argv[i],"-help") == 0){
            g2v_usage();
            return 1;
        }else{
            cerr<<"Error: undefined parameter"<<endl;
            g2v_usage();
            return 1;
        }
    }
    //
    if(gfaFile == nullptr){
        if(outDir != nullptr && gIndex){
            cout<<"Only index the graph"<<endl;
        }else{
            cerr<<"Error: lack of parameters"<<endl;
            g2v_usage();
            return 1;
        }
    }
    
    if(outDir == nullptr){
        cerr<<"Error: there is no output directory (--outDir)"<<endl;
        return 1;
    }
    gfa2view(rfChrFile,gfaFile,refName,csep,range,ex,gIndex,nocross,nthread,outDir);
}



