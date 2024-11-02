
#include <iostream>
#include <fstream>
//#include <unistd.h>
#include <set>
#include <stdlib.h>
#include <sstream>
#include <map>
#include <unordered_set>
#include <string.h>

using namespace std;

typedef struct{
    long long offset;
    int len;
} SeqDx;

string revComSeq(string &seq,map<char,char> &trans){
    string rvstr = "";
    int n = seq.size() - 1;
    for(int pos = n; pos >= 0; pos--){
        if(trans.find(seq[pos]) == trans.end()){
            cerr<<"Warning: invalid character in the sequence. "<<seq[pos]<<endl;
            //exit(1);
        }else{
            rvstr += trans[seq[pos]];
        }
    }
    return rvstr;
}

void makeTrans(map<char,char> &trans){
    trans.emplace('A','T');
    trans.emplace('T','A');
    trans.emplace('C','G');
    trans.emplace('G','C');
    trans.emplace('N','N');
    trans.emplace('a','t');
    trans.emplace('t','a');
    trans.emplace('c','g');
    trans.emplace('g','c');
    trans.emplace('n','n');
}

void node2seq(char *gfaFile,string &upDir){
    unordered_set<int> flipSet;
    int intSize = sizeof(int);
    string flipFile = upDir + "/flip.bw";
    ifstream in(flipFile.c_str());
    if(in){
        int num,node;
        in.read((char *)&num,intSize);
        for(int i = 0; i < num; ++i){
            in.read((char *)&node,intSize);
            flipSet.insert(node);
        }
        in.close();
    }
    ifstream gfh(gfaFile);
    if(! gfh){
        cerr<<"Error: file open failed. "<<gfaFile<<endl;
        exit(1);
    }
    string outFile = upDir + "/graphSeq.fa";
    string dxFile = upDir + "/graphSeq.fa.bdx";
    string snFile = upDir + "/node.sort.bw";
    string addFile = upDir + "/add.node.info";
    
    ifstream nfh(snFile.c_str());
    if(! nfh){
        cerr<<"Error: file open failed. "<<snFile<<endl;
        exit(1);
    }
    int total;
    nfh.read((char *)&total,intSize);
    nfh.close();
    
    ofstream ofh(outFile.c_str());
    if(! ofh){
        cerr<<"Error: file open failed. "<<outFile<<endl;
        exit(1);
    }
    
    string line;
    string tag,snode,seq;
    int tgLen = 0;
    long long cumLen = 0;
    stringstream strStream;
    SeqDx *sdx = new SeqDx[total];
    map<char,char> trans;
    makeTrans(trans);
    while(getline(gfh,line)){
        if(line[0] == 'S'){
            strStream << line;
            strStream >> tag >> snode >> seq;
            strStream.clear();
            strStream.str("");
            //
            int tnode;
            if(snode[0] == 's'){
                string msNode = snode.substr(1);
                ofh << ">" << msNode << endl;
                tnode = atoi(msNode.c_str());
                tgLen = snode.size() - 1;
            }else{
                ofh << ">" << snode << endl;
                tnode = atoi(snode.c_str());
                tgLen = snode.size();
            }
            if(flipSet.empty()){
                ofh << seq <<endl;
            }else{
                if(flipSet.find(tnode) != flipSet.end()){
                    string revSeq = revComSeq(seq,trans);
                    ofh << revSeq <<endl;
                }else{
                    ofh << seq <<endl;
                }
            }
            
            cumLen += tgLen + 2;
            int seqLen = seq.size();
            SeqDx tx = {cumLen,seqLen};
            sdx[tnode-1] = tx;
            //
            cumLen += seqLen + 1;
        }
    }
    //
    ifstream adfh(addFile.c_str());
    int sceNode,nwNode;
    //
    if(adfh){
        while(getline(adfh,line)){
            strStream << line;
            strStream >> tag >> sceNode >> nwNode;
            strStream.clear();
            strStream.str("");
            if(tag == "S"){
                sdx[nwNode-1] = sdx[sceNode - 1];
            }else{
                SeqDx tx = {0,sceNode};
                sdx[nwNode-1] = tx;
            }
            
        }
        
        adfh.close();
    }
    //
    ofstream xfh(dxFile.c_str());
    if(! xfh){
        cerr<<"Error: file open failed. "<<dxFile<<endl;
        exit(1);
    }
    int uSize = sizeof(SeqDx);
    for(int k = 0; k < total; k++){
        xfh.write((char *)(sdx+k),uSize);
    }
    delete []sdx;
    //
    gfh.close();
    ofh.close();
    xfh.close();
}

//
void ndsq_usage(){
    cout<<"--gfaFile     <File>    input GFA or rGFA file."<<endl;
    cout<<"--upDir       <Dir>     'upload' directory which including files generated by 'gfa2view' or 'vrpg_preprocess.py'."<<endl;
}

int main(int argc,char **argv){
    char *gfaFile = nullptr;
    string upDir = "";
    for(int i = 1; i < argc; ++i){
        if(strcmp(argv[i],"--gfaFile") == 0 || strcmp(argv[i],"-gfaFile") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --gfaFile is used but there is no value."<<endl;
                return 1;
            }
            gfaFile = argv[i];
        }else if(strcmp(argv[i],"--upDir") == 0 || strcmp(argv[i],"-upDir") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --upDir is used but there is no value."<<endl;
                return 1;
            }
            upDir = argv[i];
        }else if(strcmp(argv[i],"--help") == 0 || strcmp(argv[i],"-help") == 0){
            ndsq_usage();
            return 1;
        }else{
            cerr<<"Error: undefined parameter."<<endl;
            ndsq_usage();
            return 1;
        }
        
    }
    if(gfaFile == nullptr){
        cerr<<"Error: --gfaFile is required."<<endl;
        return 1;
    }
    if(upDir == ""){
        cerr<<"Error: --upDir is required."<<endl;
        return 1;
    }
    node2seq(gfaFile,upDir);
}



