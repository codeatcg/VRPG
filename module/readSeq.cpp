
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

typedef struct{
    long long offset;
    int len;
} SeqDx;

void readSeq(int node,string &seqFile,string &dxFile){
    ifstream in(seqFile.c_str());
    ifstream dxfh(dxFile.c_str());
    int uSize = sizeof(SeqDx);
    long long dxOffset = (long long)(node - 1) * uSize;
    dxfh.seekg(dxOffset);
    SeqDx tdx;
    dxfh.read((char *)&tdx,uSize);
    if(tdx.offset == 0){
        string outStr(tdx.len,'N');
    }else{
        in.seekg(tdx.offset);
        char *seq = new char[tdx.len+1];
        in.read(seq,tdx.len);
        seq[tdx.len] = '\0';
        
        string outStr = seq;
        //cout<<seq<<endl;
        delete []seq;
        cout<<outStr<<endl;
    }
    in.close();
    dxfh.close();
}

int main(int argc,char **argv){
    string seqFile = argv[1];
    string dxFile = argv[2];
    int node = atoi(argv[3]);
    readSeq(node,seqFile,dxFile);
}




