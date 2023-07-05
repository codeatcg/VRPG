
#include <fstream>
#include "vgraph.h"

using namespace std;

typedef struct{
    int neoID;
    char ori;
} Jnode;

inline int parseAlign(std::string &align){
    if(align == "*"){
        return 0;
    }
    
    int value = 0;
    int tvalue = 0;
    for(char x : align){
        if(x >= '0' && x <= '9'){
            tvalue = tvalue * 10 + (x - '0');
        }else{
            if(x == 'M' || x == 'I' || x == 'S' || x == 'X' || x == '='){
                value += tvalue;
            }
            tvalue = 0;
        }
    }
    return value;
}

inline void markSign(char mark,bool flip1,bool flip2,char &sign1,char &sign2){
    switch(mark){
        case '2':
            sign1 = '+';
            sign2 = '+';
            break;
        case '3':
            sign1 = '+';
            sign2 = '-';
            break;
        case '4':
            sign1 = '-';
            sign2 = '+';
            break;
        case '5':
            sign1 = '-';
            sign2 = '-';
    }
    if(flip1){
        if(sign1 == '+'){
            sign1 = '-';   
        }else{
            sign1 = '+';   
        }
    }
    if(flip2){
        if(sign2 == '+'){
            sign2 = '-';   
        }else{
            sign2 = '+';   
        }
    }
    
}

int readGFA(std::string &tmpFolder,std::ifstream &in,std::string &refStr,std::string &sep,std::unordered_map<int,int> &mNodeLen,std::map<NEdge,int> &edgeMap,std::map<NEdge,int> &jumpMap,std::ofstream &afh,std::ofstream &acfh);

int psRchrWalk(std::string &refPath,std::string &fullName,int refStart,std::unordered_map<int,int> &mNodeLen,int &neoID,std::unordered_set<int> &refNodeSet,
               std::unordered_set<int> &flipSet,std::unordered_map<int,int> &rCovMap,std::ofstream &nfh,std::ofstream &efh,std::ofstream &pfh);
               
void psNWalk(std::string &nrPath,std::string &ass,int walkStart,std::unordered_set<int> &flipSet,std::unordered_map<int,int> &tmap,std::unordered_map<int,int> &mNodeLen,std::unordered_set<int> &noutNode,std::ofstream &pfh,std::ofstream &nfh);

//

int psRchrPath(std::string &refPath,std::string &fullName,std::unordered_map<int,int> &mNodeLen,std::map<NEdge,int> &edgeMap,std::map<NEdge,int> &jumpMap,int &neoID,std::unordered_set<int> &refNodeSet,
               std::unordered_set<int> &flipSet,std::set<NEdge> &rEset,std::set<NEdge> &rJset,std::map<NEdge,Jnode> &jNeoMap,std::unordered_map<int,int> &rCovMap,std::ofstream &nfh,std::ofstream &efh,std::ofstream &pfh
                );
                
void psNPath(std::string &nrPath,std::string &ass,std::unordered_set<int> &flipSet,std::map<NEdge,Jnode> &jNeoMap,std::unordered_map<int,int> &tmap,std::unordered_map<int,int> &mNodeLen,std::unordered_set<int> &noutSet,std::ofstream &pfh,std::ofstream &nfh);

void psAllPath(std::string &tmpFolder,std::string &assFile,std::string &sepStr,int &neoID,std::unordered_map<int,int> &mNodeLen,std::map<NEdge,int> &edgeMap,std::map<NEdge,int> &jumpMap,
                std::ofstream &nfh,std::ofstream &efh,std::ofstream &covfh,std::ofstream &xcovfh,std::ofstream &cfh,std::string &pathDir
               );
               
void gfa2view(char *rfChrFile,char *gfaFile,char *refName,char *sep,int range,int ex,bool index,int nocross,int nthread,int storeDep,char *outDir);
void g2v_usage();








