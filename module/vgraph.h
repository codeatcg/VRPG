
#include <vector>
#include <map>
#include <set>
#include <cstring>
#include <list>
#include <unordered_map>
#include <unordered_set>

#define FIELDSIZE 64

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

typedef struct{
    int node;
    int reStart;
    int len;
    char name[FIELDSIZE];
    uint8_t layer;
    char strand;
} NodeGene;

typedef struct{
    float start;
    float end;
    int layer;
    int strand;    
} FigGene;

inline void assSplit(const std::string &r_chr, const std::string &sep,std::string &tName,std::string &t_hap,std::string &tchr){
    auto t_pos = r_chr.find(sep);
    if(t_pos != std::string::npos){
        tName = r_chr.substr(0,t_pos);
        auto h_pos = t_pos + sep.length();
        t_pos = r_chr.find(sep,t_pos+1);
        if(t_pos != std::string::npos){
            t_hap = r_chr.substr(h_pos,t_pos-h_pos);
            tchr = r_chr.substr(t_pos + sep.length());
        }else{
            t_hap = "0";
            tchr = r_chr.substr(h_pos);
        }
    }else{
        tName = "REF";
        t_hap = "HAP";
        tchr = r_chr;
    }
}

inline char getMark(char ori_1,char ori_2, char sign='+'){
    char mark = '0';
    if(ori_1 == sign){
        if(ori_2 == sign){
            mark = '2';
        }else{
            mark = '3';
        }
    }else{
        if(ori_2 == sign){
            mark = '4';
        }else{
            mark = '5';
        }
    }
    return mark;
}

typedef int NodeType;

inline char revMark(char mark){
    char revVal = '0';
    if(mark == '2'){
        revVal = '5';
    }else if(mark == '5'){
        revVal = '2';
    }else{
        revVal = mark;
    }
    return revVal;
}

struct LEdge{
    NodeType node1;
    NodeType node2;
    char mark;
    
    bool operator < (const LEdge& tedge)const{
        if(node1 < tedge.node1){
            return true;
        }
        if(node1 > tedge.node1){
            return false;
        }

        if(node2 < tedge.node2){
            return true;
        }
        if(node2 > tedge.node2){
            return false;
        }

        if(mark < tedge.mark){
            return true;
        }
        if(mark > tedge.mark){
            return false;
        }
        return false;
    }
};

typedef struct LEdge NEdge;

typedef struct{
    NodeType node1;
    NodeType node2;
    char mark;
} CEdge;

typedef struct{
    NodeType node;
    char mark;
} ENode;

typedef struct{
    int s_nid;
    int e_nid;
    int gNum;
} Nid;

typedef std::map<std::string,int> Ndic;

typedef struct{
    int rByte;
    int ranNum;
} ChrRange;

typedef struct{
    int ranStart;
    int ranEnd;
    long long offByte;
    int ranNum;
    
} OneRange;

typedef struct{
    int len;
    int ass;
} LenAss;

struct INode{
    int node;
    int start;
    int pend;
    int achr;
    
    bool operator < (const INode &tnode)const{
        return node < tnode.node;
    }
};
typedef struct INode SANode;

typedef struct{
    int start;
    int pend;
    int achr;
} ANode;

typedef struct{
    int node;
    int start;
    int pend;
} RNode;

typedef struct{
    char diff;
    char ori;
} LagNode;

struct PaRa{
    int frag;
    int firNode;
    char firOri;
    std::list<LagNode> lag;
    bool operator < (const PaRa &pr)const{
        return frag < pr.frag;
    }
};

typedef PaRa PathRang;

typedef struct {
    long long rOffsize;
    int rCount;
} EdRang;

typedef struct{
    int start;
    int pend;
} RanPos;

std::string getSep(std::string &sepFile);

class GraphRange{

public:
    GraphRange(std::string &t_upDir,int index);
    std::vector<Ndic> draw_node;
    std::vector<float> draw_pos;
    std::vector<float> dnode_len;
    std::vector<std::map<std::string,int> > draw_edge;
    std::vector<int> genome;
    std::vector<NodeType> nnames;
    std::vector<int> hnGroup;
    std::vector<std::string> hLinks;
    std::vector<char> hDir;
    
    std::vector<float> ndGenePos;
    std::vector<std::string> geneVec;
    std::vector<int> layerVec;
    std::vector<char> strandVec;
    
    void formatGraph(std::string &ass,std::string &sChr,int sStart,int sEnd,int ex,int wStart,int wWidth,int wCut,int wY,int queryDep,bool sim,bool refSim);
    void edgeWrite(std::string &spChrFile,int rangeSize,int ex,int nocross,int nthread,int storeDep);
    
private:
    std::string nodeFile,edgeFile,pathDir,assFile,chrFile,comChrFile,sepFile;
    std::string upDir;
    std::string sep;
    
    int indexFlag;
    void conformEdge(NodeType &node1,NodeType &node2,char mark,std::unordered_map<NodeType,std::vector<ENode> > &iedge,std::unordered_map<NodeType,std::vector<ENode> > &oedge);
    void parseEdge(std::unordered_map<NodeType,std::vector<ENode> > &iedge,std::unordered_map<NodeType,std::vector<ENode> > &oedge);
    void parseNode(std::string &sChr,int sStart,int sEnd,int ex,std::vector<NodeType> &rangeNode,std::unordered_set<NodeType> &exNode,std::unordered_map<NodeType,LenAss> &info,int &realLen);
    void hAssNode(std::string &ass,int assNum,std::map<NEdge,int> &r_edge_dict,std::unordered_map<NodeType,Nid> &nid_dict);
    //
    void eAssFind(std::vector<char> &orient,std::vector<NodeType> &nodes,std::map<NEdge,int> &r_edge_dict,std::unordered_map<NodeType,Nid> &nid_dict,std::unordered_set<int> &ndGroup);
    //
    void parseIndex(int chrNum,int sStart,int sEnd,std::unordered_map<NodeType,std::vector<ENode> > &iedge,std::unordered_map<NodeType,std::vector<ENode> > &oedge);
    int getChrNum(std::string &sChr);
    int getAssNum(std::string &ass);
    void getExNode(int chrNum,int sStart,int sEnd,int ex,std::vector<NodeType> &rangeNode,std::unordered_set<NodeType> &exNode,std::unordered_map<NodeType,LenAss> &info,int &realLen);
    void queryDbPath(int assNum,int chrNum,int sStart,int sEnd,std::vector<std::vector<char> > &oriMulti,std::vector<std::vector<int> > &nodeMulti);
    void dxAssNode(int assNum,int chrNum,int sStart,int sEnd,std::map<NEdge,int> &r_edge_dict,std::unordered_map<NodeType,Nid> &nid_dict);
    
    //
    void splitRange(int rangeNum,std::unordered_map<std::string,int> &chrMap,std::unordered_map<std::string,int> &refChrMap,std::string &rndDxFile,std::string &rndFile,std::string &nspecFile,std::string &snFile);
    void getNrefEdge(std::string &rndFile,std::string &nspecFile,std::vector<NEdge> &resEdge);
    void getChrRmEdge(std::unordered_set<int> &ntNode,std::vector<NEdge> &chrRmEdge);
    
    void parseRange(std::vector<RNode> &chrRnode,std::vector<OneRange> &arcVec,int sStart,int sEnd,int ex,std::vector<NodeType> &rangeNode,std::unordered_set<NodeType> &exNode);
    void edgeRange(std::vector<RNode> &chrRnode,std::vector<OneRange> &arcVec,int sStart,int sEnd,int ex,int nocross,int storeDep,std::vector<NEdge> &chrRmEdge,std::unordered_map<NodeType,std::vector<ENode> > &iedge,std::unordered_map<NodeType,std::vector<ENode> > &oedge,std::set<NEdge> &r_edge_dict,std::unordered_set<NodeType> &nRefNode);
    
    void oneTask(std::unordered_map<NodeType,std::vector<ENode> > &iedge,std::unordered_map<NodeType,std::vector<ENode> > &oedge,std::vector<RNode> &chrRnode,std::vector<OneRange> &acrVec,std::vector<NEdge> &chrRmEdge,
             int ex,int nocross,int frStart,int frEnd,int storeDep,std::ofstream &tndfh,std::ofstream &tbfh,int *frNrefNum,int *frEdgeNum
    );
    
    void pthTask(std::unordered_map<NodeType,std::vector<int> > &ndCutMap,std::vector<RanPos> &allpos,char *header,int dxByte,int frStart,int frEnd,std::vector<std::ifstream> &pthVec,std::vector<std::ofstream> &xpthVec,std::vector<std::ofstream> &wpthVec);
    
    void fillNode(std::string &comChrFile,std::string &ndAssLFile,std::string &nrNodeFile,std::string &nrNumFile,std::string &snFile,std::string &nrdFile);
    void mergeDx(std::string &rndDxFile,std::string &nrNumFile,std::string &mgDxFile);
    void rangePath(std::vector<char> &orient,std::vector<NodeType> &nodes,std::unordered_map<NodeType,std::vector<int> > &ndCutMap,std::list<PathRang> &allPaRa);
    
    void indexPath(std::string &assFile,std::string &eIndexFile,std::string &bEdgeFile,int nthread);
    //
    void readRefGene(std::string &ovFile,std::string &gDxFile,int chrNum,int sStart,int sEnd,std::unordered_set<NodeType> &retainID,std::vector<NodeGene> &refNodeGene);
    void getFigGene(std::string &bwGeneFile,std::string &gDxFile,int chrNum,int sStart,int sEnd,float wPerK);
};


class QueryNode{

public:
    QueryNode(std::string &t_dbDir);
    std::vector<int> ndCov;
    std::string nodeAss,nodeChr;
    int nodeStart,nodeEnd;
    std::vector<std::vector<std::string> > geneList;
    
    void queryDbNode(int node);
    void queryDbCov(int node);
    void queryGene(int node,std::string &nodeAss);
    void queryAssCov(std::vector<int> &nodeVec,std::string &ass);
private:
    int node;
    std::string dbDir,comChrFile,assFile,sepFile,sep;
    std::vector<std::string> header;
    void getHeader();
    int countHeader();
};







