
#include <vector>
#include <cstdint>
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
    int node1;
    int node2;
    int reStart1;
    int reStart2;
    char name[FIELDSIZE];
    uint8_t layer;
    char strand;
    char type;
    uint8_t num;
} GeneNode;

typedef struct{
    int node1;
    int node2;
    int reStart1;
    int reStart2;
    char name[FIELDSIZE];
    uint16_t layer;
    char strand;
    float score;
} BedNode;

typedef struct{
    float start;
    float end;
    int layer;
    char strand;
    char margin;
} FigGene;

typedef struct{
    float start;
    float end;
    float score;
    //int layer;
    char strand;
    std::string name;
} FigBed;

typedef struct{
    float start;
    float end;
    int layer;
    char strand;
    uint8_t num; 
} FigEle;

typedef struct{
    long long offset;
    int len;
} SeqDx;

inline void asmSplit(const std::string &r_chr, const std::string &sep,std::string &tName,std::string &t_hap,std::string &tchr){
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
    int asmb;
} LenAsm;

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
    int pathNum;
    int firNode;
    char firOri;
    std::list<LagNode> lag;
    std::list<int> posList;
    std::list<std::string> cigarList;
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

typedef struct{
    int pnum;
    int start;
    int end;
    int loc;
} PathPos;

bool numSort(PathPos &a, PathPos &b);

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
    std::vector<int> hEdgeAsm;
    
    //annotation gene
    std::vector<float> ndGenePos;
    std::vector<std::string> geneVec;
    std::vector<int> layerVec;
    std::vector<char> strandVec;
    std::vector<char> mgFlagVec;
    
    //annotation exon
    std::vector<float> ndExonPos;
    std::vector<std::string> rnaVec;
    std::vector<int> eLayerVec;
    std::vector<char> eStrandVec;
    std::vector<char> eFlagVec;
    std::vector<int> eNumVec;
    
    //annotation cds
    std::vector<float> ndCDSPos;
    std::vector<std::string> cdsVec;
    std::vector<int> cLayerVec;
    std::vector<int> cNumVec;
    
    //annotation track in bed format
    std::vector<std::string> tkNameVec;
    std::vector<std::string> tkDesVec;
    std::vector<int> tkColVec;
    std::vector<int> tkCumVec;
    std::vector<int> tkItem;
    std::vector<float> rBedPos;
    std::vector<std::string> rBedName;
    std::vector<int> rBedLayer;
    std::vector<float> rBedScore;
    std::vector<char> rBedStrand;
    
    std::vector<int> tickValue;
    std::vector<float> tickPos;
    
    float figScale;
    
    //alignment
    std::vector<std::string> qChr;
    std::vector<int> qStart;
    std::vector<int> qEnd;
    std::vector<std::string> qPath;
    std::vector<std::string> qCigar;
    
    void formatGraph(std::string &asmb,std::string &sChr,int sStart,int sEnd,int ex,int wStart,int wWidth,int wCut,int wY,int queryDep,int varLen,bool sim,bool refSim);
    void edgeWrite(std::string &spChrFile,int rangeSize,int ex,int nocross,int nthread,int storeDep);
    
private:
    std::string nodeFile,edgeFile,pathDir,asmFile,chrFile,comChrFile,sepFile;
    std::string upDir;
    std::string sep;
    
    int indexFlag;
    // 
    void conformEdge(NodeType &node1,NodeType &node2,char mark,std::unordered_map<NodeType,std::vector<ENode> > &iedge,std::unordered_map<NodeType,std::vector<ENode> > &oedge);
    //---------- read unindexed dada -------------------
    void parseEdge(std::unordered_map<NodeType,std::vector<ENode> > &iedge,std::unordered_map<NodeType,std::vector<ENode> > &oedge);
    void parseNode(std::string &sChr,int sStart,int sEnd,int ex,std::vector<NodeType> &rangeNode,std::unordered_map<NodeType,char> &exNode,std::unordered_map<NodeType,LenAsm> &info,int &realLen);
    //
    void hAsmNode(std::string &taskDir,bool refSim,std::unordered_map<NodeType,LenAsm> &info,int asmNum,std::map<NEdge,int> &r_edge_dict,std::unordered_map<NodeType,Nid> &nid_dict);
    void hEdge2fig(std::unordered_map<NodeType,Nid> &nid_dict,std::map<NEdge,int> &h_edge_dict,std::map<NEdge,int> &r_edge_dict);
    void hAsmNode2(bool refSim,std::vector<int> &asmNumVec,std::map<NEdge,int> &r_edge_dict,std::unordered_map<NodeType,Nid> &nid_dict);
    
    //---------- for sequence alignment -------------------
    bool f_cigar2pos(int qStart,int rStart,std::unordered_map<NodeType,LenAsm> &info,std::string &cigar,std::vector<int> &nodes,std::vector<std::string> &qCigarVec,std::vector<int> &qPosVec);
    bool f_path2pos(int rStart,std::unordered_map<NodeType,LenAsm> &info,std::vector<int> &nodes,std::vector<int> &qPosVec);
    
    //---------- intermediate functions used by functions to deal with both indexed and unindexed data -------------------
    std::string subCigar(std::string &rgCigar,int ndStart,int ndEnd);
    /*  
        Highlight assembly path:
        visCigar is true.
            
        Highlight query path (sequence alignment):
        visCigar is true if all nodes that composing the query sequence are in the indexed graph, 
        or visCigar is false.
    */
    void eAsmFind(bool visCigar,bool refSim,std::vector<char> &orient,std::vector<NodeType> &nodes,std::vector<int> &qPosVec,std::string &rgCigar,std::string &rgName,std::map<NEdge,int> &r_edge_dict,std::unordered_map<NodeType,Nid> &nid_dict,std::unordered_set<int> &ndGroup);
    void eAsmFind2(bool refSim,int asmCode,std::vector<char> &orient,std::vector<NodeType> &nodes,std::unordered_map<NodeType,Nid> &nid_dict,std::unordered_set<int> &ndGroup,std::map<NEdge,int> &h_edge_dict);
    int getChrNum(std::string &sChr);
    int getAsmNum(std::string &asmb);
    // get code number for a list of assemblies
    void getAsmNum2(std::set<std::string> &nameSet,std::vector<int> &asmNumVec);
    
    //---------- read indexed data -------------------    
    // read indexed edge data
    void parseIndex(int chrNum,int sStart,int sEnd,std::unordered_map<NodeType,std::vector<ENode> > &iedge,std::unordered_map<NodeType,std::vector<ENode> > &oedge);    
    // read indexed node data
    void getExNode(int chrNum,int sStart,int sEnd,int ex,std::vector<NodeType> &rangeNode,std::unordered_map<NodeType,char> &exNode,std::unordered_map<NodeType,LenAsm> &info,int &realLen,int &realStart);
    // read index path dada
    void queryDbPath(bool formR,int asmNum,int chrNum,int sStart,int sEnd,std::unordered_map<NodeType,LenAsm> &info,std::vector<std::vector<char> > &oriMulti,std::vector<std::vector<int> > &nodeMulti,std::vector<std::vector<int> > &qPosMulti,std::vector<std::string> &cigarMulti,std::vector<std::string> &nameMulti);    
    
    void dxAsmNode(bool refSim,int asmNum,int chrNum,int sStart,int sEnd,std::unordered_map<NodeType,LenAsm> &info,std::map<NEdge,int> &r_edge_dict,std::unordered_map<NodeType,Nid> &nid_dict);    
    void dxAsmNode2(bool refSim,std::vector<int> &asmVec,int chrNum,int sStart,int sEnd,std::unordered_map<NodeType,LenAsm> &info,std::map<NEdge,int> &r_edge_dict,std::unordered_map<NodeType,Nid> &nid_dict);
    
    //---------- create indexes for node, edge and path data -------------------
    void splitRange(int rangeNum,std::unordered_map<std::string,int> &chrMap,std::unordered_map<std::string,int> &refChrMap,std::string &rndDxFile,std::string &rndFile,std::string &nspecFile,std::string &snFile);
    //
    void getNrefEdge(std::string &rndFile,std::string &nspecFile,std::vector<NEdge> &resEdge);
    void getChrRmEdge(std::unordered_set<int> &ntNode,std::vector<NEdge> &chrRmEdge);
    
    void parseRange(std::vector<RNode> &chrRnode,std::vector<OneRange> &arcVec,int sStart,int sEnd,int ex,std::vector<NodeType> &rangeNode,std::unordered_set<NodeType> &exNode);
    void edgeRange(std::vector<RNode> &chrRnode,std::vector<OneRange> &arcVec,int sStart,int sEnd,int ex,int nocross,int storeDep,std::vector<NEdge> &chrRmEdge,std::unordered_map<NodeType,std::vector<ENode> > &iedge,std::unordered_map<NodeType,std::vector<ENode> > &oedge,std::set<NEdge> &r_edge_dict,std::unordered_set<NodeType> &nRefNode);
    
    void oneTask(std::unordered_map<NodeType,std::vector<ENode> > &iedge,std::unordered_map<NodeType,std::vector<ENode> > &oedge,std::vector<RNode> &chrRnode,std::vector<OneRange> &acrVec,std::vector<NEdge> &chrRmEdge,
             int ex,int nocross,int frStart,int frEnd,int storeDep,std::ofstream &tndfh,std::ofstream &tbfh,int *frNrefNum,int *frEdgeNum
    );
    //
    void fillNode(std::string &comChrFile,std::string &ndAsmLFile,std::string &nrNodeFile,std::string &nrNumFile,std::string &snFile,std::string &nrdFile);
    void mergeDx(std::string &rndDxFile,std::string &nrNumFile,std::string &mgDxFile);
    //
    void pthTask(bool formR,std::vector<int> &allLen,std::unordered_map<NodeType,std::vector<int> > &ndCutMap,std::vector<RanPos> &allpos,char *header,int dxByte,int frStart,int frEnd,std::vector<std::ifstream> &pthVec,std::vector<std::ofstream> &xpthVec,std::vector<std::ofstream> &wpthVec);
    void rangePath(bool formR,int num,std::vector<std::string> &qCigarVec,std::vector<int> &qPosVec,std::vector<char> &orient,std::vector<NodeType> &nodes,std::unordered_map<NodeType,std::vector<int> > &ndCutMap,std::list<PathRang> &allPaRa);
    void indexPath(bool formR,std::string &asmFile,std::string &eIndexFile,std::string &bEdgeFile,std::string &snFile,int nthread);
    
    void getAllLen(std::string &snFile,std::vector<int> &allLen);
    void cigar2pos(int qStart,int rStart,std::vector<int> &allLen,std::string &cigar,std::vector<int> &nodes,std::vector<std::string> &qCigarVec,std::vector<int> &qPosVec);
    void path2pos(int rStart,std::vector<int> &allLen,std::vector<int> &nodes,std::vector<int> &qPosVec);
    //---------- read indexed annotation data -------------------
    void readRefGene(std::string &ovFile,std::string &gDxFile,int chrNum,int sStart,int sEnd,std::unordered_set<NodeType> &retainID,std::unordered_map<NodeType,char> &exNode,std::vector<GeneNode> &refNodeGene,std::unordered_set<std::string> midExon);
    void getFigGene(std::string &bwGeneFile,std::string &gDxFile,int chrNum,int sStart,int sEnd,float wPerK,std::unordered_map<NodeType,char> &exNode);
    
    void readRefBed(std::string &ovFile,std::string &gDxFile,int chrNum,int sStart,int sEnd,std::unordered_set<NodeType> &retainID,std::unordered_map<NodeType,char> &exNode,std::vector<BedNode> &refBedNode);
    void getFigBed(std::string &tkDesFile,std::string &bwGeneFile,std::string &gDxFile,int chrNum,int sStart,int sEnd,float wPerK,std::unordered_map<NodeType,char> &exNode);
};

//
class QueryNode{

public:
    QueryNode(std::string &t_dbDir);
    std::vector<int> ndCov;
    std::string nodeAsm,nodeChr;
    int nodeStart,nodeEnd;
    std::vector<std::vector<std::string> > geneList;
    
    std::string nodeSeq;
    
    void fetchNdSeq(int node);
    void queryDbNode(int node);
    void queryDbCov(int node);
    void queryGene(int node,std::string &nodeAsm);
    void queryAsmCov(std::vector<int> &nodeVec,std::string &asmb);
private:
    int node;
    std::string dbDir,comChrFile,asmFile,sepFile,sep;
    std::vector<std::string> header;
    void getHeader();
    int countHeader();
};







