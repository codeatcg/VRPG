
#include <vector>
#include <map>
#include <set>
#include <cstring>
#define NODESIZE 12
inline void assSplit(const std::string &r_chr, const std::string &sep,std::string &tName,std::string &t_hap,std::string &tchr){
    auto t_pos = r_chr.find(sep);
    if(t_pos != std::string::npos){
        tName = r_chr.substr(0,t_pos);
        auto h_pos = t_pos + sep.length();
        t_pos = r_chr.find(sep,t_pos+1);
        if(t_pos != std::string::npos){
            t_hap = r_chr.substr(h_pos,t_pos-h_pos);
            tchr = r_chr.substr(t_pos + sep.length());
        } 
    }else{
        tName = "REF";
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
    std::string node1;
    std::string node2;
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
    char node1[NODESIZE];
    char node2[NODESIZE];
    char mark;
} CEdge;

typedef struct{
    std::string node;
    char mark;
} ENode;

typedef struct{
    int s_nid;
    int e_nid;
    int gNum;
} Nid;

typedef std::map<std::string,int> Ndic;


class GraphRange{

public:
    GraphRange(std::string &t_nodeFile,std::string &t_edgeFile,std::string &t_pathFile,std::string &t_sepFile);
    GraphRange(std::string &t_nodeFile,std::string &t_edgeFile,std::string &t_pathFile,std::string &t_sepFile,std::string &t_chrFile,std::string &t_bEdgeFile,std::string &t_eIndexFile);
    
    std::vector<Ndic> draw_node;
    std::vector<float> draw_pos;
    std::vector<float> dnode_len;
    std::vector<std::map<std::string,int> > draw_edge;
    std::vector<std::string> genome;
    std::vector<std::string> nnames;
    std::vector<int> hnGroup;
    std::vector<std::string> hLinks;
    
    void formatGraph(std::string &ass,std::string &sChr,int sStart,int sEnd,int ex,int wStart,int wWidth,int wCut,int wY);
    void edgeWrite(int rangeSize,int ex);
private:
    std::string nodeFile,edgeFile,pathFile,sepFile,sep,chrFile,bEdgeFile,eIndexFile;
    int indexFlag;
    //
    std::map<std::string, std::map<std::string,std::vector<ENode> > > edge;
    void getSep();
    void conformEdge(std::string &node1,std::string &node2,char mark);
    void parseEdge();
    void parseIndex(std::string &sChr,int sStart,int sEnd,int ex);
    void parseNode(std::string &sChr,int sStart,int sEnd,int ex,std::vector<std::string> &rangeNode,std::set<std::string> &exNode,std::map<std::string,int> &info,std::map<std::string,std::string> &assDict,int &realLen);
    void hAssNode(std::string &ass,std::map<NEdge,int> &r_edge_dict,std::map<std::string,Nid> &nid_dict);
    
    void edgeRange(std::string &sChr,int sStart,int sEnd,int ex,std::map<NEdge,int> &r_edge_dict);
    
    
};

//






