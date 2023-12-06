
#include<cstring>
#include<iostream>
int addRef_main(int argc,char **argv);
int ndg_main(int argc,char **argv);

void usage(){
    std::cout<<"Usage: GraphAnno command options"<<std::endl;
    std::cout<<"Command:"<<std::endl;
    std::cout<<"addRef       Create files for reference gene track plot."<<std::endl;
    std::cout<<"nodeGene     Create files for interactive view of genes with which a node overlaps."<<std::endl;
}
int main(int argc,char **argv){
    if(argc < 2){
        usage();
        return 1;
    }

    if(strcmp(argv[1],"addRef") == 0){
        addRef_main(argc - 1,argv + 1);
    }else if(strcmp(argv[1],"nodeGene") == 0){
        ndg_main(argc - 1,argv + 1);
    }else{
        std::cerr<<"Error: undefined command"<<std::endl;
        usage();
        return 1;
    }
}

