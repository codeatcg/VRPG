
#include "gz.h"

using namespace std;

iogzbuf* iogzbuf::open(const char *inFile,const char *mode){
    iogz = gzopen(inFile,mode);
    if(iogz == NULL){
        return (iogzbuf *)0;
    }else{
        return this;
    }
}

void iogzbuf::close(){
    gzclose(iogz);
}

streambuf::int_type iogzbuf::underflow(){
    if(gptr() < egptr()){
        return traits_type::to_int_type(*gptr());
    }
    //
    int num = gptr() - eback();
    if(num > 4){
        num = 4;
    }
    memmove(buff+(4-num),gptr()-num,num);
    int rNum = gzread(iogz,buff+4,imem-4);
    if(rNum <= 0){
        return EOF;
    }
    setg(buff+(4-num),buff+4,buff+4+rNum);
    return traits_type::to_int_type(*gptr());
}

igzstream::igzstream(const char *inFile):istream(0){
    open(inFile);
}

igzstream* igzstream::open(const char *inFile){
    if(buff.open(inFile,"rb")){
        rdbuf(&buff);
        return this;
    }else{
        return (igzstream *)0;
    }
}

void igzstream::close(){
    buff.close();
}



