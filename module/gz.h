

#include <iostream>
#include <streambuf>
#include <cstdio>
#include <zlib.h>
#include <cstring>

class iogzbuf : public std::streambuf{
    gzFile iogz;
    static const int imem = 512;
    char buff[imem];  
public:
    iogzbuf(){
        setg(buff+4,buff+4,buff+4);
    }
    iogzbuf* open(const char *inFile,const char *mode);
    void close();

    virtual int_type underflow();
};

class igzstream : public std::istream{
    iogzbuf buff;
public:
    igzstream(const char *inFile);
    igzstream* open(const char *inFile);
    void close();
};





