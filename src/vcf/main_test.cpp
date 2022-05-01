

#include <iostream>

#include "GVCFWriter.h"

int main(int argc, char* argv[]){
    if (argc!=3){
        std::cout<<"Usage: ./main-test ./outpath ./bampath"<<std::endl;
    }
    std::string outpath=argv[1];
    std::string bampath=argv[2];

    GVCFWriter* writer =new GVCFWriter(outpath, bampath);

    writer->write_vcf_record();

    delete writer;
    return 0;
}

