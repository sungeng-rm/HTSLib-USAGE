#ifndef GVCFWRITER_H
#define GVCFWRITER_H

#include <string>

#include "htslib/vcf.h"
#include "htslib/sam.h"
#include "sample.h"


class GVCFWriter
{
private:
    std::string                             m_output_path;
    std::string                             m_bam_path;
    htsFile*                                m_hts_out;
    htsThreadPool                           m_thread_pool;
    bcf_hdr_t*                              m_out_hdr;

private:
    int write_header();
    int add_sample_header(sam_hdr_t*      h_tmp);

public:
    GVCFWriter(std::string &out_path, std::string &bam_path);
    int write_vcf_record();
    ~GVCFWriter();
};



#endif //GVCFWRITER_H