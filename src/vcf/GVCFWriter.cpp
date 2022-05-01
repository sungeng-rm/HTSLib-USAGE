#include <iostream>

#include "GVCFWriter.h"
#include "htslib/thread_pool.h"

GVCFWriter::GVCFWriter(std::string &out_path, std::string &bam_path)
:m_output_path(out_path)
,m_bam_path(bam_path)
{
    m_hts_out = hts_open(m_output_path.c_str(), "w");
    if (!m_hts_out) {
        std::cout << "failed to open gvcf output path" << strerror(errno);
        exit(1);
    }

    // pool
    m_thread_pool = {};
    if (!(m_thread_pool.pool = hts_tpool_init(32))) {
        std::cout << "No Pool Mem";
        exit(1);
    }
    hts_set_opt(m_hts_out, HTS_OPT_THREAD_POOL, &m_thread_pool);

    write_header();

}

GVCFWriter::~GVCFWriter(){
    hts_close(m_hts_out);
    hts_tpool_destroy(m_thread_pool.pool);
    bcf_hdr_destroy(m_out_hdr);
}

int GVCFWriter::write_header(){
    htsFile*        bam_fp;
    sam_hdr_t*      h_tmp;
    
    //create an empty VCF header.
    m_out_hdr = bcf_hdr_init("w");

    if ( !(bam_fp = hts_open(m_bam_path.c_str(), "r")) || !(h_tmp = sam_hdr_read(bam_fp))) {
        std::cout << "No Bam path :" << strerror(errno);
        exit(1);
    }
    
    if ( add_sample_header(h_tmp)){
        std::cout<<"Add Sample Error :"<< strerror(errno);
    }

    //FORMAT
    bcf_hdr_append(m_out_hdr, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">");
    bcf_hdr_append(m_out_hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">");
    bcf_hdr_append(m_out_hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
    bcf_hdr_append(m_out_hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(m_out_hdr, "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP observed within the GVCF block\">");
    bcf_hdr_append(m_out_hdr, "##FORMAT=<ID=PGT,Number=1,Type=String,Description=\"Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another; will always be heterozygous and is not intended to describe called alleles\">");
    bcf_hdr_append(m_out_hdr, "##FORMAT=<ID=PID,Number=1,Type=String,Description=\"Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group\">");
    bcf_hdr_append(m_out_hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">");
    bcf_hdr_append(m_out_hdr, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phasing set (typically the position of the first variant in the set)\">");
    bcf_hdr_append(m_out_hdr, "##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description=\"Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)\">");
    bcf_hdr_append(m_out_hdr, "##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.\">");
    
    for(int i=0; i<h_tmp->n_targets; ++i){
        std::string chr_name=*(h_tmp->target_name + i);
        std::string len=std::to_string(*(h_tmp->target_len + i));
        std::string contig="##contig=<Id=" + chr_name + "," + "length=" + len + ">";
        bcf_hdr_append(m_out_hdr, contig.c_str());
    }
    //INFO
    bcf_hdr_append(m_out_hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">");
    bcf_hdr_append(m_out_hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">");
    bcf_hdr_append(m_out_hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
    bcf_hdr_append(m_out_hdr, "##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities\">");
    bcf_hdr_append(m_out_hdr, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">");
    bcf_hdr_append(m_out_hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">");
    bcf_hdr_append(m_out_hdr, "##INFO=<ID=ExcessHet,Number=1,Type=Float,Description=\"Phred-scaled p-value for exact test of excess heterozygosity\">");
    bcf_hdr_append(m_out_hdr, "##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher\'s exact test to detect strand bias\">");
    bcf_hdr_append(m_out_hdr, "##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description=\"Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation\">");
    bcf_hdr_append(m_out_hdr, "##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=\"Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed\">");
    bcf_hdr_append(m_out_hdr, "##INFO=<ID=MLEAF,Number=A,Type=Float,Description=\"Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed\">");
    bcf_hdr_append(m_out_hdr, "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">");
    bcf_hdr_append(m_out_hdr, "##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities\">");
    bcf_hdr_append(m_out_hdr, "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence/Quality by Depth\">");
    bcf_hdr_append(m_out_hdr, "##INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description=\"Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.\">");
    bcf_hdr_append(m_out_hdr, "##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias\">");
    bcf_hdr_append(m_out_hdr, "##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">");

    // write head
    if (bcf_hdr_write(m_hts_out, m_out_hdr) != 0) {
        std::cout << "write vcf Header Error :"<< strerror(errno);
        exit(1);
    }
    
    sam_hdr_destroy(h_tmp);
    hts_close(bam_fp);
    return 0;
}

int GVCFWriter::write_vcf_record(){
    bcf1_t *rec    = bcf_init1();
    if (!rec) {
        std::cout << "Init bcf1_t Error :"<< strerror(errno);
        exit(1);
    }
    // .. CHROM
    rec->rid = bcf_hdr_name2id(m_out_hdr, "20");
    // .. POS
    // rec->pos = 14369;
    // .. ID
    bcf_update_id(m_out_hdr, rec, "rs6054257");
    #if 0
    // .. REF and ALT
    test_update_alleles(m_out_hdr, rec);
    const char *alleles[2] = { "G", "A" };
    bcf_update_alleles_str(m_out_hdr, rec, "G,A");
    check_alleles(rec, alleles, 2);
    // .. QUAL
    rec->qual = 29;
    // .. FILTER
    int32_t tmpi = bcf_hdr_id2int(m_out_hdr, BCF_DT_ID, "PASS");
    bcf_update_filter(m_out_hdr, rec, &tmpi, 1);
    // .. INFO
    tmpi = 3;
    bcf_update_info_int32(m_out_hdr, rec, "NS", &tmpi, 1);
    tmpi = 500;
    bcf_update_info_int32(m_out_hdr, rec, "DP", &tmpi, 1);
    tmpi = 100000;
    bcf_update_info_int32(m_out_hdr, rec, "DP", &tmpi, 1);
    tmpi = 14;
    bcf_update_info_int32(m_out_hdr, rec, "DP", &tmpi, 1);
    tmpi = -127;
    bcf_update_info_int32(m_out_hdr, rec, "NEG", &tmpi, 1);
    float tmpf = 0.5;
    bcf_update_info_float(m_out_hdr, rec, "AF", &tmpf, 1);
    bcf_update_info_flag(m_out_hdr, rec, "DB", NULL, 1);
    bcf_update_info_flag(m_out_hdr, rec, "H2", NULL, 1);
    // .. FORMAT
    int32_t *tmpia = (int*)malloc(bcf_hdr_nsamples(m_out_hdr)*2*sizeof(int));
    tmpia[0] = bcf_gt_phased(0);
    tmpia[1] = bcf_gt_phased(0);
    tmpia[2] = bcf_gt_phased(1);
    tmpia[3] = bcf_gt_phased(0);
    tmpia[4] = bcf_gt_unphased(1);
    tmpia[5] = bcf_gt_unphased(1);
    bcf_update_genotypes(m_out_hdr, rec, tmpia, bcf_hdr_nsamples(m_out_hdr)*2);
    tmpia[0] = 48;
    tmpia[1] = 48;
    tmpia[2] = 43;
    bcf_update_format_int32(m_out_hdr, rec, "GQ", tmpia, bcf_hdr_nsamples(m_out_hdr));
    tmpia[0] = 0;
    tmpia[1] = 0;
    tmpia[2] = 1;
    bcf_update_format_int32(m_out_hdr, rec, "DP", tmpia, bcf_hdr_nsamples(m_out_hdr));
    tmpia[0] = 1;
    tmpia[1] = 100000;
    tmpia[2] = 1;
    bcf_update_format_int32(m_out_hdr, rec, "DP", tmpia, bcf_hdr_nsamples(m_out_hdr));
    tmpia[0] = 1;
    tmpia[1] = 8;
    tmpia[2] = 5;
    bcf_update_format_int32(m_out_hdr, rec, "DP", tmpia, bcf_hdr_nsamples(m_out_hdr));
    tmpia[0] = 51;
    tmpia[1] = 51;
    tmpia[2] = 51;
    tmpia[3] = 51;
    tmpia[4] = bcf_int32_missing;
    tmpia[5] = bcf_int32_missing;
    bcf_update_format_int32(m_out_hdr, rec, "HQ", tmpia, bcf_hdr_nsamples(m_out_hdr)*2);
    char *tmp_str[] = {"String1","SomeOtherString2","YetAnotherString3"};
    bcf_update_format_string(m_out_hdr, rec, "TS", (const char**)tmp_str, 3);
    tmp_str[0] = "LongerStringRequiringBufferReallocation";
    bcf_update_format_string(m_out_hdr, rec, "TS", (const char**)tmp_str, 3);
    tmp_str[0] = "String1";
    bcf_update_format_string(m_out_hdr, rec, "TS", (const char**)tmp_str, 3);
    if ( bcf_write1(m_hts_out, m_out_hdr, rec)!=0 ) {};

    bcf_update_genotypes(m_out_hdr, rec, tmpia, bcf_hdr_nsamples(m_out_hdr)*2);
    if ( bcf_write1(m_hts_out, m_out_hdr, rec)!=0 ) {

    }

    free(tmpia);
    free(tmpfa);
#endif
    bcf_destroy1(rec);
    return bcf_write1(m_hts_out, m_out_hdr, rec) != 0;
}

int GVCFWriter::add_sample_header(sam_hdr_t* h_tmp){
    bam_sample_t*   sm = NULL;
    int             i;
    sm = bam_smpl_init();
    bam_smpl_add(sm, m_bam_path.c_str(), sam_hdr_str(h_tmp));
    for (i = 0; i < sm->n; i++) {
        bcf_hdr_add_sample(m_out_hdr, sm->smpl[i]);
    }
    bcf_hdr_add_sample(m_out_hdr, NULL);
    bam_smpl_destroy(sm);
    return 0;

    // bcf_hdr_add_sample();
    // bcf_hdr_sync();
}

// int bcf_get_format_values(const bcf_hdr_t *m_out_hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type)