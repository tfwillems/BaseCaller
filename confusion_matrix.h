#ifndef CONFUSION_MATRIX_H_
#define CONFUSION_MATRIX_H_

#include "bamtools/include/api/BamReader.h"

void process_reads(BamTools::BamReader& bam_reader, int32_t max_read_length, std::string& ref_seq,
                   int32_t* matrix_counts, int32_t* total_counts);


void print_confusion_matrix(int32_t* matrix_counts, int32_t* total_counts, int32_t max_read_length);

#endif
