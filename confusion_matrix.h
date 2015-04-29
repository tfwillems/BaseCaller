#ifndef CONFUSION_MATRIX_H_
#define CONFUSION_MATRIX_H_

#include <iostream>

#include "bamtools/include/api/BamReader.h"

void process_reads(BamTools::BamReader& bam_reader, int32_t max_read_length, int32_t ref_id, std::string& ref_seq, std::string& fasta_dir, bool skip_soft_clipped,
                   int32_t* matrix_counts, int32_t* total_counts, int32_t& forward, int32_t& backward);


void print_confusion_matrix(int32_t* matrix_counts, int32_t* total_counts, int32_t max_read_length, std::ostream& out);

#endif
