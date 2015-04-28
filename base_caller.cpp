#include <iostream>

#include "bamtools/include/api/BamReader.h"

#include "confusion_matrix.h"
#include "error.h"
#include "io.h"
#include "seqio.h"


void compute_confusion_matrix(int32_t max_read_length, std::string bam_file, std::string fasta_file, std::string fasta_dir, bool skip_soft_clipped, std::ostream& out){
  BamTools::BamReader bam_reader;
  if (!bam_reader.Open(bam_file)) printErrorAndDie("Failed to open BAM file");

  std::string ref_seq;
  readFasta (fasta_file, fasta_dir, ref_seq);

  int32_t* matrix_counts = new int32_t [25*max_read_length]();
  int32_t* total_counts  = new int32_t [5*max_read_length]();
  process_reads(bam_reader, max_read_length, ref_seq, matrix_counts, total_counts, skip_soft_clipped);


  print_confusion_matrix(matrix_counts, total_counts, max_read_length, out);

  delete [] matrix_counts;
  delete [] total_counts;
}


int main(int argc, char* argv[]) {
  /*
 CIF* cif = read_cif(std::string(argv[1]));
 std::cout << *cif;
 delete cif;
  */

 int32_t max_read_length = 250;
 std::string bam_file    = "/data/twillems/base_calling_datasets/PhiX/PhiX_L002_R1.aligned.sorted.bam";
 std::string fasta_dir   = "/data/twillems/base_calling_datasets/PhiX/";
 std::string fasta_file  = "phiX.extended.fa";
 bool skip_soft_clipped  = false;
 std::ofstream output;
 output.open("phix_confusion_mat.inc_softclips.txt");
 compute_confusion_matrix(max_read_length, bam_file, fasta_file, fasta_dir, skip_soft_clipped, output);
 output.close();
}
