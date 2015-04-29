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
  int32_t ref_id;
  if (fasta_file.compare("N/A") == 0)
    ref_id = -2;
  else {
    readFasta(fasta_file, fasta_dir, ref_seq);
    ref_id = 0;
  }

  int32_t* matrix_counts = new int32_t [25*max_read_length]();
  int32_t* total_counts  = new int32_t [5*max_read_length]();
  int32_t forward = 0, backward = 0;
  process_reads(bam_reader, max_read_length, ref_id, ref_seq, fasta_dir, skip_soft_clipped, matrix_counts, total_counts, forward, backward);

  out << forward  << "\n"
      << backward << std::endl;
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

 std::string bam_file    = argv[1];
 std::string fasta_dir   = argv[2];
 std::string fasta_file  = argv[3];
 bool skip_soft_clipped  = std::string(argv[4]).compare("true") == 0;
 int32_t max_read_length = atoi(argv[5]);
 std::ofstream output;
 output.open(argv[6], std::ofstream::out);
 compute_confusion_matrix(max_read_length, bam_file, fasta_file, fasta_dir, skip_soft_clipped, output);
 output.close();
}
