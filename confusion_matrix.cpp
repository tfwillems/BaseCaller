#include "confusion_matrix.h"

#include <iostream>
#include "bamtools/include/api/BamAlignment.h"

#include "error.h"

inline char get_base(int base_index){
  switch(base_index){
  case 0:
    return 'A';
  case 1:
    return 'C';
  case 2:
    return 'G';
  case 3:
    return 'T';
  case 4:
    return 'N';
  default:
    printErrorAndDie("Invalid base index " + std::to_string(base_index));
  }
  return '0';
}

inline int get_base_index(char base){
  switch (base){
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'T':
    return 3;
  case 'N':
    return 4;
  default:
    printErrorAndDie("Invalid base " + std::to_string(base));
  }
  return -1;
} 

bool is_hard_clipped(BamTools::BamAlignment& aln){
  for (auto iter = aln.CigarData.begin(); iter != aln.CigarData.end(); iter++)
    if (iter->Type == 'H')
      return true;
  return false;
}

void walk_alignment_forward(BamTools::BamAlignment& aln, int32_t max_read_length, std::string& ref_seq, int32_t* matrix_counts, int32_t* total_counts){
  auto aln_seq_iter = aln.QueryBases.begin();
  auto ref_seq_iter = ref_seq.begin() + aln.Position; 
  int32_t read_pos  = 0;
  for (auto iter = aln.CigarData.begin(); iter != aln.CigarData.end(); iter++){
    switch(iter->Type){
    case 'M': case 'S':
      for (int j = 0; j < iter->Length; j++){
	if (read_pos >= max_read_length)
	  return;
	int ref_base_index  = get_base_index(*ref_seq_iter);
	int emit_base_index = get_base_index(*aln_seq_iter);
	total_counts[ref_base_index]++;
	matrix_counts[5*ref_base_index + emit_base_index]++;
	aln_seq_iter++;
	ref_seq_iter++;
	total_counts  += 5;
	matrix_counts += 25;
	read_pos++;
      }
      break;
    case 'I':
      aln_seq_iter  += iter->Length;
      total_counts  += 5*iter->Length;
      matrix_counts += 25*iter->Length;
      read_pos      += iter->Length;
      break;
    case 'D':
      ref_seq_iter += iter->Length;
      break;
    default:
      printErrorAndDie("Invalid CIGAR char " + std::to_string(iter->Type));
      break;
    }
  }
}

void walk_alignment_reverse(BamTools::BamAlignment& aln, int32_t max_read_length, std::string& ref_seq, int32_t* matrix_counts, int32_t* total_counts){
  auto aln_seq_iter = aln.QueryBases.begin() + aln.QueryBases.size()-1;
  auto ref_seq_iter = ref_seq.begin() + aln.GetEndPosition();
  int32_t read_pos  = 0;
  for (auto iter = aln.CigarData.rbegin(); iter != aln.CigarData.rend(); iter++){
    switch(iter->Type){
    case 'M': case 'S':
      for (int j = 0; j < iter->Length; j++){
	if (read_pos >= max_read_length)
	  return;
	int ref_base_index  = get_base_index(*ref_seq_iter);
	int emit_base_index = get_base_index(*aln_seq_iter);
	total_counts[ref_base_index]++;
	matrix_counts[5*ref_base_index + emit_base_index]++;
	aln_seq_iter--;
	ref_seq_iter--;
	total_counts  += 5;
	matrix_counts += 25;
	read_pos++;
      }
      break;
    case 'I':
      aln_seq_iter  -= iter->Length;
      total_counts  += 5*iter->Length;
      matrix_counts += 25*iter->Length;
      read_pos      += iter->Length;
      break;
    case 'D':
      ref_seq_iter -= iter->Length;
      break;
    default:
      printErrorAndDie("Invalid CIGAR char " + std::to_string(iter->Type));
      break;
    }
  }
}

void process_reads(BamTools::BamReader& bam_reader, int32_t max_read_length, std::string& ref_seq,
		   int32_t* matrix_counts, int32_t* total_counts){
  BamTools::BamAlignment aln;
  int32_t hard_clip_count      = 0;
  int32_t too_close_to_start   = 0;
  int32_t too_close_to_end     = 0;
  int32_t read_count           = 0;
  int32_t forward = 0, reverse = 0;
  while (bam_reader.GetNextAlignment(aln)){
    if (++read_count % 1000 == 0)
      std::cerr << "Processing read # " << read_count << std::endl;

    if (is_hard_clipped(aln)){
      hard_clip_count++;
      continue;
    }
    if (aln.Position < 100){
      too_close_to_start++;
      continue;
    } 
    if (aln.Position+max_read_length > ref_seq.size()-100){
      too_close_to_end++;
      continue;
    }
  
    if (aln.IsReverseStrand()){
      reverse++;
      walk_alignment_reverse(aln, max_read_length, ref_seq, matrix_counts, total_counts);
    }
    else {
      forward++;
      walk_alignment_forward(aln, max_read_length, ref_seq, matrix_counts, total_counts);
    }
  }
}


void print_confusion_matrix(int32_t* matrix_counts, int32_t* total_counts, int32_t max_read_length){
  int32_t* mat_count_ptr;
  int32_t* tot_count_ptr;

  for (int ref_base = 0; ref_base < 4; ref_base++){
    std::cout << get_base(ref_base) << "\n";
    for (int emit_base = 0; emit_base < 5; emit_base++){
      tot_count_ptr = total_counts  + ref_base;
      mat_count_ptr = matrix_counts + 5*ref_base + emit_base;
      std::cout << "\t" << "-> " << get_base(emit_base);
      for (int cycle = 0; cycle < max_read_length; cycle++){
	std::cout << " " << 1.0*(*mat_count_ptr)/(*tot_count_ptr);
	tot_count_ptr += 5;
	mat_count_ptr += 25;
      }
      std::cout << std::endl;
    }
  }
}

