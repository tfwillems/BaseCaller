#include "error.h"
#include "io.h"

void CIF::read_intensities(std::istream& input){
  if (data_size_ == 1){
    intensities_8_   = new int8_t[4*num_clusters_*num_cycles_];
    input.read((char*)intensities_8_, 4*num_clusters_*num_cycles_*sizeof(int8_t));
  }
  else if (data_size_ == 2){
    intensities_16_ = new int16_t[4*num_clusters_*num_cycles_];
    input.read((char*)intensities_16_, 4*num_clusters_*num_cycles_*sizeof(int16_t));
  }
  else {
    intensities_32_ = new int32_t[4*num_clusters_*num_cycles_];
    input.read((char*)intensities_32_, 4*num_clusters_*num_cycles_*sizeof(int32_t));
  }
  if (!input)
    printErrorAndDie("Failed to read CIF intensities");
}

CIF::CIF(std::istream& input){
  intensities_8_  = NULL;
  intensities_16_ = NULL;
  intensities_32_ = NULL;
  
  input.read((char*)&version_,      1);
  input.read((char*)&data_size_,    1);
  input.read((char*)&first_cycle_,  2);
  input.read((char*)&num_cycles_,   2);
  input.read((char*)&num_clusters_, 4);
  if (!input)
    printErrorAndDie("Failed to read CIF header");
  if ((1 != data_size_) && (2 != data_size_) && (4 != data_size_))
    printErrorAndDie("Invalid CIF data size: " + std::to_string(data_size_));
  read_intensities(input);
}

inline float CIF::get_value (int cluster, int cycle, int channel) const{
  switch(data_size_){
  case 1:
    return (float)intensities_8_ [cycle*num_clusters_*4 + channel*num_clusters_ + cluster];
  case 2:
    return (float)intensities_16_[cycle*num_clusters_*4 + channel*num_clusters_ + cluster];
  case 4:
    return (float)intensities_32_[cycle*num_clusters_*4 + channel*num_clusters_ + cluster];
  default:
    printErrorAndDie("Logical error in get_value");
    break;
  }
  return -1;
}

std::ostream& operator<< (std::ostream& stream, const CIF& cif){
  stream << "CIF" << "\n"
	 << "\t"  << "Version:     " << static_cast<int16_t>(cif.version_)   << "\n"
	 << "\t"  << "Data size:   " << static_cast<int16_t>(cif.data_size_) << "\n"
	 << "\t"  << "First cycle: " << cif.first_cycle_  << "\n"
	 << "\t"  << "# cycles:    " << cif.num_cycles_   << "\n"
	 << "\t"  << "# clusters:  " << cif.num_clusters_ << std::endl;
  return stream;
}

void CIF::print_intensities(std::ostream& out){
  for (int i = 0; i < num_clusters_; ++i){
    out << i << "\n";
    for (int j = 0; j < 4; ++j){
      out << "\t" << j << "\t";
      for (int k = 0; k < num_cycles_; k++)
	out << get_value(i, k, j) << " ";
      out << "\n";
    }
  }
  out << std::endl;
}


CIF* read_cif(std::string input_file){
  std::ifstream input(input_file, std::ios::binary);
  char format[4];
  input.get(format, 4);
  format[3]='\0';
  if (std::string(format).compare("CIF") != 0)
    printErrorAndDie("CIF file lacking CIF characters in first three bytes");

  CIF* cif = new CIF(input);
  input.close();
  return cif;
}
