#ifndef IO_H_
#define IO_H_

#include <fstream>
#include <iostream>
#include <string>

class CIF {
private:
  int8_t version_, data_size_;
  int16_t first_cycle_, num_cycles_;
  int32_t num_clusters_;

  // Intensity array for each possible CIF encoding
  int8_t*  intensities_8_;
  int16_t* intensities_16_;
  int32_t* intensities_32_;

  void read_intensities(std::istream& input);

public:
  CIF(std::istream& input);

  ~CIF(){
    delete [] intensities_8_;
    delete [] intensities_16_;
    delete [] intensities_32_;
  }

  float get_value (int cluster, int cycle, int channel) const;

  friend std::ostream& operator<< (std::ostream& stream, const CIF& cif);

  void print_intensities(std::ostream& out);
};

CIF* read_cif(std::string input_file);

#endif
