#include <iostream>

#include "io.h"

int main(int argc, char* argv[]) {
 CIF* cif = read_cif(std::string(argv[1]));
 std::cout << *cif;
 delete cif;
}
