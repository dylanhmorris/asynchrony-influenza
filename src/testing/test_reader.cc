#define CATCH_CONFIG_MAIN  
#include "catch.hpp"

#include "../ini/INIReader.h"
#include <string>
#include <iostream>
#include <iomanip>

TEST_CASE( "reader can parse vectors",
           "[parse vector]" ) {

  std::string input_path = "src/testing/TestINI.ini";
  INIReader reader(input_path);
  if (reader.ParseError() < 0) {
    std::cout << "\nCan't load default ini file " <<
      input_path << std::endl;
  };
  std::vector<double> vec =
    reader.GetRealVector("test_section", "vector_vals", {1.0, 25});

  REQUIRE(vec.at(0) == Approx(5.02));
  REQUIRE(vec.at(1) == Approx(3.53223));
  REQUIRE(vec.size() == 5);
}


// int main(int argc, char* argv[]){
//   if(argc < 2){
//     std::cout << "\nUSAGE: ./test_reader <path to ini file>\n" << std::endl;
//     return 0;
//   }
  
//   std::string input_path = argv[1];
//   INIReader reader(input_path);
//   if (reader.ParseError() < 0) {
//     std::cout << "\nCan't load default ini file " <<
//       input_path << std::endl;
//     return 1;
//   };
//   std::vector<double> vec =
//     reader.GetRealVector("test_section", "vector_vals", {1.0, 25});

//   for(auto val : vec)
//     std::cout << val << std::endl;

//   return 0;
// };
