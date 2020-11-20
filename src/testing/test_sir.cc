#define CATCH_CONFIG_MAIN  
#include "catch.hpp"

#include "../multisir/sir.h"
#include "../random/xoroshiro128plus.h"
#include <chrono>

using namespace std::chrono;


TEST_CASE( "Basic SIR model works as expected",
           "[basic SIR]" ) {

  xoroshiro128plus myrng(5);
  
  SIR mysir(3, 0.000002, 0.0, 0.0, &myrng, 1e10);
  
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  mysir.simulate(4, 3000, 20, (int)1e4);
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  high_resolution_clock::rep duration =
    duration_cast<microseconds>( t2 - t1 ).count();
  
  printf("simulation ran in %15.8f seconds \n", duration * 1.0e-6);

  REQUIRE(mysir.state.infected_hosts.size() == 0);
  REQUIRE(mysir.state.living_hosts.size() == 3000);
}
