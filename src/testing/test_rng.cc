#define CATCH_CONFIG_MAIN  
#include "catch.hpp"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <math.h>

#include "../random/xoroshiro128plus.h"
#include "../random/rng.h"

TEST_CASE( "xoroshiro128+ RNG can be instatiated and deleted" ) {
  xoroshiro128plus my_rng(1);
}

TEST_CASE( "test randchoice functionality for choosing from"
           "vectors of various types") {
  xoroshiro128plus my_rng(1);

  vector<int> intvec;
  vector<double> doublevec;

  intvec.push_back(1);
  intvec.push_back(1);

  doublevec.push_back(2);
  doublevec.push_back(2.5);

  int intchoice = my_rng.randchoice(intvec);
  REQUIRE(intchoice == 1);

  double doublechoice = my_rng.randchoice(doublevec);
  REQUIRE(doublechoice == 2.5);
  my_rng.runif();
  doublechoice = my_rng.randchoice(doublevec);
  REQUIRE(doublechoice == 2.0);    
}

TEST_CASE( "trying to draw from an empty vector throws error") {
  xoroshiro128plus my_rng(1);

  vector<int> emptyvec;
  emptyvec.clear();

  REQUIRE_THROWS(my_rng.randchoice(emptyvec));

}

TEST_CASE( "check that seed resets work") {
  int init_seed = 1;
  
  xoroshiro128plus my_rng(init_seed);

  uint64_t my_first_int = my_rng.next();

  my_rng.jump();

  uint64_t my_second_int = my_rng.next();

  REQUIRE(my_first_int != my_second_int);

  my_rng.set_seed(init_seed);
  uint64_t my_third_int = my_rng.next();
  REQUIRE(my_first_int == my_third_int);
}


TEST_CASE( "check that weighted_choice works") {
  int init_seed = 1;
  
  xoroshiro128plus my_rng(init_seed);

  vector<double> weights_one = {1.0, 0.0, 0.0};
  vector<int> choices;
  
  for(int k = 0; k < 100; k++){
    choices.push_back(my_rng.weighted_choice(weights_one));
  }

  REQUIRE(choices.at(0) == 0);
  REQUIRE(std::equal(choices.begin() + 1, choices.end(), choices.begin()));

  vector<double> weights_two = {0, 0.5, 0.0};
  choices.clear();
  for(int k = 0; k < 100; k++){
    choices.push_back(my_rng.weighted_choice(weights_two));
  }

  REQUIRE(choices.at(0) == 1);
  REQUIRE(std::equal(choices.begin() + 1, choices.end(), choices.begin()));

  vector<double> weights_three = {0, 0.5, 0.5};
  choices.clear();
  int choice;
  for(int k = 0; k < 5; k++){
    choice = my_rng.weighted_choice(weights_three);
    REQUIRE((choice == 1 || choice == 2));
  }
}

TEST_CASE( "test poisson with very large numbers") {
  int init_seed = 1;
  
  xoroshiro128plus my_rng(init_seed);

  long long big_mean = 1e15;
  long long choice = 0;
  
  for(int k = 0; k < 5; k++){
    choice = my_rng.rpois(big_mean);
    REQUIRE(log(choice * 1.0) - log(big_mean) < 1);
  }
}

TEST_CASE( "test latin hypercube truly is") {
  int init_seed = 1;
  
  xoroshiro128plus my_rng(init_seed);
  int n_cutpoints = 1000;
  vector<double> sample = my_rng.runif_latin(n_cutpoints);

  for(int k = 1; k < n_cutpoints + 1; k++){
    int sum_less = 0;
    for(auto &x : sample){
      if(x < 1.0 * k / n_cutpoints)
        sum_less++;
    }
    REQUIRE(sum_less == k);
  }
}
