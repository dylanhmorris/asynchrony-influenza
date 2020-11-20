#define CATCH_CONFIG_MAIN  
#include "catch.hpp"

#include "../withinhost/MinimalParamset.h"
#include "../random/xoroshiro128plus.h"


TEST_CASE( "MinimalParamset class" ) {

  string input_path = "src/testing/TestParams.mk";

  xoroshiro128plus model_rng(0);

  SECTION("Can be initialized")
    MinimalParamset test_paramset(input_path,
                                  "SUCCEED");
  
  SECTION("Throws error if given bad values"){
    MinimalParamset test_paramset(input_path,
                                  "SUCCEED");
    //test_paramset.print_paramset(std::cout);
    REQUIRE_THROWS(test_paramset.parse_input("FAIL"));
  }

  SECTION("Hypercube random paramsets work as expected"){
    MinimalParamset test_params(input_path,
                                "SUCCEED");

    vector<string> hypercube_params = {"R_W", "VB_RATIO"};

    // do this a bunch of times to make sure
    // not a seed quirk
    int n_hyper_iter = 5;
    
    for(int j = 0; j < n_hyper_iter; j++){
      int n_hyper = 1000;

      vector<double> upper_half = {};
      vector<double> lower_half = {};
    
      test_params.setup_hypercube(hypercube_params,
                                  n_hyper,
                                  model_rng);
      double r_w;
      int lowest_count = 0;

      for(int k = 0; k < n_hyper; k++){
        test_params.hypercube_next();
        r_w = test_params.get_parameter("R_W");
        if(r_w > 6)
          upper_half.push_back(r_w);
        else
          lower_half.push_back(r_w);
        if(r_w < 2)
          lowest_count++;
      }

      // check that we have a true hypercube
      REQUIRE(upper_half.size() == 500);
      REQUIRE(lower_half.size() == 500);
      REQUIRE(lowest_count == 100);
    }
  }
}
