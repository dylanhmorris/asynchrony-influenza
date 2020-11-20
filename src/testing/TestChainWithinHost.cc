#define CATCH_CONFIG_MAIN  
#include "catch.hpp"

#include "../withinhost/ChainWithinHost.h"
#include "../random/xoroshiro128plus.h"


TEST_CASE( "ChainWithinHost functor works as expected" ) {

  xoroshiro128plus model_rng(0);

  
  double beta = 1e-10;
  double r = 100;
  double mu = 0.33e-5;
  double d_v = 4;
  double delta = 6;

  
  double t_M = 2;
  double t_N = 6;

  double C_max = 4e8;
  double p_M = 0;
  double f = 1;
  double p_C = 0;


  ChainWithinHost model_functor(beta,
                                r,
                                mu,
                                d_v,
                                delta,
                                t_M,
                                t_N,
                                C_max,
                                p_M,
                                f,
                                p_C);
  vector<double> test_state;
  vector<double> test_rates;
  double test_time;
  
  SECTION("Target cells behave as expected"){
    test_state = {C_max * 1.0, 1.0, 0.0};
    test_rates = {0, 0, 0, 0, 0, 0};
    test_time = 0.5;
    model_functor(test_rates, test_state, test_time);

    REQUIRE(test_rates.at(C_PLUS) == Approx(0));
    REQUIRE(test_rates.at(C_MINUS) == Approx(beta * C_max));
  }
  
  SECTION("immunity behaves as expected"){

    // no response
    test_state = {C_max * 1.0, 1.0, 0.0};
    test_rates = {0, 0, 0, 0, 0, 0};
    test_time = 2.5;
    model_functor(test_rates, test_state, test_time);

    REQUIRE(test_rates.at(VMAJ_MINUS) == Approx(d_v));
    REQUIRE(test_rates.at(VMIN_MINUS) == Approx(0.0));

    // homotypic response
    test_time = 6.5;
    test_rates = {0, 0, 0, 0, 0, 0};
    model_functor.set_history({0, 0.5});
    model_functor.set_maj_phen(0.5);
    model_functor.set_min_phen(1.5);
    
    model_functor(test_rates, test_state, test_time);
    REQUIRE(test_rates.at(VMAJ_MINUS) == Approx(d_v + delta));
    REQUIRE(test_rates.at(VMIN_MINUS) == Approx(0));

    // cross-immune response
    test_time = 2.5;
    model_functor.set_history({0, 0.5});
    model_functor.set_maj_phen(0.75);
    model_functor.set_min_phen(10);
    test_rates = {0, 0, 0, 0, 0, 0};
    model_functor(test_rates, test_state, test_time);
    REQUIRE(test_rates.at(VMAJ_MINUS) == Approx(d_v + 0.75 * delta));
    REQUIRE(test_rates.at(VMIN_MINUS) == Approx(0));
  }
}

TEST_CASE( "ChainWithinHost mucosal sampling works as expected" ) {

  xoroshiro128plus model_rng(0);

  
  double beta = 1e-10;
  double r = 100;
  double mu = 0.33e-5;
  double d_v = 4;
  double delta = 6;

  
  double t_M = 2;
  double t_N = 6;

  double C_max = 4e8;
  double p_M = 0;
  double f = 1;
  double p_C = 0;


  ChainWithinHost model_functor(beta,
                                r,
                                mu,
                                d_v,
                                delta,
                                t_M,
                                t_N,
                                C_max,
                                p_M,
                                f,
                                p_C);
  vector<double> test_state;
  vector<double> test_rates;
  double test_time;

  model_functor.set_history({1});
  
  SECTION("no mutant if no mutant in donor host"){
    vector<double> result = model_functor.get_inoculum(0,
                                                       1000.0,
                                                       500.0,
                                                       100.0,
                                                       101.0,
                                                       model_rng,
                                                       false);
    REQUIRE(result.at(0) == Approx(500));
    REQUIRE(result.at(1) == Approx(0.0));
  }
    
  
  SECTION("inoculuation selection works"){
    vector<double> result = model_functor.get_inoculum(0.008,
                                                       1000.0,
                                                       4.0,
                                                       1.0,
                                                       5.0,
                                                       model_rng,
                                                       true);
    REQUIRE(result.at(1) == Approx(2));
    REQUIRE(result.at(0) == Approx(0.0));

    
    result = model_functor.get_inoculum(0.008,
                                        1000.0,
                                        4.0,
                                        3.0,
                                        5.0,
                                        model_rng,
                                        true);
    REQUIRE(result.at(1) == Approx(0.0));
    REQUIRE(result.at(0) == Approx(4));
  }
  
  SECTION("hypergeometric sampling works"){

    vector<double> result = model_functor.get_inoculum(0.5,
                                                       1000.0,
                                                       4.0,
                                                       1.0,
                                                       5.0,
                                                       model_rng,
                                                       true);
    REQUIRE(result.at(1) + result.at(0) == Approx(4));
    result = model_functor.get_inoculum(0.5,
                                        1000.0,
                                        4.0,
                                        1.0,
                                        5.0,
                                        model_rng,
                                        true);
    REQUIRE(result.at(1) + result.at(0) == Approx(4));
  }
}

