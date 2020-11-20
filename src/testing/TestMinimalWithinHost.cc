#define CATCH_CONFIG_MAIN  
#include "catch.hpp"

#include "../withinhost/MinimalWithinHost.h"
#include "../withinhost/MinimalParamset.h"
#include "../random/xoroshiro128plus.h"


TEST_CASE( "MinimalWithinHost works as expected" ) {

  xoroshiro128plus model_rng(0);

  
  double beta = 1e-10;
  double r_w = 100;
  double r_m = 100;
  double mu = 0.33e-5;
  double d_v = 4;
  double delta = 6;

  
  double t_M = 2;
  double t_N = 6;

  double C_max = 4e8;
  double c = 0;
  double p_M = 0;
  double f = 1;
  double p_C = 0;
  bool E_w = 0;
  bool E_m = 0;
  
  MinimalWithinHost model_functor(beta,
                                  r_w,
                                  r_m,
                                  mu,
                                  d_v,
                                  delta,
                                  t_M,
                                  t_N,
                                  C_max,
                                  c,
                                  p_M,
                                  f,
                                  p_C,
                                  E_w,
                                  E_m);
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

    REQUIRE(test_rates.at(VW_MINUS) == Approx(d_v));
    REQUIRE(test_rates.at(VM_MINUS) == Approx(0.0));

    // novel response
    test_time = 6.5;
    test_rates = {0, 0, 0, 0, 0, 0};
    model_functor(test_rates, test_state, test_time);
    REQUIRE(test_rates.at(VW_MINUS) == Approx(d_v + delta));
    REQUIRE(test_rates.at(VM_MINUS) == Approx(0));

    // recall response
    test_time = 2.5;
    model_functor.E_w_ = true;
    test_rates = {0, 0, 0, 0, 0, 0};
    model_functor(test_rates, test_state, test_time);
    REQUIRE(test_rates.at(VW_MINUS) == Approx(d_v + delta));
    REQUIRE(test_rates.at(VM_MINUS) == Approx(0));
  }
}

