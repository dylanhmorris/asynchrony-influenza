#define CATCH_CONFIG_MAIN  
#include "catch.hpp"

#include "../withinhost/MinimalWithinHost.h"
#include "../gillespie/midpointleap.h"
#include "../gillespie/firstreaction.h"
#include "../gillespie/AdaptivePRC.h"
#include "../gillespie/FunctionModel.h"
#include "../random/xoroshiro128plus.h"

class HomeworkFunctor : public GillespieRateFunctor{
public:
  HomeworkFunctor(){
    int n_var = 2;
    // initialize event matrix
    for(int k = 0; k < (2 * n_var); k++)
      event_deltas_.push_back(vector<double>(n_var, 0));

    event_deltas_[0][0] = 1;
    event_deltas_[1][0] = -1;
    event_deltas_[2][1] = 1;
    event_deltas_[3][1] = -1;
    
  };

  ~HomeworkFunctor(){};

  void operator()(vector<double> &rate_array,
                  const vector<double> &state_array,
                  double time) {
    rate_array[0] = 3.0;
    rate_array[1] = state_array[0] / 3.0;
    rate_array[2] = state_array[0] * 5.0;
    rate_array[3] = state_array[1] / 30.0;
  }

};



TEST_CASE( "Model can be initialized and runs as expected" ) {
  HomeworkFunctor hwfunctor;
  
  FunctionModel testmodel(hwfunctor.n_vars(), &hwfunctor);
  testmodel.set_var_name(0, "a");
  testmodel.set_var_name(1, "b");

  REQUIRE(testmodel.n_events() == 4);
  REQUIRE(testmodel.n_vars() == 2);

  xoroshiro128plus test_rng(502);
  xoroshiro128plus test_rng2(502);
  xoroshiro128plus test_rng3(502);
  FirstReaction teststate(&testmodel, &test_rng);
  MidpointLeap teststate2(&testmodel, &test_rng2, 0.005);
  MidpointLeap teststate3(&testmodel, &test_rng3, 0.005);
  AdaptivePRC teststate4(&testmodel, &test_rng3, 0.1);

  SECTION("Model runs as expected given FirstReaction and RNG"){
    teststate.simulate(5000, 500000, false);
    REQUIRE(teststate.state_array.at(0) == 4);
    REQUIRE(teststate.state_array.at(1) == 1400);
    REQUIRE(teststate.state_time == Approx(5000.002));
    printf("BREAKPOINT \n\n");

  }

  SECTION("FirstReaction can be re-zeroed"){
    teststate2.reset_all();
    REQUIRE(teststate2.state_array[0] == 0);
    REQUIRE(teststate2.state_array[1] == 0);
    REQUIRE(teststate2.state_time == 0.0);
    printf("BREAKPOINT \n\n");

  }
  
  SECTION("MidpointLeap simulates as expected"){
    printf("simulating with midpoint leap...\n");
    teststate2.reset_all();
    teststate2.simulate(5000.321, 1200000, false);
    REQUIRE(teststate2.state_array[0] == 11);
    REQUIRE(teststate2.state_array[1] == 1525);
    REQUIRE(teststate2.state_time == 5000.321);
  }
  
  SECTION("Simulation picks up where it left off"
          "and goes to exact tfinal, and prints"){
    teststate3.reset_all();
    teststate3.simulate(5000, 5000000, false);
    teststate3.simulate(5000.321, 1200000, false);
    REQUIRE(teststate3.state_time == 5000.321);
    REQUIRE(teststate3.state_array[0] == 11);
    REQUIRE(teststate3.state_array[1] == 1525);
 
  }

  SECTION("MidpointLeap can be properly re-zeroed"){
    printf("re-zeroing...\n");
    teststate3.reset_all();
    REQUIRE(teststate3.state_array[0] == 0);
    REQUIRE(teststate3.state_array[1] == 0);
    REQUIRE(teststate3.state_time == 0.0);

  }

  SECTION("AdaptivePRC runs"){
    printf("re-zeroing...\n");
    teststate4.reset_all();
    REQUIRE(teststate4.state_array[0] == 0);
    REQUIRE(teststate4.state_array[1] == 0);
    REQUIRE(teststate4.state_time == 0.0);
    teststate4.simulate(5000, 5000000, false);
    REQUIRE(teststate4.state_array[0] == 7);
    REQUIRE(teststate4.state_array[1] == 1294);

  }

}

