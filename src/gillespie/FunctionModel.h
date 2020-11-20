/** 
 *  @file    FunctionModel.h
 *  @author  Dylan Morris (dhmorris@princeton.edu)
 *  @date    2017-03-01  
 *  @version 0.1 
 *  
 *  @brief Class FunctionModel, which holds user-specified models 
 *  of stochastic systems from which realizations are to
 *  be simulated, where rates are specified by a single function
 *  
 */

#ifndef FUNCTION_MODEL_H_
#define FUNCTION_MODEL_H_
#include <vector>
#include <stdint.h>
#include <string>
#include <cmath>
#include "EventRecord.h"
#include "GillespieRateFunctor.h"

using namespace std;

/**  
 *  @brief Class GillespieModel, which holds user-specified models 
 *  of stochastic systems from which realizations are to
 *  be simulated, where rates are specified by a single function. 
 *  A model may have variable parameters; these
 *  should be passed as a vector of doubles to a
 *  realization of the model of class GillespieInt
 */  

class FunctionModel {

 public:

  // Default constructor for Model
  FunctionModel(int n_vars, GillespieRateFunctor *rate_fn);        

  // Destructor of Model 
  ~FunctionModel();

  // get the number of events
  int n_events();

  // get the number of variables
  int n_vars();
    
  // Update rate for all Events in GillespieModel's Event list
  void update_rates(vector<double> &rate_array,
                    const vector<double> &state_array,
                    double time);

  void  update_eta_rates(vector<double> &rate_array,
                         const vector<double> &state_array,
                         double time);
  
  // update a given state_array by performing the ith event
  // in the event list
  void update_state(int i_event,
                    vector<double> &state_array);

  // Update a given state_array by performing an
  // event a certain number of times
  void rep_event_i(int i_event, long long n_reps,
                   vector<double> &state_array);

  // calculate approximate midpoints for tau leaping
  void calc_midpoints(vector<double> &midpoint_array,
                      const vector<double> &curr_rates,
                      const vector<double> &state_array,
                      double time,
                      double timestep_size);

  // set how the ith event changes the state
  void set_delta_var(int i_event, int i_var, int new_delta);

  //possibly record most recent event
  void record(const vector<double> &old_state,
              const vector<double> &new_state,
              const double new_time,
              const int sim_id,
              vector<EventRecord> &record_vec);

  // Set the name of the ith variable
  void set_var_name(int i, string name);
  
  // Returns the name of the specified variable
  string get_var_name(int which_var);

  // Returns the magnitude of state array
  // change for the ith event and kth variable
  double get_delta(int i_event, int k_var);
 private:
  int n_vars_;
  vector<string> var_names_;
  GillespieRateFunctor *rate_fn_;
};

#endif
