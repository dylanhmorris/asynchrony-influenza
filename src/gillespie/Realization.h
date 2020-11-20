/** 
 *  @file    Realization.h
 *  @author  Dylan Morris (dhmorris@princeton.edu)
 *  
 *  @brief Class Realization holds realizations of gillespie-algorithm-
 *  simulatable rate equation models
 *  
 */

#ifndef REALIZATION_H_
#define REALIZATION_H_
#include "EventRecord.h"
#include "../random/rng.h"
#include "FunctionModel.h"
#include <vector>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <iostream>
#include <fstream>

#include <cfloat>  // for DBL_MIN
#include <cstdint> // for uint64_t
#include <string>
/**  
 *  @brief Class Realization holds gillespie systems with
 *  state arrays consisting of integers
 */  

class Realization {
 public:

  // Default constructor for Realization
  Realization(FunctionModel *model, RNG *rng);

  // Destructor of Realization
  virtual ~Realization();

  const int n_vars;
  const int n_events;
  vector<double> state_array;
  vector<double> rates;
  double state_time;


  
  
  // simulates the realization from current state_time to t_final
  virtual int simulate(double t_final, uint64_t max_iter, bool output);

  /// takes one simulation step according to the chosen algorithm 
  virtual int step() = 0;

  // checks whether all rates are zero
  bool rates_are_zero();

  // prints a header current state of the simulation to cout
  void output_header(std::ostream &outstream);

  // prints the current state of the simulation to cout
  void output_state(const double &time,
                    const vector<double> &state,
                    std::ostream& stream);

  // prints the current state of the simulation to cout
  void output_all_to_file(const string output_path);
  
  // prints the current rates of the simulation to cout
  void print_rates();

  // sets all state array variables to 0
  int zero_vars();

  // sets time to zero
  int zero_time();

  // sets rate array to zero
  int zero_rates();

  // sets time to new values
  void set_time(double new_time);
  
  // sets state to new values
  int set_state(const vector<double> &new_values);

  // completely reset realization to how it appeared on construction
  void reset_all();

  // clear saved states/events/times
  void clear_saved();

  // set a non-default vector to store event records
  void set_event_repository(vector<EventRecord> *event_repo_ptr);

  // set the simulation id (metadata)
  void set_sim_id(int new_sim_id);

  FunctionModel *model;
  RNG *rng;

  vector<vector<double>> array_ts; 
  vector<int> event_ts; 
  vector<double> time_ts; 
  
  vector<EventRecord> saved_events;

  vector<EventRecord> *event_repository;

  int sim_id;

  
};

#endif
