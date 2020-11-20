#ifndef GILLESPIE_RATE_FUNCTOR_H_
#define GILLESPIE_RATE_FUNCTOR_H_

#include <vector>
#include <stdexcept>
#include "EventRecord.h"

using namespace std;

class GillespieRateFunctor {
 public:

  virtual void operator()(vector<double> &rate_array,
                          const vector<double> &state_array,
                          double time) = 0;

  // do logic to record event
  virtual void record(const vector<double> &old_state,
                      const vector<double> &new_state,
                      const double new_time,
                      const int sim_id,
                      vector<EventRecord> &record_vec);
  
  // set the event_deltas_ value for the i_event-th
  // event and j_var-th variable
  void set_delta_var(int i_event, int j_var, double new_delta);


  // get the event_deltas_ value for the i_event-th
  // event and j_var-th variable
  double get_delta(int i_event, int j_var);

  // get the number of events
  int n_events();

  // get the number of variables
  int n_vars();
  
  /* convenience function setting event_deltas array
   * to the simple +1/-1 format for models of n variables
   * in which each variable x_i can be changed by one of  two 
   * possible events: x_i += 1 and x_i -= 1, and no event
   * changes more than one variable
   */
  void ones_event_matrix(int var_count);

  // array of arrays; ith array specifies
  // how much the ith event changes the state
 protected:
  vector<vector<double>> event_deltas_;
};

#endif

