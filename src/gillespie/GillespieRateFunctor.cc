#include "GillespieRateFunctor.h"

using namespace std;

// get the event_deltas_ value for the i_event-th
// event and j_var-th variable
double GillespieRateFunctor::get_delta(int i_event, int j_var){
  return event_deltas_.at(i_event).at(j_var);
}

// get the number of events
int GillespieRateFunctor::n_events(){
  return event_deltas_.size();
}

// get the number of variables
int GillespieRateFunctor::n_vars(){
  return event_deltas_[0].size();
}
  


// set the delta for the i_event-th event and j_var-th
// variable
void GillespieRateFunctor::set_delta_var(int i_event,
                                         int j_var,
                                         double new_delta){
  event_deltas_[i_event][j_var] = new_delta;
};

void GillespieRateFunctor::ones_event_matrix(int var_count){

  if(var_count != n_vars())
    throw runtime_error("number of variables passed to ones_event_matrix "
                        "not equal to number of variables inferred "
                        "from delta_vars_ vector size");

  if(2 * var_count != n_events())
    throw runtime_error("number of variables passed to ones_event_matrix "
                        "clashes with number of events given by "
                        "delta_vars_ vector size.");

  
  for(int i_event = 0; i_event < var_count; i_event++)
    event_deltas_[i_event][i_event] = 1;
  
  for(int i_event = var_count; i_event < 2 * var_count; i_event++)
    event_deltas_[i_event][i_event - var_count] = -1;
  
};


// Determine whether the currently occurring event should
// be logged (default always false)
void GillespieRateFunctor::record(const vector<double> &old_state,
                                  const vector<double> &new_state,
                                  const double new_time,
                                  const int sim_id,
                                  vector<EventRecord> &record_vec){}
