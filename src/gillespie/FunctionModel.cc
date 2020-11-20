/** 
 *  @file    FunctionModel.cc
 *  @author  Dylan Morris (dhmorris@princeton.edu)
 *  @date    2017-03-01  
 *  @version 0.1
 *  
 *  @brief Class FunctionModel, which holds user-specified models 
 *  of stochastic systems from which realizations are to
 *  be simulated.
 *  
 */



#include "FunctionModel.h"

using namespace std;

/**
 *   Default constructor for FunctionModel  
 *  
 *   @return nothing 
 */ 

FunctionModel::FunctionModel(int n_vars, GillespieRateFunctor *rate_fn) :
  n_vars_(n_vars),
  rate_fn_(rate_fn)
{
  // default variable names are variable indecies
  for(int k = 0; k < n_vars_; k++)
    var_names_.push_back(to_string(k));
}

/**
 *   Destructor of FunctionModel  
 *  
 *   @return nothing 
 */ 
FunctionModel::~FunctionModel() {}


/**
 *   @brief Returns number of events in model
 *  
 *   @return int number of events
 */ 

int FunctionModel::n_events(){
  return rate_fn_->n_events();
}

int FunctionModel::n_vars(){
  return n_vars_;
}




/**
 *   @brief  Updates a given state_array by performing i_event-th
 *   even in event_list_
 *
 *   @param  i_event is an int specifying which event to perform
 *
 *   @param  state_array is a double array to alter according to
 *   the event 
 *
 *   @param  time is a double representing the current 
 *   system time 
 * 
 *   @param params is a double array representing a set of
 *   values for variable model parameters

 *   @return nothing
 */
void FunctionModel::update_state(int i_event,
                                 vector<double> &state_array) {
  // increment variables by the event's
  // corresponding event_deltas_ values
  for (int k_var = 0; k_var < n_vars(); k_var++)
    state_array[k_var] += rate_fn_->get_delta(i_event, k_var);
}

/**
 *   @brief  Update rate for all Events in Model's Event list
 *  
 *   @param  state_array is an double array specifiying variable 
 *   values of a function
 * 
 *   @param  rate_array is a double array hold the
 *   new rate values
 *
 *   @param  time is a double representing the current 
 *   system time 
 * 
 *   @param params is a double array representing a set of
 *   values for variable model parameters
 *
 *   @return void
 */ 
void FunctionModel::update_rates(vector<double> &rate_array,
                                 const vector<double> &state_array,
                                 double time) {
  (*rate_fn_)(rate_array, state_array, time);

  for(int i = 0; i < n_events(); i++){
  // avoid issues with rounding of very, very small rates
    if(rate_array.at(i) < 0){
      printf("Negative rate: event %d, rate %15.8f \n",
             i, rate_array.at(i));

      for(int k = 0; k < state_array.size(); k++)
        printf("%28.8f", state_array[k]);
      
      printf("\n");
      
      throw runtime_error("Negative rate generated");

    }
    else if(rate_array.at(i) < 1.0e-50)
      rate_array.at(i) = 0.0;
  }

  
}

void FunctionModel::update_eta_rates(vector<double> &rate_array,
                                     const vector<double> &state_array,
                                     double time) {
  (*rate_fn_)(rate_array, state_array, time);
}


// calculate approximate midpoints for tau leaping
void FunctionModel::calc_midpoints(vector<double> &midpoint_array,
                                   const vector<double> &curr_rates,
                                   const vector<double> &state_array,
                                   double time,
                                   double timestep_size){

  int event_count = n_events();
  
  // calc midpoint for each variable
  for(int ivar = 0; ivar < n_vars(); ivar++){
    // calculate derivative
    double deriv = 0;
    for (int i_event = 0; i_event < event_count; i_event++)
      deriv += (curr_rates[i_event] *
                rate_fn_->get_delta(i_event, ivar));

    // use derivative to calculate midpoint
    midpoint_array[ivar] = state_array[ivar] + 0.5 * timestep_size * deriv;
  }
}


// Determine whether the currently occurring event should
// be logged (default always false)
void FunctionModel::record(const vector<double> &old_state,
                           const vector<double> &new_state,
                           const double new_time,
                           const int sim_id,
                           vector<EventRecord> &record_vec)
{
  return rate_fn_->record(old_state,
                          new_state,
                          new_time,
                          sim_id,
                          record_vec);
}


/**
 *   @brief  Updates a given state_array by performing an
 *   event a certain number of times
 *
 *   @param  i_event is an int specifying which event to perform
 *
 *   @param  state_array is a double array to alter repeatedly
 *   according to the event 
 * 
 *   @return nothing
 */
void FunctionModel::rep_event_i(int i_event, long long n_reps,
                                 vector<double> &state_array) {
  // increment variables by the event's
  // corresponding event_deltas_ values
  for (int k_var = 0; k_var < n_vars(); k_var++)
    state_array[k_var] += n_reps * rate_fn_->get_delta(i_event, k_var);
}

// Set the name of the ith variable
void FunctionModel::set_var_name(int i, string name){
  var_names_.at(i) = name;
}


// Returns the name of the specified variable
string FunctionModel::get_var_name(int which_var){
  return var_names_.at(which_var);
}

double FunctionModel::get_delta(int i_event, int k_var){
  return rate_fn_->get_delta(i_event, k_var);
}
