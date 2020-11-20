/** 
 *  @file    Realization.cc
 *  @author  Dylan Morris (dhmorris@princeton.edu)
 *  @date    2017-02-27  
 *  @version 0.1 
 *  
 *  @brief Class Realization holds realizations of gillespie-algorithm-
 *  simulatable rate equation models
 *  
 */

#include "Realization.h"
#include <stdint.h>

using std::vector;
using std::ofstream;

/**
 *   @brief  Default constructor for Realization  
 *
 *   @param model is a FunctionModel object defining the 
 *   stochastic system  
 *   
 *   @param  rng is a random number generator of class RNG
 *
 *   @return nothing 
 */ 
Realization::Realization(FunctionModel *model, RNG *rng) :
  model(model),
  rng(rng),
  n_vars(model->n_vars()),
  n_events(model->n_events()),
  state_array(n_vars, 0.0),
  rates(n_events, 0.0)
{
  reset_all();
  event_repository = &saved_events;
}

/**
 *   @brief  Destructor of Realization
 *  
 *   @return nothing 
 */ 
Realization::~Realization(){
}

/**
 *   @brief Sets state_array and state_time to new values
 *  
 *   @return int
 */ 


/**
 *   @brief Simulates the system from t_inital to t_final
 *  
 *   @return int
 */ 
int Realization::simulate(double t_final, uint64_t max_iter,
                          bool output = false){
  
  uint64_t iter_count = 0;
  bool done = 0;
  array_ts.push_back(state_array);
  time_ts.push_back(state_time);

  // initialize rate array
  model->update_rates(rates, state_array, state_time);

  /* check that stop conditions haven't already been reached
     (maybe wrap this in a function) */ 
    if(state_time > t_final){
      done = 1;
    }
    else if (iter_count > max_iter) {
      done = 1;
      throw runtime_error("Max iterations exceeded");
    }
    else if(rates_are_zero()){
      done = 1;
    }

  while(done==0){

    // print the rates (for debugging)
    //    print_rates();
    
    // take step according to method
    step();

    // increment iteration count
    iter_count++;
    
    /* check that stop conditions haven't been reached
       (maybe wrap this in a function) */ 
    if(state_time > t_final){
      done = 1;
    }
    else if (iter_count > max_iter) {
      done = 1;
      throw runtime_error("Max iterations exceeded");
    }
    else if(rates_are_zero()){
      done = 1;
    }

    if(output)
      output_state(state_time, state_array, std::cout);

  }

  if(output)
    output_state(state_time, state_array, std::cout);

  return 0;
}

/**
 *   @brief  Prints a header for the simulation output to
 *   a designated oustream
 *  
 *   @return nothing 
 */ 
void Realization::output_header(std::ostream &outstream){
  outstream << left << setw(15) << "time";

  for(int k = 0; k < n_vars; k++){
    outstream << left << setw(11) << model->get_var_name(k);
  }
  
  outstream << "\n";
}

/**
 *   @brief  Prints the current state of the simulation
 *  
 *   @return int 
 */ 

void Realization::output_state(const double &time,
                               const vector<double> &state,
                               std::ostream& stream){
  // to be modified depending on ultimately
  // chosen output format

  stream << left << setprecision(8) << setw(15) << time;
  for(int i = 0; i < state.size(); i++){
    stream << left << setprecision(5) << setw(11) << state.at(i);
  }
  stream << "\n";
}

/**
 *   @brief  Prints all saved states of the simulation
 *   to a specified file
 *  
 *   @return int 
 */ 

void Realization::output_all_to_file(const string output_path){
  ofstream outstream;  
  outstream.open(output_path);
  output_header(outstream);
  for(int i_state = 0; i_state < array_ts.size(); i_state++){
    output_state(time_ts.at(i_state),
                 array_ts.at(i_state),
                 outstream);
  }
  outstream.close();
}




/**
 *   @brief  Prints the current rates to cout
 *  
 *   @return int 
 */ 

void Realization::print_rates(){
  cout << "\n";

  cout << "Rates: ";

  cout << left << setprecision(5) << setw(15) << state_time;
  for(int i = 0; i < n_events; i++){
    cout << left << setprecision(5) << setw(15) << rates[i];
  }
  cout << "\n";
  cout << "\n";


}



/**
 *   @brief checks whether all rates are zero
 *
 *   @return bool
 */
bool Realization::rates_are_zero(){
  bool val = 1;    
  
  for (int i = 0; i < n_events; i++)    
    if(rates[i] > DBL_MIN)
      val = 0;

  return val;
}

int Realization::zero_time(){
  state_time = 0.0;
  return 0;
}

int Realization::zero_vars(){
  for(int k = 0; k < n_vars; k++)
    state_array[k] = 0;
  return 0;
}

int Realization::zero_rates(){
  for(int k = 0; k < n_events; k++)
    rates[k] = 0;
  return 0;
}


int Realization::set_state(const vector<double> &new_values){
  for(int k = 0; k < n_vars; k++)
    state_array[k] = new_values[k];
  
  return 0;
}

void Realization::set_time(double new_time){
    state_time = new_time;
}


// completely reset realization to how it appeared on construction
void Realization::reset_all(){
  zero_vars();
  zero_rates();
  zero_time();
  clear_saved();
  sim_id = -1;
}


// clear saved states/events/times
void Realization::clear_saved(){
  saved_events.clear();
  time_ts.clear();
  array_ts.clear();
  event_ts.clear();
}


void Realization::set_event_repository(vector<EventRecord> *event_repo_ptr){
  event_repository = event_repo_ptr;
};


void Realization::set_sim_id(int new_sim_id){
  sim_id = new_sim_id;
};
