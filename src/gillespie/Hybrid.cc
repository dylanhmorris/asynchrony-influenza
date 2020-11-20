/** 
 *  @file    Hybrid.cc
 *  @author  Dylan Morris (dhmorris@princeton.edu)
 *  @version 1.0 
 *  
 *  @brief  Class Hybrid implements Realization step() function 
 *          using a janky hybrid of the
 *          random corrections approximate algorithm of Hu and Li
 *          and the SSA (Gillespie, 1977)
 */
 

#include "Hybrid.h"

using namespace std;

/**
 *   @brief  Default constructor for Hybrid
 *  
 *   @param  model a GillespieModel object
 *   @param  rng is a random number generator of class RNG
 *   @param  timestep is the timestep to use for tau leaping
 * 
 *   @return nothing 
 */ 

Hybrid::Hybrid(FunctionModel *model,
               RNG *rng,
               double max_tau) :
  AdaptivePRC(model, rng, max_tau),
  waiting_times_(n_events, 0),
  timestep_(max_tau),
  verbose(false)
{}

/**
 *   @brief Destructor for Hybrid
 *  
 *   @return nothing 
 */ 

Hybrid::~Hybrid(){
}


/**
 *   @brief Take a step
 *  
 *   @return int
 */  

int Hybrid::step(){
  model->update_rates(rates, state_array, state_time);
  update_eta_rates_prod();
  timestep_ = prc_timestep();

  if(use_ssa()){
    ssa_step();
  }
  else
    prc_step();

  array_ts.push_back(state_array);
  time_ts.push_back(state_time);

  return 0;
}

void Hybrid::prc_step(){
  double timestep = timestep_;
  long long n_occurences = 0;
  new_state_ = state_array;
  
  for(int event = 0; event < n_events; event++){
    double mean = rates.at(event) * timestep +
      (0.5 * timestep * timestep * eta_rates_prod.at(event));

    // only call RNG if rate is positive
    if(mean > 0.0){
          n_occurences = rng->rpois(mean);

          // do the event n_occurences times
          if(n_occurences > 0)
            model->rep_event_i(event, n_occurences, new_state_);
    }
    
    else if(mean < -0.0001)
      throw runtime_error("Negative poisson rate geneated "
                          "beyond tolerance for double rounding");
  }

  // very hacky handling of 0 overshoots
  for(int var = 0; var < n_vars; var++){
    if( new_state_.at(var) < +0.0 )
      new_state_.at(var) = +0.0;
  }
  
  
  // update state and time
  state_array = new_state_;
  state_time += timestep;

  array_ts.push_back(state_array);
  time_ts.push_back(state_time);
};


void Hybrid::ssa_step(){
  int min_ind = 0;
  for(int i = 0; i < n_events; i++){
    waiting_times_.at(i) = rng->rexp(rates[i]);
    if(waiting_times_.at(i) < waiting_times_.at(min_ind)){
      min_ind = i;
    }
  }
  // update state and time
  state_time += waiting_times_.at(min_ind);
  model->update_state(min_ind, state_array);
}


int Hybrid::simulate(double t_final, uint64_t max_iter, bool output){

  //update rates
  model->update_rates(rates,
                      state_array,
                      state_time);


  // prep for loop
  uint64_t n_full_steps = max_iter;
  uint64_t i_step = 0;

  // take steps until rates are zero or n_full_steps is reached
  while( (i_step < n_full_steps) && (state_time < t_final) &&
         (!rates_are_zero()) ){

    step();
        
    i_step++;
  }


  if(output){
    for(int k = 0; k < time_ts.size(); k++){
      output_state(time_ts.at(k), array_ts.at(k), std::cout);
    }
  }
  return 0;
};

bool Hybrid::use_ssa(){
  for(int i = 0; i < state_array.size(); i++){
    if((state_array.at(i) < 1) &&
       (state_array.at(i) > 0))
      return true;
  }
  return false;
}
