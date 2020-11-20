/** 
 *  @file    firstreaction.cc
 *  @author  Dylan Morris (dhmorris@princeton.edu)
 *  @version 1.0 
 *  
 *  @brief   Class FirstReaction implements Realization step() function
 *           using the exact First Reaction algorithm of Gillespie (1971)
 *  
 */

#include "firstreaction.h"

/**
 *   @brief  Default constructor for FirstReaction
 *  
 *   @param  model is a GillespieModel object
 *   @param  rng is a random number generator of class RNG
 * 
 *   @return nothing 
 */ 
FirstReaction::FirstReaction(FunctionModel *model,
                             RNG *rng) :
  Realization(model, rng),
  waiting_times(n_events, 0.0)
{}

/**
 *   @brief Destructor for FirstReaction
 *  
 *   @return nothing 
 */ 
FirstReaction::~FirstReaction(){
}


/**
 *   @brief Update waiting times
 *  
 *   @return int
 */  
int FirstReaction::step(){

  int min_ind = 0;
  int i;
  for(i = 0; i < n_events; i++){
    waiting_times[i] = rng->rexp(rates[i]);
    if(waiting_times[i] < waiting_times[min_ind]){
      min_ind = i;
    }
  }
  
  // update state, time, rates


  state_time += waiting_times[min_ind];
  model->update_state(min_ind, state_array);
  model->update_rates(rates, state_array, state_time);
  array_ts.push_back(state_array);
  time_ts.push_back(state_time);
  event_ts.push_back(min_ind);

  return 0;
}

