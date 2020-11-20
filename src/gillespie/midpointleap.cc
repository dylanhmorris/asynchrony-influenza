/** 
 *  @file    midpointleap.cc
 *  @author  Dylan Morris (dhmorris@princeton.edu)
 *  @version 1.0 
 *  
 *  @brief  Class MidpointLeap implements Realization step() function 
 *          using the midpoint tau leap approximate algorithm of Gillepsie
 *          (2001). The method is analogous to the deterministic midpoint 
 *          (2nd-order Runge-Kutta) method for the numerical solution of 
 *          ordinary differential equations.
 *  
 */

#include "midpointleap.h"
using namespace std;

/**
 *   @brief  Default constructor for MidpointLeap
 *  
 *   @param  model a GillespieModel object
 *   @param  rng is a random number generator of class RNG
 *   @param  timestep is the timestep to use for tau leaping
 * 
 *   @return nothing 
 */ 

MidpointLeap::MidpointLeap(FunctionModel *model,
                           RNG *rng,
                           double timestep) :
  Realization(model, rng),
  midpoint_array_(n_vars, 0),
  timestep_(timestep),
  verbose(false)
{}

/**
 *   @brief Destructor for MidpointLeap
 *  
 *   @return nothing 
 */ 

MidpointLeap::~MidpointLeap(){}


/**
 *   @brief Update waiting times
 *  
 *   @return int
 */  

int MidpointLeap::step(){
  take_step(timestep_);
  return 0;
}

void MidpointLeap::take_step(double stepsize){
  
  // populate midpoint array with current midpoints for each variable
  model->calc_midpoints(midpoint_array_,
                        rates,
                        state_array,
                        state_time,
                        stepsize);

  // hacky negative midpoint handling -- fix if time
  for(int var = 0; var < n_vars; var++){
    if(midpoint_array_.at(var) < +0.0){
      if(!(midpoint_array_.at(var) >= -0.0) && verbose){
        std::cout << state_time << ": negative midpoint at var ";
        std::cout << model->get_var_name(var) << " : ";
        std::cout << midpoint_array_.at(var);
        std::cout << " setting array to zero at that var";
        std::cout << std::endl;
        
      }
      midpoint_array_.at(var) = +0.0;
    }
  }

  // get event rates from midpoints
  model->update_rates(rates,
                      midpoint_array_,
                      state_time + 0.5 * stepsize);
  
  /* do each event a number of
     times that is poisson distributed
     according to its midpoint-determined rate */  
  
  int n_occurences = 0;
  new_state_ = state_array;
  
  for(int event = 0; event < n_events; event++){
    double timestep_rate = rates.at(event) * stepsize;

    // only call RNG if rate is positive
    if(timestep_rate > 0.0){
          n_occurences = rng->rpois(timestep_rate);

          // do the event n_occurences times
          if(n_occurences > 0)
            model->rep_event_i(event, n_occurences, new_state_);
    }
    
    else if(timestep_rate < -0.0001)
      throw runtime_error("Negative poisson rate geneated "
                          "beyond tolerance for double rounding");
  }

  // very hacky handling of 0 overshoots and negative midpoints
  for(int var = 0; var < n_vars; var++){
    if( new_state_.at(var) < +0.0 )
      new_state_.at(var) = +0.0;
    if( midpoint_array_.at(var) < 1e-20)
      new_state_.at(var) = +0.0;
  }

  model->record(state_array,
                new_state_,
                state_time + stepsize,
                sim_id,
                *event_repository);
  
  // update state and time
  state_array = new_state_;
  state_time += stepsize;

  array_ts.push_back(state_array);
  time_ts.push_back(state_time);
};


int MidpointLeap::simulate(double t_final, uint64_t max_iter, bool output){

  //update rates
  model->update_rates(rates,
                      state_array,
                      state_time);


  // prep for loop
  uint64_t n_full_steps =
    (uint64_t)floor((t_final - state_time)/timestep_);

  if(n_full_steps < 0)
    throw runtime_error("Tried to simulate to a time in the past");

  if(n_full_steps > max_iter)
    throw runtime_error("Attempt to take more than max_iter steps");

  uint64_t i_step = 0;

  // take steps until rates are zero or n_full_steps is reached
  while( (i_step < n_full_steps) &&
         (!rates_are_zero()) ){

    step();
        
    i_step++;
  }

  // take final (smaller) step
  if( (!rates_are_zero()) &&
      (t_final > state_time + timestep_ * FINAL_STEP_TOLERANCE) )
    take_step(t_final - state_time);

  if(output){
    for(int k = 0; k < time_ts.size(); k++){
      output_state(time_ts.at(k), array_ts.at(k), std::cout);
    }
  }
  return 0;
};

