/** 
 *  @file    AdaptivePRC.cc
 *  @author  Dylan Morris (dhmorris@princeton.edu)
 *  @version 1.0 
 *  
 *  @brief   Class AdaptivePRC implements Realization step() function
 *           using the adaptive approximate Poisson Random Corrections (PRC) 
 *           method of Hu and Li
 *  
 */

#include "AdaptivePRC.h"

/**
 *   @brief  Default constructor for AdaptivePRC
 *  
 *   @param  model is a GillespieModel object
 *   @param  rng is a random number generator of class RNG
 * 
 *   @return nothing 
 */ 
AdaptivePRC::AdaptivePRC(FunctionModel *model,
                         RNG *rng, double max_tau) :
  Realization(model, rng),
  max_tau(max_tau),
  eta_rates_prod(n_events),
  max_occurences_(std::numeric_limits<long long>::max() / 100.0)
{
  etas.clear();
  for(int k = 0; k < n_events; k++){
    etas.push_back(vector<double>(n_events, 0));
  }
}

/**
 *   @brief Destructor for AdaptivePRC
 *  
 *   @return nothing 
 */ 
AdaptivePRC::~AdaptivePRC(){
}


/**
 *   @brief Update waiting times
 *  
 *   @return int
 */  
int AdaptivePRC::step(){
  model->update_rates(rates, state_array, state_time);

  update_eta_rates_prod();

  double timestep = prc_timestep();
  // update state, time, rates


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

  // very hacky handling of 0 overshoots and negative midpoints
  for(int var = 0; var < n_vars; var++){
    if( new_state_.at(var) < +0.0 )
      new_state_.at(var) = +0.0;
  }
  
  
  // update state and time
  state_array = new_state_;
  state_time += timestep;

  array_ts.push_back(state_array);
  time_ts.push_back(state_time);

  return 0;
};


void AdaptivePRC::update_etas(){
  vector<double> a_js(n_events, 0);
  vector<double> projection(n_vars, 0);
  for(int column = 0; column < n_events; column++){
    for(int ivar = 0; ivar < n_vars; ivar++){
      projection.at(ivar) = state_array.at(ivar) + model->get_delta(column, ivar);
    }
    
    model->update_eta_rates(a_js, projection, state_time);

    for(int row = 0; row < n_events; row++){
      etas.at(row).at(column) = a_js.at(row) - rates.at(row);
    }
  }
};

void AdaptivePRC::update_eta_rates_prod(){

  update_etas();
  for(int row = 0; row < n_events; row++){
    eta_rates_prod.at(row) = 0;

    for(int col = 0; col < n_events; col++)
      eta_rates_prod.at(row) += etas.at(row).at(col) * rates.at(col);
  }
}

double AdaptivePRC::prc_timestep(){
  double eta_cond = -1000;
  double max_rate = 0;

  for(int i_rate = 0; i_rate < n_events; i_rate++){
    double ev_rate = rates.at(i_rate);
    double candidate = ev_rate / (10 * (eta_rates_prod.at(i_rate) + 1e-10));
    if( (candidate < 0) && (candidate > eta_cond))
      eta_cond = candidate;

    if(ev_rate > max_rate)
      max_rate = ev_rate;
  }

  return min(min(abs(eta_cond), max_tau), max_occurences_ / max_rate);
}
