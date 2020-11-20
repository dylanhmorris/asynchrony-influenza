/** 
 *  @file    TransmissionFunctions.h
 *  @author  Dylan Morris (dhmorris@princeton.edu)
 *  
 *  @brief TranmissionFunction.h defines a 
 *  set of functions for modeling transmission 
 *  between hosts, interfacing with the defined
 *  within-host models. Some are shared across
 *  within-host models. 
 *  
 */


#ifndef TRANSMISSION_FUNCTIONS_H_
#define TRANSMISSION_FUNCTIONS_H_

#include "../gillespie/GillespieRateFunctor.h"
#include "../random/rng.h"
#include "../gillespie/Realization.h"
#include "WithinhostMacros.h" // needed for backward compatibility

#include <cmath>

/* @brief continuous probability of transmission (to a naive)  
 * given a half prob value 
 * 
 * @param virions current number of virions in the transmitting host
 *
 * @param TD50_virions number of virions at which the transmission
 * probability to a naive is 50% 
 */
double trans_prob_proportional(const double virions,
                               const double TD50_virions){
  return 1 - exp(-log(2) * virions / TD50_virions);
}


/* @brief binary probability of transmission (to a naive) given a threshold 
 * 
 * @param virions current number of virions in the transmitting host
 * 
 * @param threshold_virions transmission threshold 
 * 
 * @return 1.0 if above threshold, else 0.0
 */
double trans_prob_threshold(const double virions,
                            const double threshold_virions){
  return 1.0 * (virions > threshold_virions);
}


/* @brief last state of simulation before a given time 
 * get the last state in the simulation before a given time
 * (needed since the exact time my never be achieved with 
 * uneven timesteps)  */
vector<double> last_state_before(const double time,
                                 const vector<vector<double>> &array_ts,
                                 const vector<double> &time_ts){
  for(int k = 0; k < time_ts.size() - 1; k++){
    if(time_ts.at(k + 1) > time)
      return array_ts.at(k);
  }
  
  throw std::runtime_error("Requested time was "
                           "never reached in simulation");
}

double get_viral_load(const vector<double> &state,
                      vector<int> virion_indices){
  double total_pop = 0;
  for(int & vir_ind : virion_indices){
    total_pop += state.at(vir_ind);
  }
  return total_pop;
}

/* was the infection ever above the transmission cutoff 
 * in terms of instantaneous transmission probability */
bool ever_transmissible_proportional(const Realization &model_state,
                                     const double trans_cutoff,
                                     const double transmission_parameter,
                                     vector<int> virion_indices){
  double total_pop = 0;

  for(int k = 0; k < model_state.array_ts.size(); k++){

    total_pop = 0;

    for(int & v_ind : virion_indices){
      total_pop += model_state.array_ts.at(k).at(v_ind);
    }
    
    if( trans_prob_proportional(total_pop,
                                transmission_parameter) >
        trans_cutoff ){
      
        return true;
    }
  }
  
  return false;
};

/* was the infection ever above the transmission cutoff 
 * in terms of instantaneous transmission probability */
bool ever_transmissible_threshold(const Realization &model_state,
                                  const double trans_cutoff,
                                  const double transmission_parameter,
                                  vector<int> virion_indices){
  double total_pop = 0;

  for(int k = 0; k < model_state.array_ts.size(); k++){

    total_pop = 0;

    for(int & v_ind : virion_indices){
      total_pop += model_state.array_ts.at(k).at(v_ind);
    }
    
    if( trans_prob_threshold(total_pop,
                             transmission_parameter) >
        trans_cutoff ){
      
        return true;
    }
  }
  
  return false;
};



double get_mutant_mean_freq(const Realization &model_state,
                            const double trans_threshold,
                            bool &was_transmissible){

  double weighted_total_mutant_freq = 0;
  double transmissible_time = 0;

  for(int k = 0; k < model_state.array_ts.size() - 1; k++){
    double total_pop = model_state.array_ts.at(k).at(VW_IND) +
      model_state.array_ts.at(k).at(VM_IND);
    if( total_pop > trans_threshold ){
      double dt = model_state.time_ts.at(k + 1) -
        model_state.time_ts.at(k);
      weighted_total_mutant_freq +=
        dt * model_state.array_ts.at(k).at(VM_IND) / total_pop;
      transmissible_time += dt;
    }
  }
  double output;
  if(transmissible_time > 0){
    was_transmissible = true;
    output = weighted_total_mutant_freq / transmissible_time;
  }
  else {
    was_transmissible = false;
    output = 0;
  }

  return output;
}


double get_mutant_peak_freq_threshold(const Realization &model_state,
                                      const double trans_cutoff,
                                      const double trans_parameter){
  double peak_freq = 0;
  
  for(int k = 0; k < model_state.array_ts.size(); k++){
    double total_pop = model_state.array_ts.at(k).at(VW_IND) +
      model_state.array_ts.at(k).at(VM_IND);
    if( total_pop > trans_parameter ){
      double mutant_freq = model_state.array_ts.at(k).at(VM_IND) / total_pop;
      if(mutant_freq > peak_freq)
        peak_freq = mutant_freq;
    }
  }

  return peak_freq;
};

double get_mutant_peak_freq_proportional(const Realization &model_state,
                                         const double trans_cutoff,
                                         const double trans_parameter){
  double peak_freq = 0;
  
  for(int k = 0; k < model_state.array_ts.size(); k++){
    double total_pop = model_state.array_ts.at(k).at(VW_IND) +
      model_state.array_ts.at(k).at(VM_IND);
    if(trans_prob_proportional(total_pop,
                               trans_parameter) >
       trans_cutoff){
      double mutant_freq =
        model_state.array_ts.at(k).at(VM_IND) / total_pop;
      
      if(mutant_freq > peak_freq)
        peak_freq = mutant_freq;
    }
  }

  return peak_freq;
};

bool did_ever_mutate(const Realization &model_state){
  for(int k = 0; k < model_state.array_ts.size(); k++){
    if(model_state.array_ts.at(k).at(VM_IND) > 0)
      return true;
  }
  return false;
};

double get_emergence_time(const Realization &model_state,
                          double emergence_threshold){
  double emergence_time = 99;

  // first find a mutant that survived to make at least
  // emergence_threshold copies
  int k = 0;
  bool found_emerged = false;
  while(k < model_state.array_ts.size() - 1 && !found_emerged){
    k++;
    found_emerged =
      model_state.array_ts.at(k).at(VM_IND) > emergence_threshold;
  }
  if(found_emerged){
    bool found_emergence = 0;
    int emerge_ind = k;
    while((emerge_ind > 1) && !found_emergence){
      found_emergence =
        model_state.array_ts.at(emerge_ind - 1).at(VM_IND) < 1;
      emerge_ind--;
    }
    if(found_emergence)
      emergence_time = model_state.time_ts.at(emerge_ind + 1);
  }
  return emergence_time;
}




#endif
