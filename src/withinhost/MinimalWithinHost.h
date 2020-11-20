/** 
 *  @file    MinimalWithinHost.h
 *  @author  Dylan Morris (dhmorris@princeton.edu)
 *  
 *  @brief MinimalWithinHost.h defines a derived 
 *  class of GillespieRateFunctor for doing gillespie 
 *  simulations of a flu within host model
 *  incorporating target-cell limitations / 
 * innate immunity and a very simple constant
 *  adaptive response 
 *  
 */


#ifndef MINIMAL_WITHIN_HOST_H_
#define MINIMAL_WITHIN_HOST_H_

#include "../gillespie/GillespieRateFunctor.h"
#include "../random/rng.h"
#include "../gillespie/Realization.h"
#include "TransmissionFunctions.h"
#include "WithinhostMacros.h"

#include <cmath>



class MinimalWithinHost : public GillespieRateFunctor {
 public:

/**
 *   @brief  Default constructor for MinimalWithinHost  
 *
 *   @param  beta governs the rate of target cell reduction events
 *   per C * V
 * 
 *   @param r_w is the growth rate of wild-type virions 
 * 
 *   @param r_m is the growth rate of mutant virions
 *
 *   @param d_v is the death/decay rate of virions in the absence
 *   of adaptive immunity
 *
 *   @param k is the change in the death/decay rate of virions
 *   in the presence of adaptive immunity
 *
 *   @param mu is the symmetric mutation rate
 *
 *   @param c is symmetric cross immunity
 *
 *   @return nothing 
 **/

 MinimalWithinHost(double beta,  
                   double r_w,
                   double r_m,
                   double mu,
                   double d_v,
                   double k,
                   double t_M,
                   double t_N,
                   double Cmax,
                   double c,
                   double p_M,
                   double f,
                   double p_C,
                   double E_w,
                   double E_m) : 
  beta_(beta),
    r_w_(r_w),
    r_m_(r_m),
    mu_(mu),
    d_v_(d_v),
    k_(k),
    t_M_(t_M),
    t_N_(t_N),
    Cmax_(Cmax),
    c_(c),
    p_M_(p_M),
    f_(f),
    p_C_(p_C),
    E_w_(E_w),
    E_m_(E_m),
  logistic_response_(false),
  threshold_model_(false)
  {

      int var_count = 3;
      
      // initialize event matrix
      event_deltas_.push_back({ 1,  0,  0});
      event_deltas_.push_back({-1,  0,  0});
      event_deltas_.push_back({ 0,  1,  0});
      event_deltas_.push_back({ 0, -1,  0});
      event_deltas_.push_back({ 0,  0,  1});
      event_deltas_.push_back({ 0,  0, -1});

      // make it into a ones event matrix
      
    }
  

  // destructor
  ~MinimalWithinHost(){};
  
  // call method
  void operator()(vector<double> &rate_array,
                  const vector<double> &state_array,
                  double time) {

    // for now, only model novel
    // responses to wild-type
    double t_N_w = t_N_;
    double t_N_m = 999999;
    
    double C = state_array.at(C_IND);
    double Vw = state_array.at(VW_IND);
    double Vm = state_array.at(VM_IND);
    double E_w_eff = max(E_w_, (time > t_N_w) * 1.0);
    double E_m_eff = max(E_m_, (time > t_N_m) * 1.0);
    double Mw_eff = max(E_w_eff, c_ * E_m_eff);
    double Mm_eff = max(E_m_eff, c_ * E_w_eff);


    // calculate and set rates
    rate_array.at(C_PLUS) = p_C_ * (C + 1) * (1 - (C / Cmax_)) * (C <= Cmax_);
    rate_array.at(C_MINUS) = f_ * beta_ * C * (Vw + Vm);

    rate_array.at(VW_PLUS) = beta_ * C * r_w_ * (1 - mu_) * Vw;
    rate_array.at(VW_MINUS) =  (d_v_ + k_ * Mw_eff * M(time)) * Vw;

    rate_array.at(VM_PLUS) = beta_ * C * (r_m_ * Vm + r_w_ * mu_ * Vw);
    rate_array.at(VM_MINUS) = (d_v_ + k_ * Mm_eff * M(time)) * Vm;
  }

  double M(double time){
    double M = 0;
    if(logistic_response_){
      double x = p_M_ * (time - t_M_);
      M = 1 / (1 + exp(-x));
    }
    else { M = time > t_M_; }
    return M;
  }
  
  void logistic_response(){
    logistic_response_ = true;
  }
  
  void binary_response(){
    logistic_response_ = false;
  }

  void threshold_model(){
    threshold_model_ = true;
  }
  
  void proportional_model(){
    threshold_model_ = false;
  }

  
  
  double trans_prob(double total_pop,
                    double transmission_parameter){
    if( threshold_model_ )
      return trans_prob_threshold(total_pop,
                                  transmission_parameter);
    else
      return trans_prob_proportional(total_pop,
                                     transmission_parameter);
  }


  double get_random_transmissible_pmut(const Realization &model_state,
                                       const double trans_param,
                                       RNG &rng){

  vector<double> transmissible_freqs;
  vector<double> weights;
  double total_pop = 0;

  // first get all the transmissible moments
  for(int k = 0; k < model_state.array_ts.size() - 1; k++){
    
    total_pop = model_state.array_ts.at(k).at(VW_IND) +
      model_state.array_ts.at(k).at(VM_IND);
    
    double dt = model_state.time_ts.at(k + 1) -
      model_state.time_ts.at(k);
    
    transmissible_freqs.push_back(model_state.array_ts.at(k).at(VM_IND) /
                                  total_pop);

    weights.push_back(dt * trans_prob(total_pop,
                                      trans_param));
  }

  // now pick one at random
  double result = -1;

  if(transmissible_freqs.size() > 0){
    
    int choice = rng.weighted_choice(weights);

    result = transmissible_freqs.at(choice);
    
  }

  return result;
  
};

  
  double ever_transmissible(const Realization &model_state,
                            const double trans_cutoff,
                            const double transmission_parameter){
    
    vector<int> virion_inds = {VW_IND, VM_IND};

    if( threshold_model_ )
      return ever_transmissible_threshold(model_state,
                                          trans_cutoff,
                                          transmission_parameter,
                                          virion_inds);
    else
      return ever_transmissible_proportional(model_state,
                                             trans_cutoff,
                                             transmission_parameter,
                                             virion_inds);

  }

  double get_mutant_peak_freq(const Realization &model_state,
                              const double trans_cutoff,
                              const double trans_parameter){
    if( threshold_model_ )
      return get_mutant_peak_freq_threshold(model_state,
                                            trans_cutoff,
                                            trans_parameter);
    else
      return get_mutant_peak_freq_proportional(model_state,
                                               trans_cutoff,
                                               trans_parameter);
  }
  
  vector<double> get_inoculum(double p_mut,
                              double n_encounters_iga,
                              double final_bottleneck,
                              double wt_neut_prob,
                              double mut_neut_prob,
                              RNG &rng){
    vector<double> result = {0.0, 0.0, 0.0, 0.0};
    
    int wt = rng.rpois(n_encounters_iga * (1 - p_mut) *
                       (1 - wt_neut_prob));
    int mut = rng.rpois(n_encounters_iga * p_mut *
                        (1 - mut_neut_prob));
    
    // save pre-bottleneck
    result.at(2) = float(wt);
    result.at(3) = float(mut);
    
    if (wt + mut <= final_bottleneck) {
      result.at(0) = float(wt);
      result.at(1) = float(mut);
    }

    else {
      for(int k = 0; k < final_bottleneck; k++){
        double prob_wt = wt / (wt + mut * 1.0); // float division
        if (rng.runif() < prob_wt){
          result.at(0) += 1;
          wt -= 1;
        }
        else{
          result.at(1) += 1;
          mut -= 1;
        }
      }
    }
    
    return result;
};

  

  double E_w_;
  double E_m_;
  double beta_;
  double r_w_;
  double r_m_;
  double mu_;
  double d_v_;
  double k_;
  double t_M_;
  double t_N_;
  double Cmax_;
  double c_;
  double p_M_;
  double f_;
  double p_C_;
  bool logistic_response_;
  bool threshold_model_;
};

#endif
