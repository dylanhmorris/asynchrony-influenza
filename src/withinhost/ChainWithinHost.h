/** 
 *  @file    ChainWithinHost.h
 *  @author  Dylan Morris (dhmorris@princeton.edu)
 *  
 *  @brief ChainWithinHost.h defines a derived class of WithinHostFunctor
 *  for doing transmission chain simulations according to the minimal
 *  within host model 
 *  
 */


#ifndef CHAIN_WITHIN_HOST_H_
#define CHAIN_WITHIN_HOST_H_

#include "../gillespie/GillespieRateFunctor.h"
#include "../random/rng.h"
#include "TransmissionFunctions.h"
#include "WithinhostMacros.h"

#include "../ini/INIReader.h"

#include <cmath>
#include <vector>


class ChainWithinHost : public GillespieRateFunctor {
public:

/**
 *   @brief  Default constructor for ChainWithinHost  
 *
 *   @param  beta governs the rate of target cell reduction events
 *   per C * V
 * 
 *   @param r is the growth rate of virions 
 *
 *   @param d_v is the death/decay rate of virions in the absence
 *   of adaptive immunity
 *
 *   @param k is the change in the death/decay rate of virions
 *   in the presence of adaptive immunity
 *
 *   @param mu is the (one-way) mutation rate
 * *
 *   @return nothing 
 **/

  ChainWithinHost(double beta,  
                  double r,
                  double mu,
                  double d_v,
                  double k,
                  double t_M,
                  double t_N,
                  double C_max,
                  double p_M,
                  double f,
                  double p_C) : 
  beta_(beta),
  r_(r),
  mu_(mu),
  d_v_(d_v),
  k_(k),
  t_M_(t_M),
  t_N_(t_N),
  C_max_(C_max),
  p_M_(p_M),
  f_(f),
  p_C_(p_C),
  logistic_response_(false),
  maj_imm_(0),
  min_imm_(0),
  homotypic_protection_(1),
  one_unit_protection_(0),
  threshold_model_(false)
  {
    int var_count = STATE_LENGTH;
    
    // initialize event matrix
    event_deltas_.push_back({ 1,  0,  0});
    event_deltas_.push_back({-1,  0,  0});
    event_deltas_.push_back({ 0,  1,  0});
    event_deltas_.push_back({ 0, -1,  0});
    event_deltas_.push_back({ 0,  0,  1});
    event_deltas_.push_back({ 0,  0, -1});
  };
  

  // destructor
  ~ChainWithinHost(){};

  // call method
  void operator()(vector<double> &rate_array,
                  const vector<double> &state_array,
                  double time) {

    double C = state_array.at(C_IND);
    double V_maj = state_array.at(VMAJ_IND);
    double V_min = state_array.at(VMIN_IND);

    // calculate and set rates
    rate_array.at(C_PLUS) = p_C_ * (C + 1) * (1 - (C / C_max_)) * (C <= C_max_);
    rate_array.at(C_MINUS) = f_ * beta_ * C * (V_maj + V_min);

    rate_array.at(VMAJ_PLUS) = r_ * beta_ * C * (1 - mu_) * V_maj;
    rate_array.at(VMAJ_MINUS) = (d_v_ + k_ * maj_imm_ * M(time)) * V_maj;
    rate_array.at(VMIN_PLUS) = r_ * beta_ * C * (V_min + mu_ * V_maj);
    rate_array.at(VMIN_MINUS) = (d_v_ + k_ * min_imm_ * M(time)) * V_min;
  }


  double implied_ln_titer(double sus){
    return (TITER_ALPHA + log(1/sus - 1) / TITER_BETA);
  }

  double sigmoid_susceptibility(double antigenic_distance){
    double ln_h_titer = implied_ln_titer(1 - homotypic_protection_);
    double ln_one_titer = implied_ln_titer(1 - one_unit_protection_);
    double ln_fold_drop = ln_h_titer - ln_one_titer;
    double mean_ht = exp(ln_h_titer);
    double fold_drop = exp(ln_fold_drop);
    double implied_titer = mean_ht / pow(fold_drop, antigenic_distance);
    return (1 / (1 + exp(TITER_BETA * (log(implied_titer) -
                                       TITER_ALPHA))));
  }

  double linear_susceptibility(double antigenic_distance){
    return min(antigenic_distance, 1.0);
  }
  
  double get_susceptibility(double phenotype){

    double sus = 1.0;
    
    for(auto& hist : history_){
      double antigenic_distance = fabs(hist - phenotype);
      double cross_sus = linear_susceptibility(antigenic_distance);
      if(cross_sus < sus)
        sus = cross_sus;
    }

    return sus;

  };

  
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

  
  void set_maj_phen(double new_phenotype){
    maj_phen_ = new_phenotype;
    maj_imm_ = 1 - get_susceptibility(new_phenotype);
  }
  
  void set_min_phen(double new_phenotype){
    min_phen_ = new_phenotype;
    double new_min_sus = get_susceptibility(new_phenotype);
    min_imm_ = 1 - new_min_sus;
  }
  
  double get_maj_phen(){
    return maj_phen_;
  }
  
  double get_min_phen(){
    return min_phen_;
  }

  void set_possible_histories(vector<vector<double>> new_possible_histories){
    possible_histories_ = new_possible_histories;
  }
  
  void set_homotypic_protection(double protection_level){
    homotypic_protection_ = protection_level;
  }
  
  void set_one_unit_protection(double protection_level){
    one_unit_protection_ = protection_level;
  }

  void add_to_history(double phenotype){
    history_.push_back(phenotype);
  }
  
  void set_history(int new_history_index){
    history_ = possible_histories_.at(new_history_index);
    history_index_ = new_history_index;
  }

  int get_history_index(){
    return history_index_;
  }

  void clear_history(){
    history_.clear();
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
  
  double get_neut_prob(double phenotype,
                       double n_encounters_iga){
    double overall_protection = 1 - get_susceptibility(phenotype);
    double protection = overall_protection * homotypic_protection_;
    double prob = 1 + (log(protection) / n_encounters_iga);
    prob = max(prob, 0.0);
    prob = min(prob, 1.0);
    return prob;
  }
  
  vector<double> get_inoculum(vector<double> contact_state,
                              double transmission_parameter,
                              double n_encounters_iga,
                              double final_bottleneck,
                              double maj_phen,
                              double min_phen,
                              RNG &rng,
                              bool killing){
    
    double total_pop = get_viral_load(contact_state,
                                      {VMAJ_IND, VMIN_IND});
    
    double p_transmit = trans_prob(total_pop,
                                   transmission_parameter);
    
    bool transmission = (rng.runif() < p_transmit);
    vector<double> result = {0.0, 0.0};

    if(transmission){
      
      double p_min = contact_state.at(VMIN_IND) / total_pop;
      
  
      int maj = rng.rpois(n_encounters_iga * (1 - p_min) *
                          (1 - get_neut_prob(maj_phen, n_encounters_iga) *
                           killing));
      int min = rng.rpois(n_encounters_iga * p_min *
                          (1 - get_neut_prob(min_phen, n_encounters_iga)
                           * killing));

      if (maj + min <= final_bottleneck) {
        result.at(0) = float(maj);
        result.at(1) = float(min);
      }

      else {
        for(int k = 0; k < final_bottleneck; k++){
          double prob_maj = maj / (maj + min * 1.0); // float division
          if (rng.runif() < prob_maj){
            result.at(0) += 1;
            maj -= 1;
          }
          else{
            result.at(1) += 1;
            min -= 1;
          }
        }
      }
    }
    
      return result;
  }

  double beta_;
  double r_;
  double mu_;
  double d_v_;
  double k_;
  double t_M_;
  double t_N_;
  double C_max_;
  double c_;
  double p_M_;
  double f_;
  double p_C_;
  bool logistic_response_;
  bool threshold_model_;
  
  private:
  double maj_imm_;
  double min_imm_;
  double maj_phen_;
  double min_phen_;
  double homotypic_protection_;
  double one_unit_protection_;
  std::vector<vector<double>> possible_histories_;
  std::vector<double> history_;
  int history_index_;
};

#endif
