/** 
 *  @file    AdaptivePRC.h
 *  @author  Dylan Morris (dhmorris@princeton.edu)
 *  @version 1.0 
 *  
 *  @brief   Class AdaptivePRC implements Realization step() function
 *           using the adaptive approximate Poisson Random Corrections (PRC) 
 *           method of Hu and Li
 */


// TODO: redo this using Eigen

#ifndef ADAPTIVE_PRC_H_
#define ADAPTIVE_PRC_H_
#include "Realization.h"
#include <cmath> // for abs
#include <limits> // for std::numeric_limits

/**
 *  @brief Class AdaptivePRC implements Realization step() function
 *           using the PRC algorithm of Hu and Li.
 */
class AdaptivePRC : public Realization {
 public:
  // Default constructor for AdaptivePRC
  AdaptivePRC(FunctionModel *model, RNG *rng,
              double max_tau);
  
  // Destructor for AdaptivePRC
  ~AdaptivePRC();
 
  int step();

  void update_etas();
  void update_eta_rates_prod();
  double prc_timestep();
  double max_tau;
  vector<vector<double>> etas;
  vector<double> eta_rates_prod;
  vector<double> new_state_;

 private:
  double max_occurences_;
};

#endif
