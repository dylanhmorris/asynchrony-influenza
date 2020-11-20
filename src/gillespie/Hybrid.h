/** 
 *  @file    Hybrid.h
 *  @author  Dylan Morris (dhmorris@princeton.edu)
 *  @version 1.0 
 *  
 *  @brief  Class Hybrid implements Realization step() function 
 *          using a janky hybrid of the
 *          random corrections approximate algorithm of Hu and Li
 *          and the SSA (Gillespie, 1977)
 */

#ifndef HYBRID_H_
#define HYBRID_H_
#include "AdaptivePRC.h"
#include <algorithm> // for std::any_of

#define WT_IND 1
#define MUT_IND 2


/**
 *  @brief  Class Hybrid implements Realization step() function 
 *          using a janky hybrid of the
 *          random corrections approximate algorithm of Hu and Li
 *          and the SSA (Gillespie, 1977)
 */

class Hybrid : public AdaptivePRC {
 public:
  // Default constructor for Hybrid
  Hybrid(FunctionModel *model, RNG *rng, double max_tau);
  
  // Destructor for Hybrid
  ~Hybrid();
 
  int step();
  void ssa_step();
  void prc_step();
  int simulate(double t_final, uint64_t max_iter, bool output);
  bool use_ssa();
  bool verbose;

 private:
  double timestep_;
  vector<double> waiting_times_;
};

#endif
