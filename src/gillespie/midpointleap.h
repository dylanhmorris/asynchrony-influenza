/** 
 *  @file    midpointleap.h
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

#ifndef MIDPOINTLEAP_H_
#define MIDPOINTLEAP_H_

/* FINAL_STEP TOLERANCE declares up to how many 
times smaller than global leaping timestep will 
be permitted for a final step (e.g. with a 
timestep of 0.01, a final step smaller than 1e-14 
will be ignored). This prevents re-doing an exact 
hit due to roundoff errorwhen the timestep 
evenly divides the remaining time */

#define FINAL_STEP_TOLERANCE 1e-12

#include <stdexcept>  

#include "Realization.h"

/**  
 *  @brief  Class MidpointLeap implements Realization step() function 
 *          using the midpoint tau leap approximate algorithm of Gillepsie
 *          (2001). The method is analogous to the deterministic midpoint 
 *          (2nd-order Runge-Kutta) method for the numerical solution of 
 *          ordinary differential equations.
 */  
class MidpointLeap : public Realization {
 public:
  // Default constructor for MidpointLeap
  MidpointLeap(FunctionModel *model,
               RNG *rng,
               double timestep);
  // Destructor for MidpointLeap
  ~MidpointLeap();

  // simulate to an exact time
  int simulate(double t_final, uint64_t max_iter, bool output);
  
  // take one simulation step using MidpointLeap
  // using the default stepsize
  int step();

  // take a step with a given stepsize 
  void take_step(double stepsize);

  bool verbose;  ///< whether or not to print errors, etc. default false;

 private:
  vector<double> midpoint_array_; ///< midpoint_array_ is a double array holding the current midpoints for the tau-leap step
  vector<double> new_state_; ///< double array holding the updated state
  double timestep_; ///< size of timestep for midpoint tau leap
};

#endif 
