/** 
 *  @file    firstreaction.h
 *  @author  Dylan Morris (dhmorris@princeton.edu)
 *  @version 1.0 
 *  
 *  @brief   Class FirstReaction implements Realization step() function
 *           using the exact First Reaction algorithm of Gillespie (1971)
 */

#ifndef FIRSTREACTION_H_
#define FIRSTREACTION_H_
#include "Realization.h"

/**
 *  @brief Class FirstReaction implements GillespieInt step() function
 *           using the exact First Reaction algorithm of Gillespie (1971)
 */
class FirstReaction : public Realization {
 public:
  // Default constructor for FirstReaction
  FirstReaction(FunctionModel *model, RNG *rng);
  
  // Destructor for FirstReaction
  ~FirstReaction();
 
  int step();
 private:
  vector<double> waiting_times;  ///< array to hold waiting times to events
};

#endif
