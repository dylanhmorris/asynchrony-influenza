/**
 *  @file    rng.h
 *  @author  Dylan Morris (dhmorris@princeton.edu)
 *  @date    12/6/16  
 *  @version 1.0 
 *
 *  @brief encapsulating class for user-defined RNGs that give a 
 *  new uint64_t with next() and a new double with runif(), and 
 *  that take an integer seed.
 *  
 *
 */

#ifndef RNG_H_
#define RNG_H_
#include <stdint.h>
#include <math.h>
#include <stdexcept>
#include <cfloat> // for DBL_MIN
#include <vector>
#include <numeric> // for std::iota

using std::vector;

inline double to_double(uint64_t x);

class RNG {
 public:

  virtual void set_seed(int seed) = 0;
  
  /// get a new random int64
  virtual uint64_t next() = 0;

  /// get a new random uniform(0, 1) RV
  virtual double runif() = 0;

  /// get a new random exponential(lambda) RV
  double rexp(double lambda);

  /// get a new random binomial (n, p) RV
  int rbinom(int n, double p);

  /// get a new random poisson(mean) RV
  long long rpois(double mean);


  /// get a random vector of samples from the unit interval
  /// with latin square structure
  vector<double> runif_latin(int n_cutpoints);
  
  /// quick 2^64 calls to next (for parallelism)
  virtual void jump() = 0;
  

  /** log factorial function modified from public domain C# implementation by
   * John D. Cook (http://www.johndcook.com/blog/csharp_log_factorial/)
   * and PTRS algorithm by Wolfgang Hoermann (1993)
   */
  double log_factorial(long long k);


  /* generates a random integer between min (inclusive)
     and max (exclusive) */
  int randint(int min, int max);

  // picks a random element of a standard vector and returns it
  int randchoice(const vector<int>& input_vector);
  double randchoice(const vector<double>& input_vector);

  int weighted_choice(const vector<double> weights);

  double rgaussian(double mu, double sd);
  
  // get number of calls
  long long call_count();

  long long call_count_; ///< number of calls made

 private:
  long long poisson_knuth(double mean); ///< random poisson
  long long poisson_ptrs(double mean);  ///< random poisson
};

#endif
