
/** 
 *  @file    xoroshiro128plus.h
 *  @author  Dylan Morris (dhmorris@princeton.edu)
 *  @date    12/6/16  
 *  @version 1.0 
 *  
 *  @brief Class xoroshorio128plus implements a random number 
 *         generator of Class rng
 *  
 */


#ifndef XOROSHIRO128_H_
#define XOROSHIRO128_H_
#include <stdint.h>
#include "rng.h"

/*  new C++ class based on xorshift128+ implementation by David Blackman and Sebastiano Vigna (vigna@acm.org), which appears with the following message:

Written in 2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */


/**
 *  @brief Class xoroshiro128plus implements a random number 
 *         generator of Class rng
 *  
 */
class xoroshiro128plus : public RNG {
 public:
  /// Default constructor 
  xoroshiro128plus(int seed);

  // Destructor
  ~xoroshiro128plus();

  // set a new seed
  void set_seed(int seed);
  
  // get a new random int64
  uint64_t next();

  // get a new random uniform(0, 1) RV
  double runif();

  // quick 2^64 calls to next (for parallelism)
  void jump();

 private:
  // state of the generator
  uint64_t s[2];

  // simulated rotate
  inline uint64_t rotl(const uint64_t x, int k);

  // splitmix64 implementation transforming seed
  uint64_t splitmixstate;
  uint64_t splitmix64();

};

#endif
