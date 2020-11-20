#include "ChainModel.h"
#include "../ini/INIReader.h"
#include <cmath>
#include <fstream>
#include <algorithm>



int main(int argc, char *argv[]){
  if(argc < 8){
    printf("\nUSAGE: ./TransmissionChain <path to inits>"
           " <model name> <n iter> <bottleneck> <mut mu> <mut sd> <init phen> <chain length> <output path> "
           "\n\n"
           "path to inits: an .ini/.mk formatted initialization file "
           "containing the default parameter values for the minimal "
           "model and possibly model-specific overrides\n\n"
           "model name: name of the model "
           "(referenced in param overrides)\n\n"
           "n iter: number of runs to do with these parameters\n\n"
           "chain length: length of each simulated chain\n\n"
           "output path: path to output results\n\n"
           "which iter: which subfile\n\n"
           "use_threshold: whether to use a threshold model "
           "(versus a proportional one) for transmission probability)\n\n"
           );

    return 0;
  }

  // parse arg variables
  string input_path = argv[1];
  string model_name = argv[2];
  int n_iter = atoi(argv[3]);
  int max_chain = atoi(argv[4]);
  string out_path = argv[5];
  int which_iter = atoi(argv[6]);
  bool use_threshold = atoi(argv[7]);

    
  // get default simulation parameters
  // within host model
  MinimalParamset params(input_path,
                         model_name);

  int current_seed = params.get_integer_parameter("rng_seed") *
    (which_iter + 1);
  std::cout << "rng seed: " << current_seed << std::endl;

  bool kill = KILLING;
  std::cout << "killing: " << kill << std::endl;

  std::cout << "setting up simulator..." << std::endl;
  // Create simulator class
  ChainModel chain_model(params,
                         out_path);

  if(use_threshold){
    std::cout << "using transmission threshold model with threshold"
              << std::endl;
    chain_model.use_threshold_transmission_model();
  } else {
    std::cout << "using transmission threshold model with probability"
      " proportional to viral load" 
              << std::endl;
    chain_model.use_proportional_transmission_model();
  }

  // loop over chains to simulate, outputing results
  for(int chain_no = 0; chain_no < n_iter; chain_no++){
    std::cout << "simulating chains..." << std::endl;
    chain_model.simulate_chain(current_seed,
                               chain_no,
                               max_chain);
    
    current_seed++;
    
  }
  // end of loop over chains to simulate
  
  return 0;
  
}
// end of program
