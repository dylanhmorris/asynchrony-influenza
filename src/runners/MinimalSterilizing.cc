#include "../withinhost/MinimalWithinHost.h"
#include "../withinhost/MinimalParamset.h"

#include "../gillespie/FunctionModel.h"
#include "../random/xoroshiro128plus.h"
#include "../gillespie/Ad.h"
#include <cmath>
#include <fstream>
#include <sstream> // for c/k values in string 
#include <algorithm>
#define I_WT 1
#define I_MUT 2
#define MAX_ITER 50000000

using std::ofstream;


void print_header(ostream& outstream){
  outstream << std::left << setw(10) << "rng_seed";
  outstream << setw(10) << "k";
  outstream << setw(10) << "t_M";
  outstream << std::left << setw(15) << "transmissible";
  outstream << std::left << setw(18) << "end_mutant_freq";
  outstream << std::left << setw(16) << "emergence_time";
  outstream << std::left << setw(16) << "ever_mutated";
  outstream << std::endl;
}


int main(int argc, char *argv[]){
  if(argc < 7){
    printf("\nUSAGE: ./MinimalPSelection <path to inits> "
           " <model name> <n iter> <bottleneck> <output path stem> "
           "\n\n"
           "path to inits: an .ini/.mk formatted initialization file "
           "containing the default parameter values for the minimal "
           "model and possibly model-specific overrides\n\n"
           "model name: name of the model "
           "(referenced in param overrides)\n\n "
           "n iter: number of runs to do with these parameters\n\n"
           "t_final: time to simulate to\n\n"
           "bottleneck: bottleneck size\n\n"
           "output path: path to output results\n\n"
           );

    return 0;
  }

  // parse arg variables
  string input_path = argv[1];
  string model_name = argv[2];
  int n_iter = atoi(argv[3]);
  int t_final = atof(argv[4]);
  int bn = atoi(argv[5]);
  string out_path = argv[6];
  
  // get default simulation parameters
  // within host model
  MinimalParamset params(input_path,
                         model_name);
  // calculated and renamed params
  int transformed_seed = int(params.get_integer_parameter("rng_seed"));

  std::cout << "transformed seed: " << transformed_seed << std::endl;
  
  xoroshiro128plus model_rng(transformed_seed);
  int current_seed = model_rng.randint(0, n_iter * 100);

    // prep for output

  ofstream outstream;
  outstream.open(out_path);
  print_header(outstream);
      
  MinimalWithinHost model_functor(params.get_parameter("BETA"),
                                  params.get_parameter("r_w"),
                                  params.get_parameter("r_m"),
                                  params.get_parameter("mu"),
                                  params.get_parameter("d_v"),
                                  params.get_parameter("k"),
                                  params.get_parameter("t_M"),
                                  params.get_parameter("t_N"),
                                  params.get_parameter("C_max"),
                                  params.get_parameter("cross_imm"),
                                  params.get_parameter("p_M"),
                                  params.get_parameter("f"),
                                  params.get_parameter("p_C"),
                                  params.get_integer_parameter("E_w"),
                                  params.get_integer_parameter("E_m"));

  std::cout << "Setting up model with the following parameters..."
            << std::endl;
  params.print_paramset(std::cout);
  
  FunctionModel model(model_functor.n_vars(),
                      &model_functor);
      
  // set var names
  model.set_var_name(0, "C");
  model.set_var_name(1, "Vw");
  model.set_var_name(2, "Vm");
  // set initial state
  vector<double> init_state =
    {params.get_parameter("C_max"),
     bn * 1.0,
     0};
  
  AdaptivePRC model_state(&model,
                          &model_rng,
                          params.get_parameter("leap_timestep"));



  
  for(int i_run = 0; i_run < n_iter; i_run++){

    // simulate a experienced host
    model_functor.E_w_ = true;
    model_state.set_state(init_state);
    model_state.set_time(0);
    model_state.simulate(t_final,
                         MAX_ITER,
                         false);

    double mean_mutant_freq;
    double peak_mutant_freq;
    double end_mutant_freq;
    bool transmissible;
    bool ever_mutated;
    bool inoc;
    double emergence_time = 99;
    
    transmissible = ever_transmissible(model_state,
                                       params.get_parameter("trans_threshold"));

    end_mutant_freq =
      model_state.state_array.at(2) /
      (model_state.state_array.at(2) +
       model_state.state_array.at(1));

    emergence_time =
      get_emergence_time(model_state,
                         params.get_parameter("emergence_threshold"));

    ever_mutated = did_ever_mutate(model_state);
    
    // print out the result
    outstream << setw(10) << current_seed;
    outstream << setw(10) << params.get_parameter("k");
    outstream << setw(10) << params.get_parameter("t_M");
    outstream << left << setw(15) << boolalpha << transmissible;
    outstream << left << setw(18) << setprecision(10) << end_mutant_freq;
    outstream << left << setw(16) << emergence_time;
    outstream << left << setw(16) << boolalpha << ever_mutated;
    outstream << std::endl;
    // reset and increment seed
    model_state.reset_all();
    current_seed++;
    model_rng.set_seed(current_seed);
    if(i_run % 10000 == 0){
      std::cout << "Iteration no.: " << i_run << " of "
                << n_iter << std::endl;
    }
    
  }

  outstream.close();
  return 0;
};
