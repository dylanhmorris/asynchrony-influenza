#include "../withinhost/MinimalWithinHost.h"
#include "../withinhost/MinimalParamset.h"

#include "../gillespie/FunctionModel.h"
#include "../random/xoroshiro128plus.h"
#include "../gillespie/AdaptivePRC.h"
#include <cmath>
#include <fstream>
#include <algorithm>
#define I_WT 1
#define I_MUT 2
#define MAX_ITER 50000000
#define ITER_PRINT_FREQ 10000

using std::ofstream;


void print_header(ostream& outstream){
  outstream << std::left << setw(10) << "rng_seed";
  outstream << std::left << setw(18) << "p_mut_inoc";
  outstream << setw(18) << "c";
  outstream << setw(10) << "Vw_iga";
  outstream << setw(10) << "Vm_iga";
  outstream << std::left << setw(10) << "Vw_inoc";
  outstream << std::left << setw(10) << "Vm_inoc";
  outstream << std::left << setw(15) << "is_inoc";
  outstream << std::left << setw(15) << "transmissible";
  outstream << std::left << setw(18) << "peak_mutant_freq";
  outstream << std::left << setw(15) << "bottleneck";
  outstream << std::left << setw(16) << "trans_threshold";
  outstream << std::left << setw(16) << "emergence_time";
  outstream << std::left << setw(16) << "ever_mutated";
  outstream << std::endl;
}


int main(int argc, char *argv[]){
  if(argc < 7){
    printf("\nUSAGE: ./MinimalPSelection <path to inits> "
           "<model name> <n iter> <bottleneck> <output path stem> "
           "<threshold>"
           "\n\n"
           "path to inits: an .ini/.mk formatted initialization file "
           "containing the default parameter values for the minimal "
           "model and possibly model-specific overrides\n\n"
           "model name: name of the model "
           "(referenced in param overrides)\n\n "
           "n iter: number of runs to do with these parameters\n\n"
           "bottleneck/vb ratio: bottleneck size or vb ratio\n\n"
           "output path stem: stem of path to output results\n\n"
           "threshold: whether to use a threshold or proportional "
           "transmission probability model\n\n"
           );

    return 0;
  }

  // parse arg variables
  string input_path = argv[1];
  string model_name = argv[2];
  int n_iter = atoi(argv[3]);
  int independent_variable = atoi(argv[4]);
  string out_stem = argv[5];
  string out_path = out_stem + ".txt";
  bool threshold = atoi(argv[6]);
  
  // get default simulation parameters
  // within host model
  MinimalParamset params(input_path,
                         model_name);

  // sterilizing during timecourse model
  // varies the bottleneck; others fix
  // the bottleneck and vary the VB ratio
  if(model_name == "minimal_sterilizing"){
      std::cout << "Setting bottleneck..." << std::endl;
      params.set_parameter("BOTTLENECK",
                           independent_variable);
      if(independent_variable < 11){
          n_iter *= 10;
          std::cout << "Using 10x iterations " <<
            "because of small bottleneck" << std::endl;
      };
  } else{
    std::cout << "Setting vb ratio" << std::endl;
    params.set_parameter("VB_RATIO",
                         float(independent_variable));
  }
  // calculated and renamed params

  std::cout << "Setting up rng" << std::endl;
  xoroshiro128plus model_rng(0);
  int current_seed = params.get_integer_parameter("rng_seed") +
    n_iter * independent_variable + 29;


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

  double trans_parameter;
  double trans_cutoff = params.get_parameter("trans_cutoff");
  
  if(threshold){
    trans_parameter = params.get_parameter("trans_threshold");
    model_functor.threshold_model();
    std::cout << "Using a threshold transmission model..."
            << std::endl;
  } else {
    trans_parameter = params.get_parameter("TD50");
    model_functor.proportional_model();
    std::cout << "Using a proportional transmission probability "
      "model..." << std::endl;
  }
  
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
  vector<double> inoculator_init_state =
    {params.get_parameter("C_max"),
     params.get_integer_parameter("bottleneck") * 1.0,
     0};
  
  AdaptivePRC inoculator_state(&model,
                               &model_rng,
                               params.get_parameter("inoculator_leap_timestep"));
  
  AdaptivePRC recipient_state(&model,
                              &model_rng,
                              params.get_parameter("leap_timestep"));


  // prep for output
  ofstream outstream;  
  outstream.open(out_path);
  print_header(outstream);
  // to check whether we've seen an example of each of the following
  bool repl_visible = false;
  bool repl_invisible = false;
  bool inoc_visible = false;
  bool inoc_invisible = false;

  // to check whether the inoculator
  // has produced a transmissible
  // infection
  bool transmissible_inoculator = false;
  
  for(int i_run = 0; i_run < n_iter; i_run++){

    // simulate a naive transmitting host

    model_functor.E_w_ = false;
    transmissible_inoculator = false;
    
    while(!transmissible_inoculator){
      inoculator_state.reset_all();
      inoculator_state.set_state(inoculator_init_state);
      inoculator_state.set_time(0);
      inoculator_state.simulate(params.get_parameter("t_final"),
                                MAX_ITER,
                                false);
      transmissible_inoculator =
        model_functor.ever_transmissible(inoculator_state,
                                         trans_cutoff,
                                         trans_parameter);
    };
    
    if(i_run < 1){
      string timecourse_path = out_stem + "_naive.txt";
      inoculator_state.output_all_to_file(timecourse_path);
    }

    double p_mut =
      model_functor.get_random_transmissible_pmut(inoculator_state,
                                                  trans_parameter,
                                                  model_rng);

    // simulate an immune recipient host
    recipient_state.reset_all();
    model_functor.E_w_ = true;
    model_functor.E_m_ = false;

    // immunity acts iff immune to wildtype
    // assumed immune to mutant
    double wt_neut = model_functor.E_w_ *
      params.get_parameter("wt_neut_prob");
    double mut_neut = model_functor.E_w_ *
      params.get_parameter("mut_neut_prob");
    
    vector<double> inoculum =
      model_functor.get_inoculum(p_mut,
                                 params.get_parameter("n_encounters_iga"),
                                 params.get_integer_parameter("bottleneck"),
                                 wt_neut,
                                 mut_neut,
                                 model_rng);

    double Vw_init = inoculum.at(0);
    double Vm_init = inoculum.at(1);
    double Vw_iga = inoculum.at(2);
    double Vm_iga = inoculum.at(3);

    double peak_mutant_freq;
    bool transmissible;
    bool inoc;
    double emergence_time = 99;
    bool ever_mutated = false;
    
    if((Vw_init + Vm_init) < 1){
      peak_mutant_freq = 0;
      transmissible = false;
      inoc = false;
    } else {
    
    inoc = Vm_init > 0;

    vector<double> init_state = {params.get_parameter("C_max"),
                                 Vw_init,
                                 Vm_init};
    // naive to mutant
    recipient_state.reset_all();
    recipient_state.set_state(init_state);
    recipient_state.set_time(0);

    recipient_state.simulate(params.get_parameter("t_final"),
                             MAX_ITER,
                             false);

    transmissible =
      model_functor.ever_transmissible(recipient_state,
                                       trans_cutoff,
                                       trans_parameter);
    
    peak_mutant_freq =
      model_functor.get_mutant_peak_freq(recipient_state,
                                         trans_cutoff,
                                         trans_parameter);

    emergence_time =
      get_emergence_time(recipient_state,
                         params.get_parameter("emergence_threshold"));

    ever_mutated = did_ever_mutate(recipient_state);
    }
    
    if(inoc && transmissible && !inoc_visible && peak_mutant_freq > 0.5){
      string timecourse_path = out_stem + "_inoc_visible.txt";
      recipient_state.output_all_to_file(timecourse_path);
      bool inoc_visible = true;
    }
    else if(inoc && !transmissible && !inoc_invisible){
      string timecourse_path = out_stem + "_inoc_invisible.txt";
      recipient_state.output_all_to_file(timecourse_path);
      bool inoc_invisible = true;
    }
    else if(!inoc && transmissible && !repl_visible){
      string timecourse_path = out_stem + "_repl_visible.txt";
      recipient_state.output_all_to_file(timecourse_path);
      bool repl_visible = true;
    }
    else if(!inoc && !transmissible && !repl_invisible){
      string timecourse_path = out_stem + "_repl_invisible.txt";
      recipient_state.output_all_to_file(timecourse_path);
      bool repl_invisible = true;
    }

    
    // print out the result
    outstream << setw(10) << current_seed;
    outstream << setw(18) << setprecision(10) << p_mut;
    outstream << setw(18) << setprecision(10) <<
      params.get_parameter("cross_imm");
    outstream << setw(10) << Vw_iga;
    outstream << setw(10) << Vm_iga;
    outstream << setw(10) << Vw_init;
    outstream << setw(10) << Vm_init;
    outstream << left << setw(15) << boolalpha << inoc;
    outstream << left << setw(15) << boolalpha << transmissible;
    outstream << left << setw(18) << setprecision(10) << peak_mutant_freq;
    outstream << left << setw(15) <<
      params.get_integer_parameter("bottleneck");
    outstream << left << setw(16) << params.get_parameter("trans_threshold");
    outstream << left << setw(16) << emergence_time;
    outstream << left << setw(16) << boolalpha << ever_mutated;
    outstream << std::endl;
    // reset and increment seed
    current_seed++;
    model_rng.set_seed(current_seed);
    if(i_run % ITER_PRINT_FREQ == 0){
      std::cout << "Iteration no.: " << i_run << " of "
                << n_iter << std::endl;
    }
  }
  
  return 0;
}

