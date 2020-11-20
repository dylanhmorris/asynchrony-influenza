#include "../withinhost/MinimalWithinHost.h"
#include "../withinhost/MinimalParamset.h"

#include "../gillespie/FunctionModel.h"
#include "../random/xoroshiro128plus.h"
#include "../gillespie/AdaptivePRC.h"
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>

#define I_WT 1
#define I_MUT 2
#define MAX_ITER 50000000
#define PRINT_FREQ 10000

#define PARAM_N_IGA params.get_parameter("n_encounters_iga")
#define PARAM_BN params.get_integer_parameter("bottleneck")
#define TRANS_CUTOFF params.get_parameter("trans_cutoff")
#define EMERGE_THRESH params.get_parameter("emergence_threshold")

using std::ofstream;
using std::ostringstream;
using std::string;


void print_header(ostream& outstream){
  outstream << std::left << setw(15) << "paramset_id";
  outstream << std::left << setw(10) << "rng_seed";
  outstream << std::left << setw(18) << "p_mut_inoc";
  outstream << std::left << setw(10) << "Vw_inoc";
  outstream << std::left << setw(10) << "Vm_inoc";
  outstream << std::left << setw(15) << "is_inoc";
  outstream << std::left << setw(15) << "transmissible";
  outstream << std::left << setw(18) << "peak_mutant_freq";
  outstream << std::left << setw(16) << "emergence_time";
  outstream << std::endl;
}

void print_param_header(ostream& outstream,
                        vector<string>& param_names){
  outstream << std::left << setw(20) << "paramset_id";
  
  for (auto &param_name : param_names)
    outstream << std::left << setw(20) << param_name;

  outstream << std::endl;
  
}


void print_param_values(ostream& outstream,
                        int id,
                        vector<string>& param_names,
                        MinimalParamset& paramset){

  outstream << std::left << setw(20) << id;
  
  for (auto &param_name : param_names)
    outstream <<
      std::left << setw(20) <<
      setprecision(9) <<
      paramset.get_parameter(param_name);

  outstream << std::endl;
}



int main(int argc, char *argv[]){
  if(argc < 8){
    printf("\nUSAGE: ./MinimalWhSensitivity "
           "<path to inits> "
           "<model name> "
           "<n iter> "
           "<n paramsets> "
           "<bottleneck> "
           "<output path stem> "
           "<threshold> "
           "\n\n"
           "path to inits: an .ini/.mk formatted initialization file "
           "containing the default parameter values for the minimal "
           "model and possibly model-specific overrides\n\n"
           "model name: name of the model "
           "(referenced in param overrides)\n\n "
           "n iter: number of runs to do with each parameter set\n\n"
           "n paramsets: number of random paramsets to generate\n\n"
           "bottleneck: bottleneck size\n\n"
           "output path stem: stem of path to output results\n\n"
           "threshold: whether to use a threshold model or a "
           "proportional transmission probability model\n\n"
           );

    return 0;
  }

  // parse arg variables
  string input_path = argv[1];
  string model_name = argv[2];
  int n_iter = atoi(argv[3]);
  int n_paramsets = atoi(argv[4]);
  int bn = atoi(argv[5]);
  string out_stem = argv[6];
  bool threshold = atoi(argv[7]);


  // get default simulation parameters
  // within host model


  MinimalParamset params(input_path,
                         model_name);
  
  params.set_parameter("BOTTLENECK",
                       bn);

  // calculate mutant neut prob from z_mut, not from sigma
  params.set_parameter("USE_Z_MUT", 1);

  xoroshiro128plus model_rng(0);
  int current_seed = 0;
  // make hypercube reproducible
  current_seed  = params.get_integer_parameter("rng_seed") +
      n_iter * bn + 73;
    model_rng.set_seed(current_seed);


  vector<string> hypercube_params =
    {"",
     "R0_WH",
     "R_W",
     "MU",
     "D_V",
     "K",
     "CROSS_IMM",
     "T_M",
     "T_N",
     "C_MAX",
     "Z_WT",
     "Z_MUT_Z_WT_RATIO",
     "VB_RATIO"};

  if(threshold){
    std::cout << "using a transmission threshold..." << std::endl;
    hypercube_params.at(0) = "TRANS_THRESHOLD";
  }
  else{
    std::cout << "using a transmission probability "
      "proportional to viral load..." << std::endl;
    hypercube_params.at(0) = "TD50";
  }
  double trans_param;
                       
  params.setup_hypercube(hypercube_params,
                         n_paramsets,
                         model_rng);
  // setup output
  ostringstream out_path;
  out_path << out_stem << ".txt";
  
  ostringstream param_path;
  param_path << out_stem << "_paramsets.txt";

  ofstream outstream;  
  outstream.open(out_path.str());
  print_header(outstream);

  ofstream param_stream;  
  param_stream.open(param_path.str());
  print_param_header(param_stream,
                     hypercube_params);

  // make individual runs reproducible
  current_seed = params.get_integer_parameter("rng_seed") +
      n_iter * bn + 74;
  model_rng.set_seed(current_seed);

  for(int hyper_ind = 0; hyper_ind < n_paramsets; hyper_ind++){
    params.hypercube_next();
    params.set_parameter("R_M",
                         params.get_parameter("R_W"));
    // here force r_w = r_m (no deleteriousness)

    if(threshold)
      trans_param = params.get_parameter("trans_threshold");
    else
      trans_param = params.get_parameter("TD50");
    
    // save paramset values
    print_param_values(param_stream,
                       hyper_ind,
                       hypercube_params,
                       params);
    
    // show parameter set values, for debugging
    params.print_paramset(std::cout);
    
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
                                    0.0,
                                    0.0);

  
  FunctionModel model(model_functor.n_vars(),
                      &model_functor);

  // set var names
  model.set_var_name(0, "C");
  model.set_var_name(1, "Vw");
  model.set_var_name(2, "Vm");
  // set initial state
  vector<double> inoculator_init_state =
    {params.get_parameter("C_max"),
     params.get_integer_parameter("bottleneck")
     * 1.0,
     0.0};
  
  AdaptivePRC model_state(&model,
                          &model_rng,
                          params.get_parameter("leap_timestep"));


  // prep for output
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
    model_functor.E_w_ = 0.0;
    transmissible_inoculator = false;

    while(!transmissible_inoculator){
      model_state.reset_all();
      model_state.set_state(inoculator_init_state);
      model_state.set_time(0);
      model_state.simulate(params.get_parameter("t_final"),
                           MAX_ITER,
                           false);
      transmissible_inoculator =
        model_functor.ever_transmissible(model_state,
                                         TRANS_CUTOFF,
                                         trans_param);
    };

    if(i_run < 1){
      string timecourse_path = out_stem + "_naive.txt";
      model_state.output_all_to_file(timecourse_path);
    }

    double p_mut =
      model_functor.get_random_transmissible_pmut(model_state,
                                                  trans_param,
                                                  model_rng);

    // simulate an immune recipient host
    model_state.reset_all();
    model_functor.E_w_ = 1.0;
    model_functor.E_m_ = 0.0;

    // immunity acts iff immune to wildtype
    // assumed immune to mutant
    double wt_neut = model_functor.E_w_ *
      params.get_parameter("wt_neut_prob");
    double mut_neut = model_functor.E_w_ *
      params.get_parameter("mut_neut_prob");
    
    vector<double> inoculum =
      model_functor.get_inoculum(p_mut,
                                 PARAM_N_IGA,
                                 PARAM_BN,
                                 wt_neut,
                                 mut_neut,
                                 model_rng);

    double Vw_init = inoculum.at(0);
    double Vm_init = inoculum.at(1);
    double peak_mutant_freq;
    bool transmissible;
    bool inoc;
    double emergence_time = 99;
    
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
      model_state.set_state(init_state);
      model_state.set_time(0);

      model_state.simulate(params.get_parameter("t_final"),
                           MAX_ITER,
                           false);

      transmissible = false;  

      peak_mutant_freq =
        model_functor.get_mutant_peak_freq(model_state,
                                           TRANS_CUTOFF,
                                           trans_param);

      emergence_time =
        get_emergence_time(model_state,
                           EMERGE_THRESH);

      transmissible = model_functor.ever_transmissible(model_state,
                                                       TRANS_CUTOFF,
                                                       trans_param);

    }
    if(inoc && transmissible && !inoc_visible && peak_mutant_freq > 0.5){
      string timecourse_path = out_stem + "_inoc_visible.txt";
      model_state.output_all_to_file(timecourse_path);
      bool inoc_visible = true;
    }
    else if(inoc && !transmissible && !inoc_invisible){
      string timecourse_path = out_stem + "_inoc_invisible.txt";
      model_state.output_all_to_file(timecourse_path);
      bool inoc_invisible = true;
    }
    else if(!inoc && transmissible && !repl_visible){
      string timecourse_path = out_stem + "_repl_visible.txt";
      model_state.output_all_to_file(timecourse_path);
      bool repl_visible = true;
    }
    else if(!inoc && !transmissible && !repl_invisible){
      string timecourse_path = out_stem + "_repl_invisible.txt";
      model_state.output_all_to_file(timecourse_path);
      bool repl_invisible = true;
    }

    
    // print out the result
    outstream << std::left << setw(15) << hyper_ind;
    outstream << setw(10) << current_seed;
    outstream << setw(18) << setprecision(6) << p_mut;
    outstream << setw(10) << Vw_init;
    outstream << setw(10) << Vm_init;
    outstream << left << setw(15) << boolalpha << inoc;
    outstream << left << setw(15) << boolalpha << transmissible;
    outstream << left << setw(18) << setprecision(6) << peak_mutant_freq;
    outstream << left << setw(16) << emergence_time;
    outstream << std::endl;
    
    // reset and increment seed
    model_state.reset_all();
    current_seed++;
    model_rng.set_seed(current_seed);
    if(i_run % PRINT_FREQ == 0){
      std::cout << "Iteration no.: " << i_run << " of "
                << n_iter << std::endl;
    }
  }
  }

  outstream.close();
  param_stream.close();
  return 0;
}
