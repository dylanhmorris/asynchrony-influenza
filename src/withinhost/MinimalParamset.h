#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <iostream> // for std::ofstream
#include "../ini/INIReader.h"
#include "../random/rng.h"

using std::ofstream;
using std::ostream;
using std::vector;
using std::string;

class MinimalParamset {
public:
  MinimalParamset(string input_path,
                  string input_model_name) :
    reader(input_path)
  {

    std::transform(input_model_name.begin(),
                   input_model_name.end(),
                   input_model_name.begin(),
                   ::toupper);

    model_name = input_model_name;
    
    if (reader.ParseError() < 0) {
      std::cout << "\nCan't load default ini file " <<
        input_path << "\n";
    };
    
    
    double_parameters =
      {"R0_WH",
       "R_W",
       "R_M",
       "MU",
       "D_V",
       "K",
       "T_M",
       "T_N",
       "F",
       "P_C",
       "Z_WT",
       "Z_MUT_Z_WT_RATIO",
       "MUT_WT_NEUT_RATIO",
       "C_MAX",
       "CROSS_IMM",
       "VB_RATIO",
       "TRANS_THRESHOLD",
       "TRANS_CUTOFF",
       "TD50",
       "T_FINAL",
       "P_M",
       "LEAP_TIMESTEP",
       "INOCULATOR_LEAP_TIMESTEP",
       "EMERGENCE_THRESHOLD",
       "CHAIN_MUT_MEAN",
       "CHAIN_MUT_SD",
       "CHAIN_INIT_PHEN",
       "CHAIN_CONTACT_RATE",
       "CHAIN_ANTIGENIC_SCALE"};

    int_parameters = {"BOTTLENECK",
                      "RNG_SEED",
                      "E_W",
                      "E_M",
                      "USE_Z_MUT"};

    vector_parameters =
      {"IMMUNE_HISTORIES",
       "IMMUNE_HISTORY_DISTRIBUTION"
      };

    calculated_parameters =
      {"BETA",
       "WT_NEUT_PROB",
       "MUT_NEUT_PROB",
       "N_ENCOUNTERS_IGA"
      };

    initialize_maps();    
    parse_input("DEFAULT");
    parse_input(model_name);
    calculate_calculated_params();

  };

  ~MinimalParamset(){};


  double get_parameter(string parameter_name){
    std::transform(parameter_name.begin(),
                   parameter_name.end(),
                   parameter_name.begin(),
                   ::toupper);
    if (current_double_values.count(parameter_name) > 0)
      return current_double_values.at(parameter_name);    
    else {
      std::cout << "requested parameter: "
                << parameter_name << std::endl;
      throw std::runtime_error("could not find requested "
                               "double parameter");      
    };
  };

  int get_integer_parameter(string parameter_name){
    std::transform(parameter_name.begin(),
                   parameter_name.end(),
                   parameter_name.begin(),
                   ::toupper);

    if (current_int_values.count(parameter_name) > 0)
      return current_int_values.at(parameter_name);
    else {
      std::cout << "requested parameter: "
                << parameter_name << std::endl;
      throw std::runtime_error("could not find requested "
                               "integer parameter");      
    };
  };
  
  vector<double> get_vector_parameter(string parameter_name){
    if (current_vector_values.count(parameter_name) > 0)
      return current_vector_values.at(parameter_name);
    else {
      std::cout << "requested parameter: "
                << parameter_name << std::endl;
      throw std::runtime_error("could not find requested "
                               "integer parameter");      
    };
  };

  
  void print_paramset(ostream& outstream){
    for (auto &parameter_name : double_parameters)
      outstream << parameter_name << ": " <<
        current_double_values.at(parameter_name) << std::endl;
    
    for (auto &parameter_name : calculated_parameters)
      outstream << parameter_name << ": " <<
        current_double_values.at(parameter_name) << std::endl;

    for (auto &parameter_name : int_parameters)
      outstream << parameter_name << ": " <<
        current_int_values.at(parameter_name) << std::endl;
    
    for (auto &parameter_name : vector_parameters){
      vector<double> vec_to_print =
        current_vector_values.at(parameter_name);
      
      outstream << parameter_name << ": ";

      for(auto &param_val : vec_to_print)
        outstream << param_val << " ";
      
      outstream << std::endl;
    };
  };

  void setup_hypercube(vector<string> parameters,
                       int n_samples,
                       RNG &rng){
    hypercube_parameters = parameters;

    int n_hypercube_parameters = parameters.size();
    hypercube_samples.clear();
    
    // for each param, generate a random vector with
    // n_samples cutpoints along the unit interval
    for (auto &parameter_name : hypercube_parameters)
      hypercube_samples[parameter_name] = rng.runif_latin(n_samples);
    
    hypercube_samples_remaining = n_samples;
  };

  void hypercube_next(){
    double min_val = 0;
    double max_val = 0;

    if(hypercube_samples_remaining < 1)
      throw std::runtime_error("no hypercube samples remaining");      
    
    for (auto &parameter_name : hypercube_parameters){
      min_val = min_and_max_values.at(parameter_name).at(0);
      max_val = min_and_max_values.at(parameter_name).at(1);
      double sample = hypercube_samples.at(parameter_name).back();
      hypercube_samples.at(parameter_name).pop_back();
      set_parameter(parameter_name,
                    min_val + (max_val - min_val) * sample);
    };

    hypercube_samples_remaining--;
  };

  void initialize_maps(){
    for (auto &parameter_name : double_parameters){
      current_double_values[parameter_name] = -99.0;
      min_and_max_values[parameter_name] = {-99, -99};
    }
    
    for (auto &parameter_name : calculated_parameters)
      current_double_values[parameter_name] = -99.0;

    for (auto &parameter_name : int_parameters)
      current_int_values[parameter_name] = -99;

    for (auto &parameter_name : vector_parameters)
      current_vector_values[parameter_name] = {-99};
  };
  
  void parse_input(string prefix){

    double min_val = 0;
    double max_val = 0;
    
    for ( auto &parameter_name : double_parameters){
      current_double_values.at(parameter_name) =
        reader.GetReal("",
                       prefix + "_" + parameter_name,
                       current_double_values.at(parameter_name));

      min_val =
        reader.GetReal("",
                       prefix + "_" + parameter_name +
                       "_MIN",
                       min_and_max_values.at(parameter_name).at(0));
      max_val =
        reader.GetReal("",
                       prefix + "_" + parameter_name +
                       "_MAX",
                       min_and_max_values.at(parameter_name).at(1));

      if (max_val < min_val) {
        std::cout << "Attempt to set max value (" <<
          max_val << ") " <<
          "greater than min value (" <<
          min_val << ") " <<
          "for parameter " <<
          parameter_name << " and prefix " <<
          prefix << std::endl;
      
        throw std::runtime_error("max value for parameter "
                                 "greater than minimum value");
      }

      min_and_max_values.at(parameter_name) =
        {min_val, max_val};
    };

    for ( auto &parameter_name : int_parameters){
      current_int_values.at(parameter_name) =
        reader.GetReal("",
                       prefix + "_" + parameter_name,
                       current_int_values.at(parameter_name));
    };
    
    for ( auto &parameter_name : vector_parameters){
      current_vector_values.at(parameter_name) =
        reader.GetRealVector("",
                             prefix + "_" + parameter_name,
                             current_vector_values.at(parameter_name));
    };
    
    calculate_calculated_params();
  };
  

  void calculate_calculated_params(){

    current_double_values.at("BETA") =
      get_parameter("R0_wh") *
      get_parameter("d_v") /
      (get_parameter("C_max") *
       get_parameter("r_w"));

    current_double_values.at("N_ENCOUNTERS_IGA") =
      get_parameter("VB_RATIO") * get_integer_parameter("BOTTLENECK");
    
    double zwt = get_parameter("Z_WT");
    double neut_ratio = get_parameter("MUT_WT_NEUT_RATIO");
    double v = get_parameter("N_ENCOUNTERS_IGA");
    double candidate = 0;

    if(zwt > 0){
      candidate = neutralize_prob_from_z(zwt);
      if (candidate > 1)
        throw std::runtime_error("impossible wt "
                                 "neutralization prob (> 1)");
      if (candidate < 0)
        throw std::runtime_error("impossible wt "
                                 "neutralization prob (< 0)");
    }
    
    current_double_values.at("WT_NEUT_PROB") = candidate;

    int use = get_integer_parameter("USE_Z_MUT");
    if(use > 0){
      candidate = 0;
      double zmut = get_parameter("Z_MUT_Z_WT_RATIO") * zwt;
      if(zmut > 0){
        candidate = neutralize_prob_from_z(zmut);
        if (candidate > 1)
          throw std::runtime_error("impossible mut "
                                   "neutralization prob (> 1)");
        if (candidate < 0)
          throw std::runtime_error("impossible mut "
                                   "neutralization prob (< 0)");
      }
      current_double_values.at("MUT_NEUT_PROB") = candidate;
    } else {
      current_double_values.at("MUT_NEUT_PROB") =
        current_double_values.at("WT_NEUT_PROB") * neut_ratio;
    }

  };
  
  void set_parameter(string parameter_name,
                     double new_value){
    std::transform(parameter_name.begin(),
                   parameter_name.end(),
                   parameter_name.begin(),
                   ::toupper);

    current_double_values.at(parameter_name) = new_value;
    calculate_calculated_params();
  };

  void set_parameter(string parameter_name,
                     int new_value){
    std::transform(parameter_name.begin(),
                   parameter_name.end(),
                   parameter_name.begin(),
                   ::toupper);

    current_int_values.at(parameter_name) = new_value;
    calculate_calculated_params();
  };

  void set_parameter(string parameter_name,
                     vector<double> new_value){
    std::transform(parameter_name.begin(),
                   parameter_name.end(),
                   parameter_name.begin(),
                   ::toupper);

    current_vector_values.at(parameter_name) = new_value;
    calculate_calculated_params();
  };

  double neutralize_prob_from_z(double desired_z){
    double v = get_parameter("VB_RATIO") *
      get_integer_parameter("BOTTLENECK");
    double prob_nonzero_naive = 1 - exp(-v);
    double target = 1 - (prob_nonzero_naive * (1 - desired_z));
    return (1 + (log(target) / v));
  }

  string model_name;

private:
  vector<string> double_parameters;
  vector<string> int_parameters;
  vector<string> vector_parameters;
  vector<string> calculated_parameters;

  vector<string> hypercube_parameters;

  std::map<string, double> current_double_values;
  std::map<string, int> current_int_values;
  std::map<string, vector<double>> current_vector_values;
  
  std::map<string, vector<double>> min_and_max_values;
  std::map<string, vector<double>> hypercube_samples;
  int hypercube_samples_remaining;
  
  INIReader reader;
  
};
