/** 
 *  @file    ChainModel.h
 *  @author  Dylan Morris (dhmorris@princeton.edu)
 *  
 *  @brief ChainModel.h defines a runner class of
 *  for running transmission chain simulations
 *  
 */

#ifndef CHAIN_MODEL_H_
#define CHAIN_MODEL_H_


#include "../withinhost/ChainWithinHost.h"
#include "../withinhost/MinimalParamset.h"

#include "../gillespie/FunctionModel.h"
#include "../random/xoroshiro128plus.h"
#include "../gillespie/AdaptivePRC.h"
#include <fstream>


#define MAX_ITER 50000000
#define PARAM_N_IGA params.get_parameter("n_encounters_iga")
#define PARAM_BN params.get_integer_parameter("bottleneck")
#define KILLING  params.get_parameter("z_wt") > 0
#define CONTACT_RATE params.get_parameter("chain_contact_rate")
using std::ofstream;
using std::ostream;



class ChainModel {
  

public:
  xoroshiro128plus model_rng_;
  ChainWithinHost model_functor_;
  FunctionModel model_;
  AdaptivePRC host_state_;

  // class attributes
  MinimalParamset params;
  vector<double> imm_hist_dist;
  vector<vector<double>> imm_hists;
  double init_phen;
  double mut_mean;
  double mut_sd;
  vector<double> default_init_state;
  vector<double> current_init_state;
  vector<double> current_inoculum;
  double transmission_parameter;

  vector<vector<double>> current_donor_states;
  vector<double> current_donor_times;
  double current_donor_clearance_time;

  int current_recipient_hist_choice;
  int current_donor_hist_choice;

  int current_seed;
  int current_chain_no;
  int i_chain; // iteration within chain
  
  bool new_maj;
  bool fix;

  ofstream outstream;

  // default constructor
  ChainModel(MinimalParamset &params,
             string outpath) :
    params(params),
    model_rng_(1),
    init_phen(0),
    mut_mean(0),
    mut_sd(0),
    new_maj(false),
    fix(false),
    current_donor_hist_choice(-1),
    current_recipient_hist_choice(-1),
    current_seed(-1),
    current_chain_no(0),
    i_chain(0),
    model_functor_(params.get_parameter("beta"),
                   params.get_parameter("r_w"),
                   params.get_parameter("mu"),
                   params.get_parameter("d_v"),
                   params.get_parameter("k"),
                   params.get_parameter("t_M"),
                   params.get_parameter("t_N"),
                   params.get_parameter("C_max"),
                   params.get_parameter("p_M"),
                   params.get_parameter("f"),
                   params.get_parameter("p_C")),
    model_(model_functor_.n_vars(),
           &model_functor_),
    host_state_(&model_,
                &model_rng_,
                params.get_parameter("leap_timestep"))
  {

    // output setup
    outstream.open(outpath);
    print_header();

  // set up model given parameters

    model_functor_.set_homotypic_protection(params.get_parameter("z_wt"));

    // by default, use threshold transmission;
    // can also set after setup to use
    // proportional transmission
    use_threshold_transmission_model();
    
    std::cout << "Setting up model with the following parameters..."
              << std::endl;
    params.print_paramset(std::cout);

    //
  
    // set var names
    model_.set_var_name(0, "C");
    model_.set_var_name(1, "V_maj");
    model_.set_var_name(2, "V_min");
    

    // phenotype setup
    init_phen = params.get_parameter("chain_init_phen");
    mut_mean = params.get_parameter("chain_mut_mean");
    mut_sd = params.get_parameter("chain_mut_sd");

    // state and inoculum setup
    reset_current_inoculum();
    setup_init_state();
    set_current_init_state_to_default();
    
    // immune history setup
    setup_immunity();
  }


  // default destructor
  ~ChainModel(){}


  // print the header for output
  void print_header(){
    outstream << std::left << setw(12) << "chain_no";
    outstream << std::left << setw(12) << "rng_seed";
    outstream << std::left << setw(12) << "link_no";
    outstream << std::left << setw(8) <<  "V_maj";
    outstream << std::left << setw(8) <<  "V_min";
    outstream << std::left << setw(16) << "maj_phen";
    outstream << std::left << setw(16) << "min_phen";
    outstream << std::left << setw(20) << "donor_hist_choice";
    outstream << std::left << setw(20) << "recip_hist_choice";
    outstream << std::endl;
  }

  // print the current state of the system
  void print_chain_result(double V_maj,
                          double V_min,
                          double maj_phen,
                          double min_phen){
    outstream << std::left << setw(12) << current_chain_no;
    outstream << std::left << setw(12) << current_seed;
    outstream << std::left << setw(12) << i_chain;
    outstream << std::left << setw(8) << setprecision(2) << V_maj;
    outstream << std::left << setw(8) << setprecision(2) << V_min;
    outstream << std::left << setw(16) << setprecision(5) << maj_phen;
    outstream << std::left << setw(16) << setprecision(5) << min_phen;
    outstream << std::left << setw(20) << setprecision(2) << current_donor_hist_choice;
    outstream << std::left << setw(20) << setprecision(2) << current_recipient_hist_choice;
    outstream << std::endl;
  }

  
  // safely create an initial state, respecting canonical
  // indices set in TransmissionMacros.h
  void setup_init_state(){
    default_init_state.resize(STATE_LENGTH, 0);
    default_init_state.at(C_IND) = params.get_parameter("C_max");
    default_init_state.at(VMAJ_IND) =
      params.get_integer_parameter("bottleneck") * 1.0;
    default_init_state.at(VMIN_IND) = 0.0;
  }

  // set up immune histories from parameters
  void setup_immunity(){
    imm_hist_dist =
      params.get_vector_parameter("IMMUNE_HISTORY_DISTRIBUTION");
    
    for(auto hist : params.get_vector_parameter("IMMUNE_HISTORIES"))
      imm_hists.push_back({hist});

    imm_hists.push_back({});

    // check that everything worked, raise error if not
    if(imm_hists.size() != imm_hist_dist.size()) {
      std::cout << "number of immune histories: " << 
        imm_hists.size() << std::endl;
      std::cout <<
        "number of immune history population proportions: "
                << imm_hist_dist.size() << std::endl;
      throw std::runtime_error("Immune history distribution "
                               "does not match number of histories "
                               "plus default (naive) history");
    }

    model_functor_.set_possible_histories(imm_hists);
  }

  // reset current inoculum to zeros (to avoid
  // unexpected behavior if it is referenced
  // when it should not be
  void reset_current_inoculum(){
    current_inoculum.resize(INOCULUM_LENGTH, 0);
  }


  void reset_donor(){
    current_donor_states.clear();
    current_donor_times.clear();
    current_donor_clearance_time = 0;
  }
  
  void set_current_init_state_to_default(){
    current_init_state = default_init_state;
  }

  // get a mutant effect size rel to current phenotype
  // with a Gaussian mutation kernel
  double get_mutant_effect(){
    return model_rng_.rgaussian(mut_mean, mut_sd);
  }

 
  // simulation of individual infections
  void simulate_infection(bool verbose = false){

    // produce infection
    host_state_.reset_all();
    host_state_.set_state(current_init_state);
    host_state_.set_time(0);
    host_state_.simulate(params.get_parameter("t_final"),
                        MAX_ITER,
                        verbose);

    // check that infection actually finished
    double total_pop = get_viral_load(host_state_.state_array,
                                      {VMAJ_IND, VMIN_IND});
    if(! total_pop < 1e-1){
      std::cout << setw(8) << host_state_.state_time;
      std::cout << setw(15) << setprecision(5)
                << host_state_.state_array.at(0);
      std::cout << setw(15) << setprecision(5)
                << host_state_.state_array.at(1);
      std::cout << setw(15) << setprecision(5)
                << host_state_.state_array.at(2);
      std::cout << std::endl;

      throw std::runtime_error("Infection ongoing at end of "
                               "simulation. Check parameters "
                               "and consider increasing "
                               "t_final and max iterations "
                               "if parameters are sensible\n");
    }

    current_donor_states = host_state_.array_ts;
    current_donor_times = host_state_.time_ts;
    current_donor_clearance_time = host_state_.state_time;

  }

  // create the founder of a fresh subchain after a subchain breaks
  void create_fresh(){
    int hist_choice = imm_hist_dist.size() - 1;

    double maj_phen = init_phen;
    double min_phen = maj_phen +
      get_mutant_effect(); 

    model_functor_.clear_history();
    model_functor_.set_maj_phen(maj_phen);
    model_functor_.set_min_phen(min_phen);
    model_functor_.set_history(hist_choice);

    current_donor_hist_choice = hist_choice;
    
    set_current_init_state_to_default();
    reset_current_inoculum();
  
  }

  // create an infection from its inoculum, setting majority/minorty
  // phenotypes accordingly
  void create_from_current_inoculum(){
    double V_maj_old = current_inoculum.at(0);
    double V_min_old = current_inoculum.at(1);

    double V_maj_new = 0;
    double V_min_new = 0;
    double new_maj_phen = 0;
    double new_min_phen = 0;
    double old_maj_phen = model_functor_.get_maj_phen();
    double old_min_phen = model_functor_.get_min_phen();
    
    // if the majority stays the majority at transmission
    // set phenotypes accordingly
    if(V_maj_old > V_min_old){
      new_maj = false;
      fix = false;
      V_maj_new = V_maj_old;
      V_min_new = V_min_old;

      new_maj_phen = old_maj_phen;
      // if no minor variant transmitted,
      // choose new minor phenotype
      if(!(V_min_new > 0)) {
        new_min_phen = new_maj_phen + get_mutant_effect();
      } else {
        new_min_phen = old_min_phen;
      }
    }
    
    // otherwise, if the minority becomes the
    // majority (or coequal) at transmission,
    // set phenotypes accordingly
    else {
      new_maj = true;
      V_maj_new = V_min_old;
      V_min_new = V_maj_old;
      
      new_maj_phen = old_min_phen;

      if(V_min_new > 0){

        new_min_phen = old_maj_phen;
        fix = false;
        
      } else {
        
        new_min_phen = new_maj_phen + get_mutant_effect();
        fix = true;
      }
    }

    if(new_maj){
      print_chain_result(V_maj_new,
                         V_min_new,
                         new_maj_phen,
                         new_min_phen);
    }

    // update model functor and initial state
    model_functor_.set_maj_phen(new_maj_phen);
    model_functor_.set_min_phen(new_min_phen);

    current_init_state.at(C_IND) = params.get_parameter("C_max");
    current_init_state.at(VMAJ_IND) = V_maj_new;
    current_init_state.at(VMIN_IND) = V_min_new;

    model_functor_.set_history(current_recipient_hist_choice);
    current_donor_hist_choice = current_recipient_hist_choice;
    current_recipient_hist_choice = -1;
  }
  

  /* simulate all attempted transmissions from a single
   * simulated infection */
  void sim_transmission_attempts()
  {

    double contact_time = 0;
    vector<double> inoculum(INOCULUM_LENGTH, 0);
    
    // host attempts transmissions until
    // infection cleared
    bool successful_transmission = false;
    contact_time += model_rng_.rexp(CONTACT_RATE);
    
    while(!successful_transmission &
          contact_time < current_donor_clearance_time){

      vector<double> contact_state =
        last_state_before(contact_time,
                          current_donor_states,
                          current_donor_times);

      current_recipient_hist_choice =
        model_rng_.weighted_choice(imm_hist_dist);

      model_functor_.set_history(current_recipient_hist_choice);

      inoculum = model_functor_.get_inoculum(contact_state,
                                             transmission_parameter,
                                             PARAM_N_IGA,
                                             PARAM_BN,
                                             model_functor_.get_maj_phen(),
                                             model_functor_.get_min_phen(),
                                             model_rng_,
                                             KILLING);

      if(inoculum.at(0) + inoculum.at(1) > 0){
        successful_transmission = true;
      }
      
      contact_time += model_rng_.rexp(CONTACT_RATE);

    }

    current_inoculum = inoculum;
  }

  // main method for simulation of chains
  // outputting to a given outstream
  void simulate_chain(int seed,
                      int chain_no,
                      int max_chain_iter){
    i_chain = 0;
    current_seed = seed;
    current_chain_no = chain_no;

    new_maj = false;
    fix = false;
    
    model_rng_.set_seed(seed);
    
    // simulate first host in chain

    reset_donor();
    create_fresh();
    simulate_infection();
    
    print_chain_result(default_init_state.at(VMAJ_IND),
                       default_init_state.at(VMIN_IND),
                       model_functor_.get_maj_phen(),
                       model_functor_.get_min_phen());


    while(i_chain < max_chain_iter & !new_maj){

      i_chain++;

      sim_transmission_attempts();

      // if onward tranmission was achieved,
      // continue the subchain, else start a fresh
      // subchain
      if(current_inoculum.at(0) + current_inoculum.at(1) > 0){
        create_from_current_inoculum();
        simulate_infection();
      } else {
        create_fresh();
        simulate_infection();  
      }
    } // end loop over links in chain
  } // end simulation


  void use_threshold_transmission_model(){
    transmission_parameter = params.get_parameter("TRANS_THRESHOLD");
    model_functor_.threshold_model();
  }
  
  void use_proportional_transmission_model(){
    transmission_parameter = params.get_parameter("TD50");
    model_functor_.proportional_model();
  }

  
};  // end of definition for class Chain Model

 
#endif
