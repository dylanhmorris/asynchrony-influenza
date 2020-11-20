/** 
 *  @file    EventRecord.h
 *  @author  Dylan Morris (dhmorris@princeton.edu)
 *  @date    2017-06-19 
 *  @version 0.1 
 *  
 *  @brief Class EventRecord holds a record of MultistrainSIR events
 *  
 */

#ifndef EVENT_RECORD_H_

#define EVENT_RECORD_H_

#define TRANSMISSION_ID 0
#define RECOVERY_ID 1
#define NEW_MUTANT_ID 2
#define NEW_MAJORITY_ID 3
#define NULL_HOST_ID -99

#include <vector>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;


class EventRecord {
 public:
 EventRecord(double time,
             int event_type,
             int host,
             int secondary_host,
             vector<int> &deltas) :
  time(time),
  event_type(event_type),
  host(host),
  secondary_host(secondary_host),
  deltas(deltas)
  {};
  
  
  ~EventRecord(){};

  void print_record(){

    cout << left << setw(15) << setprecision(12) << time;
      cout << left << setw(15) << get_event_type_name();
      cout << left << setw(15) << host;
      cout << left << setw(15) << secondary_host;
      for(auto &delta : deltas)
        cout << left << setw(15) << delta;
      cout << endl;
  };

  string get_event_type_name(){
  if(event_type == TRANSMISSION_ID)
    return "transmission";
  else if(event_type == RECOVERY_ID)
    return "recovery";
  else if(event_type == NEW_MUTANT_ID)
    return "new_mutant";
  else if(event_type == NEW_MAJORITY_ID)
    return "new_majority";
  else return "NO_ID";
  };

  double time;
  int event_type;
  int host;
  int secondary_host;
  vector<int> deltas;
    
};

#endif
