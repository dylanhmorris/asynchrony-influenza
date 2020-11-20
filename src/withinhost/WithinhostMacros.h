/** 
 *  @file    WithinhostMacros.h
 *  @author  Dylan Morris (dhmorris@princeton.edu)
 *  
 *  @brief WithinhostMacros.h defines macros to 
 *  control array indexing, etc, throughout the 
 *  within-host models 
 *  
 */


#ifndef WITHIN_HOST_MACROS_H_
#define WITHIN_HOST_MACROS_H_

// Minimal Model
#define STATE_LENGTH 3
#define INOCULUM_LENGTH 2

#define C_PLUS 0
#define C_MINUS 1
#define VW_PLUS 2
#define VW_MINUS 3
#define VM_PLUS 4
#define VM_MINUS 5

#define C_IND 0
#define VW_IND 1
#define VM_IND 2


// Transmission Chain model
#define VMAJ_PLUS 2
#define VMAJ_MINUS 3
#define VMIN_PLUS 4
#define VMIN_MINUS 5

#define VMAJ_IND 1
#define VMIN_IND 2

#define TITER_ALPHA 2.844
#define TITER_BETA 1.299


#endif
