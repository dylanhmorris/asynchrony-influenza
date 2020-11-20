#!/usr/bin/env python3

from distutils.sysconfig import parse_makefile


def get_params(parameter_filepath):        
    filepath = parameter_filepath
    print("Reading ", filepath)
    return parse_makefile(parameter_filepath)

def get_param(param_name=None,
              model_name=None,
              param_dict=None,
              param_key = None,
              default=None,
              verbose=False):
    if default is None:
        default = "DEFAULT_" + param_name
    model_param = model_name + "_" + param_name
    
    default_key = default.upper()
    if param_key is not None:
        model_param_key = param_key
    else:
        model_param_key = model_param.upper()

    if verbose:
        print("Getting parameter {}...".format(model_param_key))
              
    param_value = param_dict.get(model_param_key,
                                 param_dict[default_key])
    return float(param_value)
