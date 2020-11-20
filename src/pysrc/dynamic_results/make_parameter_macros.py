#!/usr/bin/env python3

# Takes Makefile, escapes, passes to Jinja2

from distutils.sysconfig import parse_makefile
import jinja2
from jinja_tools import get_template
import os
import sys
import csv

def underscores_to_caps(word):
    """
    Replace underscores with capitalization
    (for taking env variables to LaTeX macros)
    """
    return ''.join(x.capitalize() for x in word.split('_'))

def digits_to_names(word):
    trans_dict = {"0": "nought",
                  "1": "one",
                  "2": "two",
                  "3": "three",
                  "4": "four",
                  "5": "five",
                  "6": "six",
                  "7": "seven",
                  "8": "eight",
                  "9": "nine"}
    return word.translate(str.maketrans(trans_dict))

def add_prefix_and_suffix(word, prefix, suffix):
    if prefix is not None:
        word = prefix + word
    if suffix is not None:
        word = word + suffix
    return word

def escape_word(word, prefix=None, suffix=None):
    word = digits_to_names(word)

    word = add_prefix_and_suffix(word,
                                 prefix,
                                 suffix)
    word = underscores_to_caps(word)
    return word

def escape_contents(input_dict, prefix=None):
    return {escape_word(key, prefix=prefix): entry
            for key, entry in input_dict.items()}

def isnum(param_value):
    """
    Function to separate out numerical 
    parameter values from others, as the
    template will render numerical values as
    LaTeX numbers (\num{value})
    """
    if type(param_value) is str:
        return param_value.replace(
            '.', '').replace(
            'e', '').replace(
            '-', '').isdigit()
    else:
        return True

def dict_to_list(input_dict):
    return [{"name": key, "value": value,
             "isnum": isnum(value)}
            for key, value
            in input_dict.items()]

def main(source, template, outpath,
         prefix=None):
    contents = parse_makefile(source)
    contents = escape_contents(contents,
                               prefix=prefix)
    contents = dict_to_list(contents)
    print(contents)
    template = get_template(template)
    with open(outpath, 'w') as output:
        output.write(
            template.render(
                envvars=contents
                )
            )
    return 0

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Supply source file, template file, output path, in that order"
              "optionally supply prefix to append to all macros")
    elif len(sys.argv) == 4:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
    elif len(sys.argv) > 4:
        main(sys.argv[1],
             sys.argv[2],
             sys.argv[3],
             sys.argv[4])
