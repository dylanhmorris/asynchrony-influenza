#!/usr/bin/env python3

# gets a template with latex tags

import jinja2
import functools
import os
import sys
import tex
import re

def get_template(template_filepath):
    """Get a jinja template with latex tags per
    http://eosrei.net/articles/2015/11/latex-templates-python-and-jinja2-generate-pdfs
    """
    path, filename = os.path.split(template_filepath)
    latex_jinja_env = jinja2.Environment(
    	block_start_string = '\BLOCK{',
    	block_end_string = '}',
    	variable_start_string = '\VAR{',
    	variable_end_string = '}',
    	comment_start_string = '\#{',
    	comment_end_string = '}',
    	line_statement_prefix = '%%',
    	line_comment_prefix = '%#',
    	trim_blocks = True,
    	autoescape = False,
    	loader = jinja2.FileSystemLoader(path)
    )
    template = latex_jinja_env.get_template(filename)
    return template

