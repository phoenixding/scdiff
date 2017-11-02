#!/usr/bin/env python
import pdb,sys,os
"""
author: Jun Ding
usage: scdiff running command example 
scdiff -i example.E -t example.tf_dna -k auto -o example_out
"""

# run with automatic config 
os.system("scdiff -i example.E -t example.tf_dna -k auto -o example_out")
