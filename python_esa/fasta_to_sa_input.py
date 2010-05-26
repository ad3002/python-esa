#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 26.05.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

import sys


from python_esa.sa import fasta_to_sa_input

if __name__ == '__main__':
    
    if len(sys.argv)<2:
        print "wrong params"
        sys.exit()
    
    fasta_to_sa_input(sys.argv[1], "input.fa")