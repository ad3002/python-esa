#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 27.05.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

import sys

from python_esa.sa import clean_suffix_array_data

if __name__ == '__main__':
    
    if len(sys.argv)<2:
        print "wrong params"
        sys.exit()
    
    clean_suffix_array_data(sys.argv[1])