#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 27.05.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

import sys


from python_esa.sa import reduce_sa

if __name__ == '__main__':
    
    if len(sys.argv)<6:
        print "wrong params"
        sys.exit()
    
    params = {}
    params["tf_cutoff"] = sys.argv[3]
    params["df_cutoff"] = sys.argv[4]
    params["length_cutoff"] = sys.argv[5]
    
    reduce_sa(sys.argv[1], sys.argv[2], params)