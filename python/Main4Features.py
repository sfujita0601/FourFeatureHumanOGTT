#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 22:04:49 2021

@author: fujita
"""

import numpy as np
import pandas as pd
import scipy.io
import itertools
import sys
import matplotlib as mpl
import lib
import os
import lib.general as ge
import lib.mainHelper as mH

if __name__ == '__main__':
    
    file_dir = './Data/'
    save_dir = './Result/'
    if not os.path.isdir(save_dir):
         os.makedirs(save_dir) 
    ### Preprocess
    Data=mH.Preprocess(file_dir)


    OptionDict={'Anal':1}

            
    if OptionDict['Anal']==1:
        OGTThuman_class=mH.OGTThuman(Data,save_dir)
        #OGTThuman_class.TwodData = Data
        #OGTThuman_class.AnalFigS3(Data,file_dir,save_dir)#空腹値とのt検定_各時点で
        #OGTThuman_class.Clustering(Data,save_dir)
        #OGTThuman_class.AnalTVIR()
        #OGTThuman_class.AnalTPSI()
        #OGTThuman_class.AnalFig2_3()
        #OGTThuman_class.AnalFig4()
        OGTThuman_class.AnalFig5()
        #OGTThuman_class.AnalFig6()
        #OGTThuman_class.AnalFig7()
        #OGTThuman_class.AnalFig8FigS9()    
        #OGTThuman_class.AnalSTable7()
        #OGTThuman_class.AnalAddFig1(file_dir )
        #OGTThuman_class.AnalFigS2(file_dir)
        
        
