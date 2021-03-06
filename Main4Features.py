#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 22:04:49 2021

@author: fujita
"""

import os
import lib.mainHelper as mH

if __name__ == '__main__':
    
    file_dir = './Data/'
    save_dir = './Result/'
    if not os.path.isdir(save_dir):
         os.makedirs(save_dir) 
    ### Preprocess
    Data=mH.Preprocess(file_dir)
    OGTThuman_class=mH.OGTThuman(Data,save_dir)
    OGTThuman_class.TwodData = Data
    ###OGTThuman_class.AnalFigS3(Data,file_dir,save_dir)

    ###OGTThuman_class.AnalFig2_3()
    
    ###OGTThuman_class.AnalFig4()
    ###OGTThuman_class.AnalFig5()
    ##OGTThuman_class.AnalFig6()
    ###OGTThuman_class.AnalFig7()
    
    OGTThuman_class.AnalFig8FigS9()    
    
    ###OGTThuman_class.AnalSTable7()
    ###OGTThuman_class.AnalAddFig1(file_dir )
    ###OGTThuman_class.AnalFigS2(file_dir)
        
        
