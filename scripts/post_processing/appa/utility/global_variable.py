#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""

import sys,os

#Global variable


#path
global_path = os.getcwd()
def path():
    if os.path.exists(global_path) :
        if os.path.isdir(global_path):
            return str(global_path)
        else:
            return str(os.path.dirname(global_path))
    else:
        return "~/"
    

#version
version = "1.0.8"
