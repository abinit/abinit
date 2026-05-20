#!/usr/bin/env python

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""

import os

#Global variable


#path
global_path = os.getcwd()
def path():
    if os.path.exists(global_path) :
        if os.path.isdir(global_path):
            return str(global_path)
        return str(os.path.dirname(global_path))
    return "~/"


#version
version = "1.0.8"
