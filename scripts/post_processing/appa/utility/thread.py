#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""

import threading

class new_thread(threading.Thread):

    def __init__(self,target = None,args=None):
        super(new_thread, self).__init__(target=target,args=args,)
        self._stop =threading.Event()

    def stop(self):
        self._stop.set()
        

    def stoped(self):
        self._stop.isSet()
