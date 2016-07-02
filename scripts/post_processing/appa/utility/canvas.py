#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvasQT
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

from matplotlib.figure import Figure
from matplotlib.ticker import FormatStrFormatter,ScalarFormatter
try:
    from PyQt4 import Qt,QtGui,QtCore
except:
    pass;

#---------------------------------------------------#
#---------------------------------------------------#
#------------------CANVAS(NOwindows)----------------#
#---------------------------------------------------#
#---------------------------------------------------#

class Canvas(FigureCanvas):

    def __init__(self, parent=None, width=6, height=4, dpi=100,x=0,y=0,pxlbl="",pylbl="",point = False,adjust=False):

        self.fig = Figure(figsize=(width, height), dpi=dpi,)

        FigureCanvas.__init__(self,self.fig)

        self.setPlot(x,y,pxlbl,pylbl,adjust=adjust)


    def setPlot(self,x,y,xlbl,ylbl,point =False, adjust=False):
        try :
            self.formater = ScalarFormatter(useOffset=True, useMathText=False, useLocale=None)
            self.formater.set_useOffset(0)
        except:
            self.formater = ScalarFormatter(useOffset=True, useMathText=False)

        self.axes = self.fig.add_subplot(111)
        self.axes.clear()

        if point :
            self.axes.plot(x,y,'.',markersize=15)
        else :
            self.axes.plot(x,y)

        try:
            if (max(y)-min(y)) <10**-10:
                self.axes.set_ylim(min(y)-1,max(y)+1)
        except:
            pass;
        
        self.axes.set_xlabel(xlbl)
        self.axes.set_ylabel(ylbl)
        self.axes.figure.set_facecolor('white')
        self.axes.grid('on')

        self.axes.yaxis.set_major_formatter(self.formater)
        self.axes.xaxis.set_major_formatter(self.formater)

        if adjust:
            self.adjust_x_lim(x,y)

        self.draw()

    def addLegend(self,plegend,markerscale=1):
        self.axes.legend(plegend,loc=1,markerscale=markerscale)
        self.draw()

    def addPlot(self,x,y,bar = False, point = False,marker='.',marker_size=25 ):
        if bar == True:
            self.axes.bar(x, y, width=0.01)
        elif point == True :
            self.axes.plot(x, y,marker,markersize=marker_size)
        else:
            self.axes.plot(x,y)

        self.draw()

    def adjust_x_lim(self, x, y):
        test = 0
        xmin = 0
        xmax = 0
        for i in range(len(x)):
            if y[i] >= 0.0002:
                xmin = x[i]
                break
        for i in range(len(x)):
            if y[i] <= 0.002 :
                test += 1
            if test > 60 :
                xmax = x[i-50]
                break
        self.axes.set_xlim(xmin,xmax)
        self.draw()

#---------------------------------------------------#
#---------------------------------------------------#
#------------------CANVAS(Windows)------------------#
#---------------------------------------------------#
#---------------------------------------------------#

class CanvasQT(FigureCanvasQT):

    def __init__(self, parent=None, width=6, height=4, dpi=100,x=0,y=0,pxlbl="",pylbl="",point = False,adjust=False,marker='.',marker_size=25):

        self.fig = Figure(figsize=(width, height), dpi=dpi,)

        FigureCanvasQT.__init__(self,self.fig)

        FigureCanvasQT.setSizePolicy(self,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding)
        FigureCanvasQT.updateGeometry(self)

        self.setPlot(x,y,pxlbl,pylbl,point =point,adjust=adjust,marker=marker,marker_size=marker_size)


    def setPlot(self,x,y,xlbl,ylbl,point =False, adjust=False,marker='.', marker_size=25):
        try :
            self.formater = ScalarFormatter(useOffset=True, useMathText=False, useLocale=None)
            self.formater.set_useOffset(0)
        except:
            self.formater = ScalarFormatter(useOffset=False, useMathText=False)

        self.axes = self.fig.add_subplot(111)
        self.axes.clear()

        if point :
            self.axes.plot(x,y,marker,markersize=marker_size)
        else :
            self.axes.plot(x,y)

        self.axes.set_xlabel(xlbl)
        self.axes.set_ylabel(ylbl)

        try :
            if (max(y)-min(y)) <10**-10:
                self.axes.set_ylim(min(y)-1,max(y)+1)
        except :
            pass;
        self.axes.figure.set_facecolor('white')
        self.axes.grid('on')

        self.axes.yaxis.set_major_formatter(self.formater)
        self.axes.xaxis.set_major_formatter(self.formater)

        if adjust:
            self.adjust_x_lim(x,y)

        self.draw()

    def addLegend(self,plegend,markerscale=1):
        self.axes.legend(plegend,loc=1,markerscale=markerscale)
        self.draw()

    def addPlot(self,x,y,bar = False, point = False,marker='k+',marker_size=25):
        if bar == True:
            self.axes.bar(x, y, width=0.01)
        elif point == True :
            self.axes.plot(x,y,marker,markersize=marker_size)
        else:
            self.axes.plot(x,y)

        self.draw()

    def adjust_x_lim(self, x, y):
        test = 0
        xmin = 0
        xmax = 0
        for i in range(len(x)):
            if y[i] >= 0.0002:
                xmin = x[i]
                break
        for i in range(len(x)):
            if y[i] <= 0.002 :
                test += 1
            if test > 60 :
                xmax = x[i-50]
                break
        self.axes.set_xlim(xmin,xmax)
        self.draw()
