#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""

import sys,os


#Gui
import gui.open_graph as OpenGraph

#Utility
import utility.write as Write
import utility.canvas as Canvas

from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar2

from numpy import linspace,random,array
try :
    from PyQt4 import Qt,QtGui,QtCore,Qwt5
except :
    pass;



#---------------------------------------------------#
#---------------------------------------------------#
#----------------GRAPHICS---------------------------#
#---------------------------------------------------#
#---------------------------------------------------#
class graphic(QtGui.QMainWindow):
    x = 0
    y = 0
    xlabel  = ""
    ylabel  = ""
    average = True
    def __init__(self,pX,pY,pxlbl,pylbl,average=True, adjust=False,\
                     name = "Graphics",point = False,marker='.',marker_size=25, parent= None):
        
        self.x = pX
        self.y = pY
        self.xlabel=pxlbl
        self.ylabel=pylbl
        self.average = average
        self.initUI(parent,adjust,name,point,marker,marker_size)
      	self.raise_()

    def initUI(self, parent,adjust,name,point,marker,marker_size):        

        QtGui.QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle(name)
        self.resize(1000, 1000)
        screen = QtGui.QDesktopWidget().screenGeometry()
        self.move(random.randint(0, screen.height()),random.randint( 0, screen.width()))
        self.open = QtGui.QAction('&Open', self)
        self.open.setShortcut('Ctrl+O')
        self.open.setStatusTip('import graphics')
        self.connect(self.open, QtCore.SIGNAL('triggered()'),self.openGraph)

        self.close = QtGui.QAction('&Exit', self)
        self.close.setShortcut('Ctrl+Q')
        self.close.setStatusTip('Exit application')
        self.connect(self.close, QtCore.SIGNAL('triggered()'), QtCore.SLOT('close()'))

        self.printG = QtGui.QAction('&Print', self)
        self.printG.setShortcut('Ctrl+P')
        self.printG.setStatusTip('Print graphics')
        self.connect(self.printG, QtCore.SIGNAL('triggered()'),self.printGraph)

        self.save = QtGui.QAction('&Save', self)
        self.save.setShortcut('Ctrl+S')
        self.save.setStatusTip('Save Graphic')
        self.connect(self.save, QtCore.SIGNAL('triggered()'), self.showDialog)
        
        self.setStatusBar(QtGui.QStatusBar())
        self.menubar = self.menuBar()
        self.fileMenu1 = self.menubar.addMenu('&File')
        #self.fileMenu1.addAction(self.open)
        self.fileMenu1.addAction(self.save)
        self.fileMenu1.addAction(self.printG)
        self.fileMenu1.addAction(self.close)

        self.main_widget  = QtGui.QWidget(self)
        self.main_lay  = QtGui.QHBoxLayout(self.main_widget)

        #Creation of the graphics widged
        self.graph_widget = QtGui.QWidget(self)        
        self.graph_lay = QtGui.QVBoxLayout(self.graph_widget)

        self.G = Canvas.CanvasQT(self.main_widget, width=6, height=4, dpi=100,\
                                    x=self.x,y=self.y,pxlbl=self.xlabel,pylbl=self.ylabel\
                                     ,point=point,marker=marker,marker_size=marker_size)
        self.navigation_toolbar = NavigationToolbar2(self.G, self)
        

        #creation of the option widget
        self.test = QtGui.QWidget()
        self.test.hide()
        self.option_layout = QtGui.QGridLayout(self.test)

        self.lbl2 = QtGui.QLabel("select file :", self)
        self.lbl2.setFixedSize(70,36)
        self.tefile = QtGui.QTextEdit()
        self.tefile.setReadOnly(True)
        self.tefile.setFixedSize(100,36)
 
        self.browse =  QtGui.QPushButton('&Browse', self)
        self.browse.setFixedSize(100,36)
#        self.connect(self.browse ,QtCore.SIGNAL("clicked()"),self.openFile)	

        
        self.tbro = QtGui.QTextBrowser()
        self.tbro.setFixedHeight(500)
        self.tbro.setFixedWidth(200)
        self.tbro.setLineWrapMode(QtGui.QTextEdit.NoWrap)

        self.lbl3 = QtGui.QLabel("x :", self)
        self.lbl3.setFixedSize(15,36)

        self.CBox1 = QtGui.QComboBox()
        self.CBox1.addItem("")
        self.CBox1.setFixedSize(65,36)
#        self.connect(self.CBox1,QtCore.SIGNAL('currentIndexChanged(const QString&)'),self.changeData)

        self.lbl4 = QtGui.QLabel("y :", self)
        self.lbl4.setFixedSize(15,36)

        self.CBox2= QtGui.QComboBox()
        self.CBox2.addItem("")
        self.CBox2.setFixedSize(65,36)
#        self.connect(self.CBox2,QtCore.SIGNAL('currentIndexChanged(const QString&)'),self.changeData)


        self.imp =  QtGui.QPushButton('&import', self)
        self.imp.setFixedSize(100,36)
 #       self.connect(self.imp ,QtCore.SIGNAL("clicked()"),self.close)	


        self.option_layout.addWidget(self.lbl2   , 0, 0, 1, 1, QtCore.Qt.AlignCenter)
        self.option_layout.addWidget(self.tefile , 0, 1, 1, 4,QtCore.Qt.AlignLeft)
        self.option_layout.addWidget(self.browse , 0 ,5, 1, 1,QtCore.Qt.AlignCenter)
        self.option_layout.addWidget(self.tbro   , 1, 0, 10, 6, QtCore.Qt.AlignCenter)
        self.option_layout.addWidget(self.lbl3   ,12, 0, 1, 1, QtCore.Qt.AlignRight)
        self.option_layout.addWidget(self.CBox1  ,12, 1, 1, 1, QtCore.Qt.AlignLeft)
        self.option_layout.addWidget(self.lbl4   ,12, 2, 1, 1, QtCore.Qt.AlignRight)
        self.option_layout.addWidget(self.CBox2  ,12, 3, 1, 1, QtCore.Qt.AlignLeft)
        self.option_layout.addWidget(self.imp    ,12, 5, 1, 1, QtCore.Qt.AlignLeft)
 

        self.graph_lay.addWidget(self.G)
        self.graph_lay.addWidget(self.navigation_toolbar, 0)
        self.main_lay.addWidget(self.graph_widget)
        self.main_lay.addWidget(self.test)

        if self.average:
            self.checkbox_average =QtGui.QCheckBox("Show average Value")
            self.connect(self.checkbox_average, QtCore.SIGNAL('clicked()'), self.addAverage)
            self.lbl1 = QtGui.QLabel()
            self.lbl1.setFixedWidth(250)
            average = self.averageValue()
            display  = self.displayData(average,self.getStandardDeviation(self.y ,average))
            self.lbl1.setText('Average : ' + display)
            self.navigation_toolbar.addWidget(self.lbl1)
            self.navigation_toolbar.addWidget(self.checkbox_average)

        if adjust:
            self.G.adjust_x_lim(self.x,self.y)

        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

    def showDialog(self):
        fname = QtGui.QFileDialog.getSaveFileName(self,"Save Graphics",os.getcwd())
        if fname !="":
           Write.SaveFile(fname).saveGraph(self.x,self.y,self.xlabel,self.ylabel)


    def openGraph(self):
        self.importGraph = OpenGraph.OpenGraph()
        self.importGraph.raise_()
        if self.importGraph.exec_():
            data = self.importGraph.importGraph()
        else:
            data = self.importGraph.importGraph()

        if data !=0 :
            self.addPlot(data[0],data[1],bar = False,point = False)

    def setPlot(self,x,y,pxlbl,pylbl,adjust=False,point =False):
        self.G.setPlot(x,y,pxlbl,pylbl,adjust=adjust,point = point)

    def addGraph(self,px,py,bar = False,point = False,marker='.',marker_size=25):
        self.G.addPlot(px,py,bar = bar,point = point, marker=marker,marker_size=marker_size)

    def addPlot(self,px,py,bar = False,point = False,marker='.',marker_size=25):
        self.G.addPlot(px,py,bar = bar,point = point, marker=marker,marker_size=marker_size)

    def addLegend(self,legend,markerscale = 1):
        self.G.addLegend(legend,markerscale=markerscale)

    def addAverage(self):
        if self.checkbox_average.isChecked():
            ytemp = self.averageValue()
            ya = []
            for i in range(len(self.y)):
               ya.insert(i,ytemp)
            self.G.addPlot(self.x,ya)
            self.addLegend([self.ylabel,self.lbl1.text()])
#beta            self.test.show()
        else:
#beta            self.test.hide()
            self.G.setPlot(self.x,self.y,self.xlabel,self.ylabel)        

    def closeEvent(self, event):
        del self

    def update(self,pX,pY,pxlbl,pylbl,adjust=False, name = '',point = False,marker='.', marker_size=25):
        self.setWindowTitle(pylbl)
        self.x = pX
        self.y = pY
        self.xlabel=pxlbl
        self.ylabel=pylbl

        if self.average:
            average = self.averageValue()
            display  = self.displayData(average,self.getStandardDeviation(self.y ,average))
            self.lbl1.setText('Average : ' + display)

        self.G.setPlot(self.x,self.y,self.xlabel,self.ylabel,point=point,marker=marker, marker_size=marker_size)

        if self.isCheckbox_averageChecked():
            self.addAverage()

        if adjust:
            self.G.adjust_x_lim(self.x,self.y)

        if name != '':
            self.setWindowTitle(name)


    def updatePos(self,pX,pY,pxlbl,pylbl,adjust=False, name = '',point = True,marker='.', marker_size=2):
        self.setWindowTitle(pylbl)
        self.x = pX
        self.y = pY
        self.xlabel=pxlbl
        self.ylabel=pylbl

        self.G.setPlot(self.x,self.y,self.xlabel,self.ylabel,point=point,marker=marker, marker_size=marker_size)

    def averageValue(self):
        ytemp = 0
        for i in range(len(self.y)):
            ytemp+=float(self.y[i])
        ytemp/=len(self.y)
        return ytemp

    def getStandardDeviation(self,data ,averageData = 0):
        res = 0
        if averageData == 0:
            averageData = self.getMoy(data)

        for i in range(len(data)):
            res += (data[i] - averageData)**2
        res /= len(data)

        res = res**0.5

        return res

    def displayData(self,data,deviation):
        if(abs(data) < 1e-04 and deviation < 1e-04):
            return str("%.4g" %data) + ' +/- ' + str("%.4g" %deviation)
        if(abs(data) > 1e-04 and deviation < 1e-04 ):
            return str(round(data,2)) + ' +/- ' + str("%.4g" %deviation)
        if(abs(data) > 1e-04 and deviation > 1e-04 ):
            return str(round(data,2)) + ' +/- ' + str(round(deviation,2))
        return ''


    def isCheckbox_averageChecked(self):
        try:
            return self.checkbox_average.isChecked()
        except:
            return False


    def printGraph(self):
        ######BETA###########
        #printer = Qt.QPrinter()
        #printer.A4
        #printer.HighResolution
        #printer.Color
        #anotherWidget= Qt.QPrintDialog(printer,self)
        #print
        #if(anotherWidget.exec_() != Qt.QDialog.Accepted):
            #return
            #print 'test' + printer.outputFileName()
        #else:
            #print 'test' + printer.outputFileName()
            #self.G.print_figure(printer)

        #print 'test' + printer.outputFileName()

        #print printer.outputFileName()
        #p = Qt.QPixmap.grabWidget(self.G)
        #printLabel = Qt.QLabel()
        #printLabel.setPixmap(p)
        #painter = Qt.QPainter(printer)
        #self.G.render(painter)
        #printLabel.render(painter)
        #painter.end()

        printer = Qt.QPrinter(Qt.QPrinter.HighResolution)

        printer.setOutputFileName('bode-example-%s.ps' % Qt.qVersion())

        printer.setCreator('Bode example')
        printer.setOrientation(Qt.QPrinter.Landscape)
        printer.setColorMode(Qt.QPrinter.Color)

        #docName = self.plot.title().text()
        #if not docName.isEmpty():
            #docName.replace(Qt.QRegExp(Qt.QString.fromLatin1('\n')), self.tr(' -- '))
            #printer.setDocName(docName)

        dialog = Qt.QPrintDialog(printer)
        if dialog.exec_():
            #filter = Qt.PrintFilter()
            #if (QPrinter.GrayScale == printer.colorMode()):
                #filter.setOptions(
                    #Qt.QwtPlotPrintFilter.PrintAll
                    #& ~QwtPlotPrintFilter.PrintBackground
                    #| QwtPlotPrintFilter.PrintFrameWithScales)
            self.G.print_(printer)#, filter)

#---------------------------------------------------#
#---------------------------------------------------#
#---------------------------------------------------#
#---------------------------------------------------#
#---------------------------------------------------#

