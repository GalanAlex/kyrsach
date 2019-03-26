#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import numpy as np

import solve as slv
import cagesolve as cgs

class CustomWidget2(pg.GraphicsWindow):
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')

    y = 1

    def __init__(self, parent=None, **kargs):
        pg.GraphicsWindow.__init__(self, **kargs)
        self.setParent(parent)
        self.p = self.addPlot(labels={'left': 'Amplitude', 'bottom': 'z'})

        ing = 150
        x1 = np.linspace(0, cgs.z, ing + 1)
        self.data = np.zeros(ing + 1)
        self.data1 = []

        x2 = np.linspace(0, cgs.L, cgs.K)
        #print(x2.shape)
        for i in range(ing + 1):
            self.data[i] = slv.main_sum(6*self.y/ cgs.J, 4*i/ ing)

        v = (cgs.mainsolve(cgs.J, cgs.K))
        #print(v.shape)
        #print("fqwefweqf")
        self.data1 = v[1,:]
        print(self.data1)
        #print(x2)

        #self.curve = self.p.plot(x1, self.data, pen=(16, 41, 89), symbolBrush=(16, 41, 89))#аналитическоя
        self.curve1 = self.p.plot(x2, (self.data1), pen=(242, 219, 8), symbolBrush=(242, 219, 8))#cage

        #self.p.setYRange(0.0,1)
        #self.p.setXRange(0,0.2)

        #self.p.setXRange(0.17,0.22)
        #self.p.setXRange(0,0.001)


        l = pg.LegendItem((100, 60), offset=(650, 30))
        l.setParentItem(self.p)
        #l.addItem(self.curve, 'Аналитическое решение')
        l.addItem(self.curve1, 'Численное решение')



        self.p.showGrid(x=True, y=True)
        self.p.setLabel("left", text='Amplitude', units ="a")
        self.p.setLabel("bottom", text='Z', units="s")

        #self.p.setYRange(20,21)
        #o = self.p.getAxis()
        #n = self.p.getAxis("bottom")


        #timer = pg.QtCore.QTimer(self)
        #timer.timeout.connect(self.update)
        #timer.start(1)

    def update(self):
        ing = 150
        x1 = np.linspace(0, cgs.L, ing + 1)
        x2 = np.linspace(0, cgs.L, cgs.K )
        self.data = np.zeros(ing + 1)
        v = cgs.mainsolve(cgs.J, cgs.K)
        self.data1 = v[self.y, :]

        for i in range(ing + 1):
            self.data[i] = slv.main_sum( self.y*6 / cgs.J, i *4/ ing)
        #print("for z plot y = ",self.y)

        self.p.clear()
        #self.curve = self.p.plot(x1, self.data, pen=(16, 41, 89), symbolBrush=(16, 41, 89))
        self.curve1 = self.p.plot(x2, (self.data1), pen=(38, 104, 5), symbolBrush=(38, 104, 5))

        #self.p.setYRange(0,1)
        #self.p.setXRange(0,0.2)


if __name__ == '__main__':
    w = CustomWidget2()
    w.show()
    QtGui.QApplication.instance().exec_()

    #import sys
    #if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
    #    QtGui.QApplication.instance().exec_()
