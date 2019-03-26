#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import numpy as np

import solve as slv
import cagesolve as cgs

class CustomWidget1(pg.GraphicsWindow):
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')

    z = 1
    def __init__(self, parent=None, **kargs):
        pg.GraphicsWindow.__init__(self, **kargs)
        self.setParent(parent)
        self.p = self.addPlot(labels={'left': 'Temperature', 'bottom': 'Radius'})

        ing = 150
        x1 = np.linspace(0, cgs.y, ing + 1)
        self.data = np.zeros(ing + 1)
        self.data1 = []

        x2 = np.linspace(0, cgs.Ly, cgs.J)
        for i in range(ing + 1):
            self.data[i] = slv.main_sum(i*6/ ing,  self.z*4/ cgs.K)
        v = cgs.mainsolve(cgs.J, cgs.K)
        self.data1 = v[:,0]
        m =0


        #print(self.data1)
        #print(x2.shape)

        #print(self.data1)
        #self.curve = self.p.plot(x1, self.data, pen=(16, 41, 89), symbolBrush=(16, 41, 89))
        self.curve1 = self.p.plot(x2, (self.data1), pen=(38, 104, 5), symbolBrush=(38, 104, 5))

        #self.p.setYRange(0, 1)

        # self.p.addLegend()


        l = pg.LegendItem((120, 60), offset=(600, 30))
        l.setParentItem(self.p)
        #l.addItem(self.curve, 'Аналитическое решение')
        l.addItem(self.curve1, 'Численное решение')

        self.p.showGrid(x=True, y=True)
        self.p.setLabel("left", text='Amplitude', units ="a" )
        self.p.setLabel("bottom", text='Y', units="m")


        #timer = pg.QtCore.QTimer(self)
        #timer.timeout.connect(self.update)
        #timer.start(1)

    def update(self):
        ing = 150
        x1 = np.linspace(0, cgs.Ly, ing + 1)
        x2 = np.linspace(0, cgs.Ly, cgs.J)
        self.data = np.zeros(ing + 1)
        v = cgs.mainsolve(cgs.J, cgs.K)
        self.data1 = v[:, self.z]
        for i in range(ing + 1):
            self.data[i] = slv.main_sum(i*6 / ing,  self.z*4 / cgs.K)
        #print("for Y plot z = ",self.z)
        self.p.clear()
        self.curve = self.p.plot(x1, self.data, pen=(16, 41, 89), symbolBrush=(16, 41, 89))
        self.curve1 = self.p.plot(x2, (self.data1), pen=(38, 104, 5), symbolBrush=(38, 104, 5))

        #self.p.setYRange(0, 1)





if __name__ == '__main__':
    w = CustomWidget1()
    w.show()
    QtGui.QApplication.instance().exec_()
