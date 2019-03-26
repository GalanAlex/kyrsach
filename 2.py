# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '1.ui'
#
# Created by: PyQt5 UI code generator 5.12
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import (QWidget, QPushButton, QLineEdit,QInputDialog, QApplication)
from y_plot import CustomWidget1
from z_plot import CustomWidget2

import solve as slv
import cagesolve as cgs


class Ui_MainWindow(QWidget):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1057, 731)
        self.centralWidget = QtWidgets.QWidget(MainWindow)
        self.centralWidget.setObjectName("centralWidget")
        self.gridLayoutWidget_3 = QtWidgets.QWidget(self.centralWidget)
        self.gridLayoutWidget_3.setGeometry(QtCore.QRect(530, 450, 411, 111))
        self.gridLayoutWidget_3.setObjectName("gridLayoutWidget_3")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.gridLayoutWidget_3)
        self.gridLayout_3.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_3.setSpacing(6)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.label_12 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_12.setObjectName("label_12")
        self.gridLayout_3.addWidget(self.label_12, 1, 0, 1, 1)

        self.horizontalSlider = QtWidgets.QSlider(self.gridLayoutWidget_3)
        self.horizontalSlider.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider.setObjectName("horizontalSlider")
        self.horizontalSlider.setMinimum(0)
        self.horizontalSlider.setMaximum(cgs.J-1)
        self.horizontalSlider.setValue(0)
        self.gridLayout_3.addWidget(self.horizontalSlider, 0, 0, 1, 1)
        self.horizontalSlider.valueChanged.connect(self.valuechange_z)

        self.horizontalSlider_2 = QtWidgets.QSlider(self.gridLayoutWidget_3)
        self.horizontalSlider_2.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider_2.setObjectName("horizontalSlider_2")
        self.horizontalSlider_2.setMinimum(0)
        self.horizontalSlider_2.setMaximum(cgs.K-1)
        self.horizontalSlider_2.setValue(0)
        self.gridLayout_3.addWidget(self.horizontalSlider_2, 2, 0, 1, 1)
        self.horizontalSlider_2.valueChanged.connect(self.valuechange_y)




        self.label_13 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_13.setObjectName("label_13")
        self.gridLayout_3.addWidget(self.label_13, 4, 0, 1, 1)



        self.gridLayoutWidget_4 = QtWidgets.QWidget(self.centralWidget)
        self.gridLayoutWidget_4.setGeometry(QtCore.QRect(530, 560, 211, 81))
        self.gridLayoutWidget_4.setObjectName("gridLayoutWidget_4")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.gridLayoutWidget_4)
        self.gridLayout_4.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_4.setSpacing(6)
        self.gridLayout_4.setObjectName("gridLayout_4")



        self.label_14 = QtWidgets.QLabel(self.gridLayoutWidget_4)
        self.label_14.setObjectName("label_14")
        self.gridLayout_4.addWidget(self.label_14, 1, 0, 1, 1)

        self.comboBox = QtWidgets.QComboBox(self.gridLayoutWidget_4)

        self.comboBox.setGeometry(QtCore.QRect(550, 600, 151, 23))
        self.comboBox.setObjectName("comboBox")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.activated[str].connect(self.onActivated)
        self.gridLayout_4.addWidget(self.comboBox, 2, 0, 1, 1)

        self.pushButton = QtWidgets.QPushButton(self.gridLayoutWidget_4)
        self.pushButton.setObjectName("pushButton")
        self.gridLayout_4.addWidget(self.pushButton, 1, 0, 1, 1)



        self.graphicsView_2 = CustomWidget1(self.centralWidget)
        self.graphicsView_2.setGeometry(QtCore.QRect(530, 16, 501, 421))
        self.graphicsView_2.setObjectName("graphicsView_2")
        self.graphicsView_3 = CustomWidget2(self.centralWidget)
        self.graphicsView_3.setGeometry(QtCore.QRect(10, 16, 501, 421))
        self.graphicsView_3.setObjectName("graphicsView_3")
        self.gridLayoutWidget = QtWidgets.QWidget(self.centralWidget)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(10, 460, 121, 79))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(11, 11, 11, 11)
        self.gridLayout.setSpacing(6)
        self.gridLayout.setObjectName("gridLayout")
        self.label = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)

        self.doubleSpinBox = QtWidgets.QDoubleSpinBox(self.gridLayoutWidget)
        self.doubleSpinBox.setObjectName("doubleSpinBox")
        self.doubleSpinBox.setSingleStep(0.25)
        self.doubleSpinBox.setValue(slv.Ly)
        self.doubleSpinBox.valueChanged.connect(self.valuechange_Ly)
        self.gridLayout.addWidget(self.doubleSpinBox, 1, 0, 1, 1)


        self.gridLayoutWidget_5 = QtWidgets.QWidget(self.centralWidget)
        self.gridLayoutWidget_5.setGeometry(QtCore.QRect(10, 540, 121, 61))
        self.gridLayoutWidget_5.setObjectName("gridLayoutWidget_5")
        self.gridLayout_5 = QtWidgets.QGridLayout(self.gridLayoutWidget_5)
        self.gridLayout_5.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_5.setSpacing(6)
        self.gridLayout_5.setObjectName("gridLayout_5")

        self.label_9 = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label_9.setObjectName("label_9")
        self.gridLayout_5.addWidget(self.label_9, 0, 0, 1, 1)
        self.gridLayoutWidget_6 = QtWidgets.QWidget(self.centralWidget)
        self.gridLayoutWidget_6.setGeometry(QtCore.QRect(740, 560, 201, 81))
        self.gridLayoutWidget_6.setObjectName("gridLayoutWidget_6")
        self.gridLayout_6 = QtWidgets.QGridLayout(self.gridLayoutWidget_6)
        self.gridLayout_6.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_6.setSpacing(6)
        self.gridLayout_6.setObjectName("gridLayout_6")

        self.label_11 = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.label_11.setObjectName("label_11")
        self.gridLayout_6.addWidget(self.label_11, 0, 0, 1, 1)
        self.label_15 = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.label_15.setObjectName("label_15")
        self.gridLayout_6.addWidget(self.label_15, 1, 0, 1, 1)

        self.gridLayoutWidget_7 = QtWidgets.QWidget(self.centralWidget)
        self.gridLayoutWidget_7.setGeometry(QtCore.QRect(130, 540, 121, 61))
        self.gridLayoutWidget_7.setObjectName("gridLayoutWidget_7")
        self.gridLayout_7 = QtWidgets.QGridLayout(self.gridLayoutWidget_7)
        self.gridLayout_7.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_7.setSpacing(6)
        self.gridLayout_7.setObjectName("gridLayout_7")
        self.label_10 = QtWidgets.QLabel(self.gridLayoutWidget_7)
        self.label_10.setObjectName("label_10")
        self.gridLayout_7.addWidget(self.label_10, 0, 0, 1, 1)

        self.gridLayoutWidget_9 = QtWidgets.QWidget(self.centralWidget)
        self.gridLayoutWidget_9.setGeometry(QtCore.QRect(250, 540, 121, 61))
        self.gridLayoutWidget_9.setObjectName("gridLayoutWidget_9")
        self.gridLayout_9 = QtWidgets.QGridLayout(self.gridLayoutWidget_9)
        self.gridLayout_9.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_9.setSpacing(6)
        self.gridLayout_9.setObjectName("gridLayout_9")
        self.label_17 = QtWidgets.QLabel(self.gridLayoutWidget_9)
        self.label_17.setObjectName("label_17")
        self.gridLayout_9.addWidget(self.label_17, 0, 1, 1, 1)
        self.doubleSpinBox_5 = QtWidgets.QDoubleSpinBox(self.gridLayoutWidget_9)
        self.doubleSpinBox_5.setValue(cgs.J)
        self.doubleSpinBox_5.setMaximum(150)
        self.doubleSpinBox_5.setMinimum(1)
        self.doubleSpinBox_5.valueChanged.connect(self.valuechange_J)
        self.doubleSpinBox_5.setObjectName("doubleSpinBox_5")
        self.gridLayout_9.addWidget(self.doubleSpinBox_5, 1, 1, 1, 1)


        self.gridLayoutWidget_2 = QtWidgets.QWidget(self.centralWidget)
        self.gridLayoutWidget_2.setGeometry(QtCore.QRect(130, 460, 121, 79))
        self.gridLayoutWidget_2.setObjectName("gridLayoutWidget_2")
        self.gridLayout_8 = QtWidgets.QGridLayout(self.gridLayoutWidget_2)
        self.gridLayout_8.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_8.setSpacing(6)
        self.gridLayout_8.setObjectName("gridLayout_8")
        self.label_2 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_2.setObjectName("label_2")
        self.gridLayout_8.addWidget(self.label_2, 0, 0, 1, 1)


        self.doubleSpinBox_2 = QtWidgets.QDoubleSpinBox(self.gridLayoutWidget_2)
        self.doubleSpinBox_2.setSingleStep(1)
        self.doubleSpinBox_2.setValue(slv.L)
        self.doubleSpinBox_2.valueChanged.connect(self.valuechange_L)
        self.gridLayout_8.addWidget(self.doubleSpinBox_2, 1, 0, 1, 1)


        self.gridLayoutWidget_8 = QtWidgets.QWidget(self.centralWidget)
        self.gridLayoutWidget_8.setGeometry(QtCore.QRect(250, 460, 121, 79))
        self.gridLayoutWidget_8.setObjectName("gridLayoutWidget_8")
        self.gridLayout_26 = QtWidgets.QGridLayout(self.gridLayoutWidget_8)
        self.gridLayout_26.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_26.setSpacing(6)
        self.gridLayout_26.setObjectName("gridLayout_26")
        self.label_5 = QtWidgets.QLabel(self.gridLayoutWidget_8)
        self.label_5.setObjectName("label_5")
        self.gridLayout_26.addWidget(self.label_5, 0, 0, 1, 1)

        self.doubleSpinBox_3 = QtWidgets.QDoubleSpinBox(self.gridLayoutWidget_8)
        self.doubleSpinBox_3.setObjectName("doubleSpinBox_3")
        self.doubleSpinBox_3.setSingleStep(0.1)
        self.doubleSpinBox_3.setValue(slv.n)
        self.doubleSpinBox_3.valueChanged.connect(self.valuechange_n)
        self.gridLayout_26.addWidget(self.doubleSpinBox_3, 1, 0, 1, 1)


        self.gridLayoutWidget_10 = QtWidgets.QWidget(self.centralWidget)
        self.gridLayoutWidget_10.setGeometry(QtCore.QRect(370, 460, 121, 79))
        self.gridLayoutWidget_10.setObjectName("gridLayoutWidget_10")
        self.gridLayout_27 = QtWidgets.QGridLayout(self.gridLayoutWidget_10)
        self.gridLayout_27.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_27.setSpacing(6)
        self.gridLayout_27.setObjectName("gridLayout_27")
        self.label_6 = QtWidgets.QLabel(self.gridLayoutWidget_10)
        self.label_6.setObjectName("label_6")
        self.gridLayout_27.addWidget(self.label_6, 0, 1, 1, 1)

        self.doubleSpinBox_4 = QtWidgets.QDoubleSpinBox(self.gridLayoutWidget_10)
        self.doubleSpinBox_4.setObjectName("doubleSpinBox_4")
        self.doubleSpinBox_4.setSingleStep(0.1)
        self.doubleSpinBox_4.setValue(slv.lambd)
        self.doubleSpinBox_4.valueChanged.connect(self.valuechange_lambd)
        self.gridLayout_27.addWidget(self.doubleSpinBox_4, 1, 1, 1, 1)


        self.gridLayoutWidget_11 = QtWidgets.QWidget(self.centralWidget)
        self.gridLayoutWidget_11.setGeometry(QtCore.QRect(370, 540, 121, 61))
        self.gridLayoutWidget_11.setObjectName("gridLayoutWidget_11")
        self.gridLayout_10 = QtWidgets.QGridLayout(self.gridLayoutWidget_11)
        self.gridLayout_10.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_10.setSpacing(6)
        self.gridLayout_10.setObjectName("gridLayout_10")
        self.label_18 = QtWidgets.QLabel(self.gridLayoutWidget_11)
        self.label_18.setObjectName("label_18")
        self.gridLayout_10.addWidget(self.label_18, 0, 0, 1, 1)

        self.doubleSpinBox_6 = QtWidgets.QDoubleSpinBox(self.gridLayoutWidget_11)
        self.doubleSpinBox_6.setValue(cgs.K)
        self.doubleSpinBox_6.setMaximum(150)
        self.doubleSpinBox_6.setMinimum(1)
        self.doubleSpinBox_6.valueChanged.connect(self.valuechange_K)
        self.doubleSpinBox_6.setObjectName("doubleSpinBox_6")
        self.gridLayout_10.addWidget(self.doubleSpinBox_6, 1, 0, 1, 1)
        MainWindow.setCentralWidget(self.centralWidget)
        self.mainToolBar = QtWidgets.QToolBar(MainWindow)
        self.mainToolBar.setObjectName("mainToolBar")
        MainWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.mainToolBar)
        self.statusBar = QtWidgets.QStatusBar(MainWindow)
        self.statusBar.setObjectName("statusBar")
        MainWindow.setStatusBar(self.statusBar)
        self.actionr = QtWidgets.QAction(MainWindow)
        self.actionr.setObjectName("actionr")

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.comboBox.setItemText(0, _translate("MainWindow", "Кранка-Николсона"))
        self.comboBox.setItemText(1, _translate("MainWindow", "Явная"))
        self.comboBox.setItemText(2, _translate("MainWindow", "Неявная"))
        self.pushButton.setText(_translate("MainWindow", "Update"))
        self.label_12.setText(_translate("MainWindow", "z for left graphic"))
        self.label_13.setText(_translate("MainWindow", "y for right graphic"))
        self.label.setText(_translate("MainWindow", "    Ly "))
        self.label_9.setText(_translate("MainWindow", " Lx = 4"))
        self.label_10.setText(_translate("MainWindow", " k = 2*pi/lambda"))
        self.label_17.setText(_translate("MainWindow", "J"))
        self.label_2.setText(_translate("MainWindow", "   L"))
        self.label_5.setText(_translate("MainWindow", "   n"))
        self.label_6.setText(_translate("MainWindow", "   lambda"))
        self.label_18.setText(_translate("MainWindow", "K"))
        self.actionr.setText(_translate("MainWindow", "r"))

    '''
    def upd(self):
        CustomWidget1.update(self)
        CustomWidget2.update(self)
    '''

    def valuechange_z(self):
        CustomWidget1.z = self.horizontalSlider.value()
        self.label_12.setText("Z for first graph: {:.2f}s".format(
            4*self.horizontalSlider.value() / cgs.K))  # UPD: показывает значение слайдера

    def valuechange_y(self):
        CustomWidget2.y = self.horizontalSlider_2.value()
        self.label_13.setText("Y for second graph: {:.2f}mcm".format(
            6*self.horizontalSlider_2.value() / cgs.J))

    def valuechange_Ly(self):
        slv.Ly = self.doubleSpinBox.value()

    def valuechange_L(self):
        slv.L = self.doubleSpinBox_2.value()

    def valuechange_n(self):
        slv.n = self.doubleSpinBox_3.value()

    def valuechange_lambd(self):
        slv.lambd = self.doubleSpinBox_4.value()

    def valuechange_J(self):
        cgs.L = int(self.doubleSpinBox_5.value())
        #self.horizontalSliderT.setMaximum(cgs.J)

    def valuechange_K(self):
        cgs.K = int(self.doubleSpinBox_6.value())
        #self.horizontalSliderT.setMaximum(cgs.K)

    def onActivated(self, text):
        if text == "Кранка-Николсона":
            cgs.MODE = 'c'
        if text == "Явная":
            cgs.MODE = 'o'
        if text == "Неявная" :
            cgs.MODE = 'im'


from y_plot import CustomWidget1
from z_plot import CustomWidget2


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())