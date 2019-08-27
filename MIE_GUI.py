"""
@author: Emre Gürdal
"""
import sys
import scipy as sc
from scipy.special import jv, hankel1
from scipy.interpolate import InterpolatedUnivariateSpline
from PySide2 import QtGui,QtWidgets
from PySide2.QtCore import Qt
from PySide2.QtWidgets import (QMainWindow, QApplication, QWidget, QLabel, QSlider, QLineEdit,
                               QPushButton, QRadioButton, QComboBox, QHBoxLayout, QVBoxLayout, QSizePolicy, QFileDialog)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FC
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np

class MplCanvas(FC):
    def __init__(self, parent=None, width=8, height=6, a=50, n_m=1, material="Gold", csection="Scattering"):
        fig=Figure(figsize=(width,height))
        fig.patch.set_facecolor('#8b8b8b')
        self.ax=fig.add_subplot(111)
        FC.__init__(self,fig)
        FC.setSizePolicy(self,QSizePolicy.Expanding,QSizePolicy.Expanding)
        FC.updateGeometry(self)
        self.clearPlot()
        #self.drawPlot(a,n_m,material)

#A().drawPlot(50,2)
    def clearPlot(self):
        self.ax.clear()
        self.ax.set_xlim(300,900)
        self.ax.set_xlabel('Wavelength [$nm$]')
        self.ax.set_ylabel("$\sigma$ $[nm^2]$")
        self.ax.set_title('MIE scattering')
        self.ax.grid(True)
        self.draw()

    def save_File(self):
        global Q_cross
        global wl
        try:
            file=np.array([wl,Q_cross])
        except:
          print("Error")
        return file
    
    def drawPlot(self,a,n_m,material,csection):
        global Q_cross
        global wl
        #self.ax.clear()
        #READ JOHNSON AND CHRISTY (GOLD,COPPER,SILVER) CSV:
        #https://refractiveindex.info/?shelf=main&book=Au&page=Johnson
        wl,n1=np.loadtxt(str(material)+'_real.csv',usecols=[0,1],unpack=True,delimiter=',', skiprows=1)
        wl,k1=np.loadtxt(str(material)+'_imag.csv',usecols=[0,1],unpack=True,delimiter=',', skiprows=1)
        wl=wl*1e3
        
        #INTERPOLATE REFRACTIVE INDEX
        xi=np.linspace(wl[0],wl[48],1001)
        int_n1=InterpolatedUnivariateSpline(wl,n1)
        int_k1=InterpolatedUnivariateSpline(wl,k1)
        n1=int_n1(xi)
        k1=int_k1(xi)
        wl=xi
        n_p=n1+k1*1j
        
        #MIE 
        Q_sca_spek=[]
        Q_ext_spek=[]
        Q_abs_spek=[]
        
        for i in range(0,1001):
            m=n_p[i]/n_m
            k=2*sc.pi*n_m/wl[i]
            u=k*a
            N=round(u+4*u**(1/3)+2)
            L=np.arange(1,(N+1))#[1,...,N]
        
        #SPHERICAL BESSEL FUNCTION 1.KIND
            jL=(sc.pi/(2*u))**0.5*(np.array(jv(np.insert(L,0,0)+0.5,u)))
            jLm=(sc.pi/(2*u*m))**0.5*(np.array(jv(np.insert(L,0,0)+0.5,u*m)))
        
        #SPHERICAL BESSEL FUNCTION 1.KIND (DERIVATION)
            jL_der=jL[:np.int64(N)]-L/u*jL[np.int64(L)]
            jLm_der=jLm[:np.int64(N)]-L/(u*m)*jLm[np.int64(L)]
        
        #SPHERICAL BESSEL FUNCTION 3.KIND
            hL=(sc.pi/(2*u))**0.5*(np.array(hankel1(np.insert(L,0,0)+0.5,u)))
        
        #SPHERICAL BESSEL FUNCTION 3.KIND (DERIVATION)
            hL_der=hL[:np.int64(N)]-L/u*hL[np.int64(L)]
         
        #COEFFICIENTS
            a_L=(m*jLm[np.int64(L)]*jL_der-jL[np.int64(L)]*jLm_der)/\
                (m*jLm[np.int64(L)]*hL_der-hL[np.int64(L)]*jLm_der)
            
            b_L=(jLm[np.int64(L)]*jL_der-m*jL[np.int64(L)]*jLm_der)/\
                (jLm[np.int64(L)]*hL_der-m*hL[np.int64(L)]*jLm_der)
                
        #SCATTERING, EXTINCTION, ABSORPTION
            Q_sca=1/(2*sc.pi)*(wl[i]/n_m)**2*np.sum([(2*L+1)*(np.abs(a_L)**2+\
                               np.abs(b_L)**2)])
            Q_ext=1/(2*sc.pi)*(wl[i]/n_m)**2*np.sum([(2*L+1)*(a_L+b_L).real])
            Q_abs=Q_ext-Q_sca
            
            Q_sca_spek.append((Q_sca))
            Q_ext_spek.append((Q_ext))
            Q_abs_spek.append((Q_abs))
            
        Q_sca_spek=np.array(Q_sca_spek)
        Q_ext_spek=np.array(Q_ext_spek)
        Q_abs_spek=np.array(Q_abs_spek)
         
        #FIGURE
        if csection=="Scattering":           
            self.ax.plot(wl,Q_sca_spek)
            Q_cross=Q_sca_spek
        elif csection=="Absorption":
            self.ax.plot(wl,Q_abs_spek)
            Q_cross=Q_abs_spek
        elif csection=="Extinction":
            self.ax.plot(wl,Q_ext_spek)
            Q_cross=Q_ext_spek

        self.ax.set_xlim(300,900)
        self.ax.set_xlabel('Wavelength [$nm$]')
        self.ax.set_ylabel("$\sigma$ $[nm^2]$")
        self.ax.set_title('MIE scattering')
        self.ax.grid(True)
        #self.ax.legend()
        self.draw()

class MainApp(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setAttribute(Qt.WA_DeleteOnClose)
        self.setWindowTitle("MIE")
        self.main_widget=QWidget(self)
        self.setStyleSheet('QMainWindow {background-color: #8b8b8b}')#cfcfcf
        self.setWindowIcon(QtGui.QIcon('mie_icon.ico'))
        
        #DIAMETER
        self.label_d=QLabel("Radius [nm|:",self)
        self.slider_d=QSlider(Qt.Horizontal)
        self.slider_d.setMinimum(5)
        self.slider_d.setMaximum(200)
        self.slider_d.setTickInterval(5)
        self.slider_d.setValue(50)
        self.edit_d=QLineEdit(self)
        self.edit_d.setMaxLength(5)
        a=self.slider_d.value()
        self.edit_d.setText(str(a))
        self.push_d_down=QPushButton("<",self)
        self.push_d_up=QPushButton(">",self)
        
        #REFRACTIVE INDEX
        self.label_n=QLabel("Refractiv index (medium):",self)
        self.slider_n=QSlider(Qt.Horizontal)
        self.slider_n.setMinimum(100)
        self.slider_n.setMaximum(300)
        self.slider_n.setTickInterval(10)
        self.slider_n.setValue(100)
        self.edit_n=QLineEdit(self)
        self.edit_n.setMaxLength(5)
        n_m=self.slider_n.value()/100
        self.edit_n.setText(str(n_m))
        self.push_n_down=QPushButton("<",self)
        self.push_n_up=QPushButton(">",self)
        
        #CROSS SECTION
        self.label_radio=QLabel("Cross section:",self)
        self.radiobutton_sca = QRadioButton("Scattering")
        self.radiobutton_sca.setChecked(True)
        self.radiobutton_abs = QRadioButton("Absorption")
        self.radiobutton_ext = QRadioButton("Extinction")
        csection="Scattering"

        #MATERIAL
        self.label_m=QLabel("Material (Nanosphere):",self)
        self.cb_m = QComboBox()
        list_m=["Gold","Silver","Copper"]
        self.cb_m.addItems(list_m)
        material="Gold"
        
        #PLOT
        self.plot_mie=QPushButton("PLOT",self)
        self.plot_mie.setStyleSheet('QPushButton {background-color: #12E316}')

        #CLEAR PLOT
        self.clear_mie=QPushButton("CLEAR",self)
        self.clear_mie.setStyleSheet('QPushButton {background-color: #FF0505}')

        #SAVE FILE
        self.save_mie=QPushButton("SAVE FILE",self)
        self.save_mie.setStyleSheet('QPushButton {background-color: #cfcfcf}')
            
        #LAYOUT CROSS SECTION
        self.layout_radio=QHBoxLayout()
        self.layout_radio.addWidget(self.radiobutton_sca)
        self.layout_radio.addWidget(self.radiobutton_abs)
        self.layout_radio.addWidget(self.radiobutton_ext)

        #LAYOUT DIAMETER
        self.layout_d1=QHBoxLayout()
        self.layout_d1.addWidget(self.push_d_down)
        self.layout_d1.addWidget(self.slider_d)
        self.layout_d1.addWidget(self.push_d_up)
        
        self.layout_d=QVBoxLayout()
        self.layout_d.addWidget(self.label_d)
        self.layout_d.addSpacing(1)
        self.layout_d.addLayout(self.layout_d1)
        self.layout_d.addSpacing(1)
        self.layout_d.addWidget(self.edit_d)
        
        #LAYOUT REFRACTIVE INDEX
        self.layout_n1=QHBoxLayout()
        self.layout_n1.addWidget(self.push_n_down)
        self.layout_n1.addWidget(self.slider_n)
        self.layout_n1.addWidget(self.push_n_up)
        
        self.layout_n=QVBoxLayout()
        self.layout_n.addWidget(self.label_n)
        self.layout_n.addSpacing(1)
        self.layout_n.addLayout(self.layout_n1)
        self.layout_n.addSpacing(1)
        self.layout_n.addWidget(self.edit_n)
        
        #LAYOUT MATERIAL
        self.layout_m=QVBoxLayout()
        self.layout_m.addWidget(self.label_m)
        self.layout_m.addSpacing(1)
        self.layout_m.addWidget(self.cb_m)
        
        #LAYOUT PARAMETER BELOW
        self.layout_parameter=QHBoxLayout()
        self.layout_parameter.addLayout(self.layout_n)
        self.layout_parameter.addLayout(self.layout_d)
        
        #LAYOUT PARAMETER RADIO
        self.layout_parameter_radio=QVBoxLayout()
        self.layout_parameter_radio.addWidget(self.label_radio)
        self.layout_parameter_radio.addSpacing(1)
        self.layout_parameter_radio.addLayout(self.layout_radio)
        
        #LAYOUT PARAMETER ABOVE
        self.layout_parameter_o=QHBoxLayout()
        self.layout_parameter_o.addLayout(self.layout_m)
        self.layout_parameter_o.addLayout(self.layout_parameter_radio)
        
        #LAYOUT PLOT
        self.layout_save=QHBoxLayout()
        self.layout_save.addWidget(self.plot_mie)
        self.layout_save.addWidget(self.clear_mie)
        self.layout_save.addWidget(self.save_mie)
        
        #CANVAS
        self.layout_canvas=MplCanvas(self.main_widget,width=8,height=6,a=a,n_m=n_m,material=material,csection=csection)
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)
        self.navi_toolbar=NavigationToolbar(self.layout_canvas,self)
        
        #LAYOUT CANVAS PARAMETER
        self.layout_cp=QVBoxLayout()
        self.layout_cp.addLayout(self.layout_save)
        self.layout_cp.addWidget(self.layout_canvas)
        self.layout_cp.addWidget(self.navi_toolbar)

        #SETTINGS LABEL
        self.settings=QLabel("----------------------------------------------------------------------------------SETTINGS----------------------------------------------------------------------------------",self)
        self.settings.setAlignment(Qt.AlignCenter)

        #OVERALL LAYOUT
        self.layout_all=QVBoxLayout(self.main_widget)
        self.layout_all.addLayout(self.layout_cp)
        self.layout_all.addWidget(self.settings)
        self.layout_all.addLayout(self.layout_parameter_o)
        self.layout_all.addLayout(self.layout_parameter)
              
        #CONNECT
        self.push_d_down.clicked.connect(self.connect_d_down)
        self.push_d_up.clicked.connect(self.connect_d_up)
        self.slider_d.valueChanged.connect(self.connect_d)
        self.push_n_down.clicked.connect(self.connect_n_down)
        self.push_n_up.clicked.connect(self.connect_n_up)
        self.slider_n.valueChanged.connect(self.connect_n)
        self.plot_mie.clicked.connect(self.connect_plot)
        self.save_mie.clicked.connect(self.connect_save)
        self.cb_m.activated.connect(self.connect_cb)
        self.clear_mie.clicked.connect(self.clear_gui)
       
    def connect_d(self):
        a=self.slider_d.value()
        self.edit_d.setText(str(a))
        return a
    
    def connect_d_down(self):
        return self.slider_d.setValue(self.slider_d.value()-5)
    
    def connect_d_up(self):
        return self.slider_d.setValue(self.slider_d.value()+5)
    
    def connect_n(self):
        n_m=self.slider_n.value()/100
        self.edit_n.setText(str(n_m))
        return n_m
    
    def connect_n_down(self):
        return self.slider_n.setValue(self.slider_n.value()-10)
    
    def connect_n_up(self):
        return self.slider_n.setValue(self.slider_n.value()+10)

    def connect_cb(self):
        material=self.cb_m.currentText()
        #print(material)
        return material

    def connect_cs(self):
        if self.radiobutton_sca.isChecked()==True:
            csection="Scattering"
            return csection
        elif self.radiobutton_abs.isChecked()==True:
            csection="Absorption"
            return csection
        elif self.radiobutton_ext.isChecked()==True:
            csection="Extinction"
            return csection
    
    def connect_plot(self):
        return self.layout_canvas.drawPlot(self.connect_d(),self.connect_n(),self.connect_cb(),self.connect_cs())

    def saveFileDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()","","All Files (*);;Text Files (*.txt)", options=options)
        f=np.savetxt(fileName+'.txt',self.layout_canvas.save_File().T,fmt='%.3f',delimiter=',',header='Wavelength [nm],Cross section [nm^2]')
         
    def connect_save(self):
        #print(self.layout_canvas.save_File())
        try:
            self.saveFileDialog()
        except:
            print("No data!")
    
    def clear_gui(self):
        return self.layout_canvas.clearPlot()
               
if __name__=="__main__":
    if not QApplication.instance():
        app=QApplication(sys.argv)
    else:
        app=QApplication.instance()
    MyApp=MainApp()
    #MyApp.showFullScreen()
    MyApp.show()
    sys.exit(app.exec_())
