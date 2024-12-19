#!/usr/bin/env python3
import copy
import math
import os,sys
import datetime
import uuid

import yaml


import warnings
import matplotlib, numpy
import matplotlib.patches as patches
from matplotlib.figure import Figure
from matplotlib.pyplot import Circle

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

from PyQt6.QtWidgets import QApplication, QMainWindow, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget, QDialog, QGridLayout, QPushButton, QComboBox, QLabel, QLineEdit, QCheckBox, QDateEdit, QTimeEdit, QDateTimeEdit, QFileDialog
from PyQt6.QtCore import Qt, QTime, QDate, QDateTime
from PyQt6.QtGui import QFont, QColor

from astropy import units
from astropy.utils.exceptions import AstropyWarning
from astropy.time import Time
from astropy.coordinates import EarthLocation, Angle, get_sun, get_moon
from astropy.coordinates import SkyCoord, AltAz


warnings.simplefilter('ignore', category=AstropyWarning)

class OM_Gui(QWidget):
    def __init__(self, args, parent=None):
        super().__init__()

        self.cwd = os.getcwd()  # curent working directory
        self.pwd = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # app location

        print(self.cwd)
        print(self.pwd)

        if os.path.exists(self.pwd+'/config.yaml'):
            with open(self.pwd+'/config.yaml', 'r') as cfg_file:
                self.cfg = yaml.safe_load(cfg_file)
        else:
            print("File not found: config.yaml")
            sys.exit()

        #print(self.cfg)

        self.i = 1
        self.mkUI()
        self.tel = self.tel_s.currentText()


    def update_table(self):
        try:
            self.table.cellChanged.disconnect(self.data_edit)
            self.table.cellClicked.disconnect(self.pocisniecie_tabelki)
        except TypeError:
            pass

        font = QFont()
        font.setPointSize(10)  # Ustawienie mniejszej czcionki
        self.table.setFont(font)

        self.table.setColumnCount(len(self.cfg["columns"]))
        for n,col_name in enumerate(self.cfg["columns"]):
            self.table.setHorizontalHeaderItem(n,QTableWidgetItem(col_name))


        for ob in self.ob:
            ob["index"] = -2
        i = -1
        self.table.setRowCount(0)
        self.table.clearContents()
        for ob in self.ob:
            if "name" in ob.keys() and ob["show"] and ob["active"]:
                i = i + 1
                ob["index"] = i
                if self.table.rowCount() <= i:
                    self.table.insertRow(i)  # Dodanie nowego wiersza
                for j,key in enumerate(self.cfg["columns"]):
                    if key in ob.keys():
                        item = QTableWidgetItem(str(ob[key]))
                    else:
                        item = QTableWidgetItem("")
                    if "editted" in ob.keys():
                        if j in ob["editted"]:
                            item.setBackground(QColor(170, 220, 240))
                    self.table.setItem(i, j, item)



        self.table.resizeColumnsToContents()
        max_column_width = 100  # Maksymalna szerokość kolumny
        for col in range(self.table.columnCount()):
            self.table.setColumnWidth(col, min(self.table.columnWidth(col), max_column_width))

        self.table.cellChanged.connect(self.data_edit)
        self.table.cellClicked.connect(self.pocisniecie_tabelki)

    def filter_changed(self):
        txt1 = self.filter_name_e.text()
        txt2 = self.filter_other_e.text()
        txt3 = self.filter_pi_e.text()
        txt4 = self.filter_sci_e.text()
        txt5 = self.filter_tag_e.text()
        for ob in self.ob:
            ob["show"] = True
            if "name" in ob.keys():
                if txt1.lower() in ob["name"].lower():
                    pass
                else:
                    ob["show"] = False
                if txt2.lower() in str(ob["line"]).lower():
                    pass
                else:
                    ob["show"] = False
                if txt3.lower() in str(ob["pi"]).lower():
                    pass
                else:
                    ob["show"] = False
                if txt4.lower() in str(ob["sciprog"]).lower():
                    pass
                else:
                    ob["show"] = False
                if txt5.lower() in str(ob["tag"]).lower():
                    pass
                else:
                    ob["show"] = False
        self.update_table()


    def load_objects(self):
        self.ob = []
        with open(self.master_file, 'r') as plik:
            for line in plik:
                tmp={}
                tmp["line"] = line
                tmp["show"] = False
                tmp["active"] = False
                if len(line.strip()) > 0:
                    l = line.split()
                    if "#" not in l[0]:
                        tmp["active"] = True
                        tmp["name"] = l[0]
                        tmp["ra"] = l[1]
                        tmp["dec"] = l[2]
                        tmp["other"] = []
                        tmp["show"] = True
                        tmp["index"] = None
                        tmp["editted"] = []
                        tmp["pi"] = ""
                        tmp["sciprog"] = ""
                        tmp["tag"] = ""
                        for i,w in enumerate(l):
                            if i >= 3:
                                if "seq=" in w:
                                    tmp["seq"] = w.split("=")[1]
                                elif "priority=" in w:
                                    tmp["priority"] = w.split("=")[1]
                                elif "cycle=" in w:
                                    tmp["cycle"] = w.split("=")[1]
                                elif "P=" in w:
                                    tmp["P"] = w.split("=")[1]
                                elif "hjd0=" in w:
                                    tmp["hjd0"] = w.split("=")[1]
                                elif "ph_mk=" in w:
                                    tmp["ph_mk"] = w.split("=")[1]
                                elif "ph_start=" in w:
                                    tmp["ph_start"] = w.split("=")[1]
                                elif "ph_end=" in w:
                                    tmp["ph_end"] = w.split("=")[1]
                                elif "t_start=" in w:
                                    tmp["t_start"] = w.split("=")[1]
                                elif "t_end=" in w:
                                    tmp["t_end"] = w.split("=")[1]
                                elif "uobi=" in w:
                                    tmp["uobi"] = w.split("=")[1]
                                elif "sciprog=" in w:
                                    tmp["sciprog"] = w.split("=")[1]
                                elif "pi=" in w:
                                    tmp["pi"] = w.split("=")[1]
                                elif "tag=" in w:
                                    tmp["tag"] = w.split("=")[1]
                                else:
                                    tmp["other"].append(w)
                self.ob.append(tmp)

    def plot_sky_map(self):
        self.sky_window = SkyWindow(self)
        self.sky_window.show()
        self.sky_window.raise_()


    def plot_data(self):
        i = int(self.table.currentRow())
        i_tab = [int(ob["index"]) for ob in self.ob]

        try:
            n = i_tab.index(i)
            target = self.ob[n]["name"]

            self.phase_window = PhaseWindow(self,target,self.cfg["tel"][self.tel]["data_file"],self.ob[n])
            self.phase_window.show()
            self.phase_window.raise_()
        except ValueError:
            pass

    def time_changed(self):
        try:
            self.sky_window.updateMap()
            self.sky_window.raise_()
        except AttributeError:
            pass
        try:
            self.phase_window.refresh()
            self.phase_window.raise_()
        except AttributeError:
            pass

    def date_changed(self):
        try:
            self.sky_window.updateMap()
            self.sky_window.raise_()
        except AttributeError:
            pass
        try:
            self.phase_window.refresh()
            self.phase_window.raise_()
        except AttributeError:
            pass

    def update_selection(self):
        self.table.selectRow(self.i)

    def pocisniecie_tabelki(self,i,j):
        self.i=i
        self.update_selection()
        try:
            self.sky_window.updateMap()
        except AttributeError:
            pass

    def update_tel(self):
        self.tel = self.tel_s.currentText()
        try:
            #self.load_objects()
            self.update_table()
        except AttributeError:
            pass

    def data_edit(self,i_selected,j_selected):
        if self.all_c.isChecked():
            key = self.cfg["columns"][j_selected]
            txt = self.table.item(i_selected, j_selected).text()
            indx = [x["index"] for x in self.ob]
            for n in range(self.table.rowCount()):
                i = indx.index(n)
                self.ob[i][key] = txt
                self.ob[i]["editted"].append(j_selected)
            self.update_table()
        else:
            indx = [x["index"] for x in self.ob]
            i = indx.index(i_selected)
            key = self.cfg["columns"][j_selected]
            self.ob[i][key] = self.table.item(i_selected, j_selected).text()
            self.ob[i]["editted"].append(j_selected)
            self.update_table()

    def copy_ob(self):
        indx = [x["index"] for x in self.ob]
        i = indx.index(self.i)
        tmp = copy.deepcopy(self.ob[i])
        self.ob.insert(i + 1, tmp)
        if "uobi" in self.ob[i+1].keys():
            if len(self.ob[i+1]["uobi"])>2:
                self.ob[i+1]["uobi"] = str(uuid.uuid4())[:8]
        for j,_ in enumerate(self.ob[i+1].keys()):
            self.ob[i+1]["editted"].append(j)
        self.update_table()

    def save_file(self):
        file_path, _ = QFileDialog.getSaveFileName(self, "Save File", self.cfg["master_file"],"Text Files (*.txt);;All Files (*)")
        if file_path:
            try:
                with open(file_path, "w", encoding="utf-8") as file:
                    txt = ""
                    for ob in self.ob:
                        if ob["active"]:
                            line = ""
                            line = line + f'{ob["name"]:20}    '
                            line = line + f'{ob["ra"]:15}    '
                            line = line + f'{ob["dec"]:15}    '
                            if "seq" in ob.keys():
                                if len(ob["seq"].strip())<10:
                                    line = line + f'seq={ob["seq"].strip():10}    '
                                elif  len(ob["seq"].strip())>9 and len(ob["seq"].strip())<30:
                                    line = line + f'seq={ob["seq"].strip():30}    '
                                else:
                                    line = line + f'seq={ob["seq"].strip():60}    '
                            if "priority" in ob.keys():
                                line = line + f'priority={str(ob["priority"]).strip():3}    '
                            if "cycle" in ob.keys():
                                if len(ob["cycle"].strip())>0:
                                    line = line + f'cycle={str(ob["cycle"]).strip():3}    '
                            if "ph_mk" in ob.keys():
                                if len(ob["ph_mk"].strip())>0:
                                    line = line + f'ph_mk={ob["ph_mk"].strip():10}    '
                            if "ph_start" in ob.keys():
                                if len(ob["ph_start"].strip())>0:
                                    line = line + f'ph_start={ob["ph_start"].strip():10}    '
                            if "ph_end" in ob.keys():
                                if len(ob["ph_end"].strip())>0:
                                    line = line + f'ph_end={ob["ph_end"].strip():10}    '
                            if "t_start" in ob.keys():
                                if len(ob["t_start"].strip())>0:
                                    line = line + f't_start={ob["t_start"].strip():10}    '
                            if "t_end" in ob.keys():
                                if len(ob["t_end"].strip())>0:
                                    line = line + f't_end={ob["t_end"].strip():10}    '
                            if "uobi" in ob.keys():
                                if len(ob["uobi"].strip())>0:
                                    line = line + f'uobi={ob["uobi"].strip():10}    '
                            if "sciprog" in ob.keys():
                                if len(ob["sciprog"].strip())>0:
                                    line = line + f'sciprog={ob["sciprog"].strip():10}    '
                            if "pi" in ob.keys():
                                if len(ob["pi"].strip())>0:
                                    line = line + f'pi={ob["pi"].strip():10}    '
                            if "tag" in ob.keys():
                                if len(ob["tag"].strip())>0:
                                    line = line + f'tag={ob["tag"].strip():10}    '
                            if "P" in ob.keys():
                                if len(ob["P"].strip())>0:
                                    line = line + f'P={ob["P"].strip():10}    '
                            if "hjd0" in ob.keys():
                                if len(ob["hjd0"].strip())>0:
                                    line = line + f'hjd0={ob["hjd0"].strip():10}    '

                            if "other" in ob.keys():
                                if len(ob["other"])>0:
                                    for a in ob["other"]:
                                        line = line + a + " "

                            txt = txt + line + "\n"
                        else:
                            txt = txt + ob["line"]


# columns:  ["uobi","sciprog","pi","tag","other"]

                    file.write(txt)
                    print(f'objects saved to {file_path}')
            except Exception as e:
                print(f"Error saving file: {e}")

    def load_file(self):
        file_path, _ = QFileDialog.getOpenFileName(None,"Select a File",self.cfg["master_file"],"All Files (*);;Text Files (*.txt);;Images (*.png *.jpg)")
        if file_path:
            self.master_file = file_path
            self.load_objects()
            self.update_table()

    def mkUI(self):
        self.setWindowTitle('OCM observing plan manager')
        self.setGeometry(50, 50, 1400, 800)

        grid = QGridLayout()

        w = 0
        self.tel_s = QComboBox()
        self.tel_s.currentIndexChanged.connect(self.update_tel)
        self.tel_s.addItems(self.cfg["tel"].keys())
        self.date_l = QLabel("UTC:")
        self.date_e = QDateEdit(self)
        self.date_e.setCalendarPopup(True)
        self.date_e.dateChanged.connect(self.date_changed)
        self.time_e = QTimeEdit(self)
        self.time_e.timeChanged.connect(self.time_changed)

        utc_now = QDateTime.currentDateTimeUtc()
        self.date_e.setDate(utc_now.date())
        self.time_e.setTime(utc_now.time())

        grid.addWidget(self.tel_s, w, 0)
        grid.addWidget(self.date_l, w, 1)
        grid.addWidget(self.date_e, w, 2)
        grid.addWidget(self.time_e, w, 3)

        w = w + 1
        self.filter_name_l = QLabel("Filter NAME")
        self.filter_name_e = QLineEdit("")
        self.filter_name_e.textChanged.connect(self.filter_changed)
        grid.addWidget(self.filter_name_l, w, 0)
        grid.addWidget(self.filter_name_e, w, 1)

        self.filter_sci_l = QLabel("Filter SCIPROG")
        self.filter_sci_e = QLineEdit("")
        self.filter_sci_e.textChanged.connect(self.filter_changed)
        grid.addWidget(self.filter_sci_l, w, 2)
        grid.addWidget(self.filter_sci_e, w, 3)

        w = w + 1
        self.filter_pi_l = QLabel("Filter PI")
        self.filter_pi_e = QLineEdit("")
        self.filter_pi_e.textChanged.connect(self.filter_changed)
        grid.addWidget(self.filter_pi_l, w, 0)
        grid.addWidget(self.filter_pi_e, w, 1)

        self.filter_tag_l = QLabel("Filter TAG")
        self.filter_tag_e = QLineEdit("")
        self.filter_tag_e.textChanged.connect(self.filter_changed)
        grid.addWidget(self.filter_tag_l, w, 2)
        grid.addWidget(self.filter_tag_e, w, 3)

        w = w + 1
        self.filter_other_l = QLabel("Filter OB")
        self.filter_other_e = QLineEdit("")
        self.filter_other_e.textChanged.connect(self.filter_changed)
        grid.addWidget(self.filter_other_l, w, 0)
        grid.addWidget(self.filter_other_e, w, 1)

        self.all_c = QCheckBox("Edit Column")
        self.all_c.setChecked(False)
        grid.addWidget(self.all_c, w, 3)


        w = w + 1
        self.table = QTableWidget()
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        #self.table.setSelectionMode(QTableWidget.SelectionBehavior.SingleSelection)
        self.table.setStyleSheet("selection-background-color: rgb(217,239,217); selection-color: black; ")

        grid.addWidget(self.table, w, 0, 1, 4)

        w = w + 1
        self.sky_p = QPushButton("Plot SkyMap")
        self.sky_p.clicked.connect(self.plot_sky_map)
        grid.addWidget(self.sky_p, w, 3)

        self.deactivate_p = QPushButton("Deactivate")
        grid.addWidget(self.deactivate_p, w, 0)

        w = w + 1
        self.data_p = QPushButton("Plot data")
        self.data_p.clicked.connect(self.plot_data)
        grid.addWidget(self.data_p, w, 3)

        self.copy_p = QPushButton("Copy")
        self.copy_p.clicked.connect(self.copy_ob)
        grid.addWidget(self.copy_p, w, 0)

        w = w + 1
        self.load_p = QPushButton("Load file")
        self.load_p.clicked.connect(self.load_file)
        self.save_p = QPushButton("Save")
        self.save_p.clicked.connect(self.save_file)
        self.config_p = QPushButton("\u2699")
        self.close_p = QPushButton("Close")
        self.close_p.clicked.connect(self.close)
        grid.addWidget(self.load_p, w, 0)
        grid.addWidget(self.config_p, w, 1)
        grid.addWidget(self.save_p, w, 2)

        w = w + 1
        grid.addWidget(self.close_p, w, 3)

        self.setLayout(grid)

        self.show()

class PhaseWindow(QWidget):
    def __init__(self, parent, target, data_dir, ob):
        super(PhaseWindow, self).__init__()
        self.parent = parent
        self.target = target
        self.data_dir = data_dir
        self.ob = ob

        self.setStyleSheet("font-size: 11pt;")
        self.setMinimumSize(1000,400)
        self.mkUI()
        try:
            self.get_object()
        except FileNotFoundError:
            pass
        self.refresh()


    def get_object(self):

        self.f_path = self.data_dir+self.target.lower()
        self.filters = os.listdir(self.f_path)
        #print(self.filters)
        self.file_s.addItems(self.filters)

    def refresh(self):

        obs_time = datetime.datetime.combine(self.parent.date_e.date().toPyDate(), self.parent.time_e.time().toPyTime())
        time = Time(obs_time, scale='utc')
        self.current_jd = time.jd

        self.axes.clear()
        try:
            filter = self.file_s.currentText()
            file = self.f_path+"/"+filter+"/light-curve/"+self.target.lower()+"_"+filter+"_diff_light_curve.txt"

            mag = []
            jd = []
            flag = []

            with open(file, "r") as plik:
                if plik != None:
                    for line in plik:
                        if len(line.strip()) > 0:
                            try:
                                mag.append(float(line.split()[1]))
                                jd.append(float(line.split()[3]))
                                try:
                                    flag.append(int(line.split()[9]))
                                except:
                                    flag.append(int(0))
                            except ValueError:
                                pass
            if len(mag) == len(jd) and len(jd)>0:
                if self.phase_c.isChecked():

                    if "P" in self.ob.keys():
                        self.P = self.ob["P"]
                        if "hjd0" in self.ob.keys():
                            self.jd0 = self.ob["hjd0"]
                        else:
                            self.jd0 = 2460000

                        jd = (numpy.array(jd) - float(self.jd0))/float(self.P)%1

                        self.current_jd = (self.current_jd - float(self.jd0))/float(self.P)%1

                        if "ph_start" in self.ob.keys() and "ph_end" in self.ob.keys():
                            t0 = float(self.ob["ph_start"])
                            t1 = float(self.ob["ph_end"])

                            if t0 < t1:
                                self.axes.axvspan(0, t0, color='red', alpha=0.05)
                                self.axes.axvspan(t1, 1, color='red', alpha=0.05)
                            else:
                                self.axes.axvspan(t1, t0, color='red', alpha=0.05)

                        if "ph_mk" in self.ob.keys():
                            bin = float(self.ob["ph_mk"].split("/")[2])
                            fitr = self.ob["ph_mk"].split("/")[1]
                            n_obs = float(self.ob["ph_mk"].split("/")[0])

                            covered = []
                            maska_quality = numpy.array(flag) < 2
                            for i,t in enumerate(jd[maska_quality]):
                                mk1 = t > jd[maska_quality] - bin/2
                                mk2 = t < jd[maska_quality] + bin/2
                                mk = numpy.array([a and b for a,b in zip(mk1,mk2)])
                                if len(mk[mk]) >= n_obs:
                                    covered.append([t-bin/2.,t+bin/2.])

                            covered.sort(key=lambda x: x[0])
                            merged = []
                            for odcinek in covered:
                                if not merged or merged[-1][1] < odcinek[0]:
                                    merged.append(odcinek)
                                else:
                                    merged[-1][1] = max(merged[-1][1], odcinek[1])

                            for x in merged:
                                self.axes.axvspan(x[0], x[1], color='red', alpha=0.05)

                        self.axes.set_xlim(-0.1,1.1)
                        self.axes.set_title(f"{self.target} P={self.P}")

                    else:
                        self.phase_c.setChecked(False)
                        self.axes.set_title(f"{self.target}")
                else:
                    self.axes.set_title(f"{self.target}")

                mk = numpy.array(flag) == 0
                x = numpy.array(jd)[mk]
                y = numpy.array(mag)[mk]
                self.axes.plot(x,y,".g",alpha=0.5)

                mk = numpy.array(flag) == 1
                x = numpy.array(jd)[mk]
                y = numpy.array(mag)[mk]
                self.axes.plot(x,y,".c",alpha=0.5)

                mk = numpy.array(flag) == 2
                x = numpy.array(jd)[mk]
                y = numpy.array(mag)[mk]
                self.axes.plot(x,y,".k",alpha=0.5)

                d = 0.1*(max(mag)-min(mag))
                self.axes.set_ylim(max(mag)+d,min(mag)-d)

                self.axes.axvline(x=self.current_jd, color="blue")

        except (FileNotFoundError,ValueError) as e:
            print(f"Phase Window Error: {e}")

        #self.fig.tight_layout()
        self.canvas.draw()
        self.show()



    def mkUI(self):
        grid = QGridLayout()

        self.file_s = QComboBox()
        self.file_s.currentIndexChanged.connect(self.refresh)
        self.phase_c = QCheckBox("Phase")
        self.phase_c.setChecked(True)
        #self.phase_c.setStyleSheet("QCheckBox::indicator:checked {image: url(./Icons/SwitchOn.png)}::indicator:unchecked {image: url(./Icons/SwitchOff.png)}")
        self.phase_c.clicked.connect(self.refresh)


        self.fig = Figure((2.0, 2.0), linewidth=-1, dpi=100)
        self.canvas = FigureCanvas(self.fig)
        self.axes = self.fig.add_subplot(111)
        grid.addWidget(self.file_s, 0, 0)
        grid.addWidget(self.phase_c, 0, 3)
        grid.addWidget(self.canvas,1,0,4,4)

        self.toolbar = NavigationToolbar(self.canvas,self)
        grid.addWidget(self.toolbar, 5, 0, 1, 4)

        self.close_p = QPushButton('Close')
        self.close_p.clicked.connect(lambda: self.close())
        grid.addWidget(self.close_p, 6, 3)

        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 0)
        grid.setRowStretch(0, 1)
        grid.setRowStretch(1, 0)
        grid.setRowStretch(2, 0)

        self.setLayout(grid)

#################################
#          SKY WINDOW           #
#################################

class SkyWindow(QWidget):
    def __init__(self, parent):
        super(SkyWindow, self).__init__()
        self.parent = parent

        self.setStyleSheet("font-size: 11pt;")
        self.setMinimumSize(900,400)
        self.mkUI()
        cid = self.canvas.mpl_connect('button_press_event', self.zaznaczenie)
        self.updateMap()
    def mkUI(self):
        grid = QGridLayout()

        self.fig = Figure((2.0, 2.0), linewidth=-1, dpi=100)
        self.canvas = FigureCanvas(self.fig)
        self.axes = self.fig.add_subplot(121,polar=True)
        self.axes2 = self.fig.add_subplot(122,polar=True)
        grid.addWidget(self.canvas,0,0,1,1)

        self.close_p = QPushButton('Close')
        self.close_p.clicked.connect(lambda: self.close())
        grid.addWidget(self.close_p, 1, 0)

        self.setLayout(grid)
    def updateMap(self):
        self.axes.clear()
        #self.axes.set_theta_direction(-1)
        self.axes.set_theta_zero_location('N')
        #self.axes.set_ylim([0, 360])
        #self.axes.set_rlim([0, 30])
        self.axes.set_xticks([0, 2 * 3.14 * 90 / 360, 2 * 3.14 * 180 / 360, 2 * 3.14 * 270 / 360])
        self.axes.set_xticklabels(["N", "E", "S", "W"])
        #self.axes.set_rmax(self.rmax)
        #self.axes.set_rticks([0, 20, 40, 60, 90])
        #self.axes.set_yticklabels(["", "", "", "", ""])

        self.axes.set_rticks([])
        self.axes.set_yticklabels([])

        #self.axes.bar(0, self.rmax - 90, width=2 * math.pi, bottom=90, color='k', alpha=0.05)  # tutaj zmienia sie pasek ponizej horyzoontu
        self.axes.set_rlim([-90,90+self.parent.cfg["obs_latitude"]])


        self.axes2.clear()
        self.axes2.set_theta_direction(-1)
        self.axes2.set_theta_zero_location('N')
        self.axes2.set_xticks([0, 2 * 3.14 * 90 / 360, 2 * 3.14 * 180 / 360, 2 * 3.14 * 270 / 360])
        self.axes2.set_xticklabels(["N", "E", "S", "W"])
        #self.axes.set_rmax(self.rmax)
        self.axes2.set_rticks([0, 20, 40, 60, 90])
        self.axes2.set_yticklabels(["", "", "", "", ""])


        obs_location = EarthLocation(lat=self.parent.cfg["obs_latitude"], lon=self.parent.cfg["obs_longitude"], height=self.parent.cfg["obs_elevation"])  # Warszawa

        obs_time = datetime.datetime.combine(self.parent.date_e.date().toPyDate(), self.parent.time_e.time().toPyTime())

        ra = [x["ra"] for x in self.parent.ob if x["show"] and "ra" in x.keys()]
        dec = [x["dec"] for x in self.parent.ob if x["show"] and "dec" in x.keys()]

        coords = SkyCoord(ra=ra, dec=dec, unit=('hourangle', 'deg'), frame='icrs')
        altaz = coords.transform_to(AltAz(obstime=Time(obs_time), location=obs_location))

        azimuth = altaz.az.deg
        self.altitude = altaz.alt.deg
        self.azimuth_rad = numpy.radians(azimuth)
        self.axes2.scatter(self.azimuth_rad, self.altitude,marker="*", c='green')
        self.axes2.plot(self.azimuth_rad[self.parent.i],self.altitude[self.parent.i],"*r")

        self.axes2.set_rlim([90, 0])


        lst = Time(obs_time,scale='utc', location=obs_location).sidereal_time('mean')
        ha = [15*(lst.hour - (float(x.split(":")[0])+float(x.split(":")[1])/60.+float(x.split(":")[2])/3600.)) for x in ra]

        self.dec_deg =  Angle(dec, unit="deg").deg
        self.ha_rad = numpy.radians(ha)
        self.axes.scatter(self.ha_rad, self.dec_deg,marker="*", c='green', s=50)
        self.axes.plot(self.ha_rad[self.parent.i],self.dec_deg[self.parent.i],"*r")

        altaz_frame = AltAz(obstime=Time(obs_time,scale='utc', location=obs_location), location=obs_location)
        azimuths = numpy.linspace(0, 360, 360) * units.deg
        altitude = 0 * units.deg
        altaz_coords = SkyCoord(alt=altitude, az=azimuths, frame=altaz_frame)
        icrs_coords = altaz_coords.transform_to("icrs")
        ra = icrs_coords.ra.deg  # Rektascensja w stopniach
        ra_c = 15.* lst.hour - ra
        dec_c = icrs_coords.dec.deg  # Deklinacja w stopniach
        self.axes.plot(numpy.radians(ra_c),dec_c,"--b")


        sun = get_sun(Time(obs_time))
        moon = get_moon(Time(obs_time))
        ha_sun = lst - sun.ra
        ha_moon = lst - moon.ra
        dec_sun = sun.dec
        dec_moon = moon.dec
        ha_sun_rad = ha_sun.to(units.rad).value
        ha_moon_rad = ha_moon.to(units.rad).value
        dec_sun_deg = dec_sun.to(units.deg).value
        dec_moon_deg = dec_moon.to(units.deg).value
        self.axes.plot(ha_sun_rad, dec_sun_deg, 'yo')
        self.axes.plot(ha_moon_rad, dec_moon_deg, 'ko',alpha=0.5)

        altaz_frame = AltAz(obstime=obs_time, location=obs_location)
        sun_altaz = sun.transform_to(altaz_frame)
        moon_altaz = moon.transform_to(altaz_frame)

        sun_altitude_deg = sun_altaz.alt.degree
        sun_azimuth = sun_altaz.az
        moon_altitude_deg = moon_altaz.alt.degree
        moon_azimuth = moon_altaz.az

        sun_az_rad = sun_azimuth.to(units.rad).value
        moon_az_rad = moon_azimuth.to(units.rad).value
        self.axes2.plot(sun_az_rad, sun_altitude_deg, 'yo')
        self.axes2.plot(moon_az_rad, moon_altitude_deg, 'ko',alpha=0.5)


        # Odświeżenie wykresu
        self.canvas.draw()

    def zaznaczenie(self, event):
        if event.xdata != None:
            az = float(event.xdata)
            alt = float(event.ydata)

            if event.button == 1:
                if event.inaxes == self.axes:

                    a1 = az
                    a2 = self.ha_rad
                    h1 = alt + 90
                    h2 = self.dec_deg + 90
                    delta = (h1 ** 2 + h2 ** 2 - 2 * h1 * h2 * numpy.cos(a1 - a2)) ** 0.5
                    min_i = numpy.argmin(delta)
                    if delta[min_i] < 15:
                        self.parent.i = min_i
                        self.updateMap()
                        self.parent.update_selection()

                elif event.inaxes == self.axes2:
                    a1 = az
                    a2 = self.azimuth_rad
                    h1 = alt + 180.
                    h2 = self.altitude
                    delta = (h1 ** 2 + h2 ** 2 - 2 * h1 * h2 * numpy.cos(a1 - a2)) ** 0.5
                    min_i = numpy.argmin(delta)
                    if delta[min_i] < 50:
                        self.parent.i = min_i
                        self.updateMap()
                        self.parent.update_selection()


