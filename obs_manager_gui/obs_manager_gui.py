#!/usr/bin/env python3
import copy
import math
import os,sys
import datetime
import uuid
import yaml
import requests

import warnings
import numpy
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

import ephem

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

from PyQt6.QtWidgets import QApplication, QMainWindow, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget, QDialog, QGridLayout, QPushButton, QComboBox, QLabel, QLineEdit, QTextEdit, QCheckBox, QDateEdit, QTimeEdit, QDateTimeEdit, QFileDialog
from PyQt6.QtCore import Qt, QTime, QDate, QDateTime
from PyQt6.QtGui import QFont, QColor

from astropy import units
from astropy.utils.exceptions import AstropyWarning
from astropy.time import Time
from astropy.coordinates import EarthLocation, Angle, get_sun, get_moon
from astropy.coordinates import SkyCoord, AltAz
from astropy.table import Table

from telescope_plan_generator import telescope_plan_generator as tpg


warnings.simplefilter('ignore', category=AstropyWarning)

class OM_Gui(QWidget):
    def __init__(self, args, parent=None):
        super().__init__()

        self.cwd = os.getcwd()  # curent working directory
        self.pwd = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # app location

        #print(self.cwd)
        #print(self.pwd)

        if os.path.exists(self.pwd+'/config.yaml'):
            with open(self.pwd+'/config.yaml', 'r') as cfg_file:
                self.cfg = yaml.safe_load(cfg_file)
        else:
            print("File not found: config.yaml")
            sys.exit()

        #print(self.cfg)

        self.tpg_window = None
        self.i = 1
        self.mkUI()
        self.tel = self.tel_s.currentText()
        self.update_almanac()


    def update_table(self):
        row_labels = []

        try:
            self.table.cellChanged.disconnect(self.data_edit)
            self.table.cellClicked.disconnect(self.pocisniecie_tabelki)
        except TypeError:
            pass

        m = "-99"
        tpg = False
        tmp_txt = ""
        if self.tpg_window:
            try:
                tpg_ob = self.tpg_window.p.ob
                tpg_ind_list = [t["index"] for t in tpg_ob]

                date = self.date_e.date().toPyDate()
                ut = self.time_e.time().toPyTime()
                dt = datetime.datetime.combine(date, ut)

                current_time = ephem.Date(dt)
                tpg_time = self.tpg_window.p.nightTime
                time_diff = numpy.abs([t - current_time for t in tpg_time])

                m = numpy.argmin(time_diff)

                #print(ephem.Date(current_time),ephem.Date(tpg_time[m]),time_diff[m])


                tpg = True
            except AttributeError:
                tpg = False

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

                # pisze ile widoczne i czy teraz widoczne
                red_znacznik = False
                row_txt = f'{i}'

                if tpg:
                    i_tmp = ob["index"]
                    q = tpg_ind_list.index(i_tmp)
                    vis = tpg_ob[q]["visibility"]["all"]

                    if not vis[m] and m != 0:
                        # if "113" in tpg_ob[q]["name"]:
                        #     print(vis)
                        #     print(m,vis[m])
                        red_znacznik = True

                    vis = numpy.array(vis)
                    t_vis = len(vis[vis])
                    h_vis = t_vis / 60
                    m_vis = t_vis - int(h_vis)*60

                    row_txt = row_txt + f' ({int(h_vis)}h {m_vis}m) '

                for j,key in enumerate(self.cfg["columns"]):
                    if key in ob.keys():
                        item = QTableWidgetItem(str(ob[key]))
                        if red_znacznik:
                            item.setForeground(QColor("red"))
                        else:
                            item.setForeground(QColor("black"))


                    else:
                        item = QTableWidgetItem("")
                    if "editted" in ob.keys():
                        if j in ob["editted"]:
                            item.setBackground(QColor(170, 220, 240))
                    if "deactivated" in ob.keys():
                        if ob["deactivated"]:
                            item.setBackground(QColor(150, 150, 150))
                    self.table.setItem(i, j, item)
                row_labels.append(row_txt)

        self.table.setVerticalHeaderLabels(row_labels)

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
                line = line.strip()
                tmp={}
                tmp["line"] = line
                tmp["show"] = False
                tmp["active"] = False
                if len(line.split()) > 0:

                    l = line.strip().split()
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
                        tmp["comment"] = ""

                        line = line.replace(tmp["name"], "")
                        line = line.replace(tmp["ra"], "")
                        line = line.replace(tmp["dec"], "")

                        for x in self.cfg["columns"]:
                            if x+"=" in line and "comment" not in x:
                                tmp[x] = line.split(x+"=")[1].split()[0]
                                line = line.replace(f'{x}={tmp[x]}', "")

                        if "comment=" in line:
                            comment = line.split("comment=")[1]
                            if len(comment)>0:
                                if "\"" == comment[0]:
                                    comment = comment.split("\"")[1]
                                    line = line.replace(f'comment=\"{comment}\"', "")
                                elif "(" == comment[0]:
                                    comment = comment.split("(")[1].split(")")[0]
                                    line = line.replace(f'comment=({comment})', "")
                                else:
                                    comment = comment.split()[0]
                                    line = line.replace(f'comment={comment}', "")
                            tmp["comment"] = comment
                        tmp["other"] = line.strip()
                    self.ob.append(tmp)
                else:
                    tmp["line"] = "\n"
                    self.ob.append(tmp)

    def plot_sky_map(self):
        try:
            self.sky_window = SkyWindow(self)
            self.sky_window.show()
            self.sky_window.raise_()
        except:
            pass


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
            if self.sky_window.isVisible():
                self.sky_window.updateMap()
                self.sky_window.raise_()
        except AttributeError:
            pass
        try:
            if self.phase_window.isVisible():
                self.phase_window.refresh()
                self.phase_window.raise_()
        except AttributeError:
            pass

    def date_changed(self):
        try:
            if self.sky_window.isVisible():
                self.sky_window.updateMap()
                self.sky_window.raise_()
        except AttributeError:
            pass
        try:
            if self.phase_window.isVisible():
                self.phase_window.refresh()
                self.phase_window.raise_()
        except AttributeError:
            pass

    def update_almanac(self):
        obs_time = datetime.datetime.combine(self.date_e.date().toPyDate(), self.time_e.time().toPyTime())
        time = Time(obs_time, scale='utc')
        self.almanac = sun_moon_ephem(time, self.cfg["obs_latitude"], self.cfg["obs_longitude"], self.cfg["obs_elevation"], horizon=0*units.deg)
        # print(time.to_datetime() )
        # print(self.almanac["next_sunrise"])
        # print( self.almanac["next_sunrise"] - time.to_datetime() )

        txt = ""
        txt = txt + f'sunset: {self.almanac["next_sunset"]}\n'
        txt = txt + f'sunrise: {self.almanac["next_sunrise"]}\n'
        txt = txt + f'moon: {self.almanac["moon_phase"]}\n'
        txt = txt + f'moonrise: {self.almanac["next_moonrise"]}\n'
        txt = txt + f'moonset: {self.almanac["next_moonset"]}\n'



        self.almanac_e.setText(txt)
        # "julian_date": jd,
        # "prev_sunrise": prev_sunrise,
        # "next_sunrise": next_sunrise,
        # "prev_sunset": prev_sunset,
        # "next_sunset": next_sunset,
        # "prev_moonrise": prev_moonrise,
        # "next_moonrise": next_moonrise,
        # "prev_moonset": prev_moonset,
        # "next_moonset": next_moonset,
        # "moon_phase": moon_ph

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
            self.all_c.setChecked(False)
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

    def deactivate(self):
        indx = [x["index"] for x in self.ob]
        i = indx.index(self.i)
        tmp = copy.deepcopy(self.ob[i])
        self.ob[i]["deactivated"] = True
        self.ob[i]["name"] = "# "+self.ob[i]["name"]
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

                            for x in self.cfg["columns"]:
                                if x in ob.keys() and x not in ["name","ra","dec","seq","comment"]:
                                    if len(ob[x].strip()) > 0:
                                        line = line + f'{x}={str(ob[x]).strip()}    '

                            if "comment" in ob.keys():
                                if len(ob["comment"].strip())>0:
                                    tmp = ob["comment"]
                                    line = line + f'comment=\"{tmp}\" '

                            txt = txt + line
                        else:
                            txt = txt + ob["line"]
                        if not txt.endswith("\n"):
                            txt = txt + "\n"


# columns:  ["uobi","sciprog","pi","tag","other"]
                    #print(txt)
                    file.write(txt)
                    print(f'objects saved to {file_path}')
            except Exception as e:
                print(f"Error saving file: {e}")

    def update_data(self):
        base="https://araucaria.camk.edu.pl/data/ocm/internal/Nzg3N2dVlZmQtNGIy/"

        #zb08 / targets / u_lep / Ic / light - curve / u_lep_Ic_diff_light_curve.txt
        tel = self.tel_s.currentText()
        if hasattr(self, "ob"):
            for ob in self.ob:
                if "obs_data" in ob.keys():
                    print(ob["obs_data"])
                elif "name" in ob.keys():
                    name = ob["name"].lower()
                    if "ph_mk" in ob.keys():
                        filtr = ob["ph_mk"].split("/")[1]
                    elif "seq" in ob.keys():
                        filtr = ob["seq"].split("/")[1]
                    else:
                        filtr = "V"
                    path = tel+"/targets/"+name+"/"+filtr+"/light-curve/"+name+"_"+filtr+"_diff_light_curve.txt"
                    url = base+path
                    response = requests.get(url)
                    if response.status_code == 200:
                        f_name = "./obs_data/"+path
                        os.makedirs(os.path.dirname(f_name), exist_ok=True)
                        with open(f_name, "wb") as f:
                            f.write(response.content)
                        print(f'file {path.split("/")[-1]} updated!')
                    else:
                        print(f'file {path.split("/")[-1]}', response.status_code)

    def tpg_show(self):
        self.update_almanac()
        self.tpg_window = TPGWindow(self)


    def load_file(self):
        file_path, _ = QFileDialog.getOpenFileName(None,"Select a File",self.cfg["master_file"],"All Files (*);;Text Files (*.txt);;Images (*.png *.jpg)")
        if file_path:
            self.master_file = file_path
            self.load_objects()
            for t in self.cfg["tel"].keys():
                if t in self.master_file:
                    self.tel_s.setCurrentText(t)
            self.update_table()

    def mkUI(self):
        self.setWindowTitle('OCM observing plan manager')
        self.setGeometry(50, 50, 1400, 800)

        grid = QGridLayout()


        self.almanac_e = QTextEdit()
        grid.addWidget(self.almanac_e, 0, 4,3,3)

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
        self.deactivate_p.clicked.connect(self.deactivate)
        grid.addWidget(self.deactivate_p, w, 0)

        self.deactivate_p = QPushButton("TPG")
        self.deactivate_p.clicked.connect(self.tpg_show)
        grid.addWidget(self.deactivate_p, w, 1)

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
        self.update_p = QPushButton("Update Data")
        self.update_p.clicked.connect(self.update_data)
        self.save_p = QPushButton("Save")
        self.save_p.clicked.connect(self.save_file)
        self.config_p = QPushButton("\u2699")
        self.close_p = QPushButton("Close")
        self.close_p.clicked.connect(self.close)
        grid.addWidget(self.load_p, w, 0)
        grid.addWidget(self.update_p, w, 1)
        grid.addWidget(self.config_p, w, 2)
        grid.addWidget(self.save_p, w, 3)

        w = w + 1
        grid.addWidget(self.close_p, w, 4)

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
        self.setMinimumSize(1200,600)
        self.mkUI()
        try:
            self.get_object()
        except FileNotFoundError:
            pass
        self.refresh()


    def get_object(self):

        self.f_path = self.data_dir+self.target.lower()
        self.filters = os.listdir(self.f_path)
        self.file_s.addItems(self.filters)

    def refresh(self):

        obs_time = datetime.datetime.combine(self.parent.date_e.date().toPyDate(), self.parent.time_e.time().toPyTime())
        time = Time(obs_time, scale='utc')
        self.current_jd = time.jd
        jd3h = self.current_jd + numpy.arange(1,7,1)/24.


        txt = f'moon phase: {self.parent.almanac["moon_phase"]:.0f} %'
        self.ephem_e.setText(txt)

        self.ephem_e.setStyleSheet("background-color: white;")
        if self.ob.get("max_moon_phase", False):
            if float(self.ob["max_moon_phase"]) < float(self.parent.almanac["moon_phase"]):
                self.ephem_e.setStyleSheet("background-color: lightcoral;")

        if self.ob.get("h_min", False):
            hmin = float(self.ob["h_min"])
        else:
            hmin = self.parent.cfg["tel"][self.parent.tel]["hmin"]

        if self.ob.get("h_min", False):
            hmax = float(self.ob["h_max"])
        else:
            hmax = self.parent.cfg["tel"][self.parent.tel]["hmax"]

        self.axes2.clear()
        self.axes2.set_ylim(-20, 90)
        self.axes2.set_xlim(int(self.current_jd), int(self.current_jd)+1)
        self.axes2.axvline(x=self.current_jd, color="blue")
        self.axes2.axhspan(-20, 0, xmin=0, xmax=1, facecolor='red', alpha=0.1)
        self.axes2.axhspan(0, hmin, xmin=0, xmax=1, facecolor='black', alpha=0.05)
        if self.parent.cfg["tel"][self.parent.tel]["mount_type"] == "az":
            self.axes2.axhspan(hmax, 90, xmin=0, xmax=1, facecolor='black', alpha=0.1)

        i = int(self.parent.table.currentRow())
        i_tab = [int(ob["index"]) for ob in self.parent.ob]
        n = i_tab.index(i)
        ra = self.parent.ob[n]["ra"]
        dec = self.parent.ob[n]["dec"]

        obs_location = EarthLocation(lat=self.parent.cfg["obs_latitude"], lon=self.parent.cfg["obs_longitude"], height=self.parent.cfg["obs_elevation"])
        object = SkyCoord(ra=ra, dec=dec, unit=('hourangle', 'deg'), frame='icrs')
        time_range = Time(numpy.linspace(int(self.current_jd), int(self.current_jd)+1, 100), format="jd")
        altaz_frame = AltAz(obstime=time_range, location=obs_location)
        object_altaz = object.transform_to(altaz_frame)
        alt = object_altaz.alt.deg
        sun = get_sun(time_range)
        sun_alt = sun.transform_to(altaz_frame).alt.degree
        moon = get_moon(time_range)
        moon_alt = moon.transform_to(altaz_frame).alt.degree
        sep = object_altaz.separation(moon).deg

        txt = f'separation: {min(sep):.0f}\u00B0'
        self.moon_sep_e.setText(txt)

        self.moon_sep_e.setStyleSheet("background-color: white;")
        if self.ob.get("min_moon_dist", False):
            if float(self.ob["min_moon_dist"]) > float(min(sep)):
                self.moon_sep_e.setStyleSheet("background-color: lightcoral;")


        self.axes2.plot(time_range.jd,alt,"-g")
        self.axes2.plot(time_range.jd, sun_alt, "--y")
        self.axes2.plot(time_range.jd, moon_alt, "--k")
        #self.axes2.add_patch(Rectangle((0, 0), 0.5, 10, facecolor='yellow', alpha=0.5))

        self.axes.clear()
        try:
            filter = self.file_s.currentText()
            file = self.f_path+"/"+filter+"/light-curve/"+self.target.lower()+"_"+filter+"_diff_light_curve.txt"

            mag = []
            jd = []
            flag = []

            lc_tab = Table.read(file, format="ascii")
            mag = lc_tab["mag"]
            jd = lc_tab["jd_obs"]
            flag = lc_tab["quality"]

            if len(mag) == len(jd) and len(jd)>0:

                jd = numpy.array(jd)
                mag = numpy.array(mag)
                recent_obs_mask = jd > self.current_jd - float(self.parent.cfg["last_nights_to_mark"])


                # obsluga rysowania cycle part 1.
                plot_cycle = False
                if self.ob.get("cycle", False):
                    mk = numpy.array(flag) != 2   # tylko dobre obserwacje
                    last_jd = max(numpy.array(jd)[mk])
                    end_cycle =  last_jd + float(self.ob["cycle"])
                    plot_cycle = True

                if self.phase_c.isChecked():

                    plot_cycle = False

                    if "P" in self.ob.keys():
                        self.P = self.ob["P"]
                        if "hjd0" in self.ob.keys():
                            self.jd0 = self.ob["hjd0"]
                        else:
                            self.jd0 = 2460000

                        jd = (numpy.array(jd) - float(self.jd0))/float(self.P)%1

                        self.current_jd = (self.current_jd - float(self.jd0))/float(self.P)%1
                        jd3h = (jd3h - float(self.jd0))/float(self.P)%1

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


                # last 7 days:
                jd_tmp = jd[recent_obs_mask]
                mag_tmp = mag[recent_obs_mask]
                flag_tmp = numpy.array(flag)[recent_obs_mask]

                mk = numpy.array(flag_tmp) == 0
                x = numpy.array(jd_tmp)[mk]
                y = numpy.array(mag_tmp)[mk]
                self.axes.plot(x,y,".g",alpha=1)

                mk = numpy.array(flag_tmp) == 1
                x = numpy.array(jd_tmp)[mk]
                y = numpy.array(mag_tmp)[mk]
                self.axes.plot(x,y,".c",alpha=1)

                mk = numpy.array(flag_tmp) == 2
                x = numpy.array(jd_tmp)[mk]
                y = numpy.array(mag_tmp)[mk]
                self.axes.plot(x,y,".k",alpha=0.5)


                # all dates
                mk = numpy.array(flag) == 0
                x = numpy.array(jd)[mk]
                y = numpy.array(mag)[mk]
                self.axes.plot(x,y,".g",alpha=0.1)

                mk = numpy.array(flag) == 1
                x = numpy.array(jd)[mk]
                y = numpy.array(mag)[mk]
                self.axes.plot(x,y,".c",alpha=0.1)

                mk = numpy.array(flag) == 2
                x = numpy.array(jd)[mk]
                y = numpy.array(mag)[mk]
                self.axes.plot(x,y,".k",alpha=0.05)

                d = 0.1*(max(mag)-min(mag))
                self.axes.set_ylim(max(mag)+d,min(mag)-d)

                self.axes.axvline(x=self.current_jd, color="blue")
                i_tmp = 1.
                for x in jd3h:
                    i_tmp += 1
                    self.axes.axvline(x, color="blue", alpha = 1/i_tmp)


                # obsluga rysowania cycle part 2.
                if self.ob.get("cycle", False):
                    if plot_cycle:
                        self.axes.fill_between([last_jd, end_cycle],min(mag), max(mag),color='red', alpha=0.1)

                    if end_cycle > int(self.current_jd):
                        self.axes2.fill_between([int(self.current_jd), end_cycle],-20, 90,color='red', alpha=0.1)


        except (FileNotFoundError,ValueError) as e:
            print(f"Phase Window Error: {e}")

        self.fig.tight_layout()
        self.canvas.draw()
        self.show()



    def mkUI(self):
        grid = QGridLayout()

        self.file_s = QComboBox()
        self.file_s.currentIndexChanged.connect(self.refresh)

        self.ephem_e = QLineEdit()
        self.ephem_e.setReadOnly(True)
        self.moon_sep_e = QLineEdit()
        self.moon_sep_e.setReadOnly(True)

        self.phase_c = QCheckBox("Phase")
        self.phase_c.setChecked(True)
        self.phase_c.clicked.connect(self.refresh)


        self.fig = Figure((2.0, 2.0), linewidth=-1, dpi=100)
        self.canvas = FigureCanvas(self.fig)
        self.axes = self.fig.add_subplot(211)
        self.axes2 = self.fig.add_subplot(212)
        grid.addWidget(self.file_s, 0, 0)
        grid.addWidget(self.ephem_e, 0, 1)
        grid.addWidget(self.moon_sep_e, 0, 2)
        grid.addWidget(self.phase_c, 0, 3)
        grid.addWidget(self.canvas,1,0,4,4)

        self.toolbar = NavigationToolbar(self.canvas,self)
        grid.addWidget(self.toolbar, 5, 0, 1, 4)

        self.close_p = QPushButton('Close')
        self.close_p.clicked.connect(lambda: self.close())
        grid.addWidget(self.close_p, 6, 3)

        #grid.setColumnStretch(0, 1)
        #grid.setColumnStretch(1, 1)
        #grid.setColumnStretch(2, 1)
        #grid.setColumnStretch(3, 1)
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


class TPGWindow(QWidget):
    def __init__(self, parent):
        super(TPGWindow, self).__init__()
        self.parent = parent

        self.setStyleSheet("font-size: 11pt;")
        self.setMinimumSize(100,200)
        self.mkUI()


        # DUPA
        #     p = tpg(args.tel,args.date,wind=args.wind,loud=args.loud,seed=args.seed,done_uobi=[])

    def load(self):
        tel = self.parent.tel_s.currentText()

        date = self.parent.date_e.date().toPyDate()
        ut =   self.parent.time_e.time().toPyTime()

        t0 = datetime.datetime.combine(date, ut)

        t1 = self.parent.almanac["next_sunset"]
        t2 = self.parent.almanac["prev_sunset"]

        if (t1 - t0) > (t2 - t0):
            dt = [str(date - datetime.timedelta(days=1))]
        else:
            dt = str(date)

        self.p = tpg(tel, dt)

        self.p.Initiate()
        self.p.LoadObjects()
        self.p.ob = []
        for n,ob in enumerate(self.parent.ob):
            if "active" in ob.keys() and "show" in ob.keys():
                if ob["active"] and ob["show"]:
                    line = ob["line"]
                    tmp_ob = self.p.parseObjects(line)
                    tmp_ob["index"] = ob["index"]

                    # if "min_moon_dist" in tmp_ob.keys():
                    #     print(tmp_ob["min_moon_dist"])

                    self.p.ob.append(tmp_ob)
        self.p.MakeTime()

        self.log_e.clear()
        self.log_e.setText(self.p.msg)



    def calc_vis(self):
        self.p.CalcObject()
        self.parent.update_table()
        # for n,ob in enumerate(self.p.ob):
        #     if "visibility" in ob.keys():
        #         print(ob["name"])
        self.log_e.clear()
        self.log_e.setText(self.p.msg)

    def mask_moon(self):
        self.p.MaskMoon()
        self.parent.update_table()
        self.log_e.clear()
        self.log_e.setText(self.p.msg)

    def mask_wind(self):
        self.p.MaskWind()
        self.parent.update_table()
        self.log_e.clear()
        self.log_e.setText(self.p.msg)

    def mask_cycle(self):
        self.p.MaskCycle()
        self.parent.update_table()
        self.log_e.clear()
        self.log_e.setText(self.p.msg)

    def mask_startend(self):
        self.p.MaskStartEnd()
        self.parent.update_table()
        self.log_e.clear()
        self.log_e.setText(self.p.msg)

    def mask_phstartend(self):
        self.p.MaskPhaseStartEnd()
        self.parent.update_table()
        self.log_e.clear()
        self.log_e.setText(self.p.msg)

    def mask_phase(self):
        self.p.MaskPhase()
        self.parent.update_table()
        self.log_e.clear()
        self.log_e.setText(self.p.msg)

    def run_tpg(self):
        self.load()
        self.calc_vis()
        self.mask_moon()
        self.mask_wind()
        self.mask_cycle()
        self.mask_startend()
        self.mask_phstartend()
        self.mask_phase()
        self.p.Waga()
        self.p.RandomizeList()
        self.p.allocate()
        self.log_e.clear()
        self.log_e.setText(self.p.msg)

        # wypisuje log
        # vis_list = ["h_min", "h_max", "min_moon_dist", "max_moon_phase", "wind","cycle", "t_start", "t_end",
        #             "ph_start","ph_end","ph_mk","all"]
        # txt = "name    "+"    ".join(vis_list) + "\n"
        # for n,ob in enumerate(self.p.ob):
        #     txt = txt + self.p.ob[n]["name"]
        #     for k in vis_list:
        #         if k in self.p.ob[n]["visibility"].keys():
        #             vis = numpy.array(self.p.ob[n]["visibility"][k])
        #             minutes = len(vis[vis])
        #             txt = txt+ "    "+str(minutes)
        #         else:
        #             txt = txt + "    --"
        #
        #     txt = txt + "\n"
        #
        # print(txt)

        self.p.export()  # export plan
        self.parent.update_table()


    def mkUI(self):
        grid = QGridLayout()

        self.log_e = QTextEdit()
        grid.addWidget(self.log_e, 0, 0,8,1)

        self.load_p = QPushButton('Load data / Init')
        self.load_p.clicked.connect(self.load)
        grid.addWidget(self.load_p, 0, 1)

        self.vis_p = QPushButton('Visibility')
        self.vis_p.clicked.connect(self.calc_vis)
        grid.addWidget(self.vis_p, 1, 1)

        self.moon_p = QPushButton('Moon')
        self.moon_p.clicked.connect(self.mask_moon)
        grid.addWidget(self.moon_p, 2, 1)

        self.wind_p = QPushButton('Wind')
        self.wind_p.clicked.connect(self.mask_wind)
        grid.addWidget(self.wind_p, 3, 1)

        self.cycle_p = QPushButton('Cycle')
        self.cycle_p.clicked.connect(self.mask_cycle)
        grid.addWidget(self.cycle_p, 4, 1)

        self.time_p = QPushButton('Time')
        self.time_p.clicked.connect(self.mask_startend)
        grid.addWidget(self.time_p, 5, 1)

        self.phlim_p = QPushButton('Phase limits')
        self.phlim_p.clicked.connect(self.mask_phstartend)
        grid.addWidget(self.phlim_p, 6, 1)

        self.phmk_p = QPushButton('Phase density')
        self.phmk_p.clicked.connect(self.mask_phase)
        grid.addWidget(self.phmk_p, 7, 1)

        self.tpg_p = QPushButton('run tpg')
        self.tpg_p.clicked.connect(self.run_tpg)
        grid.addWidget(self.tpg_p, 8, 1)

        self.close_p = QPushButton('Close')
        self.close_p.clicked.connect(lambda: self.close())
        grid.addWidget(self.close_p, 9, 0)

        self.setLayout(grid)
        self.show()



def sun_moon_ephem(obs_time, lat, lon, altitude, horizon=0*units.deg):

    loc = EarthLocation(lat=lat*units.deg, lon=lon*units.deg, height=altitude*units.m)
    t0 = Time(obs_time)
    jd = t0.jd

    delta_min = 1 * units.min
    times = t0 + delta_min * numpy.arange(-24*60, 24*60)

    # --- sun ---
    sun_alt = get_sun(times).transform_to(AltAz(obstime=times, location=loc)).alt - horizon
    sun_crossings = numpy.where(numpy.diff(numpy.sign(sun_alt.value)))[0]

    sunrise_times, sunset_times = [], []
    for c in sun_crossings:
        t_cross = times[c].to_datetime()
        if sun_alt[c] < 0 and sun_alt[c+1] > 0:
            sunrise_times.append(t_cross)
        else:
            sunset_times.append(t_cross)

    sunrise_times = sorted(sunrise_times)
    sunset_times = sorted(sunset_times)

    prev_sunrise = max([t for t in sunrise_times if t <= obs_time], default=None)
    next_sunrise = min([t for t in sunrise_times if t > obs_time], default=None)
    prev_sunset  = max([t for t in sunset_times if t <= obs_time], default=None)
    next_sunset  = min([t for t in sunset_times if t > obs_time], default=None)

    # --- moon ---
    moon_alt = get_moon(times).transform_to(AltAz(obstime=times, location=loc)).alt
    moon_crossings = numpy.where(numpy.diff(numpy.sign(moon_alt.value)))[0]

    moonrise_times, moonset_times = [], []
    for c in moon_crossings:
        t_cross = times[c].to_datetime()
        if moon_alt[c] < 0 and moon_alt[c+1] > 0:
            moonrise_times.append(t_cross)
        else:
            moonset_times.append(t_cross)

    moonrise_times = sorted(moonrise_times)
    moonset_times = sorted(moonset_times)

    prev_moonrise = max([t for t in moonrise_times if t <= obs_time], default=None)
    next_moonrise = min([t for t in moonrise_times if t > obs_time], default=None)
    prev_moonset  = max([t for t in moonset_times if t <= obs_time], default=None)
    next_moonset  = min([t for t in moonset_times if t > obs_time], default=None)

    # --- moon phase z ephem (astropy nie ma) ---

    obs = ephem.Observer()
    obs.lat = lat
    obs.lon = lon
    obs.date = str(obs_time)

    moon_ph = ephem.Moon(obs).phase


    return {
        "julian_date": jd,
        "prev_sunrise": prev_sunrise,
        "next_sunrise": next_sunrise,
        "prev_sunset": prev_sunset,
        "next_sunset": next_sunset,
        "prev_moonrise": prev_moonrise,
        "next_moonrise": next_moonrise,
        "prev_moonset": prev_moonset,
        "next_moonset": next_moonset,
        "moon_phase": moon_ph
    }
