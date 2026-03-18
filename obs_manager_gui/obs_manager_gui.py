#!/usr/bin/env python3
import copy
import math
import os,sys
import datetime
import uuid
from operator import index

import yaml
import requests

import warnings
import numpy
import matplotlib
import matplotlib.ticker as mticker
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

import ephem

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

from PyQt6.QtWidgets import QApplication, QMainWindow, QTableWidget, QAbstractItemView, QTableWidgetItem, QVBoxLayout, QWidget, QDialog, QGridLayout, QPushButton, QComboBox, QLabel, QLineEdit, QTextEdit, QCheckBox, QDateEdit, QTimeEdit, QDateTimeEdit, QFileDialog, QListWidget, QListWidgetItem, QVBoxLayout, QHBoxLayout
from PyQt6.QtCore import Qt, QTime, QDate, QDateTime
from PyQt6.QtGui import QFont, QColor




from astropy import units
from astropy.utils.exceptions import AstropyWarning
from astropy.time import Time
from astropy.coordinates import EarthLocation, Angle, get_sun, get_moon
from astropy.coordinates import SkyCoord, AltAz
from astropy.table import Table

from pyaraucaria.obs_plan.obs_plan_parser import ObsPlanParser
from pyaraucaria.ob_validator import ObsValidator

from tpg.telescope_plan_generator import TelescopePlanGenerator as tpg
from .obs_manager_lib import *


warnings.simplefilter('ignore', category=AstropyWarning)

class OM_Gui(QWidget):
    def __init__(self, args, parent=None):
        super().__init__()

        self.inactive_statuses = ["deactivated","inactive"]

        self.cwd = os.getcwd()  # curent working directory
        self.pwd = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # app location

        if os.path.exists(self.pwd+'/config.yaml'):
            with open(self.pwd+'/config.yaml', 'r') as cfg_file:
                self.cfg = yaml.safe_load(cfg_file)
        else:
            print("File not found: config.yaml")
            sys.exit()

        #print(self.cfg)

        self.schema_columns = ObsValidator.load_schema("tpg_schema")["properties"].keys()
        self.columns = ["ok_ob","tpg_vis"] + self.cfg["columns"]

        self.tpg_window = None
        self.i = -1
        self.mkUI()
        self.tel = self.tel_s.currentText()
        self.update_almanac()

    def update_table(self):

        show_all = self.showAll_p.isChecked()

        txt1 = self.filter_name_e.text()
        txt2 = self.filter_other_e.text()
        txt3 = self.filter_pi_e.text()
        txt4 = self.filter_sci_e.text()
        txt5 = self.filter_tag_e.text()

        if len(self.master_data) == 0:
            return

        row_labels = []

        try:
            self.table.cellChanged.disconnect(self.data_edit)
            self.table.cellClicked.disconnect(self.pocisniecie_tabelki)
        except TypeError:
            pass

        font = QFont()
        font.setPointSize(10)  # Ustawienie mniejszej czcionki
        self.table.setFont(font)

        self.table.setColumnCount(len(self.columns))
        for n,col_name in enumerate(self.columns):
            self.table.setHorizontalHeaderItem(n,QTableWidgetItem(col_name))

        i = -1
        self.table.setRowCount(0)
        self.table.clearContents()
        for data in self.master_data:

            if not show_all and not data["show"]:
                data["index"] = -2
            else:

                # filtrowanie wyswietlania - lokalnie, bo to tylko czesc wyswietlania
                show_txt = True

                if not (txt2.lower() in data["line"].lower()):
                    show_txt = False

                if data["ob"]:
                    if not (txt1.lower() in data["ob"].get("name", "").lower()):
                        show_txt = False
                    if not (txt3.lower() in data["ob"].get("pi", "").lower()):
                        show_txt = False
                    if not (txt4.lower() in data["ob"].get("sciprog", "").lower()):
                        show_txt = False
                    if not (txt5.lower() in data["ob"].get("tag", "").lower()):
                        show_txt = False
                    data["ob"] = self.clean_empty(data["ob"])

                if show_txt:

                    i = i + 1
                    data["index"] = i
                    if self.table.rowCount() <= i:
                        self.table.insertRow(i)  # Dodanie nowego wiersza

                    if data["ob"]:
                        item = QTableWidgetItem("")
                        item.setBackground(QColor("white"))
                        item.setForeground(QColor("black"))

                        for j,key in enumerate(self.columns):

                            if key in data["ob"].keys():
                                item = QTableWidgetItem(str(data["ob"][key]))
                                item.setBackground(QColor("white"))

                                red_bkg = False
                                if "validator" in data.keys():
                                    if key in data["validator"]["result"]:
                                        if not data["validator"]["result"][key]:
                                            red_bkg = True

                                if red_bkg:
                                    item.setBackground(QColor("darkRed"))  # jasnoszare tło
                                else:
                                    item.setBackground(QColor("white"))

                            else:
                                if key == "ok_ob":
                                    if "validator" in data.keys():
                                        if len(data["edited"])>0:
                                            txt = "\u2699"
                                            item = QTableWidgetItem(txt)
                                            item.setForeground(QColor("blue"))
                                            font = QFont()
                                            font.setPointSize(20)
                                            item.setFont(font)
                                        elif data["validator"]["valid"]:
                                            txt = "\u2705"
                                            item = QTableWidgetItem(txt)
                                            item.setForeground(QColor("green"))
                                        else:
                                            txt = "\u274C"
                                            item = QTableWidgetItem(txt)
                                            item.setForeground(QColor("red"))

                                    item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled)
                                    item.setBackground(QColor(200, 200, 200))
                                elif key == "tpg_vis":

                                    vis = data["tpg"].get("visibility",None)
                                    if vis:
                                        vis = numpy.array(vis["all"])
                                        t_vis = len(vis[vis])
                                        h_vis = t_vis / 60
                                        m_vis = t_vis - int(h_vis) * 60
                                        txt = f'{int(h_vis)}h {m_vis}m'

                                        item = QTableWidgetItem(txt)
                                        item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled)
                                        item.setBackground(QColor("lightGreen"))
                                    else:
                                        item = QTableWidgetItem("")
                                        item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled)
                                        item.setBackground(QColor(200, 200, 200))
                                else:
                                    item = QTableWidgetItem("")
                                    item.setBackground(QColor("white"))
                                    item.setForeground(QColor("black"))

                            if key in data["edited"]:
                                item.setBackground(QColor(170, 220, 240))

                            if data["ob"].get("status",None):
                                if data["ob"]["status"] in self.inactive_statuses:
                                    item.setBackground(QColor(180, 180, 180))

                            self.table.setItem(i, j, item)

                    else:

                        item = QTableWidgetItem(data["line"])
                        if data["parser_error"]:
                            item.setForeground(QColor("red"))
                        else:
                            item.setForeground(QColor("gray"))

                        item.setBackground(QColor("white"))
                        self.table.setItem(i, 0, item)
                        self.table.setSpan(i, 0, 1, self.table.columnCount())

                        row_labels.append(f"{i}")
                else:
                    data["index"] = -2

        self.table.setVerticalHeaderLabels(row_labels)
        self.table.resizeColumnsToContents()
        max_column_width = 100  # Maksymalna szerokość kolumny
        for col in range(self.table.columnCount()):
            self.table.setColumnWidth(col, min(self.table.columnWidth(col), max_column_width))

        self.table.cellChanged.connect(self.data_edit)
        self.table.cellClicked.connect(self.pocisniecie_tabelki)


    def load_objects(self):

        self.master_data = []

        with open(self.master_file, 'r') as plik:
            for line in plik:
                tmp = self.parse_line(line)
                self.master_data.append(tmp)

    def parse_line(self, line):
        tmp = {"ob": None, "line": None, "show": None, "parser_error": None, "index": -2, "edited": [], "tpg": {}}

        tmp["line"] = line

        if len(line.strip()) > 0:
            if len(line.split()) > 0:
                if "#" not in line.split()[0]:
                    txt = f'OBJECT {line}'
                    ob_tmp = ObsPlanParser.convert_from_string(txt)

                    if ob_tmp is None:
                        tmp["show"] = True
                        tmp["parser_error"] = True
                        print("*** parser error ***")
                    else:
                        ob = ObsValidator.convert_to_obdict(ob_tmp)
                        if ob:
                            tmp["ob"] = ob
                            tmp["show"] = True
                        else:
                            tmp["show"] = False
        return tmp

    def plot_sky_map(self):
        #try:
        self.sky_window = SkyWindow(self)
        self.sky_window.show()
        self.sky_window.raise_()
        #except:
        #    pass


    def plot_data(self):
        i = int(self.table.currentRow())
        i_tab = [int(data["index"]) for data in self.master_data]

        #try:
        n = i_tab.index(i)
        target = self.master_data[n]["ob"]["name"]

        self.phase_window = PhaseWindow(self,target,self.cfg["tel"][self.tel]["data_file"],self.master_data[n])
        self.phase_window.show()
        self.phase_window.raise_()

        #except ValueError:
        #    pass

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

        txt = ""
        txt = txt + f'sunset: {self.almanac["next_sunset"]}\n'
        txt = txt + f'sunrise: {self.almanac["next_sunrise"]}\n'
        txt = txt + f'moon: {self.almanac["moon_phase"]}\n'
        txt = txt + f'moonrise: {self.almanac["next_moonrise"]}\n'
        txt = txt + f'moonset: {self.almanac["next_moonset"]}\n'

        self.almanac_e.setText(txt)

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
        indx = [x["index"] for x in self.master_data]
        key = self.columns[j_selected]
        txt = self.table.item(i_selected, j_selected).text()

        if self.all_c.isChecked():
            for n in range(self.table.rowCount()):
                i = indx.index(n)
                if self.master_data[i].get("ob",None):
                    self.master_data[i]["ob"][key] = txt
                    self.master_data[i].setdefault("edited", []).append(key)

            self.all_c.setChecked(False)
        else:
            i = indx.index(i_selected)
            if self.master_data[i].get("ob",None):
                self.master_data[i]["ob"][key] = txt
                self.master_data[i].setdefault("edited", []).append(key)
            else:
                self.master_data[i] = self.parse_line(txt)

        self.update_table()

    def fill_uobi(self):
        indx = [x["index"] for x in self.master_data]
        if self.all_c.isChecked():
            for n in range(self.table.rowCount()):
                i = indx.index(n)
                if self.master_data[i].get("ob",None):
                    if self.master_data[i]["ob"].get("uobi",None):
                        pass
                    else:
                        self.master_data[i]["ob"]["uobi"] = str(uuid.uuid4())[:8]
                        self.master_data[i].setdefault("edited", []).append("uobi")
            self.all_c.setChecked(False)
        else:
            i = indx.index(self.i)
            if self.master_data[i].get("ob",None):
                self.master_data[i]["ob"]["uobi"] = str(uuid.uuid4())[:8]
                self.master_data[i].setdefault("edited", []).append("uobi")

        self.update_table()
        self.all_c.setChecked(False)

    def copy_ob(self):
        indx = [x["index"] for x in self.master_data]
        i = indx.index(self.i)
        tmp = copy.deepcopy(self.master_data[i])
        self.master_data.insert(i + 1, tmp)
        if self.master_data[i+1]["ob"].get("uobi",None):
            if len(self.master_data[i+1]["ob"]["uobi"])>2:
                self.master_data[i+1]["ob"]["uobi"] = str(uuid.uuid4())[:8]
        self.update_table()

    def delete_line(self):
        indx = [x["index"] for x in self.master_data]
        if self.i > -1:
            i = indx.index(self.i)
            del self.master_data[i]
        self.update_table()

    def validate_ob(self):
        BASE_SCHEMA = ObsValidator.load_schema("tpg_schema.yaml")
        COMMAND_RULES = ObsValidator.load_schema("command_rules.yaml")

        for i,data in enumerate(self.master_data):
            data["edited"] = []
            if data["ob"]:
                ob = data["ob"]
                validator = ObsValidator(BASE_SCHEMA, COMMAND_RULES)
                result = validator.validate_ob(ob)

                if "validator" not in data:
                    data["validator"] = {}

                data["validator"]["valid"] = result["valid"]
                data["validator"]["result"] = result["result"]
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


                    #print(txt)
                    file.write(txt)
                    print(f'objects saved to {file_path}')
            except Exception as e:
                print(f"Error saving file: {e}")

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

    def open_config(self):
        dialog = ColumnConfigDialog(self.schema_columns, self.columns)

        if dialog.exec():
            new_columns = dialog.get_columns()
            self.columns = new_columns
            self.update_table()

    def clean_empty(self, obs: dict) -> dict:
        return {k: v for k, v in obs.items() if v not in (None, "")}

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
        self.filter_name_e.textChanged.connect(self.update_table)
        grid.addWidget(self.filter_name_l, w, 0)
        grid.addWidget(self.filter_name_e, w, 1)

        self.filter_sci_l = QLabel("Filter SCIPROG")
        self.filter_sci_e = QLineEdit("")
        self.filter_sci_e.textChanged.connect(self.update_table)
        grid.addWidget(self.filter_sci_l, w, 2)
        grid.addWidget(self.filter_sci_e, w, 3)

        w = w + 1
        self.filter_pi_l = QLabel("Filter PI")
        self.filter_pi_e = QLineEdit("")
        self.filter_pi_e.textChanged.connect(self.update_table)
        grid.addWidget(self.filter_pi_l, w, 0)
        grid.addWidget(self.filter_pi_e, w, 1)

        self.filter_tag_l = QLabel("Filter TAG")
        self.filter_tag_e = QLineEdit("")
        self.filter_tag_e.textChanged.connect(self.update_table)
        grid.addWidget(self.filter_tag_l, w, 2)
        grid.addWidget(self.filter_tag_e, w, 3)

        w = w + 1
        self.filter_other_l = QLabel("Filter TXT")
        self.filter_other_e = QLineEdit("")
        self.filter_other_e.textChanged.connect(self.update_table)
        grid.addWidget(self.filter_other_l, w, 0)
        grid.addWidget(self.filter_other_e, w, 1)

        self.fill_uobi_p = QPushButton("Fill UOBI")
        self.fill_uobi_p.clicked.connect(self.fill_uobi)
        grid.addWidget(self.fill_uobi_p, w, 2)

        self.all_c = QCheckBox("Edit Column")
        self.all_c.setChecked(False)
        grid.addWidget(self.all_c, w, 3)

        self.showAll_p = QCheckBox("Show All")
        self.showAll_p.setChecked(True)
        grid.addWidget(self.showAll_p, w, 4)
        self.showAll_p.stateChanged.connect(self.update_table)

        w = w + 1
        self.table = QTableWidget()
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)  # zaznaczenie całego wiersza
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        #self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.setStyleSheet("selection-background-color: rgb(217,239,217); selection-color: black; ")

        grid.addWidget(self.table, w, 0, 1, 7)

        w = w + 1
        self.sky_p = QPushButton("Plot SkyMap")
        self.sky_p.clicked.connect(self.plot_sky_map)
        grid.addWidget(self.sky_p, w, 3)

        self.deleteOB_p = QPushButton("Delete OB")
        self.deleteOB_p.clicked.connect(self.delete_line)
        grid.addWidget(self.deleteOB_p, w, 0)

        self.validate_p = QPushButton("Validate OB")
        self.validate_p.clicked.connect(self.validate_ob)
        grid.addWidget(self.validate_p, w, 1)

        w = w + 1

        self.tpg_p = QPushButton("TPG")
        self.tpg_p.clicked.connect(self.tpg_show)
        grid.addWidget(self.tpg_p, w, 1)

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
        self.config_p.clicked.connect(self.open_config)
        self.close_p = QPushButton("Close")
        self.close_p.clicked.connect(self.close)
        grid.addWidget(self.load_p, w, 0)
        grid.addWidget(self.config_p, w, 2)
        grid.addWidget(self.save_p, w, 3)

        w = w + 1
        grid.addWidget(self.close_p, w, 4, 1, 3)

        self.setLayout(grid)

        self.show()


# #################################
#           Phase Window
# #################################

class PhaseWindow(QWidget):
    def __init__(self, parent, target, data_dir, data):
        super(PhaseWindow, self).__init__()
        self.parent = parent
        self.target = target
        self.data_dir = data_dir
        self.ob = data["ob"]
        self.data = data

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

    def parse_time(self,t):
        if not t:
            return None

        if "T" in t:
            t = t.replace("T", " ")

        if "/" not in t:  # tylko HH:MM
            date = str(self.obs_time).split()[0]

            start_hour = int(str(self.obs_time).split()[1].split(":")[0])
            t_hour = int(t.split(":")[0])

            if start_hour > 16 and t_hour < 16:
                date = ephem.Date(date) + 1
                date = str(ephem.Date(date)).split()[0]

            t = date + " " + t

        tmp = str(ephem.Date(t)).replace("/", "-")
        return Time(tmp, format="iso").jd

    def refresh(self):
        self._prepare_time()
        self._update_ephemeris()

        self._plot_lightcurve()
        self._plot_visibility()
        self._plot_tpg()

        self.canvas.draw()
        self.fig.subplots_adjust(hspace=0.3)
        self.show()

    # ========================
    # TIME + EPHEMERIS
    # ========================

    def _prepare_time(self):
        self.obs_time = datetime.datetime.combine(
            self.parent.date_e.date().toPyDate(),
            self.parent.time_e.time().toPyTime()
        )
        t = Time(self.obs_time, scale='utc')

        self.current_jd = t.jd
        self.jd3h = self.current_jd + numpy.arange(1, 7) / 24.

    def _update_ephemeris(self):
        moon_phase = float(self.parent.almanac["moon_phase"])

        self.ephem_e.setText(f"moon phase: {moon_phase:.0f} %")
        self.ephem_e.setStyleSheet("background-color: white;")

        if self.ob.get("max_moon_phase") and float(self.ob["max_moon_phase"]) < moon_phase:
            self.ephem_e.setStyleSheet("background-color: lightcoral;")

    # ========================
    # LIGHT CURVE
    # ========================

    def _plot_lightcurve(self):
        self.axes.clear()

        try:
            jd, mag, flag = self._load_lightcurve()
            if len(jd) == 0:
                return

            jd = numpy.array(jd)
            mag = numpy.array(mag)
            flag = numpy.array(flag)

            self.now_t = self.current_jd

            if self.phase_c.isChecked():
                jd = self._convert_to_phase(jd)
            else:
                self.axes.set_title(f"{self.target}")

            self._plot_recent_and_all(jd, mag, flag)
            self._format_lightcurve_axes(mag)

            self.axes.axvline(self.now_t, color="blue")
            self._plot_time_markers()

        except (FileNotFoundError, ValueError) as e:
            print(f"Lightcurve error: {e}")

    def _load_lightcurve(self):
        filter_name = self.file_s.currentText()

        file = self.ob.get(
            "obs_data",
            f"{self.f_path}/{filter_name}/light-curve/{self.target.lower()}_{filter_name}_diff_light_curve.txt"
        )

        tab = Table.read(file, format="ascii")
        return tab["jd_obs"], tab["mag"], tab["quality"]

    def _convert_to_phase(self, jd):
        if "P" not in self.ob:
            self.phase_c.setChecked(False)
            return jd

        P = float(self.ob["P"])
        jd0 = float(self.ob.get("hjd0", 2460000))

        self.now_t = (self.current_jd - jd0) / P % 1
        self.jd3h = (self.jd3h - jd0) / P % 1

        jd = (jd - jd0) / P % 1

        self.axes.set_xlim(-0.1, 1.1)
        self.axes.set_title(f"{self.target} P={P}")

        return jd

    def _plot_recent_and_all(self, jd, mag, flag):
        recent_mask = jd > self.current_jd - float(self.parent.cfg["last_nights_to_mark"])

        for q, color in [(0, "g"), (1, "c"), (2, "k")]:
            self.axes.plot(jd[flag == q], mag[flag == q], f".{color}", alpha=0.1)
            self.axes.plot(jd[recent_mask & (flag == q)],
                           mag[recent_mask & (flag == q)],
                           f".{color}", alpha=1)

    def _format_lightcurve_axes(self, mag):
        d = 0.1 * (max(mag) - min(mag))
        self.axes.set_ylim(max(mag) + d, min(mag) - d)

    def _plot_time_markers(self):
        alpha = 1.0
        for i, x in enumerate(self.jd3h):
            alpha /= (i + 2)
            self.axes.axvline(x, color="blue", alpha=alpha)

    # ========================
    # VISIBILITY
    # ========================

    def _plot_visibility(self):
        self.axes2.clear()

        hmin = float(self.ob.get("h_min", self.parent.cfg["tel"][self.parent.tel]["hmin"]))
        hmax = float(self.ob.get("h_max", self.parent.cfg["tel"][self.parent.tel]["hmax"]))

        t = numpy.linspace(int(self.current_jd), int(self.current_jd) + 1, 100)
        time_range = Time(t, format="jd")

        alt, sun_alt, moon_alt, sep = self._compute_altaz(time_range)

        self._format_visibility_axes(hmin, hmax)
        self._plot_visibility_lines(time_range.jd, alt, sun_alt, moon_alt)

        self._update_moon_sep(sep)
        self._plot_time_constraints(time_range)
        self._set_time_ticks(self.axes2, time_range)

    def _compute_altaz(self, time_range):
        loc = EarthLocation(
            lat=self.parent.cfg["obs_latitude"],
            lon=self.parent.cfg["obs_longitude"],
            height=self.parent.cfg["obs_elevation"]
        )

        i = int(self.parent.table.currentRow())
        i_tab = [int(data["index"]) for data in self.parent.master_data]
        n = i_tab.index(i)

        ob = self.parent.master_data[n]["ob"]

        coord = SkyCoord(ra=ob["ra"], dec=ob["dec"], unit=('hourangle', 'deg'))
        frame = AltAz(obstime=time_range, location=loc)

        alt = coord.transform_to(frame).alt.deg
        sun_alt = get_sun(time_range).transform_to(frame).alt.deg
        moon = get_moon(time_range)
        moon_alt = moon.transform_to(frame).alt.deg
        sep = coord.transform_to(frame).separation(moon).deg

        return alt, sun_alt, moon_alt, sep

    def _format_visibility_axes(self, hmin, hmax):
        self.axes2.set_ylim(-20, 90)
        self.axes2.set_xlim(int(self.current_jd), int(self.current_jd) + 1)

        self.axes2.axvline(self.current_jd, color="blue")
        self.axes2.axhspan(-20, 0, facecolor='red', alpha=0.1)
        self.axes2.axhspan(0, hmin, facecolor='black', alpha=0.05)
        self.axes2.axhspan(hmax, 90, facecolor='black', alpha=0.05)

    def _plot_visibility_lines(self, t, alt, sun_alt, moon_alt):
        self.axes2.plot(t, alt, "-g")
        self.axes2.plot(t, sun_alt, "--y")
        self.axes2.plot(t, moon_alt, ":k")

    def _update_moon_sep(self, sep):
        self.moon_sep_e.setText(f"separation: {min(sep):.0f}\u00B0")
        self.moon_sep_e.setStyleSheet("background-color: white;")

        if self.ob.get("min_moon_dist") and float(self.ob["min_moon_dist"]) > float(min(sep)):
            self.moon_sep_e.setStyleSheet("background-color: lightcoral;")

    def _plot_time_constraints(self, time_range):
        t_start = self.parse_time(self.ob.get("t_start"))
        t_end = self.parse_time(self.ob.get("t_end"))

        xmin = time_range.jd[0]
        xmax = time_range.jd[-1]

        if t_start:
            self.axes2.axvspan(xmin, t_start, color="gray", alpha=0.15)
        if t_end:
            self.axes2.axvspan(t_end, xmax, color="gray", alpha=0.15)

    # ========================
    # TPG
    # ========================

    def _plot_tpg(self):
        self.axes3.clear()

        tpg = self.data.get("tpg")
        if not tpg:
            return

        nt = tpg.get("nightTime")
        vis = tpg.get("visibility")

        if not nt or not vis:
            return

        nt = Time([ephem.Date(t).datetime() for t in nt], scale='utc').jd

        for i, (label, values) in enumerate(vis.items()):
            segments = self._split_segments(nt, values)

            for color, seg in segments.items():
                if seg:
                    self.axes3.broken_barh(seg, (i - 0.4, 0.8),
                                           facecolors=color, alpha=0.3)

        self.axes3.axvline(self.current_jd, color="blue")
        self.axes3.set_yticks(range(len(vis)))
        self.axes3.set_yticklabels(vis.keys())
        self.axes3.set_xlim(self.axes2.get_xlim())

        self._set_time_ticks(self.axes3, Time(nt, format="jd"))

    def _split_segments(self, nt, values):
        green, red = [], []

        for i in range(len(values) - 1):
            seg = (nt[i], nt[i + 1] - nt[i])
            (green if values[i] else red).append(seg)

        return {"green": green, "red": red}

    # Ticksy na wykresach

    def _format_jd_tick(self, jd, pos=None):
        dt = Time(jd, format="jd").to_datetime()
        return dt.strftime("%H:%M")

    # Funkcja ustawiająca inteligentne ticki
    def _set_time_ticks(self, ax, time_range):
        jd_start = time_range.jd[0]
        jd_end = time_range.jd[-1]

        # 1. Najważniejsze ticki
        ticks_priority = []

        # current time
        ticks_priority.append((self.current_jd, 0))  # 0 = najwyższy priorytet

        # wschód/zachód słońca
        if self.parent.almanac.get("next_sunset"):
            ticks_priority.append((Time(self.parent.almanac["next_sunset"]).jd, 1))
        if self.parent.almanac.get("next_sunrise"):
            ticks_priority.append((Time(self.parent.almanac["next_sunrise"]).jd, 1))

        # wschód/zachód księżyca
        if self.parent.almanac.get("next_moonrise"):
            ticks_priority.append((Time(self.parent.almanac["next_moonrise"]).jd, 2))
        if self.parent.almanac.get("next_moonset"):
            ticks_priority.append((Time(self.parent.almanac["next_moonset"]).jd, 2))

        # co 2 godziny
        ticks_2h = list(numpy.arange(jd_start, jd_end, 2 / 24.))
        for t in ticks_2h:
            ticks_priority.append((t, 3))

        # sortowanie po priorytecie
        ticks_priority.sort(key=lambda x: x[1])

        # filtr odległości minimalnej między tickami (np. 30 min = 0.0208 JD)
        min_distance = 60 / 60 / 24  # 30 minut w JD
        final_ticks = []

        for jd_val, prio in ticks_priority:
            if all(abs(jd_val - t) > min_distance for t in final_ticks):
                final_ticks.append(jd_val)

        final_ticks.sort()  # dla estetyki od lewej do prawej
        ax.set_xticks(final_ticks)
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(self._format_jd_tick))

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

        gs = self.fig.add_gridspec(3, 1, height_ratios=[2, 2, 1])

        self.axes = self.fig.add_subplot(gs[0])
        self.axes2 = self.fig.add_subplot(gs[1])
        self.axes3 = self.fig.add_subplot(gs[2])

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

        ra = []
        dec = []
        self.indx = []

        tabel_elements = range(self.parent.table.rowCount())
        ob_ind = [d["index"] for d in self.parent.master_data]
        for te in tabel_elements:
            ti = ob_ind.index(te)
            data = self.parent.master_data[ti]
            if data.get("ob",None):
                r = data["ob"].get("ra",None)
                d = data["ob"].get("dec",None)
                if r and d:
                    i = data["index"]
                    ra.append(r)
                    dec.append(d)
                    self.indx.append(i)

        coords = SkyCoord(ra=ra, dec=dec, unit=('hourangle', 'deg'), frame='icrs')
        altaz = coords.transform_to(AltAz(obstime=Time(obs_time), location=obs_location))

        azimuth = altaz.az.deg
        self.altitude = altaz.alt.deg
        self.azimuth_rad = numpy.radians(azimuth)
        self.axes2.scatter(self.azimuth_rad, self.altitude,marker="*", c='green')

        if self.parent.i in self.indx:
            n = self.indx.index(self.parent.i)
        else:
            n = -2

        if n > -1:
            self.axes2.plot(self.azimuth_rad[n],self.altitude[n],"*r")

        self.axes2.set_rlim([90, 0])


        lst = Time(obs_time,scale='utc', location=obs_location).sidereal_time('mean')
        ha = [15*(lst.hour - (float(x.split(":")[0])+float(x.split(":")[1])/60.+float(x.split(":")[2])/3600.)) for x in ra]

        self.dec_deg =  Angle(dec, unit="deg").deg
        self.ha_rad = numpy.radians(ha)
        self.axes.scatter(self.ha_rad, self.dec_deg,marker="*", c='green', s=50)
        if n > -1:
            self.axes.plot(self.ha_rad[n],self.dec_deg[n],"*r")

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
                        self.parent.i = self.indx[min_i]
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
                        self.parent.i = self.indx[min_i]
                        self.updateMap()
                        self.parent.update_selection()

# ######################
#        TPG
# ######################


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

        t_next = self.parent.almanac["next_sunset"]
        t_prev = self.parent.almanac["prev_sunset"]

        if t_prev <= t0 < t_next:
            # dzień → bierzemy najbliższy zachód (dzisiejszy)
            night_start = t_next
        else:
            # noc → bierzemy poprzedni zachód
            night_start = t_prev

        dt = [night_start.strftime("%Y/%m/%d")]

        self.p = tpg(tel, dt, loud=True)

        self.p.Initiate()
        self.p.LoadObjects()
        self.p.ob = []
        for n,data in enumerate(self.parent.master_data):
            if data.get("ob",None):
                line = ObsValidator.convert_from_obdict(data.get("ob"))
                line = line.split(" ", 1)[1]
                tmp_ob = self.p.parseObjects(line)
                tmp_ob["index"] = n
                self.p.ob.append(tmp_ob)
        self.p.MakeTime()

        self.log_e.clear()
        self.log_e.setText(self.p.msg)

    def calc_vis(self):
        self.p.CalcObject()
        for n,ob in enumerate(self.p.ob):
            if "visibility" in ob.keys():
                self.parent.master_data[ob["index"]]["tpg"]["visibility"] = ob["visibility"]
                self.parent.master_data[ob["index"]]["tpg"]["nightTime"] = self.p.nightTime
        self.parent.update_table()
        self.log_e.clear()
        self.log_e.setText(self.p.msg)

    def mask_moon(self):
        self.p.MaskMoon()
        for n,ob in enumerate(self.p.ob):
            if "visibility" in ob.keys():
                self.parent.master_data[ob["index"]]["tpg"]["visibility"] = ob["visibility"]
        self.parent.update_table()
        self.log_e.clear()
        self.log_e.setText(self.p.msg)
        self.log_e.setText("Moon masked")

    def mask_wind(self):
        self.p.MaskWind()
        for n,ob in enumerate(self.p.ob):
            if "visibility" in ob.keys():
                self.parent.master_data[ob["index"]]["tpg"]["visibility"] = ob["visibility"]
        self.parent.update_table()
        self.log_e.clear()
        self.log_e.setText(self.p.msg)
        self.log_e.setText("Wind masked")

    def mask_cycle(self):
        self.p.MaskCycle()
        for n,ob in enumerate(self.p.ob):
            if "visibility" in ob.keys():
                self.parent.master_data[ob["index"]]["tpg"]["visibility"] = ob["visibility"]
        self.parent.update_table()
        self.log_e.clear()
        self.log_e.setText(self.p.msg)
        self.log_e.setText("Cycle masked")

    def mask_startend(self):
        self.p.MaskStartEnd()
        for n,ob in enumerate(self.p.ob):
            if "visibility" in ob.keys():
                self.parent.master_data[ob["index"]]["tpg"]["visibility"] = ob["visibility"]
        self.parent.update_table()
        self.log_e.clear()
        self.log_e.setText(self.p.msg)
        self.log_e.setText("Time masked")

    def mask_phstartend(self):
        self.p.MaskPhaseStartEnd()
        for n,ob in enumerate(self.p.ob):
            if "visibility" in ob.keys():
                self.parent.master_data[ob["index"]]["tpg"]["visibility"] = ob["visibility"]
                print(self.parent.master_data[ob["index"]]["tpg"]["visibility"].keys())

        self.parent.update_table()
        self.log_e.clear()
        self.log_e.setText(self.p.msg)
        self.log_e.setText("Phase masked")

    def mask_phase(self):
        self.p.MaskPhase()
        for n,ob in enumerate(self.p.ob):
            if "visibility" in ob.keys():
                self.parent.master_data[ob["index"]]["tpg"]["visibility"] = ob["visibility"]
                print(self.parent.master_data[ob["index"]]["tpg"]["visibility"].keys())

        self.parent.update_table()
        self.log_e.clear()
        self.log_e.setText(self.p.msg)
        self.log_e.setText("Phase density masked")

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





class ColumnConfigDialog(QDialog):

    def __init__(self, all_columns, visible_columns=None, parent=None):
        super().__init__(parent)

        self.setWindowTitle("Configure Columns")

        if visible_columns is None:
            visible_columns = all_columns

        layout = QVBoxLayout(self)

        self.list = QListWidget()
        self.list.setDragDropMode(QListWidget.DragDropMode.InternalMove)

        layout.addWidget(self.list)

        # najpierw widoczne kolumny w zapisanej kolejności
        for name in visible_columns:
            item = QListWidgetItem(name)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Checked)
            self.list.addItem(item)

        # potem niewidoczne
        for name in all_columns:
            if name not in visible_columns:
                item = QListWidgetItem(name)
                item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                item.setCheckState(Qt.CheckState.Unchecked)
                self.list.addItem(item)

        buttons = QHBoxLayout()

        ok = QPushButton("OK")
        cancel = QPushButton("Cancel")

        ok.clicked.connect(self.accept)
        cancel.clicked.connect(self.reject)

        buttons.addWidget(ok)
        buttons.addWidget(cancel)

        layout.addLayout(buttons)

    def get_columns(self):

        result = []

        for i in range(self.list.count()):
            item = self.list.item(i)

            if item.checkState() == Qt.CheckState.Checked:
                result.append(item.text())

        return result
