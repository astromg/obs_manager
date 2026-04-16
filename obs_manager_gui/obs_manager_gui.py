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

from PyQt6.QtWidgets import QApplication, QMainWindow, QTableWidget, QAbstractItemView, QTableWidgetItem, QVBoxLayout, \
    QWidget, QDialog, QGridLayout, QPushButton, QComboBox, QLabel, QLineEdit, QTextEdit, QCheckBox, QDateEdit, \
    QTimeEdit, QDateTimeEdit, QFileDialog, QListWidget, QListWidgetItem, QVBoxLayout, QHBoxLayout, QFrame
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
from .plan_gui import Plan_Gui

warnings.simplefilter('ignore', category=AstropyWarning)

class OM_Gui(QWidget):
    def __init__(self, args, parent=None):
        super().__init__()
        self.plan_gui = None
        self.phase_window = None
        self.sky_window = None

        self.inactive_statuses = ["deactivated","inactive"]

        self.cwd = os.getcwd()  # curent working directory
        self.pwd = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # app location

        if os.path.exists(self.pwd+'/config.yaml'):
            with open(self.pwd+'/config.yaml', 'r') as cfg_file:
                self.cfg = yaml.safe_load(cfg_file)
        else:
            print("File not found: config.yaml")
            sys.exit()

        t = tpg("test",["2022/10/10"])
        self.tpg_cfg = t.cfg

        self.schema_columns = ObsValidator.load_schema("tpg_schema")["properties"].keys()
        self.extra_columns = ["ok_ob","alt","last_obs","tpg_vis"]
        self.columns = self.extra_columns + self.cfg["columns"]

        self.tpg_window = None
        self.i = -1
        self.mkUI()
        self.tel = self.tel_s.currentText()
        self.update_almanac()

        self.plan_gui = Plan_Gui(self)
        self.plan_gui.show()

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
                                            txt = "\u2714"
                                            item = QTableWidgetItem(txt)
                                            item.setForeground(QColor("green"))
                                        else:
                                            txt = "\u274C"
                                            item = QTableWidgetItem(txt)
                                            item.setForeground(QColor("red"))

                                    item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled)
                                    if j % 2 == 0:
                                        color = QColor(210, 210, 210)  # jasny szary
                                    else:
                                        color = QColor(220, 220, 220)
                                    item.setBackground(color)
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
                                        if j % 2 == 0:
                                            color = QColor(210, 210, 210)  # jasny szary
                                        else:
                                            color = QColor(220, 220, 220)
                                    item.setBackground(color)
                                elif key == "ctc":

                                    #ctc = data["tpg"].get("ob_time",None)

                                    item = QTableWidgetItem("")
                                    item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled)
                                    if j % 2 == 0:
                                        color = QColor(210, 210, 210)  # jasny szary
                                    else:
                                        color = QColor(220, 220, 220)
                                    item.setBackground(color)

                                elif key == "last_obs":

                                    last_jd = data["ob"].get("last_jd",None)
                                    if last_jd:
                                        dt = self.almanac["julian_date"] - last_jd
                                        item = QTableWidgetItem(f'{dt:.1f}')
                                    else:
                                        item = QTableWidgetItem("")

                                    item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled)
                                    if j % 2 == 0:
                                        color = QColor(210, 210, 210)  # jasny szary
                                    else:
                                        color = QColor(220, 220, 220)
                                    item.setBackground(color)


                                elif key == "alt":

                                    alt = data.get("alt_now",None)

                                    if alt:
                                        item = QTableWidgetItem(f'{alt:.0f}')
                                    else:
                                        item = QTableWidgetItem("")

                                    item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled)
                                    if j % 2 == 0:
                                        color = QColor(210, 210, 210)  # jasny szary
                                    else:
                                        color = QColor(220, 220, 220)
                                    item.setBackground(color)

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
                    # wywalic jak juz bedzie wszystko konsystentnie
                    if "OBJECT" in line.split()[0]:
                        txt = line
                    else:
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

        self.phase_window = PhaseWindow(self,target,self.tpg_cfg[self.tel]["data_file"],self.master_data[n])
        #self.phase_window = PhaseWindow(self, target, self.cfg["tel"][self.tel]["data_file"], self.master_data[n])
        self.phase_window.show()
        self.phase_window.raise_()

        #except ValueError:
        #    pass

    def time_changed(self):
        self.update_almanac()
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

        if self.plan_gui:
            self.plan_gui.update_table()



    def date_changed(self):
        self.update_almanac()
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

        if self.plan_gui:
            self.plan_gui.update_table()

    def update_almanac(self):
        self.obs_time = datetime.datetime.combine(self.date_e.date().toPyDate(), self.time_e.time().toPyTime())
        time = Time(self.obs_time, scale='utc')
        self.almanac = sun_moon_ephem(time, self.tpg_cfg["obs_lat"], self.tpg_cfg["obs_lon"], self.tpg_cfg["obs_elev"], horizon=0)
        self.almanac["julian_date"] = time.jd

        txt = f"""

        <table cellspacing="4">
        <tr><td><b>Julian date:</b></td><td>{time.jd:.5f}</td> </tr>
        <tr><td><b>Sunset:</b></td><td>{format_dt(self.almanac["next_sunset"])}</td></tr>
        <tr><td><b>Sunrise:</b></td><td>{format_dt(self.almanac["next_sunrise"])}</td></tr>
        <tr><td><b>Moon phase:</b></td><td>{self.almanac["moon_phase"]:.1f} %</td></tr>
        </table>


        """

        # <tr><td><b>Moonrise:</b></td><td>{format_dt(self.almanac["next_moonrise"])}</td></tr>
        # <tr><td><b>Moonset:</b></td><td>{format_dt(self.almanac["next_moonset"])}</td></tr>

        self.almanac_e.setHtml(txt)

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

    def last_obs(self):
        indx = [x["index"] for x in self.master_data]

        for n in range(self.table.rowCount()):
            i = indx.index(n)
            ob = self.master_data[i].get("ob",None)
            if ob:
                try:
                    if ob.get("obs_data",None):
                        fname = self.master_data[i]["ob"]["obs_data"]
                        file = self.tpg_cfg[self.tel]["data_file"]+"/"+fname
                    else:
                        filtr = ob["seq"].split("/")[1]
                        file = (self.tpg_cfg[self.tel]["data_file"] + ob["name"].lower() + "/" + filtr + "/light-curve/" + ob["name"].lower() + "_" + filtr + "_diff_light_curve.txt")
                    print(file)
                    lc_tab = Table.read(file, format="ascii")
                    jd = numpy.array(lc_tab["jd_obs"])
                    if len(jd) == 0:
                        continue
                    last_jd = max(jd)
                    ob["last_jd"] = last_jd
                except (FileNotFoundError, ValueError, KeyError):
                    pass
        self.update_table()



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

    def calc_visibility(self):
        oca = ephem.Observer()
        oca.lon = str(self.tpg_cfg["obs_lon"])
        oca.lat = str(self.tpg_cfg["obs_lat"])
        oca.elev = float(self.tpg_cfg["obs_elev"])
        oca.horizon = "0"

        obs_time = datetime.datetime.combine(self.date_e.date().toPyDate(), self.time_e.time().toPyTime())
        oca.date = str(Time(obs_time, scale='utc'))

        for i,data in enumerate(self.master_data):
            if data["ob"]:
                ob = data["ob"]
                if ob.get("ra", None) and ob.get("dec", None):
                    s = ephem.FixedBody()
                    s._ra = ob["ra"]
                    s._dec = ob["dec"]
                    s.compute(oca)
                    alt = float(s.alt) * 180.0 / ephem.pi
                    az = float(s.az) * 180.0 / ephem.pi
                    self.master_data[i]["alt_now"] = alt

        self.update_table()

    def validate_ob(self):
        BASE_SCHEMA = ObsValidator.load_schema("base_schema.yaml")
        TPG_SCHEMA = ObsValidator.load_schema("tpg_schema.yaml")

        SCHEMA = merge_schemas(BASE_SCHEMA, TPG_SCHEMA)

        COMMAND_RULES = ObsValidator.load_schema("command_rules.yaml")

        for i,data in enumerate(self.master_data):
            data["edited"] = []
            if data["ob"]:
                ob = data["ob"]
                validator = ObsValidator(SCHEMA, COMMAND_RULES)
                result = validator.validate_ob(ob)

                if "validator" not in data:
                    data["validator"] = {}

                data["validator"]["valid"] = result["valid"]
                data["validator"]["result"] = result["result"]
        self.update_table()

    def add_to_plan(self):
        indx = [x["index"] for x in self.master_data]
        if self.i > -1:
            i = indx.index(self.i)
            if self.master_data[i]["ob"]:
                ob = self.master_data[i]["ob"]
                self.plan_gui.add(ob)

    def save_file(self):
        file_path, _ = QFileDialog.getSaveFileName(self, "Save File", self.tpg_cfg[self.tel]["master_file"],"Text Files (*.txt);;All Files (*)")
        if file_path:
            try:
                with open(file_path, "w", encoding="utf-8") as file:
                    txt = ""
                    for data in self.master_data:
                        if data["ob"]:
                            txt = txt + ObsValidator.convert_from_obdict(data.get("ob"))+"\n"
                        else:
                            txt = txt + data["line"]
                    file.write(txt)
                    print(f'objects saved to {file_path}')
            except Exception as e:
                print(f"Error saving file: {e}")

    def tpg_show(self):
        self.update_almanac()
        self.tpg_window = TPGWindow(self)



    def load_file(self):
        file_path, _ = QFileDialog.getOpenFileName(None,"Select a File",self.tpg_cfg[self.tel]["master_file"],"All Files (*);;Text Files (*.txt);;Images (*.png *.jpg)")
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


            self.cfg["columns"] = [c for c in self.columns if c not in self.extra_columns]
            if os.path.exists(self.pwd+'/config.yaml'):
                with open(self.pwd+'/config.yaml', 'w') as cfg_file:
                    yaml.safe_dump(self.cfg, cfg_file)

    def clean_empty(self, obs: dict) -> dict:
        return {k: v for k, v in obs.items() if v not in (None, "")}

    def mkUI(self):
        self.setWindowTitle('OCM observing plan manager')
        #self.setGeometry(50, 50, 1400, 800)
        self.resize(1400, 800)

        grid = QGridLayout()

        self.almanac_e = QTextEdit()
        self.almanac_e.setReadOnly(True)
        self.almanac_e.setStyleSheet("background-color: rgb(235,235,235);")
        grid.addWidget(self.almanac_e, 0, 6,4,2)

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
        grid.addWidget(self.date_l, w, 3)
        grid.addWidget(self.date_e, w, 4)
        grid.addWidget(self.time_e, w, 5)

        w = w + 1
        self.filter_name_l = QLabel("Filter NAME")
        self.filter_name_e = QLineEdit("")
        self.filter_name_e.textChanged.connect(self.update_table)
        grid.addWidget(self.filter_name_l, w, 0)
        grid.addWidget(self.filter_name_e, w, 1, 1,2)

        self.filter_sci_l = QLabel("Filter SCIPROG")
        self.filter_sci_e = QLineEdit("")
        self.filter_sci_e.textChanged.connect(self.update_table)
        grid.addWidget(self.filter_sci_l, w, 3)
        grid.addWidget(self.filter_sci_e, w, 4,1,2)

        w = w + 1
        self.filter_pi_l = QLabel("Filter PI")
        self.filter_pi_e = QLineEdit("")
        self.filter_pi_e.textChanged.connect(self.update_table)
        grid.addWidget(self.filter_pi_l, w, 0)
        grid.addWidget(self.filter_pi_e, w, 1,1,2)

        self.filter_tag_l = QLabel("Filter TAG")
        self.filter_tag_e = QLineEdit("")
        self.filter_tag_e.textChanged.connect(self.update_table)
        grid.addWidget(self.filter_tag_l, w, 3)
        grid.addWidget(self.filter_tag_e, w, 4,1,2)

        w = w + 1
        self.filter_other_l = QLabel("Filter TXT")
        self.filter_other_e = QLineEdit("")
        self.filter_other_e.textChanged.connect(self.update_table)
        grid.addWidget(self.filter_other_l, w, 0)
        grid.addWidget(self.filter_other_e, w, 1,1,2)

        self.fill_uobi_p = QPushButton("Fill UOBI")
        self.fill_uobi_p.clicked.connect(self.fill_uobi)
        grid.addWidget(self.fill_uobi_p, w, 3)

        self.all_c = QCheckBox("Edit Column")
        self.all_c.setChecked(False)
        grid.addWidget(self.all_c, w, 4)

        self.showAll_p = QCheckBox("Show All")
        self.showAll_p.setChecked(True)
        grid.addWidget(self.showAll_p, w, 5)
        self.showAll_p.stateChanged.connect(self.update_table)

        w = w + 1
        self.table = QTableWidget()
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)  # zaznaczenie całego wiersza
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        #self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.setStyleSheet("selection-background-color: rgb(217,239,217); selection-color: black; ")

        grid.addWidget(self.table, w, 0, 1, 8)

        w = w + 1

        self.deleteOB_p = QPushButton("Delete OB")
        self.deleteOB_p.clicked.connect(self.delete_line)
        grid.addWidget(self.deleteOB_p, w, 0)

        self.validate_p = QPushButton("Validate OB")
        self.validate_p.clicked.connect(self.validate_ob)
        grid.addWidget(self.validate_p, w, 2)

        self.sky_p = QPushButton("Plot SkyMap")
        self.sky_p.clicked.connect(self.plot_sky_map)
        grid.addWidget(self.sky_p, w, 4,1,2)

        self.add_p = QPushButton("Add to Plan")
        self.add_p.clicked.connect(self.add_to_plan)
        grid.addWidget(self.add_p, w, 6 ,1 ,2)

        w = w + 1

        self.visibility_p = QPushButton("Visibility")
        self.visibility_p.clicked.connect(self.calc_visibility)
        grid.addWidget(self.visibility_p, w, 2)

        self.data_p = QPushButton("Plot data")
        self.data_p.clicked.connect(self.plot_data)
        grid.addWidget(self.data_p, w, 4,1,2)

        self.tpg_p = QPushButton("TPG")
        self.tpg_p.clicked.connect(self.tpg_show)
        grid.addWidget(self.tpg_p, w, 6 ,1 ,2)


        w = w + 1

        self.copy_p = QPushButton("Copy")
        self.copy_p.clicked.connect(self.copy_ob)
        grid.addWidget(self.copy_p, w, 0)

        self.last_p = QPushButton("Last obs")
        self.last_p.clicked.connect(self.last_obs)
        grid.addWidget(self.last_p, w, 2)

        self.save_p = QPushButton("Save")
        self.save_p.clicked.connect(self.save_file)
        grid.addWidget(self.save_p, w, 6,1,2)

        w = w + 1

        self.line_l = QFrame()
        self.line_l.setFrameShape(QFrame.Shape.HLine)
        self.line_l.setFrameShadow(QFrame.Shadow.Raised)
        grid.addWidget(self.line_l, w, 0, 1, 8)

        w = w + 1

        self.load_p = QPushButton("Load file")
        self.load_p.clicked.connect(self.load_file)
        grid.addWidget(self.load_p, w, 0, 1,2)

        self.config_p = QPushButton("\u2699")
        self.config_p.clicked.connect(self.open_config)
        grid.addWidget(self.config_p, w, 2)

        self.close_p = QPushButton("Close")
        self.close_p.clicked.connect(QApplication.quit)
        grid.addWidget(self.close_p, w, 6, 1, 2)

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
        self.zaznaczenie_i = -1
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

        cid = self.canvas.mpl_connect('key_press_event', self.zaznaczenie)


    def zaznaczenie(self,event):
        if event.key == "f":
            if event.xdata != None:
                x = float(event.xdata)
                y = float(event.ydata)

                dx = self.jd - x
                dy = self.mag - y
                r = dx**2+dy**2
                i = numpy.argmin(r)

                self.zaznaczenie_i = i
                print(self.fits_file[i])
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

        self.fig.subplots_adjust(hspace=0.35)
        self.canvas.draw()

        self.canvas.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.canvas.setFocus()

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
            fits_file, jd, mag, flag = self._load_lightcurve()
            if len(jd) == 0:
                return

            self.jd = numpy.array(jd)
            self.mag = numpy.array(mag)
            self.flag = numpy.array(flag)
            self.fits_file = fits_file

            self.now_t = self.current_jd

            if self.phase_c.isChecked():
                self.jd = self._convert_to_phase(self.jd)
                self._plot_phase_constraints(self.jd, self.flag)
            else:
                self.axes.set_title(f"{self.target}")

            self._plot_recent_and_all(self.jd, self.mag, self.flag)
            self._plot_zaznaczenie(self.jd, self.mag, self.flag)
            self._format_lightcurve_axes(self.mag)

            # thicker current time
            self.axes.axvline(self.now_t, color="blue", lw=2, alpha=0.9, label="now")

            plot_cycle, last_jd, end_cycle = self._handle_cycle_logic(self.jd, self.flag)
            self._plot_cycle_overlay(plot_cycle, last_jd, end_cycle, self.mag)

            self._plot_time_markers()

            # legend
            from matplotlib.lines import Line2D
            handles = [
                Line2D([0], [0], marker='o', color='w', markerfacecolor='g',
                       markersize=6, label='quality 0'),
                Line2D([0], [0], marker='o', color='w', markerfacecolor='c',
                       markersize=6, label='quality 1'),
                Line2D([0], [0], marker='o', color='w', markerfacecolor='k',
                       markersize=6, label='quality 2'),
                Line2D([0], [0], color='blue', lw=2, label='now')
            ]
            self.axes.legend(handles=handles, loc="best", fontsize=8)

        except (FileNotFoundError, ValueError) as e:
            print(f"Lightcurve error: {e}")

    def _plot_cycle_overlay(self, plot_cycle, last_jd, end_cycle, mag):
        if not plot_cycle or last_jd is None or end_cycle is None:
            return

        # lightcurve
        self.axes.fill_between(
            [last_jd, end_cycle],
            min(mag),
            max(mag),
            color='red',
            alpha=0.1
        )

        # visibility (axes2)
        if end_cycle > int(self.now_t):
            self.axes2.fill_between(
                [int(self.now_t), end_cycle],
                -20,
                90,
                color='red',
                alpha=0.1
            )

    def _handle_cycle_logic(self, jd, flag):
        plot_cycle = False
        last_jd = None
        end_cycle = None

        if self.ob.get("cycle"):
            mk = numpy.array(flag) != 2  # dobre obserwacje
            if numpy.any(mk):
                last_jd = max(numpy.array(jd)[mk])
                end_cycle = last_jd + float(self.ob["cycle"])
                plot_cycle = True

        return plot_cycle, last_jd, end_cycle

    def _plot_phase_constraints(self, jd, flag):
        # --- zakresy zabronione ---
        if "ph_start" in self.ob and "ph_end" in self.ob:
            t0 = float(self.ob["ph_start"])
            t1 = float(self.ob["ph_end"])

            if t0 < t1:
                self.axes.axvspan(0, t0, color='red', alpha=0.05)
                self.axes.axvspan(t1, 1, color='red', alpha=0.05)
            else:
                self.axes.axvspan(t1, t0, color='red', alpha=0.05)

        # --- pokrycie fazy ---
        if "ph_mk" in self.ob:
            n_obs, filt, bin_size = self.ob["ph_mk"].split("/")
            n_obs = float(n_obs)
            bin_size = float(bin_size)

            covered = []
            mask_quality = numpy.array(flag) < 2

            jd_good = jd[mask_quality]

            for t in jd_good:
                mk1 = t > jd_good - bin_size / 2
                mk2 = t < jd_good + bin_size / 2
                mk = mk1 & mk2

                if numpy.sum(mk) >= n_obs:
                    covered.append([t - bin_size / 2, t + bin_size / 2])

            # merge
            covered.sort(key=lambda x: x[0])
            merged = []

            for seg in covered:
                if not merged or merged[-1][1] < seg[0]:
                    merged.append(seg)
                else:
                    merged[-1][1] = max(merged[-1][1], seg[1])

            for x0, x1 in merged:
                self.axes.axvspan(x0, x1, color='red', alpha=0.05)

    def _load_lightcurve(self):
        filter_name = self.file_s.currentText()

        file = self.ob.get("obs_data",f"{self.target.lower()}/{filter_name}/light-curve/{self.target.lower()}_{filter_name}_diff_light_curve.txt")
        fpath = self.data_dir+"/"+file
        tab = Table.read(fpath, format="ascii")
        return tab["file"], tab["jd_obs"], tab["mag"], tab["quality"]

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

    def _plot_zaznaczenie(self, jd, mag, flag):
        if self.zaznaczenie_i >= 0:
            self.axes.plot(self.jd[self.zaznaczenie_i],self.mag[self.zaznaczenie_i],"ro",markersize=8,zorder=10)

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

        hmin = float(self.ob.get("h_min", self.parent.tpg_cfg[self.parent.tel]["hmin"]))
        hmax = float(self.ob.get("h_max", self.parent.tpg_cfg[self.parent.tel]["hmax"]))

        t = numpy.linspace((self.current_jd), (self.current_jd) + 1, 240)
        time_range = Time(t, format="jd")

        alt, sun_alt, moon_alt, sep = self._compute_altaz(time_range)

        self._format_visibility_axes(hmin, hmax)

        # better lines
        self.axes2.plot(time_range.jd, alt, color="green", lw=2.2, label="target")
        self.axes2.plot(time_range.jd, sun_alt, "--", color="goldenrod", lw=1.4, label="sun")
        self.axes2.plot(time_range.jd, moon_alt, ":", color="black", lw=1.2, label="moon")

        self._update_moon_sep(sep)

        # ---------------------------
        # TIME LIMITS as vertical lines
        # ---------------------------
        t_start = self.parse_time(self.ob.get("t_start"))
        t_end = self.parse_time(self.ob.get("t_end"))

        if t_start:
            self.axes2.axvline(t_start, color="gray", ls="--", lw=1.2, label="t_start")

        if t_end:
            self.axes2.axvline(t_end, color="gray", ls="--", lw=1.2, label="t_end")

        # ---------------------------
        # TWILIGHT OFFSETS
        # ---------------------------
        sunset = Time(self.parent.almanac["next_sunset"]).jd
        sunrise = Time(self.parent.almanac["next_sunrise"]).jd

        sunset_dt = float(self.ob.get("sunset_dt", 0.0)) / 24.0
        sunrise_dt = float(self.ob.get("sunrise_dt", 0.0)) / 24.0

        if sunset_dt > 0:
            x = sunset + sunset_dt
            self.axes2.axvline(x, color="navy", ls=":", lw=1.5, label="after sunset")

        if sunrise_dt > 0:
            x = sunrise - sunrise_dt
            self.axes2.axvline(x, color="navy", ls=":", lw=1.5, label="before sunrise")



    def _compute_altaz(self, time_range):
        loc = EarthLocation(
            lat=self.parent.tpg_cfg["obs_lat"],
            lon=self.parent.tpg_cfg["obs_lon"],
            height=self.parent.tpg_cfg["obs_elev"]
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
        self.axes2.set_xlim((self.current_jd), (self.current_jd) + 1)

        self.axes2.axvline(self.current_jd, color="blue", lw=2)

        # horizon zones
        self.axes2.axhspan(-20, 0, facecolor='lightcoral', alpha=0.10)
        self.axes2.axhspan(0, hmin, facecolor='gray', alpha=0.08)
        self.axes2.axhspan(hmax, 90, facecolor='gray', alpha=0.08)

        self.axes2.set_ylabel("Alt [deg]")

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

        labels = list(vis.keys())

        for i, label in enumerate(labels):
            values = vis[label]
            segments = self._split_segments(nt, values)

            if segments["green"]:
                self.axes3.broken_barh(
                    segments["green"],
                    (i - 0.42, 0.84),
                    facecolors="green",
                    alpha=0.35
                )

            if segments["red"]:
                self.axes3.broken_barh(
                    segments["red"],
                    (i - 0.42, 0.84),
                    facecolors="red",
                    alpha=0.15
                )

        self.axes3.axvline(self.current_jd, color="blue", lw=2)

        self.axes3.set_yticks(range(len(labels)))
        self.axes3.set_yticklabels(labels, fontsize=9)
        self.axes3.set_xlim(self.axes2.get_xlim())
        self.axes3.set_ylim(-0.8, len(labels) - 0.2)

        #self._set_time_ticks(self.axes3, Time(nt, format="jd"))
        self.axes3.set_xlabel("UT")

    def _split_segments(self, nt, values):
        green, red = [], []

        for i in range(len(values) - 1):
            seg = (nt[i], nt[i + 1] - nt[i])
            (green if values[i] else red).append(seg)

        return {"green": green, "red": red}

    # Ticksy na wykresach

    def _format_jd_tick(self, jd, pos=None):
        try:
            jd = float(jd)
            dt = Time(jd, format="jd",scale="utc").to_datetime()
            return dt.strftime("%H:%M")
        except Exception:
            return ""


    # Funkcja ustawiająca inteligentne ticki
    def _set_time_ticks(self, ax, time_range):
        jd_start = float(time_range.jd[0])
        jd_end = float(time_range.jd[-1])

        ticks = []

        # now
        ticks.append(self.current_jd)

        # sunset / sunrise
        for key in ["next_sunset", "next_sunrise",
                    "next_moonrise", "next_moonset"]:
            val = self.parent.almanac.get(key)
            if val:
                t = Time(val).jd
                if jd_start <= t <= jd_end:
                    ticks.append(t)

        # regular every 2h
        step = 2 / 24.0
        t0 = jd_start
        while t0 <= jd_end:
            ticks.append(t0)
            t0 += step

        # sort + unique
        ticks = sorted(set([round(x, 6) for x in ticks]))

        # remove too close ticks (30 min)
        filtered = []
        min_sep = 30 / 60 / 24.0

        for t in ticks:
            if not filtered or abs(t - filtered[-1]) > min_sep:
                filtered.append(t)

        ax.set_xticks(filtered)
        ax.set_xticklabels(
            [Time(t, format="jd").datetime.strftime("%H:%M") for t in filtered],
            rotation=0
        )

        ax.set_xlim(jd_start, jd_end)


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

        # bigger TPG panel
        gs = self.fig.add_gridspec(3, 1, height_ratios=[2, 2, 1.8])

        self.axes = self.fig.add_subplot(gs[0])
        self.axes2 = self.fig.add_subplot(gs[1])
        self.axes3 = self.fig.add_subplot(gs[2], sharex=self.axes2)
        self.axes2.tick_params(axis="x", labelbottom=False)

        grid.addWidget(self.file_s, 0, 0)
        grid.addWidget(self.ephem_e, 0, 1)
        grid.addWidget(self.moon_sep_e, 0, 2)
        grid.addWidget(self.phase_c, 0, 3)
        grid.addWidget(self.canvas, 1, 0, 4, 4)

        self.toolbar = NavigationToolbar(self.canvas, self)
        grid.addWidget(self.toolbar, 5, 0, 1, 4)

        self.close_p = QPushButton("Close")
        self.close_p.clicked.connect(lambda: self.close())
        grid.addWidget(self.close_p, 6, 3)

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
        self.axes.set_rlim([-90,90+self.parent.tpg_cfg["obs_lat"]])

        self.axes2.clear()
        self.axes2.set_theta_direction(-1)
        self.axes2.set_theta_zero_location('N')
        self.axes2.set_xticks([0, 2 * 3.14 * 90 / 360, 2 * 3.14 * 180 / 360, 2 * 3.14 * 270 / 360])
        self.axes2.set_xticklabels(["N", "E", "S", "W"])
        #self.axes.set_rmax(self.rmax)
        self.axes2.set_rticks([0, 20, 40, 60, 90])
        self.axes2.set_yticklabels(["", "", "", "", ""])

        obs_location = EarthLocation(lat=self.parent.tpg_cfg["obs_lat"], lon=self.parent.tpg_cfg["obs_lon"], height=self.parent.tpg_cfg["obs_elev"])  # Warszawa

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
#     TPG WINDOW
# ######################

class TPGWindow(QWidget):
    def __init__(self, parent):
        super().__init__()

        self.parent = parent
        self.p = None

        self.setWindowTitle("TPG Planner")
        self.setMinimumSize(510, 600)
        self.setStyleSheet("font-size: 11pt;")
        self.mkUI()
        self.start_changed()

    # =====================================================
    # HELPERS
    # =====================================================

    def start_changed(self):
        if self.sunset_start_c.isChecked():
            txt = f'{self.parent.date_e.date().toPyDate()}'
            self.sunset_start_e.setText(txt)
        else:
            txt = f'{self.parent.date_e.date().toPyDate()} {self.parent.time_e.time().toPyTime().replace(microsecond=0)}'
            self.sunset_start_e.setText(txt)

    def log(self, text):
        self.log_e.append(text)

    def ok(self, text):
        self.log(f"{text} <span style='color:green;'>✔</span>")

    def err(self, text):
        self.log(f"<span style='color:red;'>{text}</span>")

    def sync_visibility(self):
        for ob in self.p.ob:
            if "visibility" in ob:
                idx = ob["index"]
                self.parent.master_data[idx]["tpg"]["visibility"] = ob["visibility"]
                self.parent.master_data[idx]["tpg"]["nightTime"] = self.p.nightTime

        self.parent.update_table()

    def ensure_loaded(self):
        if self.p is None:
            self.err("Load first!")
            return False
        return True

    # =====================================================
    # CORE
    # =====================================================

    def load(self):
        try:
            self.log_e.clear()

            tel = self.parent.tel_s.currentText()

            dt = self.sunset_start_e.text()
            dt = dt.split()

            #dt = [night_start.strftime("%Y/%m/%d")]

            wind = None
            if self.avoid_wind_c.isChecked():
                wind = float(self.avoid_wind_e.text())

            fwhm = None
            if self.fwhm_c.isChecked():
                fwhm = float(self.fwhm_e.text())

            seed = None
            if self.seed_c.isChecked():
                seed = int(self.seed_e.text())


            self.p = tpg(tel,dt,loud=True,wind=wind,fwhm=fwhm,seed=seed,save_plan=self.save_c.isChecked(),add_start_makro=self.start_macro_c.isChecked(),add_end_makro=self.end_macro_c.isChecked(),)

            self.p.Initiate()
            self.p.init_ctc()

            self.p.ob = []

            for n, data in enumerate(self.parent.master_data):
                if data.get("ob", None):
                    line = ObsValidator.convert_from_obdict(data["ob"])
                    line = line.split(" ", 1)[1]
                    tmp = self.p.parseObjects(line)
                    tmp["index"] = n
                    self.p.ob.append(tmp)

            self.p.MakeTime()
            self.p.ObjectMask()

            self.ok("Loaded data")
            self.ok(f"Night start: {self.p.start_time}")

        except Exception as e:
            self.err(str(e))

    # =====================================================
    # MASK STEPS
    # =====================================================

    def calc_ob(self):
        if not self.ensure_loaded(): return
        self.p.CalcObject()
        self.ok("Calculated objects")

    def mask_vis(self):
        if not self.ensure_loaded(): return
        self.p.MaskVisibility()
        self.sync_visibility()
        self.ok("Visibility masked")

    def mask_moon(self):
        if not self.ensure_loaded(): return
        self.p.MaskMoon()
        self.sync_visibility()
        self.ok("Moon masked")

    def mask_wind(self):
        if not self.ensure_loaded(): return
        self.p.MaskWind()
        self.sync_visibility()
        self.ok("Wind masked")

    def mask_twilight(self):
        if not self.ensure_loaded(): return
        self.p.MaskTwilight()
        self.sync_visibility()
        self.ok("Twilight masked")

    def mask_cycle(self):
        if not self.ensure_loaded(): return
        self.p.MaskCycle()
        self.sync_visibility()
        self.ok("Cycle masked")

    def mask_startend(self):
        if not self.ensure_loaded(): return
        self.p.MaskStartEnd()
        self.sync_visibility()
        self.ok("Time masked")

    def mask_phstartend(self):
        if not self.ensure_loaded(): return
        self.p.MaskPhaseStartEnd()
        self.sync_visibility()
        self.ok("Phase limits masked")

    def mask_phase(self):
        if not self.ensure_loaded(): return
        self.p.MaskPhase()
        self.sync_visibility()
        self.ok("Phase density masked")

    # =====================================================
    # FULL RUN
    # =====================================================

    def run_tpg(self):
        self.load()
        if self.p is None:
            return

        self.load()
        self.calc_ob()
        self.mask_vis()
        self.mask_moon()
        self.mask_wind()
        self.mask_twilight()
        self.mask_cycle()
        self.mask_startend()
        self.mask_phstartend()
        self.mask_phase()

        self.p.Waga()
        self.p.RandomizeList()
        self.ok(f'randomization with seed:{self.p.seed}')
        self.p.allocate()
        self.p.export()
        self.p.SavePlan()

        self.parent.update_table()

        if self.p.plan_saved:
            self.ok(f'plan saved to :{self.p.plan_filename}')

        for line in self.p.plan:
            ob_tmp = ObsPlanParser.convert_from_string(line)
            ob = ObsValidator.convert_to_obdict(ob_tmp)
            self.parent.plan_gui.add(ob)

        self.ok("TPG FINISHED")

    # =====================================================
    # UI
    # =====================================================

    def mkUI(self):
        grid = QGridLayout(self)

        r = 0

        # ---------------- OPTIONS ----------------
        self.sunset_start_c = QCheckBox("Start at sunset")
        self.sunset_start_c.setChecked(True)
        self.sunset_start_c.stateChanged.connect(self.start_changed)
        self.sunset_start_e = QLineEdit()
        grid.addWidget(self.sunset_start_c, r, 0)
        grid.addWidget(self.sunset_start_e, r, 1)
        r += 1

        self.avoid_wind_c = QCheckBox("Avoid wind")
        self.avoid_wind_e = QLineEdit()
        grid.addWidget(self.avoid_wind_c, r, 0)
        grid.addWidget(self.avoid_wind_e, r, 1)
        r += 1

        self.fwhm_c = QCheckBox("Limit FWHM")
        self.fwhm_e = QLineEdit()
        grid.addWidget(self.fwhm_c, r, 0)
        grid.addWidget(self.fwhm_e, r, 1)
        r += 1

        self.seed_c = QCheckBox("Use seed")
        self.seed_e = QLineEdit()
        grid.addWidget(self.seed_c, r, 0)
        grid.addWidget(self.seed_e, r, 1)
        r += 1

        self.save_c = QCheckBox("Save plan")
        self.save_c.setChecked(True)
        grid.addWidget(self.save_c, r, 0, 1, 2)
        r += 1

        self.start_macro_c = QCheckBox("Add start macro")
        self.start_macro_c.setChecked(True)
        grid.addWidget(self.start_macro_c, r, 0, 1, 2)
        r += 1

        self.end_macro_c = QCheckBox("Add end macro")
        self.end_macro_c.setChecked(True)
        grid.addWidget(self.end_macro_c, r, 0, 1, 2)
        r += 1

        # separator
        line = QFrame()
        line.setFrameShape(QFrame.Shape.VLine)
        grid.addWidget(line, 0, 2, 10, 1)

        # ---------------- BUTTONS ----------------
        r = 0

        self.load_p = QPushButton("1. Load / Init")
        self.load_p.clicked.connect(self.load)
        grid.addWidget(self.load_p, r, 3)
        r += 1

        self.calc_p = QPushButton("2. Calc Object")
        self.calc_p.clicked.connect(self.calc_ob)
        grid.addWidget(self.calc_p, r, 3)
        r += 1

        self.vis_p = QPushButton("3. Visibility")
        self.vis_p.clicked.connect(self.mask_vis)
        grid.addWidget(self.vis_p, r, 3)
        r += 1

        self.moon_p = QPushButton("4. Moon")
        self.moon_p.clicked.connect(self.mask_moon)
        grid.addWidget(self.moon_p, r, 3)
        r += 1

        self.wind_p = QPushButton("5. Wind")
        self.wind_p.clicked.connect(self.mask_wind)
        grid.addWidget(self.wind_p, r, 3)
        r += 1

        self.twilight_p = QPushButton("6. Twilight")
        self.twilight_p.clicked.connect(self.mask_twilight)
        grid.addWidget(self.twilight_p, r, 3)
        r += 1

        self.cycle_p = QPushButton("7. Cycle")
        self.cycle_p.clicked.connect(self.mask_cycle)
        grid.addWidget(self.cycle_p, r, 3)
        r += 1

        self.time_p = QPushButton("8. Time")
        self.time_p.clicked.connect(self.mask_startend)
        grid.addWidget(self.time_p, r, 3)
        r += 1

        self.phlim_p = QPushButton("9. Phase limits")
        self.phlim_p.clicked.connect(self.mask_phstartend)
        grid.addWidget(self.phlim_p, r, 3)
        r += 1

        self.phmk_p = QPushButton("10. Phase density")
        self.phmk_p.clicked.connect(self.mask_phase)
        grid.addWidget(self.phmk_p, r, 3)
        r += 1


        self.run_p = QPushButton("RUN FULL TPG")
        self.run_p.clicked.connect(self.run_tpg)
        grid.addWidget(self.run_p, r, 0, 1, 4)

        r += 1
        # ---------------- LOG ----------------
        self.log_e = QTextEdit()
        self.log_e.setReadOnly(True)
        self.log_e.setStyleSheet("background-color: rgb(235,235,235);")
        grid.addWidget(self.log_e, r, 0, 1, 4)

        r = r + 3
        self.close_p = QPushButton("Close")
        self.close_p.clicked.connect(self.close)
        grid.addWidget(self.close_p, r, 3)

        self.show()


def sun_moon_ephem(obs_time, lat, lon, altitude=0, horizon=0):

    obs = ephem.Observer()
    obs.lat = str(lat)
    obs.lon = str(lon)
    obs.elevation = altitude
    obs.horizon = str(horizon)  # np. '0' albo '-0:34' dla refrakcji
    obs.date = str(obs_time)

    sun = ephem.Sun()
    moon = ephem.Moon()

    # --- Sun ---
    try:
        prev_sunrise = obs.previous_rising(sun).datetime()
    except:
        prev_sunrise = None

    try:
        next_sunrise = obs.next_rising(sun).datetime()
    except:
        next_sunrise = None

    try:
        prev_sunset = obs.previous_setting(sun).datetime()
    except:
        prev_sunset = None

    try:
        next_sunset = obs.next_setting(sun).datetime()
    except:
        next_sunset = None

    # --- Moon ---
    try:
        prev_moonrise = obs.previous_rising(moon).datetime()
    except:
        prev_moonrise = None

    try:
        next_moonrise = obs.next_rising(moon).datetime()
    except:
        next_moonrise = None

    try:
        prev_moonset = obs.previous_setting(moon).datetime()
    except:
        prev_moonset = None

    try:
        next_moonset = obs.next_setting(moon).datetime()
    except:
        next_moonset = None

    # --- Moon phase ---
    moon.compute(obs)
    moon_phase = moon.phase  # %

    # --- Julian Date ---
    julian_date = ephem.julian_date(obs.date)

    return {
        "julian_date": julian_date,
        "prev_sunrise": prev_sunrise,
        "next_sunrise": next_sunrise,
        "prev_sunset": prev_sunset,
        "next_sunset": next_sunset,
        "prev_moonrise": prev_moonrise,
        "next_moonrise": next_moonrise,
        "prev_moonset": prev_moonset,
        "next_moonset": next_moonset,
        "moon_phase": moon_phase
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


def merge_schemas(base: dict, extra: dict) -> dict:
    merged = base.copy()

    # properties
    merged.setdefault("properties", {})
    merged["properties"].update(extra.get("properties", {}))

    # required (jeśli masz)
    if "required" in base or "required" in extra:
        merged["required"] = list(set(base.get("required", []) + extra.get("required", [])))

    # additionalProperties – ostrożnie (tu przykład: AND)
    if "additionalProperties" in base or "additionalProperties" in extra:
        merged["additionalProperties"] = (
            base.get("additionalProperties", True)
            and extra.get("additionalProperties", True)
        )

    return merged

def format_dt(dt):
    if dt is None:
        return "—"
    return dt.strftime("%Y-%m-%d %H:%M:%S")