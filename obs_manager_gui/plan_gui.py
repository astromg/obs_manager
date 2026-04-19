
import datetime
import ephem

from PyQt6 import QtCore
from PyQt6.QtWidgets import QTableWidget, QAbstractItemView, QTableWidgetItem, QWidget, QGridLayout, QPushButton, \
    QFrame, QFileDialog, QLineEdit, QLabel, QComboBox, QHeaderView
from PyQt6.QtGui import QFont, QColor

from astropy.time import Time


from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

from pyaraucaria.obs_plan.obs_plan_parser import ObsPlanParser
from pyaraucaria.ob_validator import ObsValidator

from .obs_manager_lib import seq_time


class Plan_Gui(QWidget):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        self.plan_plot_window = None
        self.table_header = ["UT", "Name", "Alt@UT", "Az@UT", "Moon dist"]
        self.plan = []
        self.i = -1

        self.obs_time = datetime.datetime.combine(self.parent.date_e.date().toPyDate(), self.parent.time_e.time().toPyTime())
        self.mkUI()

    def add(self,ob):
        if self.i < 0:
            self.plan.append({"ob":ob})
        elif self.i == len(self.plan) - 1:
            self.plan.append({"ob":ob})
        else:
            self.plan.insert(self.i + 1, {"ob":ob})
            self.i += 1
        self.update_table()


    def update_ephem(self):

        oca = ephem.Observer()
        oca.lon = str(self.parent.tpg_cfg["obs_lon"])
        oca.lat = str(self.parent.tpg_cfg["obs_lat"])
        oca.elev = float(self.parent.tpg_cfg["obs_elev"])
        oca.horizon = "0"

        moon = ephem.Moon()

        self.obs_time = datetime.datetime.combine(self.parent.date_e.date().toPyDate(), self.parent.time_e.time().toPyTime())
        t = Time(self.obs_time, scale='utc')

        for data in self.plan:
            ob = data["ob"]
            oca.date = str(t)
            if ob.get("ra",None) and ob.get("dec",None):
                s = ephem.FixedBody()
                s._ra = ob["ra"]
                s._dec = ob["dec"]
                s.compute(oca)
                moon.compute(oca)
                alt = float(s.alt) * 180.0 / ephem.pi
                az = float(s.az) * 180.0 / ephem.pi
                moon_sep = float(ephem.separation(s, moon)) * 180.0 / ephem.pi
                data["alt"] = alt
                data["az"] = az
                data["moon_sep"] = moon_sep

            data["ut"] = str(t)

            if ob.get("ob_time", None):
                slotTime = ob["ob_time"]
            elif ob.get("seq", None):
                slotTime = seq_time(ob["seq"])
            elif ob.get("sec", None):
                slotTime = float(ob["sec"])
            else:                           # to nie przeszkodzi pozniej ut, sunrise, sunset
                slotTime = 0

            data["slotTime"] = slotTime
            t = t + ephem.second * slotTime

            if ob.get("ut", None):
                ut = ob["ut"]
                now = t.datetime
                ut_time = datetime.datetime.strptime(ut, "%H:%M:%S").time()
                ut_dt = datetime.datetime.combine(t.datetime.date(), ut_time)

                if now.hour >= 12 and ut_time.hour < 12:
                    ut_dt = ut_dt + datetime.timedelta(days=1)

                if now < ut_dt:
                    t = Time(ut_dt, scale='utc')

            elif ob.get("sunset", None):
                oca.horizon = ob["sunset"]
                sunrise_ut = oca.next_rising(ephem.Sun(), use_center=True)
                sunset_ut = oca.next_setting(ephem.Sun(), use_center=True)
                if sunset_ut < sunrise_ut:
                    t = Time(sunset_ut.datetime(), scale='utc')

            elif ob.get("sunrise", None):
                oca.horizon = ob["sunrise"]
                sunrise_ut = oca.next_rising(ephem.Sun(), use_center=True)
                sunset_ut = oca.next_setting(ephem.Sun(), use_center=True)
                if sunrise_ut < sunset_ut:
                    t = Time(sunrise_ut.datetime(), scale='utc')



    def update_table(self):
        self.update_ephem()

        for n,col_name in enumerate(self.table_header):
            self.table_t.setHorizontalHeaderItem(n,QTableWidgetItem(col_name))

        font = QFont()
        font.setPointSize(10)  # Ustawienie mniejszej czcionki
        self.table_t.setFont(font)

        i = -1
        self.table_t.setRowCount(0)
        self.table_t.clearContents()
        self.table_t.setColumnCount(len(self.table_header))
        for data in self.plan:
            ob = data["ob"]
            i += 1
            if self.table_t.rowCount() <= i:
                self.table_t.insertRow(i)

                for j,key in enumerate(self.table_header):
                    if key == "Name":
                        txt = ""
                        if ob["command_name"] == "OBJECT":
                            if ob.get("name",None):
                                txt = ob["name"]
                                item = QTableWidgetItem(txt)
                        elif ob["command_name"] == "FOCUS":
                            if ob.get("name",None):
                                txt = ob["command_name"] + " " + ob["name"]
                                item = QTableWidgetItem(txt)
                        elif ob["command_name"] == "SKYFLAT":
                            if ob.get("name",None):
                                txt = ob["command_name"] + " " + ob["name"]
                                item = QTableWidgetItem(txt)
                        elif ob["command_name"] == "DOMEFLAT":
                            txt = ob["command_name"]
                            item = QTableWidgetItem(txt)
                        elif ob["command_name"] == "STOP":
                            txt = ob["command_name"]
                            item = QTableWidgetItem(txt)
                            item.setForeground(QColor("brown"))
                        elif ob["command_name"] == "BELL":
                            txt = ob["command_name"]
                            item = QTableWidgetItem(txt)
                            item.setForeground(QColor("darkorchid"))
                        elif ob["command_name"] == "WAIT":
                            txt = ""
                            for k in ["sec", "ut", "sunrise", "sunset"]:
                                if ob.get(k,None):
                                    txt = txt + f'{k}={ob[k]}'
                            item = QTableWidgetItem(txt)
                            item.setForeground(QColor("dodgerblue"))
                        elif ob["command_name"] == "ZERO":
                            txt = ob["command_name"]
                            item = QTableWidgetItem(txt)
                        elif ob["command_name"] == "DARK":
                            txt = ob["command_name"]
                            item = QTableWidgetItem(txt)
                        else:
                            txt = ob["command_name"]
                            item = QTableWidgetItem(txt)

                        self.table_t.setItem(i, j, item)

                    elif key == "UT":
                        txt = ""
                        if data.get("ut",None):
                            txt = data["ut"].split()[1].split(":")[0] + ":" + data["ut"].split()[1].split(":")[1]
                            #txt = data["ut"]
                        item = QTableWidgetItem(txt)
                        item.setBackground(QColor(235, 235, 235))
                        self.table_t.setItem(i, j, item)

                    elif key == "Alt@UT":
                        txt = ""
                        if data.get("alt",None):
                            txt = f'{data["alt"]:.0f}'
                        item = QTableWidgetItem(txt)
                        item.setBackground(QColor(235, 235, 235))
                        self.table_t.setItem(i, j, item)

                    elif key == "Az@UT":
                        txt = ""
                        if data.get("az",None):
                            txt = f'{data["az"]:.0f}'
                        item = QTableWidgetItem(txt)
                        item.setBackground(QColor(235, 235, 235))
                        self.table_t.setItem(i, j, item)

                    elif key == "Moon dist":
                        txt = ""
                        if data.get("moon_sep",None):
                            txt = f'{data["moon_sep"]:.0f}'
                        item = QTableWidgetItem(txt)
                        item.setBackground(QColor(235, 235, 235))
                        self.table_t.setItem(i, j, item)

                    else:
                        txt = "--"
                        item = QTableWidgetItem(txt)
                        self.table_t.setItem(i, j, item)



        self.table_t.resizeColumnsToContents()

        if self.i is not None and self.i < self.table_t.rowCount():
            self.table_t.selectRow(self.i)

        if self.plan_plot_window:
            self.plan_plot_window.refresh()

    def pocisniecie_tabelki(self,i,j):
        self.i = i

    def pocisniecie_delAll(self):
        self.plan = []
        self.update_table()


    def pocisniecie_del(self):
        if len(self.plan) == 0:
            return
        if self.i < 0:
            return
        self.plan.pop(self.i)
        self.update_table()

    def pocisniecie_first(self):
        if len(self.plan) == 0:
            return
        if self.i < 0:
            return
        item = self.plan.pop(self.i)
        self.plan.insert(0, item)
        self.i = 0
        self.update_table()

    def pocisniecie_last(self):
        if len(self.plan) == 0:
            return
        if self.i < 0:
            return
        self.plan.append(self.plan[self.i])
        self.plan.pop(self.i)
        self.i = len(self.plan) - 1
        self.update_table()

    def pocisniecie_up(self):
        if len(self.plan) == 0:
            return
        if self.i <= 0:
            return
        self.plan[self.i - 1], self.plan[self.i] = self.plan[self.i], self.plan[self.i - 1]
        self.i -= 1
        self.update_table()

    def pocisniecie_down(self):
        if len(self.plan) == 0:
            return
        if self.i < 0:
            return
        if self.i == len(self.plan) - 1:
            return
        self.plan[self.i + 1], self.plan[self.i] = self.plan[self.i], self.plan[self.i + 1]
        self.i += 1
        self.update_table()

    def plot_plan(self):
        self.plan_plot_window = PlotWindow(self)

    def pocisniecie_edit(self):
        self.edit_window = EditWindow(self)

    def pocisniecie_add(self):
        self.add({"command_name":"OBJECT"})
        self.update_table()
        self.add_window = EditWindow(self)

    def pocisniecie_copy(self):
        if len(self.plan) == 0:
            return
        if self.i < 0:
            return
        self.add(self.plan[self.i]["ob"])
        self.update_table()

    def pocisniecie_load(self):
        file_path, _ = QFileDialog.getOpenFileName(None,"Select a File",self.parent.tpg_cfg["plan_catalog"],"All Files (*);;Text Files (*.txt);;Images (*.png *.jpg)")
        if file_path:
            try:
                with open(file_path, 'r') as plik:
                    for line in plik:
                        if "BELL" in line.split()[0]:
                            ob = {"command_name":"BELL"}
                        else:
                            ob_tmp = ObsPlanParser.convert_from_string(line)
                            ob = ObsValidator.convert_to_obdict(ob_tmp)
                        if ob:
                            self.add(ob)
            except Exception as e:
                print(f"Error loading Plan file {file_path}: {e}")
        self.update_table()


    def pocisniecie_save(self):
        txt = ""
        for data in self.plan:
            ob = data["ob"]
            txt = txt + ObsValidator.convert_from_obdict(ob) + "\n"

        file_path, _ = QFileDialog.getSaveFileName(self, "Save File", self.parent.tpg_cfg["plan_catalog"], "Text Files (*.txt);;All Files (*)")
        if file_path:
            try:
                with open(file_path, "w", encoding="utf-8") as file:
                    file.write(txt)
                    print(f'Plan saved to {file_path}')
            except Exception as e:
                print(f"Error saving Plan file {file_path}: {e}")

    def mkUI(self):
        self.setWindowTitle("Plan")
        self.resize(500, 800)

        self.grid = QGridLayout()

        w = 0
        self.table_t = QTableWidget()
        self.table_t.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table_t.setHorizontalHeaderLabels(self.table_header)
        self.table_t.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)  # zaznaczenie całego wiersza
        self.table_t.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table_t.setStyleSheet("selection-background-color: rgb(217,239,217); selection-color: black; ")
        self.table_t.verticalHeader().hide()
        self.table_t.verticalHeader().setDefaultSectionSize(15)

        self.grid.addWidget(self.table_t, w, 0, 7, 5)

        w = w + 7
        self.plotPlan_p = QPushButton('Plot Plan')
        self.grid.addWidget(self.plotPlan_p, w, 4, 1, 1)

        w = w + 1
        self.line_l = QFrame()
        self.line_l.setFrameShape(QFrame.Shape.HLine)
        self.line_l.setFrameShadow(QFrame.Shadow.Raised)
        self.grid.addWidget(self.line_l, w, 0, 1, 5)

        w = w + 1
        # self.addStop_p = QPushButton("STOP \u2B23")
        # self.addBell_p = QPushButton('BELL \u266A')
        # self.addWait_p = QPushButton('WAIT \u266A')
        #
        # self.grid.addWidget(self.addStop_p, w, 2)
        # self.grid.addWidget(self.addBell_p, w, 3)

        w = w + 1
        self.edit_p = QPushButton('Edit OB')
        self.edit_p.clicked.connect(self.pocisniecie_edit)

        self.add_p = QPushButton('Add OB')
        self.add_p.clicked.connect(self.pocisniecie_add)

        self.copy_p = QPushButton('Copy OB')
        self.copy_p.clicked.connect(self.pocisniecie_copy)

        self.grid.addWidget(self.add_p, w, 4)
        self.grid.addWidget(self.edit_p, w, 2)
        self.grid.addWidget(self.copy_p, w, 0)

        w = w + 1
        self.line_l = QFrame()
        self.line_l.setFrameShape(QFrame.Shape.HLine)
        self.line_l.setFrameShadow(QFrame.Shadow.Raised)
        self.grid.addWidget(self.line_l, w, 0, 1, 5)

        w = w + 1
        self.del_p = QPushButton('Del')
        self.up_p = QPushButton('Up')
        #self.swap_p = QPushButton('Swap')
        self.first_p = QPushButton('First')

        self.grid.addWidget(self.del_p, w, 0)
        self.grid.addWidget(self.up_p, w, 2)
        #self.grid.addWidget(self.swap_p, w, 3)
        self.grid.addWidget(self.first_p, w, 4)

        w = w + 1
        self.delAll_p = QPushButton('Del All')
        self.down_p = QPushButton('Down')
        self.last_p = QPushButton('Last')

        self.grid.addWidget(self.delAll_p, w, 0)
        self.grid.addWidget(self.down_p, w, 2)
        self.grid.addWidget(self.last_p, w, 4)

        w = w + 1
        self.line_l = QFrame()
        self.line_l.setFrameShape(QFrame.Shape.HLine)
        self.line_l.setFrameShadow(QFrame.Shadow.Raised)
        self.grid.addWidget(self.line_l, w, 0, 1, 5)

        w = w + 1
        self.load_p = QPushButton('Load Plan')
        self.load_p.clicked.connect(self.pocisniecie_load)
        self.save_p = QPushButton('Save Plan')
        self.save_p.clicked.connect(self.pocisniecie_save)

        self.grid.addWidget(self.load_p, w, 0, 1, 2)
        self.grid.addWidget(self.save_p, w, 3, 1, 2)


        self.table_t.cellClicked.connect(self.pocisniecie_tabelki)


        # self.save_p.clicked.connect(self.savePlan)
        # self.plan_t.cellClicked.connect(self.pocisniecie_tabelki)
        # self.plan_t.horizontalHeader().sectionClicked.connect(self.pocisniecie_headera)
        #
        self.plotPlan_p.clicked.connect(self.plot_plan)
        # self.next_p.clicked.connect(self.setNext)
        # self.skip_p.clicked.connect(self.setSkip)
        self.up_p.clicked.connect(self.pocisniecie_up)
        self.down_p.clicked.connect(self.pocisniecie_down)
        self.del_p.clicked.connect(self.pocisniecie_del)
        self.delAll_p.clicked.connect(self.pocisniecie_delAll)
        self.first_p.clicked.connect(self.pocisniecie_first)
        self.last_p.clicked.connect(self.pocisniecie_last)
        #self.swap_p.clicked.connect(self.pocisniecie_swap)
        # self.addStop_p.clicked.connect(self.pocisniecie_addStop)
        # self.addBell_p.clicked.connect(self.pocisniecie_addBell)

        self.setLayout(self.grid)
        self.table_t.setColumnWidth(0, 30)



# #############################################
# ######### OKNO WYKRESU (PLOT PLAN) ##########
# #############################################

class PlotWindow(QWidget):
    def __init__(self, parent):
        super(PlotWindow, self).__init__()
        self.parent = parent
        self.obs_time = self.parent.obs_time

        self.setStyleSheet("font-size: 11pt;")
        #self.set_initial_geometry(100,100,1800,600)
        self.setMinimumSize(1800,600)
        self.mkUI()
        self.refresh()
        self.close_p.clicked.connect(lambda: self.close())

    def refresh(self):
        self.axes.clear()

        self.oca = ephem.Observer()
        self.oca.lon = str(self.parent.parent.tpg_cfg["obs_lon"])
        self.oca.lat = str(self.parent.parent.tpg_cfg["obs_lat"])
        self.oca.elev = float(self.parent.parent.tpg_cfg["obs_elev"])
        self.oca.horizon = "0"

        self.obs_time = self.parent.obs_time
        self.oca.date = str(Time(self.obs_time, scale='utc'))
        self.t_now = self.oca.date

        # liczenie wschodu slonca i zachodu
        t1 = self.oca.next_setting(ephem.Sun(), use_center=True) - self.oca.date
        t2 = self.oca.next_rising(ephem.Sun(), use_center=True) - self.oca.date
        if t1 < t2 :
            self.t0 = self.oca.next_setting(ephem.Sun(),use_center=True)
        else:
            self.t0 = self.oca.previous_setting(ephem.Sun(),use_center=True)
        self.oca.date = self.t0
        self.t_end = self.oca.next_rising(ephem.Sun(),use_center=True)

        # liczenie zmierzchu
        self.oca.horizon = "-18"
        self.t0_dusk = self.oca.next_setting(ephem.Sun(),use_center=True)
        self.t_end_dusk = self.oca.next_rising(ephem.Sun(),use_center=True)

        if len(self.parent.plan) > 0:

            colors = ["c", "m"]
            j = 0

            for i, data in enumerate(self.parent.plan):

                ob = data["ob"]
                slotTime = data.get("slotTime", 0)
                start_t = ephem.Date(data["ut"])

                if i < len(self.parent.plan) - 1:
                    next_t = ephem.Date(self.parent.plan[i + 1]["ut"])
                else:
                    next_t = start_t + ephem.second * slotTime

                color = colors[j % len(colors)]

                if ob.get("ra") and ob.get("dec"):
                    j += 1

                if slotTime < 60:
                    fontsize = 3
                elif slotTime < 5 * 60:
                    fontsize = 5
                elif slotTime < 10 * 60:
                    fontsize = 7
                else:
                    fontsize = 9


                if ob.get("command_name", None):
                    if ob.get("command_name") == "STOP":
                        self.axes.axvline(x=start_t, color="red")
                        self.axes.text(start_t + 2*ephem.minute, 7, "STOP", rotation=90, fontsize=8, color="red",va="bottom")
                    elif ob.get("command_name") == "BELL":
                        self.axes.axvline(x=start_t, color="skyblue")
                        self.axes.text(start_t+ 2*ephem.minute, 7, "BELL", rotation=90, fontsize=8, color="skyblue",va="bottom")
                    elif ob.get("command_name") == "DOMEFLAT":
                        self.axes.fill_betweenx([0, 4], start_t, next_t, color="paleturquoise", alpha=0.5)
                        self.axes.text((start_t+next_t)/2, 7, "DOMEFLAT", rotation=90, fontsize=8, color="paleturquoise",va="bottom")

                if ob.get("ra") and ob.get("dec"):

                    ra = ob["ra"]
                    dec = ob["dec"]

                    end_t = start_t + ephem.second * slotTime

                    t_tab = []
                    alt_tab = []

                    t = start_t
                    while t <= end_t:
                        self.oca.date = t

                        star = ephem.FixedBody()
                        star._ra = str(ra)
                        star._dec = str(dec)
                        star.compute(self.oca)

                        alt = float(star.alt) * 180.0 / ephem.pi

                        t_tab.append(t)
                        alt_tab.append(alt)

                        t += ephem.minute

                    self.axes.plot(t_tab, alt_tab, color=color, linewidth=2)
                    self.axes.text(start_t, 93, ob.get("name", "target"), rotation=90, fontsize=fontsize, color=color,va="bottom", ha="left")

                    # if end_t < next_t:
                    #     self.axes.fill_betweenx([0, 4], end_t, next_t, color="red", alpha=0.9)

                elif ob.get("ut"):
                    if start_t < next_t:
                        self.axes.fill_betweenx([0, 4], start_t, next_t, color="purple", alpha=0.5)
                        self.axes.text((start_t+next_t)/2, 7, f"WAIT UT {ob['ut']}", rotation=90, fontsize=8, color="purple", va="bottom")

                elif ob.get("sunset"):
                    if start_t < next_t:
                        self.axes.fill_betweenx([0, 4], start_t, next_t, color="darkorange", alpha=0.5)
                        self.axes.text((start_t+next_t)/2, 7, "WAIT SUNSET", rotation=90, fontsize=8, color="red", va="bottom")

                elif ob.get("sunrise"):
                    if start_t < next_t:
                        self.axes.fill_betweenx([0, 4], start_t, next_t, color="darkorange", alpha=0.5)
                        self.axes.text((start_t+next_t)/2, 7, "WAIT SUNRISE", rotation=90, fontsize=8, color="red", va="bottom")

                elif ob.get("sec"):
                    self.axes.fill_betweenx([0, 4], start_t, next_t, color="blue", alpha=0.5)
                    self.axes.text((start_t+next_t)/2, 7, f"WAIT {int(slotTime)}s", rotation=90, fontsize=8, color="blue",va="bottom")


        self.axes.grid(True, alpha=0.25)
        self.axes.set_ylim(0, 90)
        #self.axes.set_xlim(self.t0-2*ephem.hour,self.t_end+2*ephem.hour)
        self.axes.fill_betweenx([0, 35], self.t0_dusk, self.t_end_dusk, color="grey", alpha=0.1)
        self.axes.fill_betweenx([80, 90], self.t0_dusk, self.t_end_dusk, color="grey", alpha=0.1)
        self.axes.fill_betweenx([0, 90], self.t0, self.t0_dusk, color="yellow", alpha=0.1)
        self.axes.fill_betweenx([0, 90], self.t_end_dusk, self.t_end, color="yellow", alpha=0.1)
        self.axes.axvline(x=self.t_now, color="blue")
        txt = str(self.t_now).split()[1].split(":")[0] + ":" + str(self.t_now).split()[1].split(":")[1]
        self.axes.text(self.t_now + 2*ephem.minute, 30, f"{txt}", rotation=90, fontsize=12)

        xtics = [self.t0, self.t0_dusk, self.t_end_dusk, self.t_end]
        t =  ephem.Date(self.t0_dusk+30*ephem.minute)
        while t < ephem.Date(self.t_end_dusk-30*ephem.minute):
            t = ephem.Date(t) + ephem.hour
            h = str(ephem.Date(t)).split()
            xtics.append( ephem.Date(h[0]+" "+h[1].split(":")[0]+":00:00"))
        xtics_labels = [str(x).split()[1].split(":")[0]+":"+str(x).split()[1].split(":")[1] for x in xtics]
        self.axes.set_xticks(xtics)
        self.axes.set_xticklabels(xtics_labels,rotation=45,minor=False)

        self.axes.set_yticks([0, 35, 80, 90])
        self.axes.set_yticklabels(["0 deg", "35 deg", "80 deg", "90 deg"])

        #self.axes.set_ylabel("altitude")
        #self.axes.set_xlabel("UT")
        self.fig.subplots_adjust(bottom=0.12,top=0.8,left=0.08,right=0.98)
        #self.fig.tight_layout()

        self.canvas.draw()



    def mkUI(self):
        grid = QGridLayout()
        self.fig = Figure((1.0, 1.0), linewidth=-1, dpi=100)
        self.canvas = FigureCanvas(self.fig)
        self.axes = self.fig.add_subplot(111)
        grid.addWidget(self.canvas,0,0,1,2)

        self.toolbar = NavigationToolbar(self.canvas,self)
        grid.addWidget(self.toolbar, 1, 0, 1, 2)

        self.close_p = QPushButton('Close')
        grid.addWidget(self.close_p, 2, 1)

        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 0)
        grid.setRowStretch(0, 1)
        grid.setRowStretch(1, 0)
        grid.setRowStretch(2, 0)

        self.setLayout(grid)
        self.show()

# ########################################
#              EDIT WINDOW
# ########################################

class EditWindow(QWidget):
    def __init__(self, parent):
        super().__init__()

        self.parent = parent

        self.base_schema = ObsValidator.load_schema("base_schema.yaml")
        self.command_rules = ObsValidator.load_schema("base_rules.yaml")
        self.validator = ObsValidator(self.base_schema, self.command_rules)

        self.setWindowTitle(" OB EDIT WINDOW")
        self.resize(600, 600)
        self.setStyleSheet("font-size: 11pt;")

        self.updating = False
        self.initial_load = True

        self.mkUI()
        self.load_initial()

    def save_ob(self):
        if self.validate_current():
            self.parent.plan[self.parent.i]["ob"] = self.collect_table_data()
            self.parent.update_table()

    def load_initial(self):
        try:
            self.ob = self.parent.plan[self.parent.i]["ob"]
            txt = ObsValidator.convert_from_obdict(self.ob)
        except Exception:
            print("ERROR: Loading OB")
            txt = "ERROR"

        self.block_e.setText(txt)
        self.block_e.setCursorPosition(0)
        self.initial_load = False


    def build_block(self):
        data = self.collect_table_data()
        txt = self.validator.convert_from_obdict(data)
        if txt:
            return txt
        return ""

    def collect_table_data(self):
        data = {"command_name": self.type_s.currentText()}

        for row in range(self.tab_t.rowCount()):
            key = self.tab_t.item(row, 0).data(QtCore.Qt.ItemDataRole.UserRole)
            val_item = self.tab_t.item(row, 1)
            if not val_item:
                continue

            val = val_item.text().strip()
            if val != "":
                data[key] = val

        return data

    def schema_prop(self, key):
        return self.base_schema.get("properties", {}).get(key, {})

    def tooltip_for(self, key):
        prop = self.schema_prop(key)

        desc = prop.get("description", "")
        typ = prop.get("type", "")
        enum = prop.get("enum", None)

        lines = []

        if desc:
            lines.append(desc)

        if typ:
            lines.append(f"Type: {typ}")

        if enum:
            lines.append("Allowed: " + ", ".join(map(str, enum)))

        return "\n".join(lines)

    def example_for(self, key):
        prop = self.schema_prop(key)
        ex = prop.get("examples", [])
        if ex:
            return str(ex[0])
        return ""


    def rebuild_table(self, command_name, values=None):
        self.updating = True

        self.tab_t.blockSignals(True)
        self.tab_t.setRowCount(0)

        allowed = self.command_rules[command_name]["allowed"]

        visible = [x for x in allowed if x != "command_name"]

        for r, key in enumerate(visible):
            self.tab_t.insertRow(r)

            # parameter
            item0 = QTableWidgetItem(key)
            item0.setBackground(QColor(235, 235, 235))
            item0.setFlags(item0.flags() & ~QtCore.Qt.ItemFlag.ItemIsEditable)
            item0.setData(QtCore.Qt.ItemDataRole.UserRole, key)
            item0.setToolTip(self.tooltip_for(key))
            self.tab_t.setItem(r, 0, item0)

            # value
            val = ""
            if values and key in values:
                val = str(values[key])

            item1 = QTableWidgetItem(val)
            self.tab_t.setItem(r, 1, item1)

            # example
            item2 = QTableWidgetItem(self.example_for(key))
            item2.setBackground(QColor(235, 235, 235))
            item2.setFlags(item2.flags() & ~QtCore.Qt.ItemFlag.ItemIsEditable)
            self.tab_t.setItem(r, 2, item2)

        self.tab_t.blockSignals(False)
        self.updating = False


    def block_changed(self):
        if not self.initial_load:
            self.status_l.setText("\u2699 OB changed")
            self.status_l.setStyleSheet("color: blue; font-weight: normal;")

        if self.updating:
            return

        txt = self.block_e.text()
        ob_tmp = ObsPlanParser.convert_from_string(txt)
        if ob_tmp:
            ob = ObsValidator.convert_to_obdict(ob_tmp)
            cmd = ob.get("command_name", "OBJECT")
        else:
            ob = {"command_name", "OBJECT"}
            cmd = "OBJECT"

        if cmd not in self.base_schema["properties"]["command_name"]["enum"]:
            ob = {"command_name", "OBJECT"}
            cmd = "OBJECT"

        self.updating = True
        self.type_s.setCurrentText(cmd)
        self.rebuild_table(cmd, ob)
        self.updating = False

    def command_changed(self):
        if self.updating:
            return

        cmd = self.type_s.currentText()
        current = self.collect_table_data()
        self.rebuild_table(cmd, current)

        self.refresh_block()

    def table_changed(self):
        self.status_l.setText("\u2699 OB changed")
        self.status_l.setStyleSheet("color: blue; font-weight: normal;")

        if self.updating:
            return
        self.refresh_block()

        for r in range(self.tab_t.rowCount()):
            it = self.tab_t.item(r, 1)
            if it:
                it.setBackground(QColor("white"))


    def refresh_block(self):
        self.updating = True
        txt = self.build_block()
        self.block_e.setText(txt)
        self.updating = False


    def validate_current(self):
        data = self.collect_table_data()

        result = self.validator.validate_ob(data)

        row_map = {}
        for r in range(self.tab_t.rowCount()):
            key = self.tab_t.item(r, 0).data(QtCore.Qt.ItemDataRole.UserRole)
            row_map[key] = r


        self.tab_t.blockSignals(True)
        self.updating = True

        # clear colors
        for r in range(self.tab_t.rowCount()):
            it = self.tab_t.item(r, 1)
            if it:
                it.setBackground(QColor("white"))

        # apply colors
        for key, state in result["result"].items():
            if key not in row_map:
                continue

            r = row_map[key]
            color = QColor(217, 239, 217) if state is True else QColor(255, 180, 80)

            it = self.tab_t.item(r, 1)
            if it:
                it.setBackground(color)

        self.tab_t.blockSignals(False)
        self.updating = False

        if result["valid"]:
            self.status_l.setText("✅ Valid OB ")
            #self.status_l.setText("\u2714 Valid OB ")
            self.status_l.setStyleSheet("color: green; font-weight: bold;")
        else:
            self.status_l.setText("❌ Validation errors")
            #self.status_l.setText("\u274C Validation errors")
            self.status_l.setStyleSheet("color: orange; font-weight: bold;")

        return result["valid"]

    def mkUI(self):
        grid = QGridLayout()

        r = 0
        # block line
        self.block_e = QLineEdit()
        self.block_e.textChanged.connect(self.block_changed)
        grid.addWidget(self.block_e, r, 0, 1, 3)

        r += 1
        # command selector
        self.type_l = QLabel("TYPE")
        self.type_s = QComboBox()
        self.type_s.addItems(list(self.command_rules.keys()))
        self.type_s.currentTextChanged.connect(self.command_changed)

        grid.addWidget(self.type_l, r, 0)
        grid.addWidget(self.type_s, r, 1, 1, 2)

        r += 1
        # table
        self.tab_t = QTableWidget()
        self.tab_t.setColumnCount(3)
        self.tab_t.setHorizontalHeaderLabels(["Parameter", "Value", "Example"])
        self.tab_t.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.tab_t.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        self.tab_t.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)

        self.tab_t.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)  # zaznaczenie całego wiersza
        self.tab_t.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.tab_t.setStyleSheet("selection-background-color: rgb(200,220,255); selection-color: black; ")
        self.tab_t.verticalHeader().hide()

        self.tab_t.itemChanged.connect(self.table_changed)

        grid.addWidget(self.tab_t, r, 0, 1, 3)

        r += 1
        # status
        self.status_l = QLabel("Not validated")
        grid.addWidget(self.status_l, r, 0, 1, 2)

        self.validate_p = QPushButton("Validate OB")
        self.validate_p.clicked.connect(self.validate_current)
        grid.addWidget(self.validate_p, r, 1, 1 ,2)

        r += 1

        self.save_p = QPushButton("Save")
        self.save_p.clicked.connect(self.save_ob)
        grid.addWidget(self.save_p, r, 2)

        self.close_p = QPushButton("Close")
        self.close_p.clicked.connect(self.close)
        grid.addWidget(self.close_p, r, 0)

        self.setLayout(grid)
        self.show()