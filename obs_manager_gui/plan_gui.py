
import datetime
import ephem

from PyQt6.QtWidgets import QTableWidget, QAbstractItemView, QTableWidgetItem, QWidget, QGridLayout, QPushButton, \
    QFrame, QFileDialog
from PyQt6.QtGui import QFont

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
        self.table_header = ["UT", "Name", "Alt@UT", "Moon dist"]
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
        oca.lon = str(self.parent.cfg["obs_longitude"])
        oca.lat = str(self.parent.cfg["obs_latitude"])
        oca.elev = float(self.parent.cfg["obs_elevation"])
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
                slotTime = ob["sec"]
            else:                           # to nie przeszkodzi pozniej ut, sunrise, sunset
                slotTime = 0

            data["slotTime"] = slotTime
            t = t + ephem.second * slotTime

            if ob.get("ut", None):
                ut = ob["ut"]
                ut_time = datetime.datetime.strptime(ut, "%H:%M:%S").time()
                ut_dt = datetime.datetime.combine(t.datetime.date(), ut_time)
                if ut_dt < t.datetime:
                    ut_dt = ut_dt + datetime.timedelta(days=1)
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
                        elif ob["command_name"] == "FOCUS":
                            if ob.get("name",None):
                                txt = ob["command_name"] + " " + ob["name"]
                        elif ob["command_name"] == "SKYFLAT":
                            if ob.get("name",None):
                                txt = ob["command_name"] + " " + ob["name"]
                        elif ob["command_name"] == "DOMEFLAT":
                            txt = ob["command_name"]
                        elif ob["command_name"] == "STOP":
                            txt = ob["command_name"]
                        elif ob["command_name"] == "WAIT":
                            txt = ""
                            for k in ["sec", "ut", "sunrise", "sunset"]:
                                if ob.get(k,None):
                                    txt = txt + f'{k}={ob[k]}'
                        elif ob["command_name"] == "ZERO":
                            txt = ob["command_name"]
                        elif ob["command_name"] == "DARK":
                            txt = ob["command_name"]

                        item = QTableWidgetItem(txt)
                        self.table_t.setItem(i, j, item)

                    elif key == "UT":
                        txt = ""
                        if data.get("ut",None):
                            txt = data["ut"].split()[1].split(":")[0] + ":" + data["ut"].split()[1].split(":")[1]
                        item = QTableWidgetItem(txt)
                        self.table_t.setItem(i, j, item)

                    elif key == "Alt@UT":
                        txt = ""
                        if data.get("alt",None):
                            txt = f'{data["alt"]:.0f}'
                        item = QTableWidgetItem(txt)
                        self.table_t.setItem(i, j, item)

                    elif key == "Moon dist":
                        txt = ""
                        if data.get("moon_sep",None):
                            txt = f'{data["moon_sep"]:.0f}'
                        item = QTableWidgetItem(txt)
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

        # max_column_width = 100  # Maksymalna szerokość kolumny
        # for col in range(self.table_t.columnCount()):
        #     self.table_t.setColumnWidth(col, min(self.table_t.columnWidth(col), max_column_width))





    # def update_selection(self):
    #     self.table_t.selectRow(self.i)
    #



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
        pass

    def pocisniecie_copy(self):
        if len(self.plan) == 0:
            return
        if self.i < 0:
            return
        self.add(self.plan[self.i]["ob"])
        self.update_table()

    def pocisniecie_load(self):
        file_path, _ = QFileDialog.getOpenFileName(None,"Select a File",self.parent.cfg["master_file"],"All Files (*);;Text Files (*.txt);;Images (*.png *.jpg)")
        if file_path:
            try:
                with open(file_path, 'r') as plik:
                    for line in plik:
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

        file_path, _ = QFileDialog.getSaveFileName(self, "Save File", self.parent.cfg["master_file"], "Text Files (*.txt);;All Files (*)")
        if file_path:
            try:
                with open(file_path, "w", encoding="utf-8") as file:
                    file.write(txt)
                    print(f'Plan saved to {file_path}')
            except Exception as e:
                print(f"Error saving Plan file {file_path}: {e}")

    def mkUI(self):
        self.setWindowTitle("Plan")
        self.resize(600, 800)

        self.grid = QGridLayout()

        w = 0
        self.table_t = QTableWidget()
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
        self.line_l = QFrame()
        self.line_l.setFrameShape(QFrame.Shape.HLine)
        self.line_l.setFrameShadow(QFrame.Shadow.Raised)
        self.grid.addWidget(self.line_l, w, 0, 1, 5)

        w = w + 1
        self.edit_p = QPushButton('Edit OB')
        self.edit_p.clicked.connect(self.pocisniecie_edit)
        self.copy_p = QPushButton('Copy OB')
        self.copy_p.clicked.connect(self.pocisniecie_copy)

        self.grid.addWidget(self.edit_p, w, 4)
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
        self.oca.lon = str(self.parent.parent.cfg["obs_longitude"])
        self.oca.lat = str(self.parent.parent.cfg["obs_latitude"])
        self.oca.elev = float(self.parent.parent.cfg["obs_elevation"])
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


        # Rysowanie

        if len(self.parent.plan)>0:
            color = ["c", "m"]
            self.t = self.t_now

            j=0
            for data in self.parent.plan:
                ob = data["ob"]
                fontsize = 9
                if j==len(color): j=0
        #         tmp_ok = False
        #         if self.parent.current_i > -1 and i >= self.parent.current_i: tmp_ok = True
        #         if i >= self.parent.next_i: tmp_ok = True
        #         if 'skip' in self.parent.plan[i].keys():
        #             if self.parent.plan[i]['skip']:
        #                 tmp_ok = False
        #         if 'skip_alt' in self.parent.plan[i].keys():
        #             if self.parent.plan[i]['skip_alt']:
        #                 tmp_ok = False
        #         if 'ok' in self.parent.plan[i].keys():
        #             if not self.parent.plan[i]['ok']:
        #                 tmp_ok = False
        #
        #         if tmp_ok:
        #             if 'type' in self.parent.plan[i].keys():
        #                 if self.parent.plan[i]["type"] == "STOP":
        #                     self.axes.axvline(x=self.t, color="red",alpha=0.5)
        #                     self.axes.text(self.t,2,"STOP",rotation=90,fontsize=fontsize)
        #
        #             if "wait" in self.parent.plan[i].keys():
        #                 if len(self.parent.plan[i]["wait"]) > 0:
        #                     slotTime = float(self.parent.plan[i]["wait"])
        #                     self.axes.fill_betweenx([0, 2], self.t, self.t+ephem.second*slotTime, color="r", alpha=0.5)
        #                     self.axes.text(self.t, 3, f"WAIT {int(slotTime)}s", rotation=90, fontsize=fontsize)
        #                     self.t = self.t + ephem.second * slotTime
        #
        #
        #             if "wait_ut" in self.parent.plan[i].keys():
        #                 if len(self.parent.plan[i]["wait_ut"]) > 0:
        #                     wait_ut = ephem.Date(str(ephem.Date(self.t)).split()[0] + " " + self.parent.plan[i]["wait_ut"])
        #                     if self.t < wait_ut:
        #                         self.axes.fill_betweenx([0, 2], self.t, wait_ut, color="r",
        #                                                 alpha=0.5)
        #                         self.axes.text(self.t, 3, f"WAIT UT {wait_ut}", rotation=90, fontsize=fontsize)
        #                         self.t = wait_ut
        #
        #             if "wait_sunset" in self.parent.plan[i].keys():
        #                 if len(self.parent.plan[i]["wait_sunset"]) > 0:
        #                     self.oca.horizon = self.parent.plan[i]["wait_sunset"]
        #                     wait_ut = self.oca.next_setting(ephem.Sun(), use_center=True)
        #                     if self.t < wait_ut:
        #                         self.axes.fill_betweenx([0, 2], self.t, wait_ut, color="r",
        #                                                 alpha=0.5)
        #                         self.axes.text(self.t, 3, f"WAIT SUNSET {wait_ut}", rotation=90, fontsize=fontsize)
        #                         self.t = wait_ut
        #
        #             if "wait_sunrise" in self.parent.plan[i].keys():
        #                 if len(self.parent.plan[i]["wait_sunrise"]) > 0:
        #                     self.oca.horizon = self.parent.plan[i]["wait_sunrise"]
        #                     wait_ut = self.oca.next_rising(ephem.Sun(), use_center=True)
        #                     if self.t < wait_ut:
        #                         self.axes.fill_betweenx([0, 2], self.t, wait_ut, color="r",
        #                                                 alpha=0.5)
        #                         self.axes.text(self.t, 3, f"WAIT SUNRISE {wait_ut}", rotation=90, fontsize=fontsize)
        #                         self.t = wait_ut
        #

                slotTime = data["slotTime"]

                if slotTime < 60:
                    fontsize = 2
                if slotTime < 60 * 5:
                    fontsize = 5
                if slotTime < 60 * 10:
                    fontsize = 7
                else:
                    fontsize = 9

                if "ra" in ob.keys():
                    ra = ob["ra"]
                    dec = ob["dec"]
                    t_tab = []
                    alt_tab = []

                    t = self.t
                    while t <= self.t + ephem.second * slotTime:
                        self.oca.date = t
                        star = ephem.FixedBody()
                        star._ra = str(ra)
                        star._dec = str(dec)
                        star.compute(self.oca)

                        alt = float(star.alt) * 180.0 / ephem.pi
                        az = float(star.az) * 180.0 / ephem.pi

                        t_tab.append(t)
                        alt_tab.append(alt)

                        t = t + ephem.minute

                    self.axes.plot(t_tab,alt_tab,color=color[j])
                    self.axes.text(self.t, 93, f"{ob['name']}", color=color[j], rotation=90, fontsize=fontsize)
                    j=j+1

                self.t = self.t + ephem.second * slotTime


        self.axes.set_ylim(0, 90)
        #self.axes.set_xlim(self.t0-2*ephem.hour,self.t_end+2*ephem.hour)
        self.axes.fill_betweenx([0, 35], self.t0_dusk, self.t_end_dusk, color="grey", alpha=0.1)
        self.axes.fill_betweenx([80, 90], self.t0_dusk, self.t_end_dusk, color="grey", alpha=0.1)
        self.axes.fill_betweenx([0, 90], self.t0, self.t0_dusk, color="yellow", alpha=0.1)
        self.axes.fill_betweenx([0, 90], self.t_end_dusk, self.t_end, color="yellow", alpha=0.1)
        self.axes.axvline(x=self.t_now, color="blue")
        txt = str(self.t_now).split()[1].split(":")[0] + ":" + str(self.t_now).split()[1].split(":")[1]
        self.axes.text(self.t_now, 82, f"{txt}", rotation=90, fontsize=12)

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
        self.show()





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
