
import datetime
import ephem

from PyQt6.QtWidgets import  QTableWidget, QAbstractItemView, QTableWidgetItem, QWidget,  QGridLayout, QPushButton, QFrame
from PyQt6.QtGui import QFont

from astropy.time import Time

from pyaraucaria.obs_plan.obs_plan_parser import ObsPlanParser
from pyaraucaria.ob_validator import ObsValidator

from .obs_manager_lib import seq_time


class Plan_Gui(QWidget):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        self.table_header = ["UT", "Name", "Alt", "Moon dist"]
        self.plan = []
        self.i = -1

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
        # oca.horizon = str(self.cfg[self.tel]["hsun"])

        moon = ephem.Moon()
        #sun = ephem.Sun()

        obs_time = datetime.datetime.combine(self.parent.date_e.date().toPyDate(), self.parent.time_e.time().toPyTime())
        t = Time(obs_time, scale='utc')

        for data in self.plan:
            ob = data["ob"]

            if ob.get("ob_time", None):
                slotTime = ob["ob_time"]
            elif ob.get("seq", None):
                slotTime = seq_time(ob["seq"])
            else:
                slotTime = 10

            t = t + ephem.second * slotTime
            oca.date = str(t)


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
                        txt = ob["name"]
                        item = QTableWidgetItem(txt)
                        self.table_t.setItem(i, j, item)
                    elif key == "UT":
                        txt = data["ut"].split()[1].split(":")[0] + ":" + data["ut"].split()[1].split(":")[1]
                        item = QTableWidgetItem(txt)
                        self.table_t.setItem(i, j, item)
                    elif key == "Alt":
                        txt = f'{data["alt"]:.0f}'
                        item = QTableWidgetItem(txt)
                        self.table_t.setItem(i, j, item)
                    elif key == "Moon dist":
                        txt = f'{data["moon_sep"]:.0f}'
                        item = QTableWidgetItem(txt)
                        self.table_t.setItem(i, j, item)
                    else:
                        txt = "--"
                        item = QTableWidgetItem(txt)
                        self.table_t.setItem(i, j, item)

                    # ["UT", "Name", "Alt", "Moon dist"]


        self.table_t.resizeColumnsToContents()

        if self.i is not None and self.i < self.table_t.rowCount():
            self.table_t.selectRow(self.i)

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
        # self.next_p = QPushButton('NEXT \u2192')
        # self.addStop_p = QPushButton("STOP \u2B23")
        # self.addBell_p = QPushButton('BELL \u266A')
        # self.skip_p = QPushButton('SKIP \u23ED')
        #
        # self.grid.addWidget(self.next_p, w, 0)
        # self.grid.addWidget(self.addStop_p, w, 2)
        # self.grid.addWidget(self.addBell_p, w, 3)
        # self.grid.addWidget(self.skip_p, w, 4)

        w = w + 1
        self.line_l = QFrame()
        self.line_l.setFrameShape(QFrame.Shape.HLine)
        self.line_l.setFrameShadow(QFrame.Shadow.Raised)
        self.grid.addWidget(self.line_l, w, 0, 1, 5)

        w = w + 1
        self.copy_p = QPushButton('Copy OB')
        self.grid.addWidget(self.copy_p, w, 2)

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
        self.save_p = QPushButton('Save Plan')


        self.grid.addWidget(self.save_p, w, 3, 1, 2)


        self.table_t.cellClicked.connect(self.pocisniecie_tabelki)


        # self.save_p.clicked.connect(self.savePlan)
        # self.plan_t.cellClicked.connect(self.pocisniecie_tabelki)
        # self.plan_t.horizontalHeader().sectionClicked.connect(self.pocisniecie_headera)
        #
        # self.plotPlan_p.clicked.connect(self.plot_plan)
        # self.next_p.clicked.connect(self.setNext)
        # self.skip_p.clicked.connect(self.setSkip)
        self.up_p.clicked.connect(self.pocisniecie_up)
        self.down_p.clicked.connect(self.pocisniecie_down)
        self.del_p.clicked.connect(self.pocisniecie_del)
        self.delAll_p.clicked.connect(self.pocisniecie_delAll)
        self.first_p.clicked.connect(self.pocisniecie_first)
        self.last_p.clicked.connect(self.pocisniecie_last)
        #self.swap_p.clicked.connect(self.pocisniecie_swap)
        # self.copy_p.clicked.connect(self.pocisniecie_copy)
        # self.addStop_p.clicked.connect(self.pocisniecie_addStop)
        # self.addBell_p.clicked.connect(self.pocisniecie_addBell)

        self.setLayout(self.grid)
        self.table_t.setColumnWidth(0, 30)
