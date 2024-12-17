#!/usr/bin/env python3


import os
import sys
from PyQt6.QtWidgets import QApplication

from obs_manager_gui.obs_manager_gui import OM_Gui

def main():
    app = QApplication(sys.argv)

    main = OM_Gui(None)

    try:
        sys.exit(app.exec())
    except Exception as e:
        print(f"Application exited with an error: {e}", file=sys.stderr)

if __name__ == '__main__':
    main()