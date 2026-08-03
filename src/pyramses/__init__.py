#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Python library for RAMSES dynamic simulator."""

__package_name__ = "pystepss-ulg"
__version__ = '0.0.52'
__author__ = "Petros Aristidou and Bastien Ewbank"
__copyright__ = "Petros Aristidou"
__license__ = "Apache-2.0"
__maintainer__ = "Bastien Ewbank"
__email__ = "bastien.ewbank@uliege.be"
__url__ = "https://stepss.sps-lab.org"
__status__ = "3 - Alpha"

import sys
from warnings import warn

from .cases import cfg
from .globals import __runTimeObs__, __which
from .simulator import sim, sim_exe 
from .extractor import extractor, curplot, cur


if sys.platform in ('win32', 'cygwin'):
    checkGnuplot = __which('gnuplot.exe')
else:
    checkGnuplot = __which('gnuplot')
if checkGnuplot is None:
    warn("RAMSES: Gnuplot executable could not be found in the system path, so the runtime observables are disabled.")
    __runTimeObs__ = False
else:
    __runTimeObs__ = True
