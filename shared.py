# coding=utf-8
from __future__ import print_function, division

import glob
import hashlib
import logging
import sys
import seaborn as sns
import numpy as np
from neuron import h

import settings

sys.path.append('../')
# noinspection PyUnresolvedReferences
from use_files import create_dir

create_dir = create_dir
# set up logging to file - see previous section for more details

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-30s %(levelname)-8s [%(filename)-20s:%(lineno)4d] \t %(message)s',
                    datefmt='%m-%d %H:%M')
logger = logging.getLogger('shared')
# to log to file add:
# filename='/temp/myapp.log',
# filemode='w'

t_vec = h.Vector()


def INIT():
    # load mod files
#     logger.info("load dll")
#     h.nrn_load_dll(settings.NRNMECH_PATH)
    logger.info("load hoc files")
    # load hoc files including usefulFns.hoc
    for hoc_file in glob.glob(settings.HOC_PATH + "/*.hoc"):
        h.load_file(hoc_file.replace("\\", "/"))
    # show GUI
    if settings.NEURON_GUI:
        # noinspection PyUnresolvedReferences
        from neuron import gui
        # h.showV()
        h.showRunControl()
        # h.topology()

    # general properties
    h.celsius = 37
    h.v_init = -65
    logger.info("celsius={} and v_init={}".format(h.celsius, h.v_init))

    np.random.seed(settings.RANDOM_SEED)

    # MATPLOTLIBPLOT CONFIG
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    logging.getLogger("matplotlib").setLevel(logging.WARNING)

    def set_figure_dpi(dpi=None, figure_formats=['png2x']):
        """Set resolution and format of figures.

        Parameters
        ----------
        dpi : int, optional
            Resolution of png output in dots per inch.
        figure_formats : list of strings
            Only concerns the IPython environment; see
            `IPython.core.display.set_matplotlib_formats` for more details. For
            setting the default format for saving figures, directly set
            `file_format_figs`.
        """
        try:
            import IPython
            IPython.core.display.set_matplotlib_formats(*figure_formats)
        except:
            pass
        global _dpi
        if dpi is not None: _dpi = dpi
        # need to set the following two lines as older Jupyter notebooks seem to use
        # 'savefig.dpi' and more rescent ones 'figure.dpi'
        rcParams['savefig.dpi'] = _dpi
        rcParams['figure.dpi'] = _dpi

    rcParams['pdf.fonttype'] = 42  # Output Type 3 (Type3) or Type 42 (TrueType)
    plt.rc('font', family='Arial', size=settings.SMALL_SIZE)  # controls default text sizes
    rcParams['text.latex.unicode'] = True
    rcParams['axes.titlesize'] = settings.BIGGER_SIZE  # fontsize of the axes title
    rcParams['axes.labelsize'] = settings.MEDIUM_SIZE  # fontsize of the x and y labels
    rcParams['xtick.labelsize'] = settings.SMALL_SIZE  # fontsize of the tick labels
    rcParams['ytick.labelsize'] = settings.SMALL_SIZE  # fontsize of the tick labels
    rcParams['legend.fontsize'] = settings.MEDIUM_SIZE  # legend fontsize
    rcParams['figure.titlesize'] = settings.BIGGEST_SIZE  # fontsize of the figure title
    set_figure_dpi(settings.FIG_RES_ADJUST * rcParams['figure.dpi'])
    # mpl.rcParams['savefig.dpi'] = settings.FIG_RES_ADJUST * mpl.rcParams['savefig.dpi']   # doesn't work

