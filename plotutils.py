from __future__ import division

import logging

import numpy as np
from matplotlib import pyplot as plt, colors
from neuron import h

import settings
from shared import t_vec

logger = logging.getLogger('plotutils')


def annotate_cols_rows(axes, cols=None, rows=None, row_pad=5):
    if rows is None:
        rows = []
    if cols is None:
        cols = []
    annotate_cols(axes, cols)
    annotate_rows(axes, rows, pad=row_pad)


def annotate_cols(axes, labels):
    """SET COLUMN TITLE
    """
    if len(axes) > 0:
        for ax, col in zip(axes[0], labels):
            ax.set_title(col)


def annotate_rows(axes, labels, pad=5):
    """SET ROW TITLE
    """
    if len(axes) > 0:
        for ax, row in zip(axes[:, 0], labels):
            ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                        xycoords=ax.yaxis.label, textcoords='offset points',
                        size='large', ha='right', va='center')


def FORMAT_FIG(fig=None, ax=None,
               remove_top=True,
               remove_bottom=False,
               remove_right=True,
               remove_left=False,
               adjust_left=settings.FIG_ADJUST_LEFT,
               adjust_top=settings.FIG_ADJUST_TOP,
               scalebar=None):
    """

    :param fig:
    :param ax:
    :param remove_top:
    :param remove_bottom:
    :param remove_right:
    :param remove_left:
    :param adjust_left:
    :param adjust_top:
    :param scalebar:

    :type ax: np.ndarray or plt.Axes
    :type scalebar: dict or bool


    """
    if type(ax) == list or type(ax) == np.ndarray:
        for sub_ax in ax:
            FORMAT_FIG(fig, sub_ax, remove_top=remove_top, remove_bottom=remove_bottom,
                       remove_right=remove_right, remove_left=remove_left,
                       adjust_left=adjust_left, adjust_top=adjust_top)
    else:
        if scalebar is not None:
            from scalebars import add_scalebar
            if type(scalebar) is dict:
                sb = add_scalebar(ax, **scalebar)
            else:
                sb = add_scalebar(ax)
        elif ax is not None:
            if remove_top:
                ax.spines['top'].set_visible(False)
            if remove_right:
                ax.spines['right'].set_visible(False)
            if remove_bottom:
                ax.spines['bottom'].set_visible(False)
                x_axis = ax.axes.get_xaxis()
                x_axis.set_visible(False)
            if remove_left:
                ax.spines['left'].set_visible(False)
                y_axis = ax.axes.get_yaxis()
                y_axis.set_visible(False)
        if fig is not None:
            fig.tight_layout()
            # tight_layout doesn't take row labels into account. We'll need
            # to make some room. These numbers are manually tweaked.
            # You could automatically calculate them, but it's a pain.
            if adjust_left:
                fig.subplots_adjust(left=adjust_left)
            if adjust_top:
                fig.subplots_adjust(top=adjust_top)


def plot_input_events(ax=None, input_events={}, cmap='coolwarm', y_offset=0, marker='v',
                      inhib_syn_type=None, exc_syn_type=None):
    cmap_coolwarm = plt.get_cmap(cmap)
    for key, event_vec in input_events.items():
        if inhib_syn_type is not None and inhib_syn_type in key:
            # first color in coolwarm is blue
            base_color = colors.rgb2hex(cmap_coolwarm(0)[:3])
            y_offset_iter = y_offset - 5
        elif exc_syn_type is not None and exc_syn_type in key:
            # first color in coolwarm is red
            base_color = colors.rgb2hex(cmap_coolwarm(255)[:3])
            y_offset_iter = y_offset + 5
        else:
            base_color = '#000000'
            y_offset_iter = y_offset
        color = opacity(50, base_color)  # base color at 50% opacity
        events = np.array(event_vec)
        ax.plot(events, np.ones(len(events)) * y_offset_iter,
                color=color, linestyle='none', marker=marker, markersize=6, label=None)


def plot_var(x, y,
             title=None,
             xlabel=settings.TIME + " " + settings.UNITS(settings.ms),
             ylabel=None,
             ax=None):
    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(x, y,label=ylabel)
    if title is not None:
        ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xlabel == 'time (ms)':
        ax.set_xlim([0, h.tstop])
    return ax


def plot_v(vec_hoc, title=None, ax=None, section=None, location=None):
    if section is None:
        section = 'soma'
    if location is None:
        location = 0.5
    y = vec_hoc['v'][section][location]
    return plot_var(t_vec, y, title=title,
                    ylabel=" ".join([settings.MEMBRANE_POTENTIAL, settings.UNITS(settings.mV)]),
                    ax=ax)


def plot_cli(vec_hoc, title=None, ax=None, section=None, location=None):
    if section is None:
        section = 'soma'
    if location is None:
        location = 0.5
    y = vec_hoc['cli'][section][location]
    return plot_var(t_vec, y, title=title,
                    ylabel=" ".join([settings.CLI, settings.UNITS(settings.mM)]),
                    ax=ax)


def fill_gaps(ax,**kwargs):
    logger.info("Filling gaps")
    handles, labels = ax.get_legend_handles_labels()
    ys = {}
    min_has_max = {}
    for line in handles * 2:
        x_values = line._x
        y_values = line._y
        min_x = min(x_values)
        max_x = max(x_values)
        ys[min_x] = y_values[x_values == min_x][0]
        ys[max_x] = y_values[x_values == max_x][0]
        if min_x in min_has_max and min_has_max[min_x] == -10:
            min_has_max[min_x] = True
        else:
            min_has_max[min_x] = -1
        if max_x in min_has_max and min_has_max[max_x] == -1:
            min_has_max[max_x] = True
        else:
            min_has_max[max_x] = -10

    # sort the keys and remove the minimum value
    sorted_keys = sorted(min_has_max)
    for i in range(len(sorted_keys)):
        x1 = sorted_keys[i - 1]
        x2 = sorted_keys[i]
        if (min_has_max[x1] == -10 or min_has_max[x1] is True) and min_has_max[x2] == -1:
            y1 = ys[x1]
            y2 = ys[x2]
            ax.plot([x1, x2], [y1, y2], label='none',**kwargs)


# COLORS
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
            cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def opacity(level, color):
    if level > 1:
        level /= 100
    _opacity = "%0.2X" % round((level * 255))  # note: not sure if round or int works better here
    if len(color) == 9:
        # already has an opacity applied, so reset it by removing last 2
        color = color[:-2]
    return color + _opacity


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    h,l,s = colorsys.rgb_to_hls(*mc.to_rgb(c))
    l = 1 - amount * (1 - l)
    l = min(1, l)
    l = max(0, l)
    return colorsys.hls_to_rgb(h,l,s)
