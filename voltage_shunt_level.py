# coding=utf-8
import logging
import time

import numpy as np
import pandas as pd
from cycler import cycler
from matplotlib import pyplot as plt, colors, colorbar
from neuron import h

import settings
from baseneuron import BaseNeuron
from morphology import SingleDend
from plotutils import opacity, plot_input_events, FORMAT_FIG, lighten_color
from pynrnutils import get_base_vm_cli, hoc_run
from shunt_level import integrate
from synapses import get_synapse_type

logger = logging.getLogger('voltage_shunt_level')


def voltage_shunt_level_plot(new_once=True, subplots=True):
    if new_once and not subplots:
        fig = plt.figure()
        # noinspection PyTypeChecker
        axes = np.ndarray(5, dtype=plt.Axes)
        size = (2, 10)
        axes[0] = plt.subplot2grid(size, (0, 0), colspan=9)
        axes[3] = axes[0].twinx()
        axes[4] = plt.subplot2grid(size, (0, 9))
        axes[1] = plt.subplot2grid(size, (1, 0), colspan=10)
        axes[2] = axes[1].twinx()  # add second axis to this axis
        # Set ax's patch invisible
        axes[1].patch.set_visible(False)
        axes[2].patch.set_visible(False)
        # place axes[1] in front of axes[2]
        axes[1].set_zorder(axes[2].get_zorder() + 1)
        return fig, axes
    elif new_once and subplots:
        fig, axes = plt.subplots(2, 1)
        axes = np.append(axes, axes[1].twinx())  # add second axis to this axis
        # Set ax's patch invisible
        axes[1].patch.set_visible(False)
        axes[2].patch.set_visible(False)
        # place axes[1] in front of axes[2]
        axes[1].set_zorder(axes[2].get_zorder() + 1)
        # note that axes[3] refers to the second legend for axes[0]
        axes = np.append(axes, axes[0].twinx())
        return fig, axes
    else:
        return None


def voltage_shunt_level(plot=None,  # Tuple of (fig,axes)
                        plot_df_sl_atten=False, # plot shunt level attenutation
                        v_init=-65,  # mV
                        cli=5.0,  # mM
                        t=200,  # ms
                        tm=0,  # integration window (ms)
                        g_i=0.001,  # uS
                        e_offset=-0.0,  # v_init + e_offset
                        g_e=0.001,
                        inhib_syn_type='GABAa',
                        inh_netstim_args=None,
                        exc_syn_type='ampa_Tian',
                        exc_netstim_args=None,
                        hotspot_n_loc_insert=None,
                        inhib_n_loc_insert=None,
                        show_input=0,  # 1 for just one 2 for both
                        v_sample_show_section=False,
                        neuron=None,
                        neuron_args=None,
                        colormap=None,
                        iter_label=None,
                        x_units = settings.um): # settings.um or 'X'
    if neuron is None:
        neuron = SingleDend(L=707, diam=1, nseg=273)
    if neuron_args is None:
        neuron_args=dict()
        
    if v_init == np.nan or v_init is None or v_init == 'auto' or cli == 'auto':
        v_init, cli = get_base_vm_cli(neuron,**neuron_args)

    if isinstance(neuron, BaseNeuron):
        neuron = neuron
    else:
        neuron = neuron(**neuron_args)
        
    if inhib_n_loc_insert is None:
        inhib_n_loc_insert = [(neuron.dend, [0.8])]
    elif type(inhib_n_loc_insert[0]) is not tuple:
        inhib_n_loc_insert = [(neuron.dend, inhib_n_loc_insert)]
    if hotspot_n_loc_insert is None:
        hotspot_n_loc_insert = [(neuron.dend, [0.5])]
    elif type(hotspot_n_loc_insert[0]) is not tuple:
        hotspot_n_loc_insert = [(neuron.dend, hotspot_n_loc_insert)]
    
    if inh_netstim_args is None:
        inh_netstim_args = dict(hz=10,
                                start=20,
                                noise=0.5,
                                duration=t,
                                weight=1,  # weight is number of channels activated at the synapse
                                n_sources=False,
                                own_stream=False)
    if exc_netstim_args is None:
        exc_netstim_args = dict(hz=10,
                                start=20,
                                noise=0.5,
                                duration=t,
                                weight=1,  # weight is number of channels activated at the synapse
                                n_sources=False,
                                own_stream=False)
        
    assert x_units == settings.um or x_units == 'X'

    # sim params
    if tm <= 0 or tm > t:
        tm = t  # time integration window
    T = t / tm  # normalised integration window

    std = 0.1

    # add inhibitory synapses
    logger.debug("inhib_n_loc_insert = {}".format(inhib_n_loc_insert))
    inhib_synapses = []
    inhib_synapse_locations = []
    inhib_n_loc_insert_actual_short = []
    if inhib_syn_type is not None:
        syn_args = get_synapse_type(inhib_syn_type, gmax=0, e=v_init + e_offset, std_i=0)
        for sections, n_or_locations in inhib_n_loc_insert:
            new_inhsyn = neuron.add_synapses(inhib_syn_type,
                                             sections=sections,
                                             locations=n_or_locations,
                                             **syn_args)
            inhib_synapses += new_inhsyn
        for i, syn in enumerate(inhib_synapses):
            x = float("{:.5f}".format(syn['object'].get_loc()))
            h.pop_section()
            sec_name = syn['sec'].hname() if not hasattr(syn['sec'], 'name') else syn['sec'].name
            inhib_n_loc_insert_actual_short.append((sec_name, x))
            inhib_synapse_locations.append(x)
        logger.debug("actual locations: {}".format(inhib_n_loc_insert_actual_short))

    # add excitatory synapses
    exc_synapses = []
    hotspot_synapse_locations = []
    hotspot_n_loc_insert_actual_short = []
    if exc_syn_type is not None:
        logger.debug("hotspot_n_loc_insert = {}".format(hotspot_n_loc_insert))
        exc_syn_args = get_synapse_type(exc_syn_type, gmax=g_e, e=0)
        for sections, n_or_locations in hotspot_n_loc_insert:
            new_excsyn = neuron.add_synapses(exc_syn_type,
                                             sections=sections,
                                             locations=n_or_locations,
                                             **exc_syn_args)
            exc_synapses += new_excsyn
        # get locations where the synapses were actually placed (as it depends on nseg)
        for i, syn in enumerate(exc_synapses):
            try:
                x = float("{:.5f}".format(syn['object'].get_loc()))
                h.pop_section()
            except RuntimeError as re:
                raise
            sec_name = syn['sec'].hname() if not hasattr(syn['sec'], 'name') else syn['sec'].name
            hotspot_n_loc_insert_actual_short.append((sec_name, x))
            hotspot_synapse_locations.append(x)
        logger.debug("actual locations: {}".format(hotspot_n_loc_insert_actual_short))

    # set random seeds for repeatable simulations (netstim method uses settings.RANDOM_SEED)
    for synapse in neuron.synapses:
        try:
            synapse['object'].new_seed(settings.RANDOM_SEED)
        except AttributeError:
            pass
    if exc_syn_type is not None and 'exfluct' not in exc_syn_type:
        neuron.netstim(synapses=exc_synapses, **exc_netstim_args)
    if inhib_syn_type is not None and 'inhfluct' not in inhib_syn_type:
        neuron.netstim(synapses=inhib_synapses, **inh_netstim_args)

    # record input events
    netstims = neuron.netstims
    input_events = {}
    netcons = h.List()
    dummy_nmda = h.nmda_Tian(0.5)  # // used as target in scheme to record spikes
    for key, stim in netstims.items():
        event_vec = h.Vector()
        # // last arg 0 keeps sim unchanged
        # // source, target, thresh, delay, weight
        dummy_netcon = h.NetCon(stim,
                                dummy_nmda, 0, 0, 0)
        dummy_netcon.record(event_vec)
        input_events[key] = event_vec
        netcons.append(dummy_netcon)

    # record voltage everywhere on the dendrites of the neuron
    neurons = [neuron]
    locations = {sec: 'all' for sec in neuron.sections}
    # locations[simple_neuron.soma] = [0.5]
    record_args = [
        {
            "record_var": "v",
            "locations":  locations
        }]

    if settings.NEURON_GUI:

        g = h.Graph(0)
        # g.size(0,5,-80,40)
        #        label      var    col brush    section
        g.addvar("simple_neuron.dend[0].v(.5)", 'v(0.5)', 2, 2, sec=neuron.dend[0])
        g.addvar("simple_neuron.dend[0].v(.1)", 'v(0.1)', 2, 3, sec=neuron.dend[0])
        g.addvar("simple_neuron.dend[0].v(.9)", 'v(0.9)', 2, 4, sec=neuron.dend[0])
        g.addvar("some.v(.5)", 'v(0.5)', 3, 2, sec=neuron.soma)

        # |     Plot view   |       Window size     |
        # x1, y1, dx, dy    | x0, y0, dx, dy
        g.view(0, -80, t, 100, 400, 400, 3000, 1000)
        g.exec_menu("Keep Lines")
        if g not in h.graphList[0]:
            h.graphList[0].append(g)

    neuron.set_cl(cli)
    hoc_run(v_init=v_init, tstop=t, quick=False, record_from=neurons, record_args=record_args)
    logger.debug("creating Vd dataframe")
    logger.debug("\t converting recordings from HocObjects to numpy arrays")
    start = time.time()
    df = neuron.convert_recordings()
    logger.debug("conversion took {:.2f}s".format(time.time() - start))
    df_v = df['v']

    logger.debug("integrating Vd with t={} tm={}".format(t, tm))
    start = time.time()
    integral_v_d = integrate(df_v, window=tm)
    logger.debug("integrating took {:.2f}s".format(time.time() - start))

    logger.info("activating inhibitory conductances with g={}".format(g_i))
    # activate conductances of inhibitory synapse(s)
    # initialise to prevent issues with parameter assignment mid-way through simulations
    h.finitialize()
    if inhib_syn_type is not None and 'inhfluct' in inhib_syn_type:
        params = {'g_i0': g_i, 'std_i': g_i * std}
    else:
        params = {'gmax': g_i}
    for inh_syn in inhib_synapses:
        for key, value in params.items():
            try:
                inh_syn['object'].__setattr__(key, value)
            except LookupError:
                if key == 'gmax':
                    # try alternative gmax param name
                    setattr(inh_syn['object'], 'g_max', value)

    # reset seeds
    for synapse in neuron.synapses:
        try:
            synapse['object'].new_seed(settings.RANDOM_SEED)
        except AttributeError:
            pass
    for key, stim in netstims.items():
        stim.seed(settings.RANDOM_SEED)

    hoc_run(v_init=v_init, tstop=t, quick=False, record_from=neurons, record_args=record_args)
    logger.debug("creating Vd* dataframe")
    start = time.time()
    df_v_star = neuron.convert_recordings()['v']
    logger.debug("conversion took {:.2f}s".format(time.time() - start))
    df_v_star = df_v_star.reindex(df_v_star.columns, axis='columns')

    logger.debug("integrating Vd* with t={} tm={}".format(t, tm))
    start = time.time()
    integral_v_d_star = integrate(df_v_star, window=tm)
    logger.debug("integrating took {:.2f}s".format(time.time() - start))

    logger.info("calculating shunt level")
    df_sl_tail = ((integral_v_d_star - integral_v_d) / integral_v_d_star)
    # re-orientate
    df_sl = df_sl_tail.T.sort_index()

    logger.info("calculating attenuation")
    if df_sl.sum().values[0] == 0.0:
        # no hotspot or no shunting
        df_sl_attenuation = df_sl
    else:
        if len(hotspot_n_loc_insert_actual_short) == 0:
            attenuation_calculation_locations = set(inhib_n_loc_insert_actual_short)
        else:
            # at hotspot
            attenuation_calculation_locations = set(hotspot_n_loc_insert_actual_short)

        column_index = pd.MultiIndex.from_product([attenuation_calculation_locations,
                                                   df_sl.columns.values],
                                                  names=['hotspot_loc', 'tm_value'])
        df_sl_attenuation = pd.DataFrame(columns=column_index,
                                         index=df_sl.index)
        for hotspot_sec_name, hotspot_loc in attenuation_calculation_locations:
            values = df_sl.loc[(hotspot_sec_name, hotspot_loc)].values
            initial_attenuation_calc = np.divide(df_sl, values)
            for tm_value in initial_attenuation_calc:
                step_2_calc = initial_attenuation_calc[tm_value]
                # adjust relative to location
                greater_than_hotspot_idx = step_2_calc > 1.0
                step_3 = step_2_calc
                # DataFrame.update only applies to non-nan values and does so in-place
                # https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.update.html
                step_3.update(1.0 / step_2_calc[greater_than_hotspot_idx])
                # assign
                df_sl_attenuation.loc[:, ((hotspot_sec_name, hotspot_loc), tm_value)] = step_3

    if plot is False:
        return df_sl, df_sl_attenuation
    # graphs
    if plot is None:
        fig, axes = voltage_shunt_level_plot()
    else:
        fig, axes = plot

    # because of the floating point mantissa, need to adjust the locations to the indexed values
    # Note this is deprecated in favour of formatting locations to 5 decimal points
    # for i, loc in enumerate(hotspot_loc):
    #     hotspot_loc[i] = loc + min(abs(df_sl.index.values - loc))
    # for i, loc in enumerate(inhib_loc):
    #     inhib_loc[i] = loc + min(abs(loc - df_sl.index.values))

    # choose samples to plot
    section_sample = inhib_n_loc_insert[0][0]
    if type(section_sample) is list:
        section_sample = section_sample[0]
    if type(section_sample) is str:
        if '[' in section_sample:
            section_sample, index = section_sample.split('[')
            index = int(index[:-1]) + 1  # ignore ']'
        else:
            index = '1'
        section_sample_name = "{}_{}".format(section_sample, index)
        exists = getattr(neuron, section_sample_name, False)
        if not exists:
            section_sample_name = section_sample
    else:
        section_sample_name = section_sample.hname() if not hasattr(section_sample, 'name') else section_sample.name
    cols = df_v[section_sample_name].columns
    mid = int(len(cols) / 2)
    quart = int(mid / 2)
    col_samples = cols[[0, quart, mid, mid + quart, -1]].values
    v_samples = [(section_sample_name, v_sample) for v_sample in col_samples]

    # also include places of input
    col_samples = np.sort(list(set(np.append(col_samples, hotspot_synapse_locations + inhib_synapse_locations))))
    # v_samples = [(section_sample_name, v_sample) for v_sample in col_samples]
    v_samples = sorted(set(v_samples + hotspot_n_loc_insert_actual_short + inhib_n_loc_insert_actual_short),
                       key=lambda tup: (tup[0], tup[1]))
    if colormap is None:
        default_color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    else:
        rgb_colors = plt.get_cmap(colormap, len(v_samples) * 2)
        # convert to hex
        default_color_cycle = [colors.rgb2hex(rgb_colors(i)[:3])
                               for i in range(int(rgb_colors.N / 2), rgb_colors.N)]
        assert len(default_color_cycle) == len(v_samples)
    for i, (sample, color) in enumerate(zip(col_samples, default_color_cycle)):
        if sample in hotspot_synapse_locations or sample in inhib_synapse_locations or sample == min(col_samples):
            default_color_cycle[i] = opacity(50, color)
        else:
            default_color_cycle[i] = opacity(50, color)
    color_cycle = cycler('color', default_color_cycle[:len(v_samples)])
    axes[0].set_prop_cycle(color_cycle)

    # plot
    df_v_to_plot = df_v.loc[:, v_samples]
    axes[0].plot(df_v_to_plot, linestyle='-')
    # plot input events
    if show_input >= 1:
        plot_input_events(axes[0], input_events, cmap='coolwarm', y_offset=0,
                          inhib_syn_type=inhib_syn_type,
                          exc_syn_type=exc_syn_type)

    axes[0].plot(df_v_star.loc[:, v_samples], linestyle='--')
    if show_input >= 2:
        plot_input_events(axes[0], input_events, cmap='coolwarm', y_offset=-25,
                          inhib_syn_type=inhib_syn_type,
                          exc_syn_type=exc_syn_type)

    logger.info("Setting distance measurement from soma")
    # set loc origin to soma
    where = 0
    if 'radial' in neuron.soma.name:
        where = 1
    h.distance(0, where, sec=neuron.soma)

    # skip the 'axon' section
    excluded = 'axon'
    if excluded in df_sl.index.levels[0]:
        indices = df_sl.index.levels[0].difference([excluded])
        # indx = pd.IndexSlice[:, indices.values]
        df_sl = df_sl.loc[indices, :]

    logger.info("distance units for shunt level is {}".format(x_units))
    df_distance = pd.DataFrame(columns=[settings.um, 'X'])
    for sec in neuron.sections:
        sec_name = sec.hname() if not hasattr(sec, 'name') else sec.name
        if sec_name == 'axon':
            continue
        d = h.distance(0.0, sec=sec)
        x = float("{:.5f}".format(0.0))
        df_temp_um = pd.DataFrame({settings.um: d, 'X': x}, index=[(sec_name, x)])
        df_distance = df_distance.append(df_temp_um)
        for seg in sec:
            d = h.distance(seg.x, sec=sec)
            x = float("{:.5f}".format(seg.x))
            df_temp_um = pd.DataFrame({settings.um: d, 'X': x}, index=[(sec_name, x)])
            df_distance = df_distance.append(df_temp_um)
        d = h.distance(1.0, sec=sec)
        x = float("{:.5f}".format(1.0))
        df_temp_um = pd.DataFrame({settings.um: d, 'X': x}, index=[(sec_name, x)])
        df_distance = df_distance.append(df_temp_um)

    # create a MultiIndex
    df_distance = df_distance.reindex(pd.MultiIndex.from_tuples(df_distance.index)).sort_index()
    alt_index = pd.MultiIndex.from_tuples(
            zip(df_distance[x_units].index.droplevel(level=1), df_distance[x_units]),
            names=['compartment_name', x_units])

    if len(hotspot_synapse_locations) > 0:
        hotspot_indices, = np.where(col_samples == hotspot_synapse_locations[-1])
        hotspot_color_index = hotspot_indices[-1]
    else:
        hotspot_color_index = 0

    # SET COLORS
    section_names = df_sl.index.levels[0]
    if colormap is not None:
        rgb_colors = plt.get_cmap(colormap, 2 + len(section_names) * 2)
        # convert to hex
        axes1_color_cycle = [colors.rgb2hex(rgb_colors(i)[:3])
                             for i in range(int(rgb_colors.N / 2), rgb_colors.N)]

    else:
        axes1_color_cycle = default_color_cycle

    grey_rgb_colors_sections = plt.get_cmap("Greys", 2 + len(section_names) * 2)
    grey_color_cycle_sections = [colors.rgb2hex(grey_rgb_colors_sections(i)[:3])
                                 for i in range(int(grey_rgb_colors_sections.N / 2), grey_rgb_colors_sections.N)]

    axes[1].set_prop_cycle(cycler('color', axes1_color_cycle))
    axes[2].set_prop_cycle(cycler('color', axes1_color_cycle))

    logger.info("plotting shunt level")
    # if x_units == 'X':
    for section_name in section_names:
        if section_name == 'axon':
            continue
        axes[1].plot(df_distance.loc[section_name, x_units],
                     df_sl.loc[section_name],
                     linestyle='-',
                     # color=default_color_cycle[hotspot_color_index],
                     label=section_name)
    # else:
    #     df_sl_to_plot = df_sl
    #     df_sl_to_plot.index = alt_index
    #     axes[1].plot(df_sl_to_plot,
    #                  linestyle='-')

    if tm < t:
        for i, tm_value in enumerate(df_sl):
            # y = df_sl.max()[tm_value]
            # x = df_sl.loc[:, tm_value][df_sl.loc[:, tm_value] == y].index.values
            x = df_sl.loc[:, tm_value].index.values[0]
            y = df_sl.loc[:, tm_value].values[0]
            axes[1].text(x, y, "tm={}$_{}$".format(tm, i), rotation=0,
                         ha='center', va='center',
                         bbox=dict(boxstyle="square",
                                   ec=(1, 1, 1, 0.8),
                                   fc=(1, 1, 1, 0.8),
                                   ))
    # DEPRECATED
    # if x_units == settings.um:
    #     # AXES[1] fill in the gaps
    #     fill_gaps(axes[1], linestyle='-', color=opacity(100, axes1_color_cycle[0]))

    logger.info("placing markers")
    if len(inhib_synapse_locations) > 0:
        for inh_section_name, inh_section_loc in inhib_n_loc_insert_actual_short:
            axes[1].plot(df_distance.loc[(inh_section_name, inh_section_loc), x_units],
                         df_sl.loc[(inh_section_name, inh_section_loc)].values,
                         linestyle='none',
                         color=opacity(100, default_color_cycle[0]),
                         marker='v', label=settings.gi)
    if len(hotspot_n_loc_insert_actual_short) > 0:
        for hp_section_name, hp_section_loc in hotspot_n_loc_insert_actual_short:
            axes[1].plot(df_distance.loc[(hp_section_name, hp_section_loc), x_units],
                         df_sl.loc[(hp_section_name, hp_section_loc)].values, linestyle='none',
                         color=opacity(100, default_color_cycle[-1]),
                         marker='o', label='Hotspot')
            dark_color = lighten_color(default_color_cycle[hotspot_color_index], 1.5)
            # plot attenuation center axis
            axes[2].axvline(df_distance.loc[(hp_section_name, hp_section_loc), x_units],
                            linestyle='--',
                            color=dark_color)

    if plot_df_sl_atten:
        logger.info("plotting shunt level attenuation")
        for sec_name in df_sl_attenuation.index.levels[0]:
            if sec_name in df_sl_attenuation.index and sec_name in df_distance.index:
                axes[2].plot(df_distance.loc[sec_name, x_units],
                             df_sl_attenuation.loc[sec_name, :],
                             linestyle=':',
                             # color=default_color_cycle[hotspot_color_index],
                             label='{}.{}'.format(sec_name, df_distance.loc[sec_name, 'X']))
        axes[2].set_ylabel(settings.SHUNT_LEVEL_ATTENUATION + "\nat 'Hotspot'")
        axes[2].set_ylim(auto=True)
        axes[2].set_ylim(bottom=min([0, axes[2].get_ylim()[0], df_sl_attenuation.min().values[0]]),
                     top=max([1.01, axes[2].get_ylim()[1], df_sl_attenuation.max().values[0]]))
    else:
        axes[2].yaxis.set_visible(False)


    # DRAW Sections
    logger.debug("adding sections to plot")
    upper_bound_secs = {}
    for sec_name in df_sl_attenuation.index.levels[0]:
        if sec_name not in df_distance.index:
            continue
        upper_bound_sec = df_distance.loc[sec_name, x_units].max()
        if upper_bound_sec in upper_bound_secs:
            upper_bound_secs[upper_bound_sec] += 1
        else:
            upper_bound_secs[upper_bound_sec] = 1
        if sec_name[-1].isdigit() and upper_bound_secs[upper_bound_sec] > 1:
            # section is same distance away from soma (default), skip plotting
            continue
        axes[2].axvline(upper_bound_sec, ymax=2, lw=0.5,
                        linestyle='-', color=opacity(10, '#000000'))
        axes[2].text(upper_bound_sec, 0.5,
                     "{}".format(sec_name), rotation=90,
                     ha='right', va='center', color=opacity(10, '#000000'),
                     fontsize=settings.SMALLEST_SIZE)
    else:
        # end of for loop
        for sec_loc, num_repeats in upper_bound_secs.items():
            if num_repeats > 1:
                axes[2].text(sec_loc, 0.5,
                             "[ x {}]".format(num_repeats), withdash=True,
                             rotation=90,
                             ha='left', va='center', color=opacity(10, '#000000'),
                             fontsize=settings.SMALLEST_SIZE)
    logger.info("customising figure legends")

    # FIGURE SETTINGS
    vd_legend = axes[0].get_legend()
    if vd_legend is not None:
        labels_Vd = [text.get_text() for text in vd_legend.texts]  # contains Text objects
        labels_Vd = filter(lambda x: not x.startswith('('), labels_Vd)
    else:
        labels_Vd = []
    # remove generated ones from dataframe
    labels_kv = {}
    for prev_label_vd in labels_Vd:
        d_start = prev_label_vd.find("{") + 1
        if prev_label_vd[d_start:d_start + 4] == 'soma':
            prev_label_loc = 0.00
        else:
            prev_label_loc = float(prev_label_vd[d_start:d_start + 4])
        labels_kv[prev_label_loc] = prev_label_vd
    # lines,labels_Vd_star = axes[0].get_legend_handles_labels()
    for i, sample in enumerate(v_samples):
        sec_name = sample[0]
        loc = sample[1]
        extra = ''
        if loc in hotspot_synapse_locations:
            extra = '(hotspot)'
        elif loc in inhib_synapse_locations:
            extra = '(g_i)'
        elif loc == min(col_samples):
            loc = np.nan  # a ".replace('nan','')" hack is required below to remove the 'nan' as .2f requires a float
            extra = 'soma'
        vd_label = settings.Vd.replace('d', "{:.2f}{}".format(loc, extra).replace('nan', ''))
        label = vd_label.replace("V_", '')
        if vd_label not in labels_Vd:
            labels_Vd.append(vd_label)
        if loc == np.nan:
            loc = 0.00
        if loc not in labels_kv or len(label) > labels_kv[loc]:
            labels_kv[loc] = label
        # labels_Vd_star.append(settings.Vdstar.replace('d', "{:.2f}{}".format(loc, extra).replace('nan', '')))

    # created ordered tuple on labels_kv keys
    ordered_loc_labels = sorted(labels_kv.items(), key=lambda (k, v): k)
    ordered_labels = [label for loc, label in ordered_loc_labels]
    ordered_locs = [loc for loc, label in ordered_loc_labels]
    if np.nan in ordered_locs:
        index = ordered_locs.index(np.nan)
        ordered_locs[index] = 0.00

    # legend of voltage trace
    # we now have all labels from all sims
    rgb_colors = plt.get_cmap('Greys', len(labels_Vd) + 1)
    custom_lines = []
    for i, label in enumerate(ordered_labels):
        if 'hotspot' in label or 'g_i' in label:
            color = opacity(50, colors.rgb2hex(rgb_colors(i + 1)[:3]))
        else:
            color = opacity(50, colors.rgb2hex(rgb_colors(i + 1)[:3]))
        from matplotlib.lines import Line2D
        custom_lines.append(Line2D([], [], color=color, linestyle='-'))
    l = axes[0].legend(custom_lines, ordered_labels,
                       bbox_to_anchor=(0., 1.02, 1., .102), loc=4, ncol=4, mode=None,
                       borderaxespad=0., fontsize=settings.SMALLEST_SIZE)
    if len(axes) >= 5:
        # colorbar as legend
        full_grey_cmap = plt.get_cmap('Greys')

        # Get the colormap colors
        my_cmap = full_grey_cmap(np.arange(full_grey_cmap.N))
        # Set alpha
        my_cmap[:, -1] = 1
        # Create new colormap
        my_cmap = colors.ListedColormap(my_cmap)
        # bounds = v_samples
        # norm = colors.BoundaryNorm(v_sample_locs, full_grey_cmap.N)
        cb2 = colorbar.ColorbarBase(axes[4], cmap=my_cmap,
                                    # norm=norm,
                                    ticks=ordered_locs,  # optional
                                    spacing='proportional',
                                    orientation='vertical')
        cb2.ax.tick_params(axis='y', which='minor', colors='blue')
        cb2.set_ticklabels(ordered_labels)
        cb2.set_label("Location, d (X)")

    # plt.colorbar(rgb_colors,ax=axes[0])
    l.set_title(u"{} \u2500  \t {} \u254C".format(settings.Vd, settings.Vdstar),
                prop={
                    'family':  'serif',
                    'stretch': 'normal',
                    'weight':  'bold',
                    'size':    settings.LARGE_SIZE})

    labels = [settings.gi, "'Hotspot'"]
    # for i in df_sl.loc[inhib_locs]:
    # for i in df_sl.loc[hotspot_locs]:
    # skip the solid line (handle[0])
    custom_lines = [Line2D([], [], color=opacity(100, grey_color_cycle_sections[0]),
                           linestyle='none', marker='v'),
                    Line2D([], [], color=opacity(100, grey_color_cycle_sections[-1]),
                           linestyle='--', marker='o')]
    # just put g_i and Hotspot in labels
    axes[1].legend(custom_lines, labels,
                   bbox_to_anchor=(0.0175, 0.95, 1., -0.7), loc=2, ncol=2, mode=None,
                   borderaxespad=0.)

    logger.info("setting xlabels and ylabels")

    axes[0].set_ylabel(" ".join([settings.MEMBRANE_POTENTIAL, settings.UNITS(settings.mV)]))
    axes[1].set_ylabel(settings.SHUNT_LEVEL)

    axes[0].set_xlabel(" ".join([settings.TIME, settings.UNITS(settings.ms)]))
    axes[1].set_xlabel("Location ({})".format(x_units) if x_units != settings.um
                       else "Distance from soma ({})".format(settings.um))
    axes[2].set_xlabel("$d$ ({})".format(axes[2].get_xlabel()))

    axes[0].set_ylim(auto=True)
    axes[1].set_ylim(auto=True)
    axes[1].set_ylim(min(0, df_sl.min().values[0]), max(axes[1].get_ylim()[1] + 0.002, df_sl.max().values[0]))
    axes[0].set_xlim(0, t)
    # axes[1].set_xlim(0, 1 if x_units)
    # axes[2].set_xlim(0, 1)

    logger.info("formatting figure")
    # FORMAT_FIG(ax=axes[1], remove_bottom=True)
    FORMAT_FIG(fig, axes, remove_right=True, adjust_left=False, adjust_top=False)

    # add extra legend on first plot for different iterations of this method
    if iter_label is not None:
        logger.info("extra legend for for voltage")
        if iter_label == 'auto':
            iter_label = "{}{}{}".format(settings.SYNAPSE_MAP[exc_syn_type] if exc_syn_type is not None else '',
                                         ' + ' if (exc_syn_type is not None and inhib_syn_type is not None) else '',
                                         settings.SYNAPSE_MAP[inhib_syn_type] if inhib_syn_type is not None else '')

        # re plot
        axes[3].plot(df_v_to_plot.iloc[:, 0], label=iter_label, color=default_color_cycle[
            hotspot_color_index])
        axes[3].set_xlim(axes[0].get_xlim())
        axes[3].set_ylim(axes[0].get_ylim())
        axes[3].xaxis.set_visible(False)
        axes[3].yaxis.set_visible(False)
        axes[3].set_frame_on(False)
        axes[3].legend(title="Synapse Type",
                       bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=1, mode=None,
                       borderaxespad=0., fontsize=settings.SMALL_SIZE)
        # FORMAT_FIG(ax=axes[3], remove_right=True)

    # if plot is None and iter_label is not None:
    #     fig.suptitle(iter_label)

    # ADD SCALEBARS

    # FORMAT_FIG(ax=axes[0], scalebar=dict(
    #         matchx=True, matchy=True, hidex=True, hidey=True,
    #         loc=5,
    #         labelx=" ".join([settings.TIME, settings.UNITS(settings.ms)]),
    #         labely=" ".join([settings.MEMBRANE_POTENTIAL, settings.UNITS(settings.mV)]),
    # ))
    # FORMAT_FIG(ax=axes[1], scalebar=dict(
    #         matchx=True, matchy=True, hidex=True, hidey=True,
    #
    #         loc=4,
    #         labelx="Location ({})".format(x_units) if x_units != settings.um
    #                 else "Distance from soma ({})".format(settings.um),
    #         labely=settings.SHUNT_LEVEL
    # ))

    return df_sl, df_sl_attenuation, fig, axes
