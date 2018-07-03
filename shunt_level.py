# coding=utf-8
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from neuron import h

import settings
from baseneuron import logger
from morphology import SingleDend
from plotutils import annotate_cols_rows, FORMAT_FIG
from pynrnutils import hoc_run

zz = h.Impedance()


def calc_z(sec, where, FREQ=0):
    logger.debug("calculating impedance in sec {} at {} ".format(sec.hname(), where))
    h.distance(0, where, sec=sec)
    zz.loc(where, sec=sec)
    zz.compute(FREQ, 1, sec=sec)
    # logger.debug("{} \t {} \t {} \t {} \t {} \t {} \t {}".format("x",
    #                                                              "distance",
    #                                                              "input",
    #                                                              "input_phase",
    #                                                              "transfer",
    #                                                              "transfer_phase",
    #                                                              "ratio"))
    # for seg in sec:
    #     x = seg.x
    #     logger.debug("{} \t {} \t {} \t {} \t {} \t {} \t {}".format(x,
    #                                                                  h.distance(x, sec=sec),
    #                                                                  zz.input(x, sec=sec),
    #                                                                  zz.input_phase(x, sec=sec),
    #                                                                  zz.transfer(x, sec=sec),
    #                                                                  zz.transfer_phase(x, sec=sec),
    #                                                                  zz.ratio(x, sec=sec)))


# noinspection PyPep8Naming
def calc_shunt_level_steady(synapses=None, sections=None, g='g'):
    """ Shunt Level Calculation from
    Gidon, A., & Segev, I. (2012). Principles Governing the Operation of Synaptic Inhibition in Dendrites. Neuron,
        75(2), 330–341. https://doi.org/10.1016/j.neuron.2012.05.015

    Method calls 'calcZ(WHERE)' method in 'usefulFns.hoc' to calculate impedance (after accessing relevant section)
    https://www.neuron.yale.edu/neuron/static/new_doc/analysis/programmatic/impedance.html

    for each synapse at location i
        Set origin for calculations to synapse location i (this is a specific synapse on a section)
        for each location d along sections
            calculate voltage attenuation (Impedance.ratio(d)) from i to d
            calculate input resistance (Impedance.input(d)) at d
        for each location d along sections
            Set origin for calculations to location d
            calculate voltage attenuation (Impedance.ratio(i)) from d to i
            calculate input resistance (Impedance.input(i)) at i
        calculate shunt level for synapse at location i

    Symbols according to table in paper
    X,(Xi)  Electrotonic distance (in units of the space constant, l) from origin (to location i); (dimensionless).
    L 		Electrotonic length (in units of l) of a dendritic branch; (dimensionless).
    V 		Steady membrane potentials, as a deviation from the resting potential; (volt).
    Ri  	Input resistance at location i;(U).
    DRi 	Change in Ri due to synaptic conductance perturbation; (U).
    gi 		Steady synaptic conductance perturbation at location i; (S).
    SL 		Shunt level; (0%SL%1; dimensionless).
    SLi 	Shunt level DRi / Ri due to activation of single or multiple conductance perturbations; (0%SL%1;
            dimensionless).
    Ri,j  	Transfer resistance between location i and location j;(U).
    SLi,j	Attenuation of SL (SLj/ SLi) for a single conductance perturbation at location i;(0%SLi,j%1;
            dimensionless).
    Ai,j	Voltage attenuation, Vj/Vi, for current perturbation at location i;(0%Ai,j%1; dimensionless).
    p 		Dendritic-to-somatic conductance ratio; (G dendrite/G soma; dimensionless).
    RN  	Input resistance at X = 0 for a semi-infinite cable; (U).
    B 		Cable boundary condition; (G dendrite/GN; dimensionless)

    In equations, i = input location and d = other location
    DRd = Rd - Rd* = (gi x (Ri,d)^2)/(1+ gi x Ri)       Equation 4
    Ri,d = Rd,i = Ri x Ai,d = Rd x Ad,i                 Equation 5
    SLd = DRd / Rd = [(gi x Ri)/(1 + gi x Ri)] x Ai,d x Ad,i
    SLi,d = SLd,i = Ai,d x Ad,i
    :rtype: pd.DataFrame

    """
    # create DataFrame
    labels = []  # list of tuples (synapse label, secname, loc) where loc in range [0:1]
    for syn in synapses:
        for sec in sections:
            labels += [(syn['label'], sec.hname(), seg.x) for seg in sec]
    index = pd.MultiIndex.from_tuples(labels, names=['synapse', 'section', 'loc'])
    Adi = settings.Adi
    Aid = settings.Aid
    Ri = settings.Ri
    Rd = settings.Rd
    gi = settings.gi
    um = settings.um
    columns = [um, "distance(id)", "distance(di)", Aid, Adi, Ri, Rd, gi]
    df = pd.DataFrame(columns=columns, index=index)

    logger.debug(
        "{:2s} {:13s} {:2s} {:10s} {:2s} {:6s} {}".format("i", "synapse_label", "x", "section.hname()", "d",
                                                          "seg.x", "section.x"))
    for syn_index, syn in enumerate(synapses):
        if 'object' not in syn or syn['object'] is None:
            # steady state
            g_i = syn['g']
            i = syn['loc']
            synapse_label = syn['label']
        else:
            # transient synapses
            g_i = getattr(syn['object'], g)
            i = syn['object'].get_loc()
            synapse_label = syn['label']

        sec_ref = h.SectionRef(sec=syn['sec'])
        # h.calcZ(sec_ref, i)
        # for sec_index, sec in enumerate(sections):
        #     for seg_index, seg in enumerate(sec):
        #         logger.debug(
        #             "{:2d} {:13s} {:2d} {:10s} {:2d} {:.5f} {:3.2f}um".format(syn_index, synapse_label, sec_index,
        #                                                                       sec.hname(), seg_index, seg.x,
        #                                                                       seg.x * sec.L))
        #         d = seg.x
        #         df.loc[(synapse_label, sec.hname(), d),
        #                     ("um", "$g_i$")] = seg.x * sec.L, g_i
        #         distance_id_h = h.distance(d, sec=sec)
        #         # retrieve impedance results
        #         attenuation_id_h = h.zz.ratio(d, sec=sec)
        #         input_resistance_d_h = h.zz.input(d, sec=sec)
        #         df.loc[(synapse_label, sec.hname(), d),  # assign column values (below) to row (this line)
        #                     ("distance(id)", "$A_{i,d}$", "$R_d$")] = distance_id_h, attenuation_id_h, input_resistance_d_h
        calc_z(syn['sec'], i)
        for sec_index, sec in enumerate(sections):
            for seg_index, seg in enumerate(sec):
                logger.debug(
                    "{:2d} {:13s} {:2d} {:10s} {:2d} {:.5f} {:3.2f}um".format(syn_index, synapse_label, sec_index,
                                                                              sec.hname(), seg_index, seg.x,
                                                                              seg.x * sec.L))
                d = seg.x
                df.loc[(synapse_label, sec.hname(), d),
                       (um, gi)] = seg.x * sec.L, g_i
                # distance_id_h, attenuation_id_h, input_resistance_d_h = df.loc[(synapse_label, sec.hname(), d),
                #                        ("distance(id)", "$A_{i,d}$", "$R_d$")]

                distance_id = h.distance(d, sec=sec)
                # retrieve impedance results
                attenuation_id = zz.ratio(d, sec=sec)
                input_resistance_d = zz.input(d, sec=sec) * 1e6  # convert from MOhm to Ohm

                # assert distance_id_h == distance_id, "{} != {}".format(distance_id_h,distance_id)
                # assert attenuation_id_h == attenuation_id, "{} != {}".format(attenuation_id_h,attenuation_id)
                # assert input_resistance_d_h == input_resistance_d, "{} != {}".format(input_resistance_d_h,input_resistance_d)
                df.loc[(synapse_label, sec.hname(), d),  # assign column values (below) to row (this line)
                       ("distance(id)", Aid, Rd)] = distance_id, attenuation_id, input_resistance_d

        for sec_index, sec in enumerate(sections):
            sec_ref = h.SectionRef(sec=sec)
            L = h.lambda_f(h.freq)
            logger.debug("Space constant L = {}".format(L))
            for seg_index, seg in enumerate(sec):
                logger.debug(
                    "{:2d} {:13s} {:2d} {:10s} {:2d} {:.5f} {:3.2f}um".format(syn_index, synapse_label, sec_index,
                                                                              sec.hname(), seg_index, seg.x,
                                                                              seg.x * sec.L))
                d = seg.x
                # new impedance calculation
                # h.calcZ(sec_ref, d)
                # distance_di_h = h.distance(i, sec=syn['sec'])
                # # retrieve impedance results
                # attenuation_di_h = h.zz.ratio(i, sec=syn['sec'])
                # input_resistance_i_h = h.zz.input(i, sec=syn['sec'])
                calc_z(sec, d)
                distance_di = h.distance(i, sec=syn['sec'])
                # retrieve impedance results
                attenuation_di = zz.ratio(i, sec=syn['sec'])
                input_resistance_i = zz.input(i, sec=syn['sec']) * 1e6  # convert from MOhm to Ohm
                df.loc[(synapse_label, sec.hname(), d),
                       ("distance(di)", Adi, Ri)] = distance_di, attenuation_di, input_resistance_i
                # assert distance_di_h == distance_di, "{} != {}".format(distance_di_h,distance_di)
                # assert round(attenuation_di_h,5) == round(attenuation_di,5), "{} != {}".format(attenuation_di_h,attenuation_di)
                # assert input_resistance_i_h == input_resistance_i, "{} != {}".format(input_resistance_i_h,input_resistance_i)
                # df.loc[(synapse_label, sec.hname(), d),
                #             ("distance(di)", "$A_{d,i}$", "$R_i$")] = distance_di_h, attenuation_di_h, input_resistance_i_h
            # does a vector/matrix operation so SL is calculated for every d
            R_i, A_id, A_di, R_d = df.loc[
                (synapse_label, sec.hname()), (Ri, Aid, Adi, Rd)].values.T
            df.loc[(synapse_label, sec.hname()), settings.Adi + " alt"] = R_i * A_id / R_d
            SLd = settings.SLd
            df.loc[(synapse_label, sec.hname()), SLd] = ((g_i * R_i) / (1 + g_i * R_i)) * A_id * A_di
    return df


def integrate(df, how='trapz', window=20, rolling=False):
    """Numerically integrate the time series.

    :type df: pd.DataFrame
    :param df:
    :param how: the method to use (trapz by default)
    :param window: the integration window, tm, in ms (20 by default)

    :return

    Available methods:
     * trapz - trapezoidal
     * cumtrapz - cumulative trapezoidal
     * simps - Simpson's rule
     * romb - Romberger's rule

    See http://docs.scipy.org/doc/scipy/reference/integrate.html for the method details.
    or the source code
    https://github.com/scipy/scipy/blob/master/scipy/integrate/quadrature.py
    """
    from scipy import integrate

    available_rules = set(['trapz', 'cumtrapz', 'simps', 'romb'])
    if how in available_rules:
        rule = integrate.__getattribute__(how)
    else:
        print('Unsupported integration rule: %s' % (how))
        print('Expecting one of these sample-based integration rules: %s' % (str(list(available_rules))))
        raise AttributeError

    # df = df-(-80.0)
    # df = df - df.iloc[0,0]
    if rolling:
        rolling_window = df.rolling(window=int(window / h.dt))
        result = rolling_window.apply(rule)  # integrate along the index (time)
    else:
        # shift dataframe to only have points at edge of window (leftmost will always be 0
        t_points = df[::int(window / h.dt)].index.values
        if len(t_points) <= 1:
            t_points = df[::(len(df) - 1)].index.values
            logger.info("window was too big, changing window from {} to [}".format(window, (len(df) - 1) * h.dt))
        result = pd.DataFrame(columns=df.columns, index=t_points)
        for t in range(1, len(t_points)):
            result.loc[t_points[t]] = df[t_points[t - 1]:t_points[t]].apply(rule)
    return result.apply(np.abs)[1::]    # skip the first series (when t=0 and all values are nan)


def create_fig_shunt_level_analytical_vs_synapses():
    # neuron with no inputs
    simple_neuron_steady = SingleDend(L=707, diam=1, nseg=81)
    # neuron to have inputs
    simple_neuron_before = SingleDend(L=707, diam=1, nseg=81)
    simple_neuron_after = SingleDend(L=707, diam=1, nseg=81)

    g_i = 0.001 * 1e-6  # convert from uS to S
    g_e = g_i * 1e6  # convert to uS for mod file
    std = 0.1
    loc = [0.6]
    # add 1 excitatory conductances
    exfluct_simple_synapses_before = simple_neuron_before.add_synapses('exfluct_simple',
                                                                       sections=simple_neuron_before.dend,
                                                                       locations=loc, g_e0=0, std_e=0)
    exfluct_simple_synapses_after = simple_neuron_after.add_synapses('exfluct_simple',
                                                                     sections=simple_neuron_after.dend,
                                                                     locations=loc, g_e0=g_e, std_e=g_e * std)
    exfluct_simple_synapses_after[0]['label'] = "exfluct_simple_after"
    steady_synapses = [
        {
            'object': None, 'label': 'steady', 'g': g_i, 'sec': simple_neuron_steady.dend[0],
            'loc':    syn['object'].get_loc()} for
        syn in exfluct_simple_synapses_before]
    hoc_run(v_init=-65, tstop=0)
    df_steady = calc_shunt_level_steady(synapses=steady_synapses,
                                        sections=simple_neuron_steady.dend, g=g_i)
    df_active_before = calc_shunt_level_steady(synapses=exfluct_simple_synapses_before,
                                               sections=simple_neuron_before.dend, g='g_e0')
    df_active_after = calc_shunt_level_steady(synapses=exfluct_simple_synapses_after,
                                              sections=simple_neuron_after.dend, g='g_e0')
    # exfluct_simple_synapses_before = simple_neuron_before.add_synapses('ampa_Tian',
    #                                                                    sections=simple_neuron_before.dend,
    #                                                                    locations=loc, g_max=0)
    # exfluct_simple_synapses_after = simple_neuron_after.add_synapses('ampa_Tian',
    #                                                                  sections=simple_neuron_after.dend,
    #                                                                  locations=loc)
    # simple_neuron_after.netstim(synapses=exfluct_simple_synapses_after, **dict(hz=10,
    #                                                                            start=0,
    #                                                                            noise=0,
    #                                                                            duration=10,
    #                                                                            weight=1,
    #                                                                            # weight is number of channels
    #                                                                            # activated at the synapse
    #                                                                            n_sources=False,
    #                                                                            own_stream=False))
    # exfluct_simple_synapses_after[0]['label'] = "ampa_Tian"
    # steady_synapses = [
    #     {
    #         'object': None, 'label': 'steady', 'g': g_i, 'sec': simple_neuron_steady.dend[0],
    #         'loc':    syn['object'].get_loc()} for
    #     syn in exfluct_simple_synapses_before]
    # hoc_run(v_init=-65, tstop=1)
    # df_steady = calc_shunt_level_steady(synapses=steady_synapses,
    #                                     sections=simple_neuron_steady.dend, g=g_i)
    # df_active_before = calc_shunt_level_steady(synapses=exfluct_simple_synapses_before,
    #                                            sections=simple_neuron_before.dend, g='g_max')
    # df_active_after = calc_shunt_level_steady(synapses=exfluct_simple_synapses_after,
    #                                           sections=simple_neuron_after.dend, g='g_max')
    column_groups = [[settings.Aid, settings.Adi], [settings.Ri, settings.Rd], [settings.SLd]]
    neurons = [simple_neuron_steady, simple_neuron_before, simple_neuron_after]
    synapses = [steady_synapses, exfluct_simple_synapses_before, exfluct_simple_synapses_after]
    dataframes = [df_steady, df_active_before, df_active_after]
    from cycler import cycler
    # plotting
    fig, axes = plt.subplots(len(neurons), len(column_groups))
    default_color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    df_all = pd.DataFrame(index=df_active_before.index.levels[2])
    for ax_i, (df, syn_group, neuron) in enumerate(zip(dataframes, synapses, neurons)):
        for syn in syn_group:
            for cg, column_group in enumerate(column_groups):
                color_cycle = cycler('color', default_color_cycle[:len(column_group) * 2])
                df_group = df.loc[(syn['label'], neuron.dend[0].hname()),
                                  column_group]
                if 'steady' in syn['label']:
                    df_group.plot(ax=axes[ax_i, cg], linestyle='-')
                elif 'after' in syn['label']:
                    df_group = df_group.rename({var: var + "*" for var in df_group.columns}, axis='columns')
                    ax_i = 1
                    axes[ax_i, cg].set_prop_cycle(color_cycle)
                    df_group.plot(ax=axes[ax_i, cg], linestyle='--')
                else:
                    axes[ax_i, cg].set_prop_cycle(color_cycle)
                    df_group.plot(ax=axes[ax_i, cg], linestyle='-')

                if 'steady' not in syn['label']:
                    # add columns
                    df_all = pd.concat([df_all, df_group], axis=1)
    # plot change too
    ax_i = 2
    for cg, column_group in enumerate(column_groups):
        # plot the change too
        for var in column_group:
            try:
                df_all.loc[:, "$\Delta$" + var] = (df_all.loc[:, var] - df_all.loc[:, var + "*"])
                df_all.loc[:, "$\Delta$" + var].plot(ax=axes[ax_i, cg], linestyle='-', label="$\Delta$" + var)
            except ZeroDivisionError:
                continue

    var = settings.Rd
    df_sl_rd = df_all.loc[:, "$\Delta$" + var] / df_all.loc[:, var]
    g_i_Ri = g_i * df_all.loc[:, settings.Ri]
    df_sl = g_i_Ri / (1 + g_i_Ri) * df_all[settings.Aid] * df_all[settings.Adi]
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    df_sl.plot(ax=axes[1, 2], label=r"$\frac{{{giri}}}{{{giri_bottom}}}$ {aid} {adi}".format(
            giri=(settings.gi + "  " + settings.Ri).replace("$", ""),
            giri_bottom="1 + " + (settings.gi + "  " + settings.Ri).replace("$", ""),
            aid=settings.Aid, adi=settings.Adi),
               color=colors[1])
    df_sl_rd.plot(ax=axes[1, 2], label=r"$\frac{{\Delta {rd}}}{{{rd}}}$".format(rd=settings.Rd.replace("$", "")),
                  color=colors[2])

    for ax in axes.flatten():
        ax.set_ylim(bottom=0)
        leg = ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4,
                        ncol=2, mode=None, borderaxespad=0., fontsize=settings.SMALLEST_SIZE)
    annotate_cols_rows(axes, cols=None, rows=["Steady", "Mod Synapse", "Changes"])
    FORMAT_FIG(fig, axes, adjust_top=0.95)
