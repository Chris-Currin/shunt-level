# coding=utf-8
from __future__ import print_function, division

import logging

import numpy as np
import pandas as pd
from nrnutils import Mechanism, Section, PROXIMAL

import settings
from plotutils import plot_v, plot_cli
from shared import h, t_vec

logger = logging.getLogger('baseneuron')


# noinspection PyAttributeOutsideInit,PyPep8Naming
class BaseNeuron(object):
    geom_nseg_called = False

    def __init__(self, name="BaseNeuron", call_geom_nseg=True, add_kcc2=False, *args, **kwargs):
        self.name = name
        self.kcc2_inserted = False
        self.mechanisms = []
        self.sections = []
        self.synapses = []
        self._synapses = {
            'all':        self.synapses,
            'exc':        [],
            'inh':        [],
            'exc_labels': [],
            'inh_labels': []}
        self.netcons = {}
        self.netstims = {}
        self.netstims_ranstreams = {}
        self.vec_hoc = {}
        self.vec_np = pd.DataFrame()
        self.apc = None
        # NOTE that the NEURON+PYTHON docs use these steps. Because nrnutils is used, there are changes
        # self.create_sections()
        # self.build_topology()
        # self.build_subsets()
        # self.define_geometry()
        # self.define_biophysics()
        self.create_mechanisms(**kwargs)
        self.soma = None
        self.axon = None
        self.build_sections(**kwargs)
        self.set_ions(**kwargs)

        # logger.debug("P sections: {}".format([h.psection(sec=sec) for sec in self.sections]))
        total_segments = 0
        for sec in self.sections:
            # logger.debug(h.psection(sec=sec))
            total_segments += sec.nseg
        logger.debug("total nseg BEFORE geom_nseg: {}".format(total_segments))
        if call_geom_nseg:
            self.geom_nseg()
            total_segments = 0
            for sec in self.sections:
                # logger.debug(h.psection(sec=sec))
                total_segments += sec.nseg
            logger.debug("total nseg AFTER geom_nseg: {}".format(total_segments))
        logger.debug("Neuron built")

    def geom_nseg(self,
                  freq=100,  # Hz, frequency at which AC length constant will be computed
                  d_lambda=0.1):
        if self.geom_nseg_called:
            logger.debug("'geom_nseg' method should only be called once")
        # h.geom_nseg()
        for sec in self.sections:
            sec.nseg = int((sec.L / (d_lambda * h.lambda_f(freq)) + 0.9) / 2) * 2 + 1
        self.geom_nseg_called = True

    def create_mechanisms(self, g_na_bar=2000, g_k_bar=5, g_im_bar=0.00012,
                          g_pas_k=0.000125, p_na=0.23, p_cl=0.4,
                          **kwargs):
        g_pas_na = g_pas_k * p_na
        g_pas_cl = g_pas_k * p_cl
        # mechanisms
        self.na = Mechanism('na', gbar=g_na_bar)
        self.kv = Mechanism('kv', gbar=g_k_bar)
        self.im = Mechanism('im', gkbar=g_im_bar)
        self.pasghk = Mechanism('pasghk', gclpbar=g_pas_cl, gnapbar=g_pas_na, gkpbar=g_pas_k)
        self.mechanisms.append(self.na)
        self.mechanisms.append(self.kv)
        self.mechanisms.append(self.im)
        self.mechanisms.append(self.pasghk)

    def build_sections(self, Ra=160, cm=1, axon=True, soma=True, *args, **kwargs):
        # membrane properties across the cell
        self.Ra = Ra  # Ohm cm        cytoplasmic resistivity
        self.cm = cm  # uF/cm^2     specific membrane capacitance
        # geometry and topology
        if soma:
            if 'soma_L' in kwargs:
                soma_l = kwargs['soma_L']
            else:
                soma_l = 15
            if 'soma_diam' in kwargs:
                soma_diam = kwargs['soma_diam']
            else:
                soma_diam = 15

            self.soma = Section(L=soma_l, diam=soma_diam, Ra=self.Ra, cm=self.cm,
                                mechanisms=self.mechanisms)
            self.soma.name = 'soma'

            self.sections.append(self.soma)

        if axon:
            self.axon = Section(L=100, diam=0.2, Ra=self.Ra, cm=self.cm,
                                mechanisms=self.mechanisms,
                                parent=self.soma,
                                connection_point=PROXIMAL)
            self.axon.name = 'axon'
            self.sections.append(self.axon)

    def insert_spike_count(self, section=None, location=0.99):
        if section is None:
            section = self.soma
        self.apc = section.add_synapse('apc', 'APCount1', locations=[location])

    @property
    def num_spikes(self):
        if self.apc is None:
            raise Exception("APCount1 not initialised. Call BaseNeuron.insert_spike_count(section,location)")
        return self.apc.n

    def set_ions(self, sections=None,
                 ki=140, ko=5,
                 nai=8, nao=145,
                 cli=5, clo=135,
                 hco3i=12, hco3o=23, **kwargs):
        if sections is None:
            sections = self.sections
        for sec in sections:
            sec.insert("k_ion")
            sec.insert("na_ion")
            sec.insert("cl_ion")
            sec.insert("hco3_ion")
            h.ion_style("hco3_ion", 1, 2, 1, 0, 1)  # last arg if for assigning from global
            for seg in sec:
                seg.ki = ki
                seg.ko = ko
                seg.nai = nai
                seg.nao = nao

        # set globals
        self.set_cl(cli=cli, clo=clo)
        h.hco3i0_hco3_ion = hco3i
        h.hco3o0_hco3_ion = hco3o

    @staticmethod
    def get_g_pas_ions(g_pas_k=0.000125, p_na=0.23, p_cl=0.4):
        return g_pas_k, g_pas_k * p_na, g_pas_k * p_cl

    def add_kcc2(self, **kwargs):
        kcc2 = Mechanism('KCC2', **kwargs)
        for section in self.sections:
            kcc2.insert_into(section)
        self.kcc2_inserted = True
        return self

    def remove_kcc2(self):
        for section in self.sections:
            # https://www.neuron.yale.edu/phpBB/viewtopic.php?t=3335
            h('uninsert KCC2', sec=section)
        self.kcc2_inserted = False

    def set_cl(self, cli=5, clo=None):
        h.cli0_cl_ion = cli
        if clo is not None:
            h.clo0_cl_ion = clo
        for sec in self.sections:
            for seg in sec:
                seg.cli = cli
                if clo is not None:
                    seg.clo = clo
        logger.debug("cli set to {}".format(cli))

    def add_synapses(self, syn_type='AlphaSynapse', sections=None, n=0, locations=None, **kwargs):
        location_is_n = False
        if n == 0 and type(locations) is int:
            # 'locations was' actually 'n'
            n = locations
            locations = None
            location_is_n = True
        assert n > 0 or (locations is not None and len(locations) > 0), "'total' or 'locations' must be specified"
        if sections is None:
            sections = self.sections
        elif type(sections) is str:
            if '[' in sections:
                section, index = sections.split('[')
                index = int(index[:-1])  # ignore ']'
                sections = getattr(self, section)[index]
            else:
                sections = getattr(self, sections)

        if type(sections) is not list:
            sections = [sections]

        if n == 0:
            n = len(locations)
        elif n < 0:
            raise BaseException("'n' incorrectly specified as '{}'".format(n))

        if locations is None:
            locations = np.linspace(0.0, 1.0, int(n / len(sections)) + 2)
            logger.debug("number of synapses: {}".format(n))
            logger.debug("number of branches: {}".format(len(sections)))
            logger.debug("number of synapses per branch: {}".format(int(n / len(sections))))
            # remove first and last locations (which would be 0 and 1)
            locations = locations[1:-1]

        logger.info("placing '{}' on '{}' at locations {}".format(syn_type, [str(self.name) +
                                                                             section.hname() for section in
                                                                             sections],
                                                                  locations))

        if len(locations) == 1 and n > 1 and not location_is_n:
            locations = list(locations) * n  # place n synapses in the same location
        num_unique_locs = len(set(locations))
        all_new_synapses = []
        for i, section in enumerate(sections):
            assert section in self.sections
            if section.nseg < num_unique_locs:
                section.nseg = num_unique_locs
            synapse_label = section.hname() + '_' + syn_type
            logger.debug(synapse_label)
            section.add_synapses(synapse_label, syn_type, locations=locations, **kwargs)
            new_synapses = getattr(section, synapse_label)
            if type(new_synapses) is list:
                new_synapses = [{'label': synapse_label + "_" + str(s), 'object': syn, 'sec': section} for s, syn in
                                enumerate(new_synapses)]
            else:
                new_synapses = [{'label': synapse_label, 'object': new_synapses, 'sec': section}]
            # logger.debug("actual locations: {}".format([syn['object'].get_loc() for syn in new_synapses]))
            self.synapses += new_synapses
            all_new_synapses += new_synapses
            if 'GABA' in syn_type or 'influct' in syn_type:
                self._synapses['inh'] += new_synapses
            else:
                self._synapses['exc'] += new_synapses
        # logger.debug("Synapses:")
        # logger.debug(pformat(self.synapses))

        return all_new_synapses

    def netstim(self, synapses=None, hz=5, start=0, duration=5, noise=0, weight=1, delay=0,
                n_sources=1,
                own_stream=False):
        """
        Stimulate the synapses with NET_RECEIVE blocks
        Uses RandomStream in "ranstream.hoc"



        :param synapses:        list of synapses to receive stimulation
        :param hz:              frequency of stimulation (Hz)
        :param start:           start of stimulation (ms)
        :param duration:        duration of stimulation (ms)
        :param noise:           noise of Poisson distribution [0-1]
        :param weight:          weight of stimulation to synapse (first argument for NET_RECEIVE block)
        :param delay:           delay between netstim activation and the target receiving the input (can be list)
        :param n_sources:       whether a single netstim should stim all synapses (default:1)
                                or each synapses has its own input (0)
        :param own_stream:      give each netstim it's own random stream. Will only have an effect is self.stim is None

        :type delay: List or int

        """
        if synapses is None:
            synapses = self.synapses
        if duration is None:
            duration = h.tstop
        elif duration > h.tstop:
            h.tstop = duration
        if type(delay) is list:
            assert len(delay) == len(synapses)
        netstim_every_n_synapses = int(len(synapses) / n_sources) + 1 if n_sources > 0 else 1
        prev_net_stim_obj = None
        for i, syn in enumerate(synapses):
            if syn['label'] not in self.netstims:
                # create a new NetStim object every n synapses and when i==0
                if i % netstim_every_n_synapses == 0:
                    net_stim_obj = h.NetStim()
                    prev_net_stim_obj = net_stim_obj
                else:
                    net_stim_obj = prev_net_stim_obj

                self.netstims[syn['label']] = net_stim_obj
                if own_stream:
                    # give each netstim it's own random stream
                    ran_stream = h.RandomStream(len(self.netstims))
                    net_stim_obj.noiseFromRandom(ran_stream.r)
                    ran_stream.r.negexp(1)  # must specify negexp distribution with mean = 1
                    ran_stream.start()
                    # store object in memory else hoc's initialisation fails as it can't find the object
                    self.netstims_ranstreams[syn['label']] = ran_stream
                net_stim_obj.seed(settings.RANDOM_SEED)
            else:
                # get existing object if this method is called again
                # Note that the seed is set again to repeat the simulation (change start time to continue the run
                # instead)
                net_stim_obj = self.netstims[syn['label']]
                net_stim_obj.seed(settings.RANDOM_SEED)

            if i % netstim_every_n_synapses == 0:
                net_stim_obj.interval = 1000 / hz
                net_stim_obj.number = hz * (
                        duration / 1000)  # number = freq[1/s] * (duration[ms] / 1000[convert ms to s])
                net_stim_obj.start = start
                net_stim_obj.noise = noise
            # NetCon for connecting NetStim to synapses
            if syn['label'] not in self.netcons:
                # create net connect object. args: source, target, thresh, delay, weight
                net_con_obj = h.NetCon(net_stim_obj, syn['object'], 0, 0, weight)
                self.netcons[syn['label']] = net_con_obj
            else:
                net_con_obj = self.netcons[syn['label']]

            if type(delay) is list:
                net_con_obj.delay = delay[i]
            else:
                net_con_obj.delay = delay
            net_con_obj.weight[0] = weight  # NetCon weight is a vector.

    def record(self, recordings=None):
        """Record variables at various locations of the neurons

        :param recordings: list of variables to record, with their locations
        :return:
        """
        if recordings is None:
            recordings = [
                {
                    "record_var": "v",
                    "locations":  {
                        self.soma: [0.5]
                    }
                }]
        for recording_dict in recordings:
            record_var, locations = recording_dict["record_var"], recording_dict["locations"]
            self.vec_hoc[record_var] = {}
            for compartment, loc in locations.items():
                if type(compartment) is str:
                    for sec in self.sections:
                        if sec.hname() == compartment or (hasattr(sec, 'name') and sec.name == compartment):
                            compartment = sec
                            break
                    else:
                        raise Exception("compartment with name {} not found in sections".format(compartment))
                if compartment not in self.sections:
                    logger.debug("Not a valid compartment {}.{}".format(self.name, compartment.hname()))
                    continue
                    # raise Exception("Not a valid compartment {}".format(compartment.hname()))

                if loc == 'all' or loc is None:
                    loc = [0] + [seg.x for seg in compartment] + [1]
                elif loc == 'middle':
                    loc = [0.5]
                elif type(loc) is not list:
                    loc = [loc]

                if hasattr(compartment, 'name'):
                    compartment_name = compartment.name
                else:
                    compartment_name = compartment.hname()
                self.vec_hoc[record_var][compartment_name] = {}

                for x in loc:
                    x = float("{:.5f}".format(x))
                    vec_ref = self.vec_hoc[record_var][compartment_name][x] = h.Vector()
                    try:
                        vec_ref.record(compartment(x).__getattribute__("_ref_{}".format(record_var)))
                    except NameError as ne:
                        logger.warning(ne)
                        raise ne

        t_vec.record(h.__getattribute__("_ref_t"))

    def convert_recordings(self):
        """
        Convert the recordings from the latest to a pd.DataFrame object for easier selecting of data

        :return: A dataframe with simulation time as the index and var/compartment/seg.x as the columns
        """
        index = np.array(t_vec)
        self.vec_np = pd.DataFrame()
        for record_var, compartments in self.vec_hoc.items():
            for compartment_name, loc in compartments.items():
                # create a new dataframe for this compartment, with time as index and
                #   each column will be in form var/compartment/x
                columns = pd.MultiIndex.from_product([[record_var], [compartment_name], sorted(loc.keys())],
                                                     names=['record_var', 'compartment_name', 'seg.x'])
                df_t = pd.DataFrame(columns=columns, index=index)
                for x, rec in loc.items():
                    # logger.debug("{:5s} {:10s} {:.3f}".format(record_var, compartment_name, x))
                    df_t.loc[:, (record_var, compartment_name, x)] = np.array(rec)
                # sort the columns
                # self.vec_np[compartment_name]=self.vec_np[compartment_name].sort_index(axis=1)
                # add the columns of the new dataframe to the main dataframe
                if self.vec_np.shape[1] > 0:
                    self.vec_np = pd.concat([self.vec_np, df_t], axis=1)
                else:
                    # assign it the first time so MultiIndex is used for the columns
                    self.vec_np = df_t
        return self.vec_np

    def plot(self, ax=None, section=None, location=None):
        assert len(ax) == 2
        plot_v(self.vec_hoc, ax=ax[0], section=section, location=location)
        plot_cli(self.vec_hoc, ax=ax[1], section=section, location=location)
