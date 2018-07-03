# coding=utf-8
import logging
import time

from neuron import h

logger = logging.getLogger('pynrnutils')


def hoc_run(v_init=None, tstop=None, quick=True, scale=False, record_from=None, record_args=None):
    if record_from is not None:
        if type(record_from) is not list:
            record_from = [record_from]
        for recording_object in record_from:
            recording_object.record(record_args)
    if v_init is None:
        v_init = h.v_init
    if tstop is not None:
        h.tstop = tstop
    if quick:
        logger.debug("using cvode")
        h.useCV()
    logger.debug("initialising to {}".format(v_init))
    h.v_init = v_init
    h.finitialize(v_init)
    if h.tstop > 0:
        logger.debug("running NEURON")
        start = time.time()
        h.run()
        logger.debug("running NEURON took {:.2f}s".format(time.time() - start))
    if scale and False:
        logger.debug("scaling all graphs using 'View = plot'")
        # scale graphs after run
        #   h.graphList[0] will be empty if there aren't any graphs
        for graph in h.graphList[0]:
            graph.exec_menu("View = plot")


def get_base_vm_cli(neuron_model=None, default_cli=4, **kwargs):
    """
    Determine steady-state internal chloride concentration by
    1) instantiating a class that extends BaseNeuron (NOTE: class should not have spiking at default values)
    2) adding KCC2 using the add_kcc2() method
    3) setting [Cl]_i to an arbitrary value (can be set in method invocation)
    4) running a fast simulation for a long time
    5) checking if chloride is at steady state (repeat 4 until it is)
    6) return steady state Vm and [Cl]_i

    :param neuron_model: class extending BaseNeuron
    :param default_cli: arbitrary value of [Cl]_i hopefully close to steady state
    :return: steady state Vm and [Cl]_i
    """
    from baseneuron import BaseNeuron
    remove_kcc2 = False
    if neuron_model is None:
        neuron_model = BaseNeuron
    if isinstance(neuron_model, BaseNeuron):
        # is instantiation of class
        base = neuron_model
        if not base.kcc2_inserted:
            base.add_kcc2()
            remove_kcc2 = True
        base.set_cl(default_cli)

    else:
        # is BaseNeuron class (or child class)
        base = neuron_model(**kwargs)
        base.add_kcc2()
        base.set_cl(default_cli)
    h.tstop = 50000
    h.useCV()
    h.finitialize(-65)

    h.run()

    def at_steady_state(continue_dt):
        """
        check if [Cl]_i is at steady state
        :param continue_dt: amount of time to run
        :return: [Cl]_i if at steady state, False otherwise
        """
        v_start = base.soma(.5)._ref_v[0]
        cli_start = base.soma(.5)._ref_cli[0]
        h.continuerun(h.tstop + continue_dt)
        h.tstop += continue_dt
        v_after = base.soma(.5)._ref_v[0]
        cli_after = base.soma(.5)._ref_cli[0]
        if v_after - v_start < 1e-6 and cli_after - cli_start < 1e-6:
            return cli_after
        else:
            return False

    num_steady_state_checks = 0
    while not at_steady_state(1):
        h.continuerun(h.tstop + 10000)
        h.tstop += 10000
        num_steady_state_checks += 1
        if num_steady_state_checks > 10:
            print("not at steady state even after {} ms".format(50000 + num_steady_state_checks * 10000))
            exit(-1)

    h.disableCV()
    vm,cli = base.soma(.5)._ref_v[0], base.soma(.5)._ref_cli[0]
    logger.info("steady state [Cl]_i {}".format(cli))
    logger.info("steady state Vm {}".format(vm))
    logger.info("took {} ms (simulation time)".format(50000 + num_steady_state_checks * 10000))
    if remove_kcc2:
        base.remove_kcc2()
    h.finitialize(vm)
    return vm,cli
