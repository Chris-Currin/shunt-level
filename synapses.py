SYNAPSE_TYPE = {
    'inhfluct':               dict(
            pcl=0.8,  # : permeability fraction of Cl
            phco3=0.2,  # : permeability fraction of HCO3
            g_i0=0.0573,  # (umho)	: average inhibitory conductance
            std_i=0.0066,  # (umho)	: standard dev of inhibitory conductance
            tau_i=10.49,  # (ms)	: time constant of inhibitory conductance
            Tdur=1.0,  # (ms)	: transmitter duration
    ), 'inhfluct_simple':     dict(
            e_fluct=-65,  # (mV) 	: reversal potential
            g_i0=0.0573,  # (umho)	: average inhibitory conductance
            std_i=0.0066,  # (umho)	: standard dev of inhibitory conductance
            tau_i=10.49,  # (ms)	: time constant of inhibitory conductance
            Tdur=1.0,  # (ms)	: transmitter duration
    ), 'exfluct':             dict(
            g_e0=0.0121,  # (umho)	: average excitatory conductance
            std_e=0.0030,  # (umho)	: standard dev of excitatory conductance
            tau_e=2.728,  # (ms)	: time constant of excitatory conductance
            Tdur=1.0,  # (ms)	: transmitter duration
    ), 'exfluct_simple':      dict(
            e_fluct=0,  # (mV)    : reversal potential
            g_e0=0.0121,  # (umho)	: average excitatory conductance
            std_e=0.0030,  # (umho)	: standard dev of excitatory conductance
            tau_e=2.728,  # (ms)	: time constant of excitatory conductance
            Tdur=1.0,  # (ms)	: transmitter duration
    ), 'GABAa':               dict(
            pcl=0.8,  # : permeability fraction of Cl
            phco3=0.2,  # : permeability fraction of HCO3
            Tdur=1.0,  # (ms)	: transmitter duration (rising phase)
            Alpha=5.0,  # (/ms)	: forward (binding) rate
            Beta=0.18,  # (/ms)	: backward (unbinding) rate
            gmax=0.00005,  # (uS)	: max conductance for 1 channel
    ), 'ProbAMPANMDA2_RATIO': dict(
            tau_r_AMPA=0.2,  # (ms)  : dual-exponential conductance profile
            tau_d_AMPA=1.7,  # (ms)  : IMPORTANT: tau_r < tau_d
            tau_r_NMDA=2.04,  # (ms) : dual-exponential conductance profile
            tau_d_NMDA=75.2,  # (ms) : IMPORTANT: tau_r < tau_d
            Use=0.0,  # default 1
            # (1)   : Utilization of synaptic efficacy (just initial values! Use, Dep and Fac are overwritten by
            # BlueBuilder assigned values)
            Dep=0,  # (ms) default 100 : relaxation time constant from depression
            Fac=0,  # (ms) default 10 :  relaxation time constant from facilitation
            e=0,  # (mV)  : AMPA and NMDA reversal potential
            mgVoltageCoeff=0.08,  # (mV) : the coefficient for the voltage dependence of the equation
            gmax=0.00005,  # (uS) : weight conversion factor (from nS to uS)
    ), 'NMDA_Poirazi':        dict(
            Alpha=4,  # (/ms /mM)	: forward (binding) rate (4)
            Beta=0.01,  # (/ms)
            e=0,  # (mV)		: reversal potential
            gmax=0,
    ), 'nmda_Tian':           dict(
            tau1=3.18,  # (ms) <1e-9,1e9>     :rise time constant
            tau2=57.14,  # (ms) <1e-9,1e9>	:decay time constant
            tau3=2000,  # (ms) <1e-9,1e9>	    :decay time constant
            g_max=0.000045,  # (uS)			: single channel conductance
            e=0,  # (mV)
    ), 'ampa_Tian':           dict(
            tau1=0.5,  # (ms) <1e-9,1e9>   : rise time constant
            tau2=2.6,  # (ms) <1e-9,1e9>	: decay time constant
            g_max=0.000010,  # (uS) <1e-9,1e9>	: single channel conductance
            e=0,  # (mV)
    )}


def get_synapse_type(syn_type, gmax=0, e=0, **kwargs):
    obj = SYNAPSE_TYPE[syn_type]
    for key in obj:
        if key.startswith('g'):
            obj[key] = gmax
        elif key.startswith('e'):
            obj[key] = e
    for keyword in kwargs:
        for key in obj:
            if key == keyword:
                obj[key] = kwargs[keyword]
    return obj
