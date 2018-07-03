# coding=utf-8
import platform

# ----------------------------------------------------------------------------------------------------------------------
# VARIABLES (specifically for rendering in plots)
# ----------------------------------------------------------------------------------------------------------------------
# IONS
CL = cl = "Cl$^-$"
CLI = CL_i = cli = r"[{}]$_i$".format(cl)
MILLIMOLAR = mM = "mM"
# CHLORIDE
STATIC_CHLORIDE_STR_ABBR = "Static {}".format(cl)
STATIC_CHLORIDE_STR_LONG = "Static Chloride"
DYNAMIC_CHLORIDE_STR_ABBR = "Dynamic {}".format(cl)
DYNAMIC_CHLORIDE_STR_LONG = "Dynamic Chloride"

# SYNAPSES
GABAA = GABAa = '$GABA_{A}$'
GABAAR = GABAaR = GABAA + 'R'

# map synapses in in synapses.py to nicely formatted names
SYNAPSE_MAP = {
            'inhfluct':'$inhfluct$',
          'inhfluct_simple': '$inhfluct_{simple}$',
          'exfluct': '$exfluct$',
          'exfluct_simple': '$exfluct_{simple}$',
          'GABAa': GABAa,
          'ProbAMPANMDA2_RATIO': '$AMPA + NMDA$',
          'NMDA_Poirazi': '$NMDA_{Poirazi}$',
          'nmda_Tian': '$NMDA_{Tian}$',
          'ampa_Tian': '$AMPA_{Tian}$',
}

# SLi 	Shunt level DRi / Ri due to activation of single or multiple conductance perturbations; (0%SL%1;
#             dimensionless).
SHUNT_LEVEL = SL = SLd = "$SL_{d}$"
# SLi,j	Attenuation of SL (SLj/ SLi) for a single conductance perturbation at location i;(0%SLi,j%1;
#             dimensionless).
SHUNT_LEVEL_ATTENUATION = SLid = "SL attenuation ($SL_{i,d}$)"
# Input resistance at location i;(U).
INPUT_RESISTANCE_i = Ri = "$R_i$"
INPUT_RESISTANCE_d = Rd = "$R_d$"
# Change in Ri due to synaptic conductance perturbation; (U).
CHANGE_IN_INPUT_RESISTANCE_i = DRi = "$\Delta R_i$"
# Voltage attenuation, Vj/Vi, for current perturbation at location i;(0%Ai,j%1; dimensionless).
VOLTAGE_ATTENUATION_i_d = Aid = "$A_{i,d}$"
VOLTAGE_ATTENUATION_d_i = Adi = "$A_{d,i}$"
# Dendritic-to-somatic conductance ratio; (G dendrite/G soma; dimensionless).
RHO = p = "r$\rho$"
# Conductance perturbation at location i
CONDUCTANCE_i = gi = "$g_{i}$"
# Membrane Potential
MEMBRANE_POTENTIAL = "Membrane Potential"
MILLIVOLTS = mV = "mV"
VOLTAGE_SYMBOL = Vm = "$Vm$"
VOLTAGE_d = Vd = "$V_{d}$"
VOLTAGE_d_star = Vdstar = "$V_{d}^*$"
# Distance
MICROMETERS = um = "$\mu m$"
DISTANCE = "Distance"

# Time
TIME = "Time"
MILLISECONDS = ms = 'ms'
SECONDS = s = 's'


def UNITS(text):
    return '(' + text + ')'


def ITALICISE(text):
    return '$' + text + '$'


# ----------------------------------------------------------------------------------------------------------------------
# NEURON
# ----------------------------------------------------------------------------------------------------------------------
NEURON_GUI = False
NEURON_RECOMPILE = False
HOC_PATH = "./"
MOD_PATH = "./"
NRNMECH_PATH = './'
if platform.system() == 'Linux' or platform.system() == 'Darwin':
    NRNMECH_PATH = MOD_PATH + "x86_64/.libs/libnrnmech.so"
elif platform.system() == 'Windows':
    NRNMECH_PATH = MOD_PATH + "nrnmech.dll"
else:
    print("unknown system")
    exit(-1)
NRNMECH_PATH = NRNMECH_PATH.replace("\\", "/")
print(NRNMECH_PATH)
# ----------------------------------------------------------------------------------------------------------------------
# RANDOM
# ----------------------------------------------------------------------------------------------------------------------
RANDOM_SEED = 0
# max # events in a NetStim's stream  (Adjacent streams will be correlated by this offset.)
# // before it begins to repeat values already generated
# // by a different stream.
# // set to 0 and all NetStims will produce identical streams
RANDOM_STREAM_OFFSET = 1000

# ----------------------------------------------------------------------------------------------------------------------
# MATPLOTLIB
# ----------------------------------------------------------------------------------------------------------------------
SMALLEST_SIZE = 6
SMALL_SIZE = 8
MEDIUM_SIZE = 10
LARGE_SIZE = 12
BIGGER_SIZE = 14
BIGGEST_SIZE = 16
FIG_SIZE = (15, 15)
FIG_ADJUST_LEFT = 0.20
FIG_ADJUST_TOP = 0.95
FIG_RES_ADJUST = 1.5  # for 4 for 4K

# Have colormaps separated into categories:
# http://matplotlib.org/examples/color/colormaps_reference.html
cmaps = {
    'Perceptually Uniform Sequential': [
        'viridis', 'plasma', 'inferno', 'magma'],
    'Sequential':                      [
        'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
        'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
        'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn'],
    'Sequential (2)':                  [
        'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
        'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
        'hot', 'afmhot', 'gist_heat', 'copper'],
    'Diverging':                       [
        'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
        'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic'],
    'Qualitative':                     [
        'Pastel1', 'Pastel2', 'Paired', 'Accent',
        'Dark2', 'Set1', 'Set2', 'Set3',
        'tab10', 'tab20', 'tab20b', 'tab20c'],
    'Miscellaneous':                   [
        'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
        'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'hsv',
        'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar']}
