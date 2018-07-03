# coding=utf-8
from __future__ import division

import logging

from nrnutils import Section, Mechanism, DISTAL, PROXIMAL

from baseneuron import BaseNeuron

logger = logging.getLogger('morphology')


# noinspection PyAttributeOutsideInit,PyPep8Naming
class SingleDend(BaseNeuron):
    def __init__(self, name="SingleDend", L=707, diam=1, nseg=81, geom_nseg=False, *args, **kwargs):
        # values to be used for dendrites
        self.L = L
        self.diam = diam
        self.nseg = nseg
        self.dend = []
        # don't auto-calculate nseg
        super(SingleDend, self).__init__(name, geom_nseg=geom_nseg, *args, **kwargs)
        # reset soma area: set L and diam to 2.8012206 , so that area = 24.651565 (from rho = 0.1)
        # self.soma.L = 2.8012206
        # self.soma.diam = 2.8012206

    def create_mechanisms(self,**kwargs):
        Rm = 20000  # Ohm cm^2 membrane resistance
        self.pas = Mechanism('pas', e=-65, g=(1.0 / Rm))
        assert self.pas.parameters['g'] > 0, "division not working correctly"
        self.mechanisms.append(self.pas)

    def build_sections(self, Ra=100, **kwargs):
        super(SingleDend, self).build_sections(Ra=Ra, **kwargs)
        self.dend.append(Section(L=self.L, diam=self.diam, nseg=self.nseg, Ra=self.Ra, cm=self.cm,
                                 mechanisms=[self.pas],
                                 parent=self.soma))
        self.dend[0].name = 'dend_1'
        self.sections += self.dend

    # def set_ions(self, **kwargs):
    #     # do nothing
    #     pass

