: file for steady state conductance
NEURON {
	POINT_PROCESS SSG
	RANGE i,g,e,tstart,tstop,active
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	e (mV)
	g (umho)
	tstart = 0
	tstop = 1e100
	active = 0
}
ASSIGNED {
	v (mV)
	i (nA)
}

INITIAL{
	net_send(tstart,1)
	net_send(tstop,0)
}

BREAKPOINT {
	i = g*(active)*(v - e)
}

NET_RECEIVE(w) {
	active = flag
}
