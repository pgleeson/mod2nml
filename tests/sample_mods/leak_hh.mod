TITLE leak_hh.mod   
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	    (S) = (siemens)
}

? interface
NEURON {
        SUFFIX leak_hh
        NONSPECIFIC_CURRENT il
        RANGE gl, el
	    THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
        gl = .0003 (S/cm2)	<0,1e9>
        el = -54.3 (mV)
}
 
 
ASSIGNED {
        v (mV)
        celsius (degC)
        il (mA/cm2)
}
 
? currents
BREAKPOINT {
        il = gl*(v - el)
}
 
