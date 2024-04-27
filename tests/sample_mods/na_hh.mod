TITLE na_hh.mod   
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	    (S) = (siemens)
}

? interface
NEURON {
        SUFFIX na_hh
        USEION na READ ena WRITE ina REPRESENTS CHEBI:29101
        RANGE gnabar, gna
        GLOBAL minf, hinf,  mtau, htau
	    THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
        gnabar = .12 (S/cm2)	<0,1e9>
}
 
STATE {
        m h 
}
 
ASSIGNED {
        v (mV)
        celsius (degC)
        ena (mV)

	    gna (S/cm2)
        ina (mA/cm2)
        minf hinf
	    mtau (ms) htau (ms)
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h
    	ina = gna*(v - ena)
}
 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
}

? states
DERIVATIVE states {  
        rates(v)
        m' =  (minf-m)/mtau
        h' = (hinf-h)/htau
}
 


? rates
PROCEDURE rates(v(mV)) {  
        LOCAL  alpha, beta, sum

UNITSOFF
        :"m" sodium activation system
        alpha = .1 * alpha_f(-(v+40),10)
        beta =  4 * exp(-(v+65)/18)
        sum = alpha + beta
	    mtau = 1/sum
        minf = alpha/sum

        :"h" sodium inactivation system
        alpha = .07 * exp(-(v+65)/20)
        beta = 1 / (exp(-(v+35)/10) + 1)
        sum = alpha + beta
	    htau = 1/sum
        hinf = alpha/sum
}

FUNCTION alpha_f(x,y) {  
        alpha_f =  x/(exp(x/y) - 1)
}
 
 
UNITSON
