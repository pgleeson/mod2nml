TITLE hh_clean.mod   squid sodium, potassium, and leak channels
 
COMMENT
  potassium, and leakage channels found in the squid giant axon membrane.
  ("A quantitative description of membrane current and its application 
  conduction and excitation in nerve" J.Physiol. (Lond.) 117:500-544 (1952).)
 Membrane voltage is in absolute mV and has been reversed in polarity
  from the original HH convention and shifted to reflect a resting potential
  of -65 mV.
 Remember to set a squid-appropriate temperature
 (e.g. in HOC: "celsius=6.3" or in Python: "h.celsius=6.3").
 See squid.hoc for an example of a simulation using this model.
 SW Jaslove  6 March, 1992
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	    (S) = (siemens)
}

? interface
NEURON {
        SUFFIX hhclean
        USEION na READ ena WRITE ina REPRESENTS CHEBI:29101
        USEION k READ ek WRITE ik REPRESENTS CHEBI:29103
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gl, el, gna, gk
        GLOBAL minf, hinf, ninf, mtau, htau, ntau
	    THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
        gnabar = .12 (S/cm2)	<0,1e9>
        gkbar = .036 (S/cm2)	<0,1e9>
        gl = .0003 (S/cm2)	<0,1e9>
        el = -54.3 (mV)
}
 
STATE {
        m h n
}
 
ASSIGNED {
        v (mV)
        celsius (degC)
        ena (mV)
        ek (mV)

	    gna (S/cm2)
	    gk (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf hinf ninf
	    mtau (ms) htau (ms) ntau (ms)
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h
    	ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
	    ik = gk*(v - ek)      
        il = gl*(v - el)
}
 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
}

? states
DERIVATIVE states {  
        rates(v)
        m' =  (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
}
 


? rates
PROCEDURE rates(v(mV)) {  
        LOCAL  alpha, beta, sum, q10

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

        :"n" potassium activation system
        alpha = .01*alpha_f(-(v+55),10) 
        beta = .125*exp(-(v+65)/80)
	    sum = alpha + beta
        ntau = 1/sum
        ninf = alpha/sum
}

FUNCTION alpha_f(x,y) {  
        alpha_f =  x/(exp(x/y) - 1)
}
 
 
UNITSON
