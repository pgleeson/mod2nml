TITLE k_hh.mod   
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	    (S) = (siemens)
}

? interface
NEURON {
        SUFFIX k_hhab
        USEION k READ ek WRITE ik REPRESENTS CHEBI:29103
        GLOBAL alpha, beta
        RANGE gkbar, gk
	    THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
        gkbar = .036 (S/cm2)	<0,1e9>
}
 
STATE {
        n
}
 
ASSIGNED {
        v (mV)
        celsius (degC)
        ek (mV)

	    gk (S/cm2)
        ik (mA/cm2)
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gk = gkbar*n*n*n*n
	    ik = gk*(v - ek)      
}
 
INITIAL {
	rates(v)
	n = alpha/(alpha+beta)
}

DERIVATIVE states {  
        rates(v)
        n' = (1-n)*alpha - n*beta
}
 

PROCEDURE rates(v(mV)) {  
        LOCAL x

        x = -(v+55)
        alpha = .01*x/(exp(x/10) - 1)
        beta = .125*exp(-(v+65)/80)
}

 
