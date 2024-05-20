from collections import OrderedDict
import sympy as sp

import nmodl_nml as m2n

#def test_hhrates():
#    #with open("sample_mods/k_hh_ab.mod") as f:
#    with open("sample_mods/hh_ab.mod") as f:
#        mod = f.read()
#    process_current_law(m2n.parse_mod(mod))

def test_find_current():
    mod = """
    NEURON {
        USEION na READ ena WRITE ina
        RANGE gna
    }
    BREAKPOINT {
        LOCAL drf
        drf = (v - ena)
        ina = gna*drf
    }
    """
    m2n.process_current_law(m2n.parse_mod(mod))

def test_match_gates():
    mod = """
    NEURON {
        USEION na READ ena WRITE ina
        RANGE gna
    }
    STATE { m h }
    BREAKPOINT {
        ina = gna*m*m*m*h*(v - ena)
    }
    """
    m2n.process_current_law(m2n.parse_mod(mod))

def test_find_gates():
    mod = """
    NEURON {
        USEION k READ ek WRITE ik REPRESENTS CHEBI:29103
        RANGE gkbar, gk
        RANGE gna
    }
    PARAMETER { gkbar = .036 }
    STATE { n }
    BREAKPOINT {
        SOLVE states METHOD cnexp
        gk = gkbar*n*n*n*n
	    ik = gk*(v - ek)
    }
    DERIVATIVE states {
        alpha = .01*-(v+55)/(exp(-(v+55)/10)-1)
        beta = .125*exp(-(v+65)/80)
        n' = alpha*(1-n)  - beta*n
    }
    """
    m2n.process_current_law(m2n.parse_mod(mod))

