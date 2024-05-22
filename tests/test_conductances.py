from collections import OrderedDict
import sympy as sp

import nmodl_nml as m2n

def test_hhrates():
    #with open("sample_mods/k_hh_ab.mod") as f:
    with open("sample_mods/hh_ab.mod") as f:
        mod = f.read()
    m2n.process_current_law(m2n.parse_mod(mod))

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
    currs = m2n.find_currents(m2n.parse_mod(mod))
    assert str(currs['ina']['ina']) == 'gna*(-ena + v)'

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
    ast = m2n.parse_mod(mod)
    currs = m2n.find_currents(ast)
    gna = m2n.conductance('ina', currs['ina'])
    assert str(gna) == 'gna*h*m**3'
    states = m2n.get_state_vars(ast)
    gates = m2n.match_cond_states(gna, states)
    gate_pow = {str(g): str(p) for g, p in gates.items()}
    assert gate_pow['m'] == 'm**3'
    assert gate_pow['h'] == 'h'

def test_gate_dynamics():
    mod = """
    NEURON {
        USEION k READ ek WRITE ik REPRESENTS CHEBI:29103
        RANGE gkbar, gk
        RANGE gna
        GLOBAL alpha, beta
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
        n' = alpha*(1-n) - beta*n
    }
    """
    ast = m2n.parse_mod(mod)
    states = m2n.get_state_vars(ast)
    dxs = m2n.get_gate_dynamics(ast)
    simp_dxs = m2n.simplify_derivatives(dxs,states)
    print(simp_dxs)


def test_gate_dynamics_ti():
    mod = """
    NEURON {
        USEION k READ ek WRITE ik REPRESENTS CHEBI:29103
        RANGE gkbar, gk
        RANGE gna
        GLOBAL inf, tau
    }
    PARAMETER { gkbar = .036 }
    STATE { n }
    BREAKPOINT {
        SOLVE states METHOD cnexp
        gk = gkbar*n*n*n*n
	    ik = gk*(v - ek)
    }
    DERIVATIVE states {
        tau = 10
        inf = 1/(1+exp(-(v+65)/80))
        n' = (inf-n)/tau
    }
    """
    ast = m2n.parse_mod(mod)
    states = m2n.get_state_vars(ast)
    dxs = m2n.get_gate_dynamics(ast)
    simp_dxs = m2n.simplify_derivatives(dxs,states)
    print(simp_dxs)

def test_subroutine():
    with open("sample_mods/k_hh.mod") as f:
        mod = f.read()
    AST = m2n.parse_mod(mod)
    m2n.nmodl.symtab.SymtabVisitor().visit_program(AST)
    #m2n.nmodl.visitor.ConstantFolderVisitor().visit_program(AST)
    m2n.nmodl.visitor.InlineVisitor().visit_program(AST)
    #m2n.nmodl.visitor.LocalVarRenameVisitor().visit_program(AST)
    m2n.process_current_law(AST)

