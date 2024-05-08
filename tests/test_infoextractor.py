import nmodl.dsl as nmodl
from nmodl import ast, visitor, symtab
from collections import OrderedDict
import sympy as sp
#from functools import lru_cache

#@lru_cache
def parse_mod(modstring):
    driver = nmodl.NmodlDriver()
    ast = driver.parse_string(modstring)
    return ast

def sympy_context_for_current(current, blk_vars, blk_assignments):
    # TODO: only do this for NAMEs relevant to the current
    # lookup_visitor = visitor.AstLookupVisitor()
    # names = {n.get_name() for n in
    #               lookup_visitor.lookup(bp, ast.AstNodeType.NAME)}
    # names_blk = {}

    # locals / globals
    ctxt = OrderedDict([(var, sp.symbols(var, real=True))
                        for var in blk_vars
                        if var not in blk_assignments.keys()])

    for var, eq in blk_assignments.items():
        ctxt[var] = sp.sympify(eq, ctxt)
        if var == current:
            break
    return ctxt

def get_global_vars(modast):
    symtab.SymtabVisitor().visit_program(modast)
    table = modast.get_symbol_table()
    st = symtab.NmodlType
    props = (st.global_var | st.local_var | st.range_var |
             st.param_assign | st.extern_var | st.prime_name |
             st.assigned_definition | st.read_ion_var |
             st.write_ion_var | st.nonspecific_cur_var |
             st.electrode_cur_var | st.constant_var |
             st.extern_neuron_variable | st.state_var)  # | st.factor_def
    global_vars = {var.get_name() for var in
                   table.get_variables_with_properties(props)}
    return global_vars

def get_local_vars(bp):
    table = bp.get_statement_block().get_symbol_table()
    props = symtab.NmodlType.local_var
    local_vars = {var.get_name() for var in
                  table.get_variables_with_properties(props)}
    return local_vars

def get_assignments(bp):
    lookup_visitor = visitor.AstLookupVisitor()
    binexprs = lookup_visitor.lookup(bp, ast.AstNodeType.BINARY_EXPRESSION)
    assignments = OrderedDict()
    for expr in binexprs:
        if expr.lhs.is_var_name() and expr.op.eval() == '=':
            assignments[expr.lhs.get_node_name()] = nmodl.to_nmodl(
                expr.rhs, {nmodl.ast.UNIT, nmodl.ast.UNIT_DEF})

    return assignments

def get_state_vars(modast):
    lookup_visitor = visitor.AstLookupVisitor()
    state_blk = lookup_visitor.lookup(modast, ast.AstNodeType.STATE_BLOCK)
    state_vars = {n.get_node_name() for n in
                  lookup_visitor.lookup(modast, ast.AstNodeType.ASSIGNED_DEFINITION)}
    return state_vars

def find_currents(modast):
    lookup_visitor = visitor.AstLookupVisitor()
    useions = lookup_visitor.lookup(modast, ast.AstNodeType.USEION)
    bps = lookup_visitor.lookup(modast, ast.AstNodeType.BREAKPOINT_BLOCK)

    currents = {}
    for useion in useions:
        # print('using ion', useion.name, 'writelist',
        #      [w.get_node_name()  for w in useion.writelist])
        for w in useion.writelist:
            currents[w.get_node_name()] = useion.get_node_name()

    globs = get_global_vars(modast)
    locs = get_local_vars(bps[0])  # TODO: assuming single BP
    asgns = get_assignments(bps[0])

    scope = locs | globs
    curr_eqs = {c: sympy_context_for_current(c, scope, asgns)
                for c in currents}

    return curr_eqs

def match_cond_states(cond, states):
    print(f'trying to match conductance {cond} to product of powers of {states}')

    #wilds = sp.numbered_symbols('_a', cls=sp.Wild)
    #pattern = sp.Wild('gbar')
    #for state in states:
    #    sym = next((s for s in cond.free_symbols if s.name == state), 1)
    #    pattern = pattern * sym**next(wilds)
    #print(cond.match(pattern))

    state_exps = {}
    state_syms = {s for s in cond.free_symbols if s.name in states} #TODO: should be obj prop
    rest = cond
    for s in state_syms:
        rest, state_exps[s] = rest.as_independent(s, as_Mul=True)
    state_exps['gbar'] = rest
    return state_exps

    print(state_exps)
def process_current_law(ast):
    currents = find_currents(ast)
    for current, ctxt in currents.items():
        print('Found current', current)
        #print('\tsymbolic context:', list(ctxt.items()))

        print('\tCalculating conductance')
        curr_eq = ctxt[current]
        v = ctxt['v']
        g = curr_eq.diff(v)
        print('\tdI/dv:', g)
        if g != 0 and g.diff(v) == 0:
            print(f'\t{current} seems ohmic!')

        states = get_state_vars(ast)
        print('\tfound state variables', states)
        print(match_cond_states(g, states))

    # ver se g é produto de state variables - serão as rates

    # criar um gatehhrate para cada
    # para a cinética, fazer otimização e tentar match
    # (minf-m)/tau:  gatehhtauinf
    # (1-n)*alpha - n*beta : forwardrate (HHexp HHSigmoid HHExplinear)


#def test_hhrates():
#    #with open("sample_mods/k_hh_ab.mod") as f:
#    with open("sample_mods/hh_ab.mod") as f:
#        mod = f.read()
#    process_current_law(parse_mod(mod))

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
    process_current_law(parse_mod(mod))

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
    process_current_law(parse_mod(mod))

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
    process_current_law(parse_mod(mod))
