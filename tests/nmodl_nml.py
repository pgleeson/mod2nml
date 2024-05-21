from collections import OrderedDict
import nmodl.dsl as nmodl
from nmodl import ast, visitor, symtab
import symbolic as sym
import nml_helpers as nml

#from functools import lru_cache
#@lru_cache
def parse_mod(modstring):
    driver = nmodl.NmodlDriver()
    ast = driver.parse_string(modstring)
    return ast


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

def get_local_vars(block):
    table = block.get_statement_block().get_symbol_table()
    props = symtab.NmodlType.local_var
    local_vars = {var.get_name() for var in
                  table.get_variables_with_properties(props)}
    return local_vars

def get_assignments(blk):
    lookup_visitor = visitor.AstLookupVisitor()
    binexprs = lookup_visitor.lookup(blk, ast.AstNodeType.BINARY_EXPRESSION)
    assignments = OrderedDict()
    for expr in binexprs:
        if expr.lhs.is_var_name() and expr.op.eval() == '=':
            n = expr.lhs.get_node_name()
            if expr.lhs.name.is_prime_name():
                n = n+'\\prime'
            assignments[n] = nmodl.to_nmodl(
                expr.rhs, {nmodl.ast.UNIT, nmodl.ast.UNIT_DEF})

    return assignments

def get_state_vars(modast):
    lookup_visitor = visitor.AstLookupVisitor()
    state_blk = lookup_visitor.lookup(modast, ast.AstNodeType.STATE_BLOCK)
    assert len(state_blk) == 1 #TODO: is >1 STATE even allowed?
    state_vars = {n.get_node_name() for n in
                  lookup_visitor.lookup(state_blk[0], ast.AstNodeType.ASSIGNED_DEFINITION)}
    return state_vars


def scope_for_block(block):
    prog = block
    while prog.is_program() is False:
        prog = prog.parent
    globs = get_global_vars(prog)
    locs = get_local_vars(block)
    return locs | globs


def find_currents(modast):
    lookup_visitor = visitor.AstLookupVisitor()
    useions = lookup_visitor.lookup(modast, ast.AstNodeType.USEION)
    bps = lookup_visitor.lookup(modast, ast.AstNodeType.BREAKPOINT_BLOCK)

    currents = {}
    for useion in useions:
        for w in useion.writelist:
            currents[w.get_node_name()] = useion.get_node_name()

    scope = scope_for_block(bps[0])# TODO: assuming single BP
    asgns = get_assignments(bps[0])

    curr_eqs = {c: sym.sympy_context_for_var(c, scope, asgns)
                for c in currents}

    return curr_eqs

def match_cond_states(cond, states):
    print(f'\t matching conductance {cond} to product of powers of {states}', end=': ')

    state_exps = {}
    state_syms = {s for s in cond.free_symbols if s.name in states} #TODO: should be obj prop
    rest = cond
    for s in state_syms:
        rest, state_exps[s] = rest.as_independent(s, as_Mul=True)
    state_exps['gbar'] = rest
    print(state_exps)
    return state_exps


def get_gate_dynamics(modast):
    lookup_visitor = visitor.AstLookupVisitor()
    deqs = lookup_visitor.lookup(modast, ast.AstNodeType.DERIVATIVE_BLOCK)
    dyn_eqs = {}
    if deqs:
        primes = lookup_visitor.lookup(deqs[0], ast.AstNodeType.PRIME_NAME)
        scope = scope_for_block(deqs[0])
        asgns = get_assignments(deqs[0])
        #dyn_eqs = {v.get_node_name(): sym.sympy_context_for_var(v.get_node_name(), scope, asgns)
        #         for v in primes}
        dyn_eqs = {v.get_node_name(): (scope, asgns) for v in primes}
    return dyn_eqs

def simplify_derivatives(exprs, states):
    res = {}
    for s in states:
        try:
            locvars, exprseq = exprs[s]
        except KeyError:
            print("Couldn't find dynamics for variable", s)
            return None
        ctxt = {s:sym.sp.Symbol(s, real=True) for s in locvars}
        #ctxt[s] = sym.sp.Symbol(s, real=True)
        replaced, replacements = nml.replace_standards_in_sequence(exprseq, ctxt)
        res[s] = nml.match_dynamics(replaced, ctxt[s], replacements)
    return res


def conductance(current, ctxt):
    print('\tCalculating conductance')
    curr_eq = ctxt[current]
    v = ctxt['v']
    g = curr_eq.diff(v)
    print('\tdI/dv:', g)
    if g != 0 and g.diff(v) == 0:
        print(f'\t{current} seems ohmic!')
    return g

def process_current_law(ast):
    for current, ctxt in find_currents(ast).items():
        print('Found current', current)
        #print('\tsymbolic context:', list(ctxt.items()))
        g = conductance(current, ctxt)

        states = get_state_vars(ast)
        print('\tfound state variables', states)

        gbar_n_gates = match_cond_states(g, states)
        print(gbar_n_gates)

        dxs = get_gate_dynamics(ast)
        simp_dxs = simplify_derivatives(dxs,states)

