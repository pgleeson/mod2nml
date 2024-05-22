import sympy as sp
from symbolic import *

class baseHHrate:
    def __init__(self, r, m, s):
        self.rate = r
        self.midpoint = m
        self.scale = s
    def __repr__(self):
        return f'{type(self).__name__}[rate:{self.rate:.3f}, '\
                f'midpoint:{self.midpoint}, scale:{self.scale}]'

class HHExpRate(baseHHrate):
    @classmethod
    def match(cls, expr):
        ret = None
        if c := extract_exp_arg(expr):
            cann, V, s = c
            (r,) = sp.solve(cann - rate * sp.exp(V), rate)
            if len(r.free_symbols) < 1:  # no spurious solutions found
                ret = cls(r, s[midpoint], s[scale])
        return ret
        # V = Wild('V', real=True)
        # r = Wild('r', real=True)
        # vv = (v-midpoint)/scale
        # ret = {}
        # if m := expr.match(r*exp(V)):
        #    sol = solve(m[V]-vv, midpoint, scale)
        #    ret['rate'] = m[r]
        #    ret['midpoint'] = sol[midpoint]
        #    ret['scale'] = sol[scale]
        # return ret

    def as_symbolic(self):
        x = (v - self.midpoint) / self.scale
        return self.rate * sp.exp(x)


class HHSigmoidRate(baseHHrate):
    @classmethod
    def match(cls, expr):
        ret = None
        if c := extract_exp_arg(expr, inv=True):
            cann, V, s = c
            (r,) = sp.solve(cann - rate / (1 + sp.exp(-V)), rate)
            if len(r.free_symbols) < 1:  # no spurious solutions found
                ret = cls(r, s[midpoint], s[scale])
        return ret
        # V = Wild('V', real=True)
        # r = Wild('r', real=True)
        # vv = -(v-midpoint)/scale
        # ret = {}
        # if m := expr.match(r/(1+exp(V))):
        #    sol = solve(m[V]-vv, midpoint, scale)
        #    ret['rate'] = m[r]
        #    ret['midpoint'] = sol[midpoint]
        #    ret['scale'] = sol[scale]
        # return ret

    def as_symbolic(self):
        x = (v - self.midpoint) / self.scale
        return self.rate / (1 + sp.exp(-x))


class HHExpLinearRate(baseHHrate):
    @classmethod
    def match(cls, expr):
        ret = None
        if c := extract_exp_arg(expr, inv=True):
            cann, V, s = c
            vv = (v - s[midpoint]) / s[scale]
            (r,) = sp.solve(cann + vv.subs(s) * rate / (sp.exp(-V) - 1), rate)
            if len(r.free_symbols) < 1:  # no spurious solutions found
                ret = cls(r, s[midpoint], s[scale])
        return ret
        # V = Wild('V', real=True)
        # lin = Wild('lin', real=True)
        # vv = (v-midpoint)/scale
        # ret = {}
        # if m := expr.match(lin/(exp(V)-1)):
        #    sol = solve(m[V]-vv, midpoint, scale)
        #    m[lin] = -m[lin]
        # elif m := expr.match(lin/(1-exp(V))):
        #    sol = solve(m[V]-vv, midpoint, scale)
        # if m:
        #    sol[scale] = -sol[scale]
        #    ret['rate'] = solve(m[lin]/rate-vv.subs(sol),rate)[0]
        #    ret['midpoint'] = sol[midpoint]
        #    ret['scale'] = sol[scale]
        # return ret

    def as_symbolic(self):
        x = (v - self.midpoint) / self.scale
        return self.rate * x / (1 - sp.exp(-x))

def match_dynamics(exprs, statevar, replacements={}):
    n = statevar
    a = sp.Wild("a", exclude=[n])
    b = sp.Wild("b", exclude=[n])
    expr = exprs[statevar.name+'\\prime']
    ddn = expr.diff(n)  # d(dn/dt)/dn = -1/tau = -1/(alpha+beta) for 1st order kinetics
    #eqs = {}
    ret = None
    if m := (
        sp.solve(expr, n)[0].match(a / (a + b))
    ):  # if that matches, we have inf = alpha/(alpha + beta)
        alpha = m[a]
        beta = m[b]
        for r,v in replacements.items():
            if exprs[r] == alpha:
                alpha = v
            if exprs[r] == beta:
                beta = v
        ret = GateHHrates(alpha, beta)
        #eqs["alpha"] = m[a]
        #eqs["beta"] = m[b]
    elif s := sp.solve(expr, n):
        tau = -1/ddn
        inf = s[0]
        for r,v in replacements.items():
            if exprs[r] == tau:
                tau = v
            if exprs[r] == inf:
                inf = v
        ret = GateHHtauInf(tau, inf)
        #eqs["tau"] = -1 / ddn
        #eqs["inf"] = s[0]

    return ret
    #return eqs

class GateHHrates:
    def __init__(self, alpha, beta):
        self.forward = alpha
        self.reverse = beta
    def __repr__(self):
        return f'{type(self).__name__}[forward:{self.forward}, '\
                f'reverse:{self.reverse}]'

class GateHHtauInf:
    def __init__(self, tau, inf):
        self.time_course = tau
        self.steady_state = inf
    def __repr__(self):
        return f'{type(self).__name__}[time_course:{self.time_course}, '\
                f'steady_state:{self.steady_state}]'


def match_standard_rates(expr, subs={}):
    m = None
    for r in [cls for cls in baseHHrate.__subclasses__()]:
        if m := r.match(expr):
            break
    return m

def replace_standards_in_sequence(seq, ctxt):
    def absorb_scalar(name, expr): # TODO: factor out
        stds_to_reps = {v[1]:v[0] for v in replacements.values()}
        stds = [r[1] for r in replacements.values()]
        a = sp.Wild('a', properties=[lambda k: k.is_Float])
        s = sp.Wild('s', properties=[lambda k: k.is_Atom])
        if expr.has(*stds):
            if m := expr.match(a*s):
                std = stds_to_reps[m[s]]
                std.rate *= m[a]
                sym = next(syms)
                replacements[name] = (std, sym)
                ctxt[name] = sym
        return
    replacements = {}
    syms = sp.numbered_symbols("std",real=True)
    for name, expr in seq.items():
        syex = sp.S(expr, ctxt)
        if m := match_standard_rates(syex):
            syex = next(syms)
            replacements[name] = (m, syex)
        ctxt[name] = syex
        absorb_scalar(name, syex) #TODO: other optimizations possible here
    return ctxt, replacements

