import sympy as sp

v, rate, midpoint, scale = sp.symbols("v,rate,midpoint,scale", real=True)


def extract_exp_arg(expr, inv=False):
    ret = ()
    exps = expr.atoms(sp.exp)
    if exps:
        argexp = exps.pop().args[0]  # TODO: think about cases with >1 exp
        vv = (v - midpoint) / scale
        if inv:
            vv = -vv
        s = sp.solve(argexp - vv, midpoint, scale)
        V = sp.Symbol("V", real=True)
        cann = expr.subs(vv.subs(s), -V if inv else V)
        ret = cann, V, s
    return ret


class baseHHrate:
    def __init__(self, r, m, s):
        self.rate = r
        self.midpoint = m
        self.scale = s


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


STD_RATES = [cls for cls in baseHHrate.__subclasses__()]


def match_standard_rates(expr):
    matched = None
    for r in STD_RATES:
        if m := r.match(expr):
            matched = m
            break
    return matched


def replace_standards_in_sequence(seq, ctxt):
    syms = sp.numbered_symbols("std")
    for name, expr in seq.items():
        syex = sp.S(expr, ctxt)
        if match_standard_rates(sp.S(expr, ctxt).doit()):
            syex = next(syms)
        ctxt[name] = syex

    return ctxt


def match_alpha_beta_tau_inf(expr, statevar):
    n = statevar
    a = sp.Wild("a", exclude=[n])
    b = sp.Wild("b", exclude=[n])

    ddn = expr.diff(n)  # d(dn/dt)/dn = -1/tau = -1/(alpha+beta) for 1st order kinetics
    eqs = {}
    if m := (
        sp.solve(expr, n)[0].match(a / (a + b))
    ):  # if that matches, we have inf = alpha/(alpha + beta)
        eqs["alpha"] = m[a]
        eqs["beta"] = m[b]
    elif s := sp.solve(expr, n):
        eqs["tau"] = -1 / ddn
        eqs["inf"] = s[0]

    return eqs
