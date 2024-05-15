from collections import OrderedDict

from sympy import (
    S,
    Symbol,
    Wild,
    exp,
    numbered_symbols,
    parse_expr,
    simplify,
    solve,
    symbols,
)

v, rate, midpoint, scale = symbols("v,rate,midpoint,scale", real=True)

#    x = (v-midpoint)/scale
#    hhexp = rate * exp(x)
#    hhsigmoid = rate / (1 + exp(-x))
#    hhexplin = rate * x / (1 - exp(-x))


def extract_exp_arg(expr, inv=False):
    ret = ()
    exps = expr.atoms(exp)
    if exps:
        argexp = exps.pop().args[0]  # TODO: think about cases with >1 exp
        vv = (v - midpoint) / scale
        if inv:
            vv = -vv
        s = solve(argexp - vv, midpoint, scale)
        V = Symbol("V", real=True)
        cann = expr.subs(vv.subs(s), -V if inv else V)
        ret = cann, V, s
    return ret


def match_HHExpRate(expr):
    ret = {}
    if c := extract_exp_arg(expr):
        cann, V, s = c
        (r,) = solve(cann - rate * exp(V), rate)
        if len(r.free_symbols) < 1:  # no spurious solutions found
            ret["midpoint"] = s[midpoint]
            ret["scale"] = s[scale]
            ret["rate"] = r
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


def match_HHSigmoidRate(expr):
    ret = {}
    if c := extract_exp_arg(expr, inv=True):
        cann, V, s = c
        (r,) = solve(cann - rate / (1 + exp(-V)), rate)
        if len(r.free_symbols) < 1:  # no spurious solutions found
            ret["midpoint"] = s[midpoint]
            ret["scale"] = s[scale]
            ret["rate"] = r
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


def match_HHExpLinearRate(expr):
    ret = {}
    if c := extract_exp_arg(expr, inv=True):
        cann, V, s = c
        vv = (v - s[midpoint]) / s[scale]
        (r,) = solve(cann + vv.subs(s) * rate / (exp(-V) - 1), rate)
        if len(r.free_symbols) < 1:  # no spurious solutions found
            ret["midpoint"] = s[midpoint]
            ret["scale"] = s[scale]
            ret["rate"] = r
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


def test_match_hhexp():
    beta = 0.125 * exp(-(v + 65) / 80)
    rms = match_HHExpRate(beta)
    assert rms["rate"] == 0.125
    assert rms["midpoint"] == -65
    assert rms["scale"] == -80


def test_match_hhexplinear():
    alpha = -0.01 * (v + 60) / (exp(-(v + 60) / 10) - 1)
    rms = match_HHExpLinearRate(alpha)
    assert rms["rate"] == 0.1
    assert rms["midpoint"] == -60
    assert rms["scale"] == 10


def test_match_hhsigmoid():
    minf = 1 / (1 + exp((-v - 48) / 10))
    rms = match_HHSigmoidRate(minf)
    assert rms["rate"] == 1
    assert rms["midpoint"] == -48
    assert rms["scale"] == 10


def match_standard_forms(expr):
    matched = None
    if m := match_HHExpRate(expr):
        print(expr, "matches HHExpRate!")
        x = (v - m["midpoint"]) / m["scale"]
        matched = m["rate"] * exp(x)
        print("\t", m)

    elif m := match_HHExpLinearRate(expr):
        print(expr, "matches HHExpLinearRate!")
        print("\t", m)
        # m['scale'] = -m['scale']
        x = (v - m["midpoint"]) / m["scale"]
        matched = m["rate"] * x / (1 - exp(-x))

    elif m := match_HHSigmoidRate(expr):
        print(expr, "matches HHSigmoidRate!")
        print("\t", m)
        x = (v - m["midpoint"]) / m["scale"]
        matched = m["rate"] / (1 + exp(-x))
    return matched


def test_match_standard_forms():
    for ex in [
        ".125*exp(-(v+65)/80)",
        "1 / ( 1 + exp( ( - v - 48 ) / 10 ) )",
        "-.01*(v+55)/(exp(-(v+55)/10)-1)",
    ]:
        expr = parse_expr(ex, local_dict={"v": v})
        matched = match_standard_forms(expr)
        assert simplify(expr / matched) == 1


# mAlpha = (0.182 * (v- -32))/(1-(exp(-(v- -32)/6)))
# mBeta  = (0.124 * (-v -32))/(1-(exp(-(-v -32)/6)))
# hAlpha = (-0.015 * (v- -60))/(1-(exp((v- -60)/6)))
# hBeta  = (-0.015 * (-v -60))/(1-(exp((-v -60)/6)))
# minf  = 1 / ( 1 + exp( ( - v - 48 ) / 10 ) )


def replace_standards_in_sequence(seq, ctxt):
    syms = numbered_symbols("std")
    for name, expr in seq.items():
        syex = S(expr, ctxt)
        if match_standard_forms(S(expr, ctxt).doit()):
            syex = next(syms)
        ctxt[name] = syex

    return ctxt


def match_alpha_beta_tau_inf(expr, statevar):
    n = statevar
    a = Wild("a", exclude=[n])
    b = Wild("b", exclude=[n])

    ddn = expr.diff(n)  # d(dn/dt)/dn = -1/tau = -1/(alpha+beta) for 1st order kinetics
    eqs = {}
    breakpoint()
    if m := (
        solve(expr, n)[0].match(a / (a + b))
    ):  # if that matches, we have inf = alpha/(alpha + beta)
        eqs["alpha"] = m[a]
        eqs["beta"] = m[b]
    elif s := solve(expr, n):
        eqs["tau"] = -1 / ddn
        eqs["inf"] = s[0]

    return eqs


def test_match_simple_odes():

    n = symbols("n", real=True)  # state var

    dnti = parse_expr("(ninf-n)/ntau", local_dict={"n": n})
    print(match_alpha_beta_tau_inf(dnti, n))

    dnab = parse_expr("(1-n)*alphan-n*betan", local_dict={"n": n})
    print(match_alpha_beta_tau_inf(dnab, n))


def test_match_multiexpr_odes():
    n = Symbol("n", real=True)  # state var
    ctxt = dict([("n", n), ("v", v)])
    exprs = OrderedDict()

    exprs["alphan"] = "-.01*(v+55)/(exp(-(v+55)/10)-1)"
    exprs["betan"] = ".125*exp(-(v+65)/80)"
    exprs["sum"] = "alphan + betan"
    exprs["taun"] = "1/sum"
    exprs["ninf"] = "alphan/sum"
    exprs["dn"] = "(ninf-n)/taun"

    replaced = replace_standards_in_sequence(exprs, ctxt)
    print(replaced)
    print(match_alpha_beta_tau_inf(replaced["dn"], n))
