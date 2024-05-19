import sympy as sp
from collections import OrderedDict

from nml_helpers import *


def test_match_hhexp():
    beta = 0.125 * sp.exp(-(v + 65) / 80)
    r = HHExpRate.match(beta)
    assert r.rate == 0.125
    assert r.midpoint == -65
    assert r.scale == -80


def test_match_hhexplinear():
    alpha = -0.01 * (v + 60) / (sp.exp(-(v + 60) / 10) - 1)
    r = HHExpLinearRate.match(alpha)
    assert r.rate == 0.1
    assert r.midpoint == -60
    assert r.scale == 10


def test_match_hhsigmoid():
    minf = 1 / (1 + sp.exp((-v - 48) / 10))
    r = HHSigmoidRate.match(minf)
    assert r.rate == 1
    assert r.midpoint == -48
    assert r.scale == 10


def test_match_standard_rates():
    for ex in [
        ".125*exp(-(v+65)/80)",
        "1 / ( 1 + exp( ( - v - 48 ) / 10 ) )",
        "-.01*(v+55)/(exp(-(v+55)/10)-1)",
    ]:
        expr = sp.parse_expr(ex, local_dict={"v": v})
        matched = match_standard_rates(expr).as_symbolic()
        assert sp.simplify(expr / matched) == 1


# mAlpha = (0.182 * (v- -32))/(1-(exp(-(v- -32)/6)))
# mBeta  = (0.124 * (-v -32))/(1-(exp(-(-v -32)/6)))
# hAlpha = (-0.015 * (v- -60))/(1-(exp((v- -60)/6)))
# hBeta  = (-0.015 * (-v -60))/(1-(exp((-v -60)/6)))
# minf  = 1 / ( 1 + exp( ( - v - 48 ) / 10 ) )


def test_match_simple_odes():
    n = sp.Symbol("n", real=True)  # state var

    dnti = sp.parse_expr("(ninf-n)/ntau", local_dict={"n": n})
    m = match_alpha_beta_tau_inf(dnti, n)
    assert m["inf"].name == "ninf"
    assert m["tau"].name == "ntau"

    dnab = sp.parse_expr("(1-n)*alphan-n*betan", local_dict={"n": n})
    m = match_alpha_beta_tau_inf(dnab, n)
    assert m["alpha"].name == "alphan"
    assert m["beta"].name == "betan"


def test_match_multiexpr_odes():
    n = sp.Symbol("n", real=True)  # state var
    ctxt = dict([("n", n), ("v", v)])
    exprs = OrderedDict()

    exprs["alphan"] = "-.01*(v+55)/(exp(-(v+55)/10)-1)"
    exprs["betan"] = ".125*exp(-(v+65)/80)"
    exprs["sum"] = "alphan + betan"
    exprs["taun"] = "1/sum"
    exprs["ninf"] = "alphan/sum"
    exprs["dn"] = "(ninf-n)/taun"

    replaced = replace_standards_in_sequence(exprs, ctxt)
    m = match_alpha_beta_tau_inf(replaced["dn"], n)
    assert m["alpha"] == replaced["alphan"]
    assert m["beta"] == replaced["betan"]
