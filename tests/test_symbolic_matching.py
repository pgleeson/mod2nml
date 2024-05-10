from sympy import symbols, Wild, solve, evaluate, parse_expr
from sympy import exp, log
v, rate, midpoint, scale, x, y, m, h, n, alpha, beta, ninf, taun = symbols(
    'v,rate,midpoint,scale,x,y,m,h,n,alpha,beta,ninf,taun', real=True)

#    x = (v-midpoint)/scale
#    hhexp = rate * exp(x)
#    hhsigmoid = rate / (1 + exp(-x))
#    hhexplin = rate * x / (1 - exp(-x))


def match_hhexprate(expr):
    V = Wild('V', real=True)
    r = Wild('r', real=True)
    vv = (v-midpoint)/scale
    ret = {}
    if m := expr.match(r*exp(V)):
        sol = solve(m[V]-vv, midpoint, scale)
        ret['rate'] = m[r]
        ret['midpoint'] = sol[midpoint]
        ret['scale'] = sol[scale]
    return ret

def test_match_hhexprate():
    beta = .125*exp(-(v+65)/80)
    rms = match_hhexprate(beta)
    assert rms['rate'] == 0.125
    assert rms['midpoint'] == -65
    assert rms['scale'] == -80

def test_match_standard_forms():
    for ex in ['.125*exp(-(v+65)/80)']:

        expr = parse_expr(ex, local_dict={'v':v})

        if  m := match_hhexprate(expr):
            print(ex, 'matches HHExpRate!')
            print('\t', m)

        if  m := match_hhexprate(expr):
            print(ex, 'matches HHExpRate!')
            print('\t', m)

