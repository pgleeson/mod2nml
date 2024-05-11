from sympy import symbols, Wild, solve, evaluate, parse_expr, simplify
from sympy import exp, log
#v, rate, midpoint, scale, x, y, m, h, n, alpha, beta, ninf, taun = symbols(
    #'v,rate,midpoint,scale,x,y,m,h,n,alpha,beta,ninf,taun', real=True)

v, rate, midpoint, scale = symbols('v,rate,midpoint,scale', real=True)

#    x = (v-midpoint)/scale
#    hhexp = rate * exp(x)
#    hhsigmoid = rate / (1 + exp(-x))
#    hhexplin = rate * x / (1 - exp(-x))


def match_HHExpRate(expr):
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

def match_HHSigmoidRate(expr):
    V = Wild('V', real=True)
    r = Wild('r', real=True)
    vv = -(v-midpoint)/scale
    ret = {}
    if m := expr.match(r/(1+exp(V))):
        sol = solve(m[V]-vv, midpoint, scale)
        ret['rate'] = m[r]
        ret['midpoint'] = sol[midpoint]
        ret['scale'] = sol[scale]
    return ret

def match_HHExpLinearRate(expr):
    V = Wild('V', real=True)
    lin = Wild('lin', real=True)
    vv = (v-midpoint)/scale
    ret = {}
    if m := expr.match(lin/(exp(V)-1)):
        sol = solve(m[V]-vv, midpoint, scale)
        m[lin] = -m[lin]
    elif m := expr.match(lin/(1-exp(V))):
        sol = solve(m[V]-vv, midpoint, scale)
    if m:
        sol[scale] = -sol[scale]
        ret['rate'] = solve(m[lin]/rate-vv.subs(sol),rate)[0]
        ret['midpoint'] = sol[midpoint]
        ret['scale'] = sol[scale]
    return ret

def test_match_hhexp():
    beta = .125*exp(-(v+65)/80)
    rms = match_HHExpRate(beta)
    assert rms['rate'] == 0.125
    assert rms['midpoint'] == -65
    assert rms['scale'] == -80

def test_match_hhexplinear():
    alpha = -.01*(v+60)/(exp(-(v+60)/10)-1)
    rms = match_HHExpLinearRate(alpha)
    assert rms['rate'] == 0.1
    assert rms['midpoint'] == -60
    assert rms['scale'] == 10

def test_match_hhsigmoid():
    minf  = 1 / ( 1 + exp( ( - v - 48 ) / 10 ) )
    rms = match_HHSigmoidRate(minf)
    assert rms['rate'] == 1
    assert rms['midpoint'] == -48
    assert rms['scale'] == 10


def test_match_standard_forms():
    for ex in ['.125*exp(-(v+65)/80)',
               '1 / ( 1 + exp( ( - v - 48 ) / 10 ) )',
               '-.01*(v+55)/(exp(-(v+55)/10)-1)']:

        expr = parse_expr(ex, local_dict={'v':v})

        if  m := match_HHExpRate(expr):
            print(ex, 'matches HHExpRate!')
            x = (v-m['midpoint'])/m['scale']
            matched = m['rate'] * exp(x)
            print('\t', m)

        if  m := match_HHExpLinearRate(expr):
            print(ex, 'matches HHExpLinearRate!')
            print('\t', m)
            #m['scale'] = -m['scale']
            x = (v-m['midpoint'])/m['scale']
            matched = m['rate'] * x / (1 - exp(-x))

        if  m := match_HHSigmoidRate(expr):
            print(ex, 'matches HHSigmoidRate!')
            print('\t', m)
            x = (v-m['midpoint'])/m['scale']
            matched = m['rate'] / (1 + exp(-x))
        assert simplify(expr / matched) == 1
# mAlpha = (0.182 * (v- -32))/(1-(exp(-(v- -32)/6)))
# mBeta  = (0.124 * (-v -32))/(1-(exp(-(-v -32)/6)))
# hAlpha = (-0.015 * (v- -60))/(1-(exp((v- -60)/6)))
# hBeta  = (-0.015 * (-v -60))/(1-(exp((-v -60)/6)))
# minf  = 1 / ( 1 + exp( ( - v - 48 ) / 10 ) )



#def test_match_odes():
#    tauinf = (ninf - n)/taun
#    alphabeta = (1 - n)*alpha - n * beta
#
#    al = -.01*(v+55)/(sp.exp(-(v+55)/10)-1)
#    be = .125*sp.exp(-(v+65)/80)
#    tau = 1/(al+be)
#    inf = alpha/tau
#    dn = (inf-n)/tau
#
#    sp.solve(dn-tauinf, ninf,taun)

