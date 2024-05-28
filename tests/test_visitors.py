import nmodl.dsl as nmodl
from mod2nml.rename_visitor import RenameReusedVisitor

mod = """
DERIVATIVE states   {
    a = 10
    b = 20
    z = 2 * a * f(a)
    a = 30
    cai' = a + b + g(a^2)
    a = 1
    cai' = a + h(a^a)
}
"""

def test_rename_reused():
    driver = nmodl.NmodlDriver()
    modast = driver.parse_string(mod)

    param_visitor = RenameReusedVisitor()
    modast.accept(param_visitor)

    print(nmodl.to_nmodl(modast))

