import nmodl.dsl as nmodl
from nmodl import ast, visitor


with open("sample_mods/nap.mod") as f:
    mod = f.read()

driver = nmodl.NmodlDriver()
modast = driver.parse_string(mod)

lookup_visitor = visitor.AstLookupVisitor()

def test_useion():
    useions = lookup_visitor.lookup(modast, ast.AstNodeType.USEION)
    for useion in useions:
        print('using ion', useion.name)

def test_states():
    states = lookup_visitor.lookup(modast, ast.AstNodeType.STATE_BLOCK)
    for state in states:
        print(nmodl.to_nmodl(state))

