import nmodl.dsl as nmodl
from nmodl import ast, visitor, symtab


with open("sample_mods/k_hh_ab.mod") as f:
    mod = f.read()

driver = nmodl.NmodlDriver()
modast = driver.parse_string(mod)

lookup_visitor = visitor.AstLookupVisitor()


def test_find_kinetics():
    useions = lookup_visitor.lookup(modast, ast.AstNodeType.USEION)
    bps = lookup_visitor.lookup(modast, ast.AstNodeType.BREAKPOINT_BLOCK)

    symv = symtab.SymtabVisitor()
    symv.visit_program(modast)
    table = modast.get_symbol_table()

    for useion in useions:
        print('using ion', useion.name)
        # locate ion current
        ix = table.lookup(f"i{useion.name}")
        print(ix)

    import pdb;pdb.set_trace()
