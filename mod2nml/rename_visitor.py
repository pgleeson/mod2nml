from nmodl import visitor
from nmodl import to_nmodl as N


class RenameReusedVisitor(visitor.AstVisitor):
    def __init__(self):
        visitor.AstVisitor.__init__(self)
        self.visited = []
        self.renamed = {}

    def visit_binary_expression(self, node):
        if node.op.eval() == '=':
            lhs = node.lhs
            n = lhs.name.value.value
            if n in self.visited:
                nn = self.renamed.get(n, n)
                new_name = nn + '_0'
                lhs.name.value.value = new_name
                self.renamed[n] = new_name
            self.visited.append(n)
            lhs.accept(self)
            node.rhs.accept(self)

    def visit_var_name(self, node):
        n = node.name.value.value
        if n in self.renamed:
            node.name.value.value = self.renamed[n]
