import argparse

from mod2nml.nmodl_nml import compile_mod

def main():
    parser = argparse.ArgumentParser(description='Compile nmodl files to NeuroML')
    parser.add_argument('modfile')
    parser.add_argument('-t','--tests', help='Generate NEURON scripts to test', required=False)
    parser.add_argument('-o','--output', help='NeuroML file to output. Default is stdout', required=False)
    args = vars(parser.parse_args())

    compiled = compile_mod(args['modfile'])
    if fn := args['output']:
        from neuroml.writers import NeroMLWriter
        NeuroMLWriter.write(compiled, fn)
    else:
        print(compiled)
