# Parsing nmodl and generating NeuroML

One of the most time consuming and error-prone parts of converting NEURON models to NeuroML is the translation of nmodl mechanisms (.mod files). We intend to automate this process as much as possible.

Here are some test using BlueBrain's [nmodl parser](https://github.com/BlueBrain/nmodl/tree/master) instead of our previous experiments which relied on a [parser built from scratch](https://github.com/borismarin/pynmodl).

**This is a proof of concept!** Do not expect it to generate correct NeuroML for NMODL files in the wild. See [tests](tests/) for examples.


## Installation

`pip install git+https://github.com/borismarin/mod2nml.git`



## Usage

After installation, a command line application will be available. 
As an example of current capabilities, try it with the sample file [hh_clean.mod](tests/sample_mods/hh_clean.mod).
```console
a@b:~$ mod2nml hh_clean.mod
```

