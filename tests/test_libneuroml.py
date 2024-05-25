import neuroml
from neuroml.writers import NeuroMLWriter
from neuroml.utils import validate_neuroml2

chan = neuroml.IonChannelHH(
    id="k_hh",
    conductance="10pS",
    species="k",
)

n_gate = neuroml.GateHHRates(id="n", instances="4")

n_gate.forward_rate = neuroml.HHRate(
    type="HHExpLinearRate", rate="0.1per_ms", midpoint="-55mV", scale="10mV"
)

n_gate.reverse_rate = neuroml.HHRate(
    type="HHExprate", rate="0.125per_ms", midpoint="-65mV", scale="-80mV"
)


chan.gate_hh_rates.append(n_gate)

doc = neuroml.NeuroMLDocument()
doc.ion_channel_hhs.append(chan)

doc.id = "k_hh"
print(doc)
NeuroMLWriter.write(doc, '/tmp/k_hh.xml')



validate_neuroml2('/tmp/k_hh.xml')
