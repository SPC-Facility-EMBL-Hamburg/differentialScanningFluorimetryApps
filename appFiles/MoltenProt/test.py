
import numpy as np
from main import  ManyDsfFitters

panta_file = "./www/panta.xlsx"
qPCR_file = "./www/qPCRdemoFile.xls"
quantStudio_file = "./www/quantStudio.txt"
prometheus_file = "./www/demo.xlsx"
uncle_file = "./www/UNCLE_multi_channel.xlsx"
tycho_file = "./www/tychoFile.xlsx"
aunty_file = "./www/AUNTY_multi_channel.xlsx"
agilent_file = "./www/MX3005P.txt"

all_files = [
    panta_file,
    qPCR_file,
    quantStudio_file,
    prometheus_file,
    uncle_file,
    tycho_file,
    aunty_file,
    agilent_file
]

fitters = ManyDsfFitters()


for i,file in enumerate(all_files):

    fitters.add_experiment(file,str(i))

    fitter_i = fitters.experiments[str(i)]

    assert fitter_i.conditions is not None

    assert len(fitter_i.conditions) > 0


fitters.set_signal("350nm")


fitters.filter_by_temperature(30,84)


all_conditions = fitters.get_experiment_properties('conditions',mode="all",flatten=True)

n_conditions = len(all_conditions)

colors = ['#FF0000'] * n_conditions
print(colors)
fitters.set_colors(colors)

