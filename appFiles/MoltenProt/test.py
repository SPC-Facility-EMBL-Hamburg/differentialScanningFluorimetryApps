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

fitters.add_experiment(prometheus_file,"test1")

fitters.set_signal(fitters.all_signals[0])
fitters.filter_by_temperature(26,39)

fitters.set_colors(["red"]*48)
fitters.select_conditions([True]*48)
fitters.estimate_fluo_derivates(3)

fitters.sort_by_conditions_name(True)
