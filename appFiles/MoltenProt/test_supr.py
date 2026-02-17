from main import  ManyDsfFitters

file = "/home/os/Downloads/Experiment 12_02_2026-ramp-plate2.supr"

fitters = ManyDsfFitters()

fitters.add_experiment(file,"test1")

fitters.set_signal(fitters.all_signals[0])