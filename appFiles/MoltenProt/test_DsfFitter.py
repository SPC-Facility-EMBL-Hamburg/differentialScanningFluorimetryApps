# To be run with pytest

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

def test_dsf_import():

    for i,file in enumerate(all_files):

        fitters.add_experiment(file,str(i))

        fitter_i = fitters.experiments[str(i)]

        assert fitter_i.conditions is not None

        assert len(fitter_i.conditions) > 0

def test_set_signal():

    fitters.set_signal("350nm")

    assert fitters.available_experiments == ['0','3','5']

def test_temperature_range():

    fitters.filter_by_temperature(30,60)

    assert np.min(fitters.experiments['0'].temps) >= 30 + 273.15
    assert np.max(fitters.experiments['0'].temps) <= 60 + 273.15

def test_set_colors():

    all_conditions = fitters.get_experiment_properties('conditions')

    n_conditions = np.sum([len(x) for x in all_conditions])

    colors = ['#FF0000'] * n_conditions

    fitters.set_colors(colors)

    assert len(fitters.experiments['0'].all_colors) == len(all_conditions[0])

def test_set_conditions():

    all_conditions = fitters.get_experiment_properties('conditions')

    n_conditions = np.sum([len(x) for x in all_conditions])

    new_conditions = ['Condition'] * n_conditions

    fitters.set_conditions(new_conditions)

    updated_conditions = fitters.get_experiment_properties('conditions')

    assert all([all([cond == 'Condition' for cond in conds]) for conds in updated_conditions])

def test_select_conditions():

    all_conditions = fitters.get_experiment_properties('conditions')

    n_conditions = np.sum([len(x) for x in all_conditions])

    boolean_mask = [True] * n_conditions

    fitters.select_conditions(boolean_mask)

    assert fitters.experiments['0'].fluo.shape[1] == 8 # The panta file has 8 conditions
    assert fitters.experiments['1'].fluo is None # The signal 350nm is not available in qPCR file

def test_median_filter():

    fitters.median_filter(8)

    assert fitters.experiments['0'].fluo.shape[1] == 8 # The panta file has 8 conditions
    assert fitters.experiments['1'].fluo is None # The signal 350nm is not available in qPCR file

def test_estimate_fluo_derivates():

    fitters.estimate_fluo_derivates()

    assert fitters.experiments['0'].derivative.shape[1] == 8 # The panta file has 8 conditions
    assert fitters.experiments['0'].derivative2.shape[1] == 8 # The panta file has 8 conditions

def test_decompose_spectra():

    fitters.decompose_spectra()

    assert 'SVD Second coefficient' in fitters.all_signals

def test_set_baseline_types():

    fitters.set_baseline_types(2,2)

    assert fitters.experiments['0'].fit_kN
    assert fitters.experiments['0'].fit_qN

def test_estimate_baselines_parameters():

    fitters.delete_experiments(['0','1','2','4','5','6','7'])

    fitters.estimate_baselines_parameters()

    np.testing.assert_allclose(fitters.experiments['3'].qUs[0],0.366,rtol=1e-2)


