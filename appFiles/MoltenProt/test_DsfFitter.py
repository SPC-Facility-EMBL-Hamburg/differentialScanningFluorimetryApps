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

    fitters.filter_by_temperature(30,84)

    assert np.min(fitters.experiments['0'].temps) >= 30 + 273.15
    assert np.max(fitters.experiments['0'].temps) <= 84 + 273.15

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

    fitters.median_filter(2)

    assert fitters.experiments['0'].fluo.shape[1] == 8 # The panta file has 8 conditions
    assert fitters.experiments['1'].fluo is None # The signal 350nm is not available in qPCR file

    # Reset back to no filter
    fitters.median_filter(0)

def test_estimate_fluo_derivates():

    fitters.estimate_fluo_derivates(temp_window_length=8)

    assert fitters.experiments['0'].derivative.shape[1] == 8 # The panta file has 8 conditions
    assert fitters.experiments['0'].derivative2.shape[1] == 8 # The panta file has 8 conditions

    expected = [348.44,345.85,342.19,337.56]

    np.testing.assert_allclose(fitters.experiments['3'].tms_from_deriv[:4],expected,rtol=1e-2)

def test_decompose_spectra():

    fitters.decompose_spectra()

    assert 'SVD Second coefficient' in fitters.all_signals

def test_set_baseline_types():

    fitters.set_baseline_types(2,2)

    assert fitters.experiments['0'].fit_kN
    assert fitters.experiments['0'].fit_qN

def test_estimate_baselines_parameters():

    fitters.delete_experiments(['0','1','2','4','5','6'])

    fitters.estimate_baselines_parameters()

    np.testing.assert_allclose(fitters.experiments['3'].qUs[0],-1.01,rtol=1e-2)

def test_equilibrium_two_state():

    boolean_mask = [True] * 4 + [False] * 44 + [True] * 3
    fitters.select_conditions(boolean_mask)
    fitters.estimate_baselines_parameters()
    fitters.equilibrium_two_state(delta_cp=0)

    np.testing.assert_allclose(fitters.experiments['3'].dG_std[1],13.99,rtol=1e-2)

def test_empirical_two_state():

    fitters.empirical_two_state()

    np.testing.assert_allclose(fitters.experiments['3'].score[1],479.21,rtol=1e-2)

def test_filter_by_relative_error():

    mask = fitters.filter_by_relative_error(100)

    assert mask == [False] + [True] * 3

def test_filter_by_fitting_std_error():

    mask = fitters.filter_by_fitting_std_error(10)

    assert mask == [True] * 4

def test_filter_by_param_values():

    mask = fitters.filter_by_param_values('Tm',50+273.15,75+273.15)

    assert mask == [True] * 4

def test_get_params():

    params_all = fitters.get_experiment_properties('params_all',flatten=True)

    assert len(params_all) == 4 # Four conditions were fitted

    params_name = fitters.get_experiment_properties('params_name',flatten=True)

    assert len(params_name) == 8

def test_remove_all_experiments():

    exps = list(fitters.experiments.keys())

    for exp in exps:

        fitters.delete_experiments(exp)

    assert len(fitters.experiments) == 0

def test_demo_and_aunty_set_signal():

    fitters.add_experiment(prometheus_file,"prometheus_demo")
    fitters.add_experiment(aunty_file,"aunty_demo")

    fitters.set_signal(fitters.all_signals[0])  

    fitters.set_conditions(['Cond1']*49)
    fitters.set_colors(['#00FF00']*49)
    fitters.select_conditions([True]*49)

    fitters.estimate_fluo_derivates(temp_window_length=8)

    fitters.set_signal("350nm")
    fitters.set_conditions(['Cond1']*49)
    fitters.set_colors(['#00FF00']*49)
    fitters.select_conditions([True]*49)

    fitters.estimate_fluo_derivates(temp_window_length=8)

    # Verify we have only 48 derivatives

    derivatives = fitters.get_experiment_properties('derivative',flatten=True)

    print(derivatives)

    assert len(derivatives[0]) == 48
