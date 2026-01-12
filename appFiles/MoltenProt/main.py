
"""
Python module to handle many DSF fitters at the same time
"""
import os
import numpy as np

from moltenprot_shiny import DsfFitter

class ManyDsfFitters:

    """
    Useful to work with many different dsf experiments
    """

    def __init__(self):


        self.experiments = {}  # Dictionary where each key value pair corresponds to one DSF experiment
        self.all_signals = None  # List of all signals available in all experiments
        self.available_experiments = None  # List of names of all experiments available, that can handle the selected signal

    def add_experiment(self, file,name=None):

        if name is None:

            # Extract file basename without extension
            name = os.path.splitext(os.path.basename(file))[0]

        if name not in self.experiments.keys():

            fitter = DsfFitter()

            fitter.import_file(file)

            self.experiments[name] = fitter

            self.set_unique_signals()

        return None

    def set_unique_signals(self):

        all_signals = self.get_experiment_properties('signals')
        all_signals = np.concatenate(all_signals)

        if len(all_signals) > 0:

            _, idx = np.unique(all_signals, return_index=True)

            self.all_signals = all_signals[np.sort(idx)].tolist()

        return None

    def delete_experiments(self, names):

        if not isinstance(names, list):

            names = [names]

        for name in names:

            if name in self.experiments.keys():

                del self.experiments[name]

                self.set_unique_signals()

            if name in self.available_experiments:

                self.available_experiments.remove(name)

        return None

    def set_signal(self,signal):

        self.available_experiments = []

        for name, fitter in self.experiments.items():

            if signal in fitter.signals:

                fitter.set_signal(signal)

                self.available_experiments.append(name)

            else:

                fitter.fluo  = None
                fitter.temps = None

        return None

    def filter_by_temperature(self,min_temp,max_temp):

        for name in self.available_experiments:

            self.experiments[name].filter_by_temperature(min_temp,max_temp)

        return None

    def sort_by_conditions_name(self):

        for _, fitter in self.experiments.items():

            fitter.sort_by_conditions_name()

        return None

    def set_colors(self,color_list):

        counter = 0

        for _, fitter in self.experiments.items():

            conditions_i_n = len(fitter.conditions_original)

            color_list_i = color_list[counter:counter+conditions_i_n]

            counter += conditions_i_n

            fitter.set_colors(color_list_i)

        return None

    def set_conditions(self,conditions_list):

        counter = 0

        for _, fitter in self.experiments.items():

            conditions_i_n = len(fitter.conditions_original)

            conditions_list_i = conditions_list[counter:counter+conditions_i_n]

            counter += conditions_i_n

            fitter.set_conditions(conditions_list_i)

        return None

    def select_conditions(self,boolean_mask):

        counter = 0

        for _, fitter in self.experiments.items():

            conditions_i_n = len(fitter.conditions_original)

            boolean_mask_i = boolean_mask[counter:counter+conditions_i_n]

            counter += conditions_i_n

            fitter.select_conditions(boolean_mask_i)

        return None

    def get_experiment_properties(self, variable):

        return [getattr(self.experiments[name], variable) for name in self.experiments.keys()]

    def median_filter(self,n_degree_window=8):

        for name in self.available_experiments:

            self.experiments[name].median_filter(n_degree_window)

        return None

    def estimate_fluo_derivates(self,temp_window_length=8):

        for name in self.available_experiments:

            self.experiments[name].estimate_fluo_derivates(temp_window_length)

        return None

    def decompose_spectra(self,min_wl=0,max_wl=800):

        for _, fitter in self.experiments.items():

            if fitter.full_spectrum:

                fitter.decompose_spectra(min_wl,max_wl)

        self.set_unique_signals()

        return None

    def set_baseline_types(self,poly_order_native=1,poly_order_unfolded=1):

        for name in self.available_experiments:

            self.experiments[name].set_baseline_types(poly_order_native,poly_order_unfolded)

        return None

    def estimate_baselines_parameters(self,baseline_degree_window=12):

        for name in self.available_experiments:

            self.experiments[name].estimate_baselines_parameters(baseline_degree_window)

        return None

    def equilibrium_two_state(self):

        for name in self.available_experiments:

            self.experiments[name].equilibrium_two_state()

        return None

