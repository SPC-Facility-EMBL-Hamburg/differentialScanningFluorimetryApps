
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

            read_status = fitter.import_file(file)

            self.experiments[name] = fitter

            self.set_unique_signals()

            return read_status

        return False

    def set_unique_signals(self):

        all_signals = self.get_experiment_properties('signals')
        all_signals = np.concatenate(all_signals)

        if len(all_signals) > 0:

            _, idx = np.unique(all_signals, return_index=True)

            self.all_signals = all_signals[np.sort(idx)].tolist()

        full_spectrum = self.get_experiment_properties('full_spectrum')
        self.full_spectrum = np.any(full_spectrum)

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

    def apply_to_available_experiments(self,method_name,**kwargs):

        for name in self.available_experiments:

            obj = self.experiments[name]
            method = getattr(obj, method_name,None)
            method(**kwargs)

        return None

    def filter_by_temperature(self,min_temp,max_temp):

        kwargs = {'min_temp':min_temp,'max_temp':max_temp}

        self.apply_to_available_experiments('filter_by_temperature',**kwargs)

        return None

    def sort_by_conditions_name(self,sort=False):

        for _, fitter in self.experiments.items():

            fitter.sort_by_conditions_name(sort)

        return None

    def set_colors(self,color_list):

        counter = 0

        for _, fitter in self.experiments.items():

            conditions_i_n = len(fitter.conditions_original)

            color_list_i = color_list[counter:counter+conditions_i_n]

            counter += conditions_i_n

            fitter.set_colors(color_list_i)

        return None

    def set_conditions(self,conditions_list,original=False):

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

    def get_experiment_properties(self, variable,flatten=False):

        list = [getattr(self.experiments[name], variable) for name in self.experiments.keys()]

        if flatten:

            list = np.concatenate(list)

        return list

    def median_filter(self,n_degree_window=8):

        kwargs = {'n_degree_window':n_degree_window}
        self.apply_to_available_experiments('median_filter',**kwargs)

        return None

    def estimate_fluo_derivates(self,temp_window_length=8):

        kwargs = {'temp_window_length':temp_window_length}
        self.apply_to_available_experiments('estimate_fluo_derivates',**kwargs)

        return None

    def decompose_spectra(self,min_wl=0,max_wl=800):

        for _, fitter in self.experiments.items():

            if fitter.full_spectrum:

                fitter.decompose_spectra(min_wl,max_wl)

        self.set_unique_signals()

        return None

    def set_baseline_types(self,poly_order_native=1,poly_order_unfolded=1):

        kwargs = {'poly_order_native':poly_order_native,'poly_order_unfolded':poly_order_unfolded}
        self.apply_to_available_experiments('set_baseline_types',**kwargs)

        return None

    def estimate_baselines_parameters(self,baseline_degree_window=12):

        kwargs = {'baseline_degree_window':baseline_degree_window}
        self.apply_to_available_experiments('estimate_baselines_parameters',**kwargs)

        return None

    def equilibrium_two_state(self,delta_cp):

        kwargs = {'delta_cp':delta_cp}
        self.apply_to_available_experiments('equilibrium_two_state',**kwargs)
        self.model_name = "EquilibriumTwoState"

        return None

    def equilibrium_three_state(self,t1min=25,t1max=90,t2min=25,t2max=90):

        kwargs = {'t1min':t1min,'t1max':t1max,'t2min':t2min,'t2max':t2max}
        self.apply_to_available_experiments('equilibrium_three_state',**kwargs)
        self.model_name = "EquilibriumThreeState"

        return None

    def empirical_two_state(self):

        self.apply_to_available_experiments('empirical_two_state')
        self.model_name = "EmpiricalTwoState"
        return None

    def empirical_three_state(self,t1min=20,t1max=70,t2min=40,t2max=90):

        kwargs = {'t1min':t1min,'t1max':t1max,'t2min':t2min,'t2max':t2max}
        self.apply_to_available_experiments('empirical_three_state',**kwargs)
        self.model_name = "EmpiricalThreeState"
        return None

    def irreversible_two_state(self,scan_rate):

        kwargs = {'scan_rate':scan_rate}
        self.apply_to_available_experiments('irreversible_two_state',**kwargs)
        self.model_name = "IrreversibleTwoState"
        return None

    def filter_by_relative_error(self,threshold_percentage):

        kwargs = {'threshold_percentage':threshold_percentage}
        self.apply_to_available_experiments('filter_by_relative_error',**kwargs)

        return None

    def filter_by_fitting_std_error(self,threshold_std_error):

        kwargs = {'threshold_std_error':threshold_std_error}
        self.apply_to_available_experiments('filter_by_fitting_std_error',**kwargs)

        return None

    def filter_by_param_values(self,param_name,low_value,high_value):

        kwargs = {'param_name':param_name,'low_value':low_value,'high_value':high_value}
        self.apply_to_available_experiments('filter_by_param_values',**kwargs)

        return None

