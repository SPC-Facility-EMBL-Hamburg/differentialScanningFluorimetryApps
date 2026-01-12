import json

import pandas as pd
import numpy  as np

from scipy.signal     import savgol_filter

from scipy.optimize   import curve_fit
from scipy.integrate  import solve_ivp
from natsort          import index_natsorted

from svd import (
    apply_svd
)

# Load custom helper functions
from helpers import (
    is_blank,
    temperature_to_kelvin,
    detect_file_type,
    detect_encoding,
    find_indexes_of_non_signal_conditions,
    find_repeated_words,
    remove_words_in_string,
    get_quantstudio_start_line,
    generate_merged_quantstudio_df,
    subset_panta_data,
    get_sheet_names_of_xlsx,
    get_barycenter,
    median_filter_from_fluo_and_temp_vectors,
    filter_fluo_by_temp,
    filter_temp_by_temp,
    fit_line_robust,
    fit_quadratic_robust,
    get_temp_at_maximum_of_derivative,
    get_temp_at_minimum_of_derivative,
    get_two_state_deltaG,
    get_two_state_tonset,
    get_eq_three_state_combined_deltaG,
    get_empirical_two_state_score,
    get_empirical_three_state_score,
    get_irrev_two_state_pkd,
    generate_2D_signal_matrix
)

from models import(

    get_fit_fx_two_state_signal_fit,
    get_fit_fx_three_state_signal_fit

)

from constants import temp_standard

class DsfFitter:

    """
    Class for nanoDSF protein stability studies
    This class was written by Osvaldo Burastero and is
    based on the work done by Vadim Kotov et al.

    If you use this script please cite:

		Kotov, Vadim, et al. "In‐depth interrogation of protein thermal unfolding data with MoltenProt."
		Protein Science 30.1 (2021): 201-217.

		Burastero, Osvaldo, et al. "eSPC: an online data-analysis platform for molecular biophysics."
		Acta Crystallographica Section D: Structural Biology 77.10 (2021).

    No warranty whatsoever
    If you have questions or would like to add new models please contact me:
    	oburastero@gmail.com, oburastero@embl-hamburg.de

    Remember -  All models are wrong but some models are useful
                or in other words, any model is at best a useful fiction!!!!!!

    """

    def __init__(self):

        """
        We initialize all attributes with None.
        """

        # Dictionary containing the raw signal data - one key per signal type, such as 350 nm or 330 nm
        self.signal_data_dictionary = None

        # Dictionary containing the raw temperature data in Kelvin units - one key per signal type, such as 350 nm or 330 nm
        self.temp_data_dictionary   = None

        # List with the possible signals, matches the keys in self.signal_data_dictionary and self.temp_data_dictionary
        self.signals = None

        # List with the conditions names, can be modified
        self.conditions = None

        # List with the original conditions, read from the file
        self.conditions_original = None

        # List of colors for plotting, one per condition
        self.all_colors = None

        # List of colors for plotting, one per selected condition
        self.colors = None

        # Minimum wavelength, in case of reading files with whole spectrum data
        self.min_wavelength = None

        # Maximum wavelength, in case of reading files with whole spectrum data
        self.max_wavelength = None

        # Selected signal type, one of self.signals
        self.signal_type = None

        # Indexes to reorder the conditions back to its original order
        self.idx_no_sorted = None

        # Average time step
        self.dt = None

        # Heating rate in degrees per second. Must be set by the User.
        self.scan_rate = None

        # Signal matrix of size n*m. n signals and m selected conditions
        self.fluo = None

        # 1D numpy array
        self.temps = None

        # Boolean to indicate if we have full spectrum data
        self.full_spectrum = False

        # Maximum temperature in the selected temperature range
        self.max_temp = None

        # Minimum temperature in the selected temperature range
        self.min_temp = None

        # First derivative of self.fluo
        self.derivative = None

        # Second derivative of self.fluo
        self.derivative2 = None

        # Empirical estimation of the Tm based on the derivative values
        self.tms_from_deriv = None

        # Initial Tms based on an heuristic algorithm
        self.tms_initial = None

        # Model name, such as EquilibriumTwoState
        self.model_name = None

        # Fitted parameters values
        self.params_all = None

        # Std error of the fitted parameters
        self.errors_abs_all = None

        # Indexes of the conditions that were successfully fitted
        self.fitted_conditions_indexes = None

        # Relative error of the fitted parameters
        self.errors_percentage_all = None

        # Predicted signal
        self.fluo_predictions_all = None

        # Std error of the fitted values
        self.std_error_estimate_all = None

        # List of booleans, one per fitted conditions, describing if the fitted parameters are 'far' from the fitting bounds.
        self.parameters_far_from_bounds = None

        # Intercepts of the native state
        self.bNs = None

        # Intercepts of the unfolded state
        self.bUs = None

        # Slopes of the native state
        self.kNs = None

        # Slopes of the unfolded state
        self.kUs = None

        # Quadratic terms of the native state
        self.qNs = None

        # Quadratic terms of the unfolded state
        self.qUs = None

        # Signal of the intermediate
        self.bIs = None

        # Boolean to indicate if the baseline are constant, linear or quadratic
        self.fit_kN = None
        self.fit_kU = None
        self.fit_qN = None
        self.fit_qU = None

        # Score provided to the user to rank the fitted conditions
        self.score = None

        # Combined DG when fitting a three-state model
        self.dG_comb_std = None

        # Temperature where 1% of the protein is unfolded. It can be an empirical parameter from the Empirical model or a thermodynamic derived quantity from the Equilibrium model
        self.T_onset = None

        # DG at the standard temperature
        self.dG_std = None

        # Delta CP value, must be set directly by the User
        self.cp = None

        # Delta CP component to correct the DG value
        self.dCp_component = None

        # Names of the fitted parameters
        self.params_name = None

        # List of conditions that were fitted
        self.fitted_conditions = None

        # Columns of the signal matrix that were fitted
        self.fitted_fluo = None

        # Score for the irreversible two-state unfolding model
        self.pkd = None

        # Score for the three-state empirical model
        self.T_eucl_comb = None


    def init_dictionary_to_store_fluo_temp_data(self):

        """
        Create 2 dictionaries that will share the same keys, i.e., "350nm", but will have different values

        signal_data_dictionary keys contain the signal data (one for each signal) - an n*m matrix where n is the number of measurements and
                                                                                    m is the number of different samples

        temp_data_dictionary   keys contain the temperature (one for each signal) - a vector of length n (same n as before)

        This function has to be called by all the functions that load the DSF data
        For input file examples check spc.embl-hamburg.de

        """

        self.signal_data_dictionary = {}
        self.temp_data_dictionary   = {}

        return None

    def import_file(self,file):

        file_type = detect_file_type(file)

        if file_type == "prometheus":

            sheet_names = get_sheet_names_of_xlsx(file)
            self.load_nano_dsf_xlsx(file,sheet_names)
            return True

        else:

            read_fx_map = {
                'tycho' : self.load_tycho_xlsx,
                'thermofluor' : self.load_thermofluor_xlsx,
                'panta' : self.load_panta_xlsx,
                'QuantStudio' : self.load_quantstudio_txt,
                'MX3005P' : self.load_agilent_mx3005p_qPCR_txt,
                'csv' : self.load_csv_file,
                'supr' : self.load_supr_dsf,
                'uncle' : self.load_uncle_multi_channel,
                'aunty' : self.load_aunty_xlsx_file
            }

            if file_type in read_fx_map.keys():

                read_fx = read_fx_map.get(file_type)
                read_fx(file)
                return True

            else:
                return False

    def load_nano_dsf_xlsx(self, processed_dsf_file, sheet_names=None):

        """
        Load the xlsx file (processed) generated by the Nanotemper Prometheus instrument that has one sheet called 'Overview' with a column called
        'Sample ID' with the names of the samples, and four sheets called 'Ratio', '330nm', '350nm' and 'Scattering'.
        The first column of the signal sheet ('Ratio', '330nm', '350nm', 'Scattering') should be called 'Time [s]'.
        The second column should have the temperature data and all subsequent columns store the fluorescence data.
        The order of the fluorescence columns should match the order of the 'Sample ID' column in the 'Overview' sheet.

        Input  - file path of the xlsx file and signal that we want to load (should match the sheet names)

        Result - self.signals has the signals we loaded i.e., ["350nm","330nm","Scattering","Ratio"]
                 self.signal_data_dictionary and self.temp_data_dictionary have the signal and temperature data
                 self.conditions (and self.conditions_original) has the sample names, if present.
        """

        if sheet_names is None:
            sheet_names = ["350nm", "330nm", "Scattering", "Ratio"]

        # Load an Excel file (this needs to be the processed file!)
        xls = pd.ExcelFile(processed_dsf_file)

        conditions_df = pd.read_excel(xls, "Overview")

        conditions    = conditions_df[['Sample ID']].values.flatten().astype(str)

        self.conditions_original = conditions
        # Add position index to avoid problems with empty names

        self.conditions = [cond + " P" + str(i+1) if is_blank(cond) else cond for i,cond in enumerate(conditions)]

        self.init_dictionary_to_store_fluo_temp_data()

        possible_signals    = ["350nm","330nm","Scattering","Ratio","Turbidity"]

        include             = []

        # Change signal name if unfolding curve is present
        for sn in sheet_names:

            include_value = any([ps in sn and "deriv" not in sn.lower() for ps in possible_signals])

            include.append(include_value)

        sheet_names_to_load = np.array([s for (i, s) in zip(include, sheet_names) if i])

        self.signals  = sheet_names_to_load

        for signal in self.signals:

            dat = pd.read_excel(xls, signal, index_col=None, header=None)

            indices   = np.argwhere(dat.iloc[:, 0].values == 'Time [s]')
            first_row = int(indices[0][0]) + 1

            # Find all the columns with the word temperature
            column_headers = dat.iloc[first_row-1, :].values

            ids_temperature = [i for i,x in enumerate(column_headers) if 'temperature' in x.lower()]
            ids_time        = [i for i,x in enumerate(column_headers) if 'time' in x.lower()]

            # If we have more than one temperature column, remove the extra temperature and time columns
            if len(ids_temperature) > 1:

                ids_temperature = ids_temperature[1:] # We keep the first temperature colummn
                ids_time        = ids_time[1:]        # We keep the first time column

                # Remove columns of dataframe by index
                dat = dat.drop(dat.columns[ids_temperature+ids_time], axis=1)

            fluo   = np.array(dat.iloc[first_row:, 2:]).astype('float')
            temp   = np.array(dat.iloc[first_row:, 1]).astype('float')
            temp   = temperature_to_kelvin(temp)

            # Reduce data so we can plot and fit the data faster.
            while len(temp) > 700:
                fluo = fluo[::2]
                temp = temp[::2]

            self.signal_data_dictionary[signal]   = fluo
            self.temp_data_dictionary[signal]     = temp

        self.full_spectrum = False

        return None

    def load_tycho_xlsx(self,file):
        
        """
        Load the xlsx file generated by the Nanotemper Prometheus Tycho instrument.
        This file has one sheet called ‘Results’ (with 6 columns named '#', 'Capillary label',..., 'Sample Brightness'),
        and one sheet called 'Profiles_raw' where the fluorescence data is stored.

        The ‘Profiles_raw’ sheet columns should have the following structure:

        One row with information about the recorded signal, e.g., ‘Ratio 350 nm / 330 nm’,  ‘Brightness @ 330 nm’, ‘Brightness @ 350 nm’.
        One row with the capillary numbers.
        One row with the time, temperature and sample names.
        The remaining rows store the temperature and signal data.

        Input  - file path of the xlsx file

        Result - self.signals has the signals we loaded i.e., ["350nm","330nm","Scattering","Ratio"]
                 self.signal_data_dictionary and self.temp_data_dictionary have the signal and temperature data
                 self.conditions (and self.conditions_original) has the sample names, if present.

        """

        xls = pd.ExcelFile(file)

        # Retrieve the conditions names - in sheet "Results"
        df = pd.read_excel(xls, "Results")

        for col_index, column in enumerate(df):

            column_values      = df[column]
            has_desired_value = column_values.isin(["Capillary label"])

            if any(has_desired_value):

                row_index_begin   = np.flatnonzero(has_desired_value)[0]+1
                initial_ratio_col = df.iloc[row_index_begin:,col_index+4]
                try:
                    row_index_end     = np.flatnonzero(initial_ratio_col.isnull())[0] + row_index_begin
                except:
                    row_index_end     = row_index_begin+6
                break

        conditions         = np.array(df.iloc[row_index_begin:row_index_end,col_index])
        number_of_conditions = len(conditions)

        self.conditions_original = conditions
        # Add position index to avoid problems with empty names
        self.conditions = [cond + " P" + str(i+1) if is_blank(cond) else cond for i,cond in enumerate(conditions)]

        self.init_dictionary_to_store_fluo_temp_data()

        # Retrieve the fluorescence signal - in sheet "Profiles_raw"
        df = pd.read_excel(xls, "Profiles_raw")

        for col_index, column in enumerate(df):

            column_values      = df[column]
            has_desired_value = column_values.isin(["Temperature [°C]"])

            if any(has_desired_value):

                row_index_begin   = np.flatnonzero(has_desired_value)[0]+1
                row_index_end     = np.flatnonzero(column_values[row_index_begin:].isnull())[0]+row_index_begin

                temperature     = np.array(df.iloc[row_index_begin:row_index_end,col_index]).astype('float')
                temperature     = temperature_to_kelvin(temperature)
                break

        column_names  = df.iloc[row_index_begin-1,:]

        start_indexes = [i for i, j in enumerate(column_names) if j == conditions[0]]

        signals      = np.array(df.iloc[row_index_begin-3,start_indexes])
        signals_clean = []
        
        for s in signals:
            if "ratio" in s.lower():
                signals_clean.append("Ratio")
            else:  
                if "330" in s.lower() :
                    signals_clean.append("330nm")
                if "350" in s.lower():
                    signals_clean.append("350nm")

        for i, signal in enumerate(signals_clean):

            col_start_index   = start_indexes[i]
            col_end_index     = col_start_index + number_of_conditions

            fluo            = np.array(df.iloc[row_index_begin:row_index_end,col_start_index:col_end_index]).astype('float')
            temp            = temperature

            while (len(temp)) > 600:

                fluo = fluo[::2]
                temp = temp[::2]

            self.signal_data_dictionary[signal]   = fluo
            self.temp_data_dictionary[signal]     = temp

        self.signals = np.array(signals_clean)

        self.full_spectrum = False

        return None

    def load_thermofluor_xlsx(self, thermofluor_file):

        """
        Load DSF Thermofluor xls file and extract data
        The xls file generated by the ThermoFluor assay in a qPCR instrument.
        This file has one sheet called 'RFU' where the first row has the sample positions (header), the first column has the temperature data,
        and all subsequent columns store the fluorescence data.

        Input  - file path of the xls file and signal that we want to load

        Result -
                 self.signal_data_dictionary and self.temp_data_dictionary have the signal and temperature data
                 self.conditions (and self.conditions_original) has the sample names.

        """

        xls  = pd.ExcelFile(thermofluor_file)
        dat = pd.read_excel(xls, "RFU",header=None)
        conditions = np.array(dat.iloc[0, 1:])

        self.conditions_original = conditions
        self.conditions = conditions

        self.init_dictionary_to_store_fluo_temp_data()

        signal = "DSF_RFU"
        self.signal_data_dictionary[signal]   = np.array(dat.iloc[1:,1:]).astype('float')
        self.temp_data_dictionary[signal]     = np.array(dat.iloc[1:, 0]).astype('float')
        self.temp_data_dictionary[signal]     = temperature_to_kelvin(self.temp_data_dictionary[signal])

        self.signals = np.array([signal])

        self.full_spectrum = False

        return None

    def load_panta_xlsx(self, panta_file):

        """
        Load the xlsx file generated by a Prometheus Panta instrument.
        This file has one sheet called ‘Overview’ with a column called 'Sample ID' with the names of the samples,
        and one sheet called ‘Data Export’ where all the data is stored. The ‘Data Export’ sheet columns should have the following order:

        Temperature capillary 1 ; Ratio capillary 1 ; … ; Temperature capillary 1 ; 350 nm capillary 1 ; … ;
        Temperature capillary 1 ; 330 nm capillary 1 ; … ; Temperature capillary 1 ; scattering capillary 1 ; … ;
        Temperature capillary 2 ; Ratio capillary 2; … ;  Temperature capillary n ; Ratio capillary n.

        Columns whose names include "Derivative" are not read.

        Input  - file path of the xlsx file and signal that we want to load (should match the sheet names)

        Result - self.signals has the signals we loaded i.e., ["350nm","330nm","Scattering","Ratio"]
                 self.signal_data_dictionary and self.temp_data_dictionary have the signal and temperature data
                 self.conditions (and self.conditions_original) has the sample names, if present.

        """

        try:

            data          = pd.read_excel(panta_file, "Data Export")

        except:

            data          = pd.read_excel(panta_file, "melting-scan") # Alternative format of PANTA

        column_names  = [str.lower(c) for c in data.columns]

        pos_350       = [i for i,x in enumerate(column_names) if "350"        in x and "deriv" not in x and "330" not in x]
        pos_330       = [i for i,x in enumerate(column_names) if "330"        in x and "deriv" not in x and "350" not in x]
        scattering    = [i for i,x in enumerate(column_names) if "scattering" in x and "deriv" not in x]
        pos_ratio     = [i for i,x in enumerate(column_names) if "ratio"      in x and "deriv" not in x]
        pos_turb      = [i for i,x in enumerate(column_names) if "turbidity"  in x and "deriv" not in x]

        possible_signals    = ["350nm","330nm","Scattering","Ratio","Turbidity"]
        signals             = []

        self.init_dictionary_to_store_fluo_temp_data()

        all_positions = [pos_350,pos_330,scattering,pos_ratio,pos_turb]

        for positions,signal in zip(all_positions,possible_signals):

            if len(positions) > 0:

                fluo, temp = subset_panta_data(data,positions)

                self.signal_data_dictionary[signal]   = fluo
                self.temp_data_dictionary[signal]     = temp 

                signals.append(signal)

        self.signals        = np.array(signals)

        try:

            conditions_df = pd.read_excel(panta_file, "Overview")
            conditions    = conditions_df[['Sample ID']].values.flatten().astype(str)

        except:

            conditions = np.repeat('',fluo.shape[1])

        self.conditions_original = conditions
        # Add position index to avoid problems with empty names
        self.conditions = [cond + " P" + str(i+1) if is_blank(cond) else cond for i,cond in enumerate(conditions)]

        self.full_spectrum = False

        return None

    def load_quantstudio_txt(self, qs_file):

        """
        Input: A txt file ('QSfile') where column 2 has the well position,
        column 3 the temperature and column 4 the fluorescence signal. Index starts at 1!!!

        --- Caution ---
        The first rows of the file are comments that start with '*' and are not read
        The temperature and signal column have commas that need to be deleted

        Columns are separated by spaces

        Input  - file path of the txt file

        Result -
                 self.signal_data_dictionary and self.temp_data_dictionary have the signal and temperature data
                 self.conditions (and self.conditions_original) has the sample names, if present.

        """

        start_row = get_quantstudio_start_line(qs_file)
        data      = pd.read_csv(qs_file, skiprows=start_row, sep=r"\s+", header=None)

        u, ind     = np.unique(data.iloc[:,1], return_index=True)
        conditions = u[np.argsort(ind)]

        self.conditions_original = conditions
        self.conditions = conditions

        fluo , temp = generate_merged_quantstudio_df(data)

        self.init_dictionary_to_store_fluo_temp_data()
        signal = "Fluorescence"

        self.signal_data_dictionary[signal]   = fluo
        self.temp_data_dictionary[signal]     = temp
        self.signals = np.array([signal])

        self.full_spectrum = False

        return None

    def load_agilent_mx3005p_qPCR_txt(self,filename):

        """
        Input: A txt file where the 2nd column has the fluorescence data, and
        the 3rd column the temperature.
        Wells are separated by rows containing a sentence like this one: 'Segment  2 Plateau  1 Well  1'
        """

        dfs      = []
        well_num = []

        with open(filename, 'r') as f:

            ls  = f.read().splitlines()

            for i,line in enumerate(ls):
                if line.startswith('Segment') and 'Well' in line:
                    # Get the well number
                    well_num .append(line.split()[-1])

                    fluorescence = []
                    temperature  = []

                    for line2 in ls[i+2:]:
                        if line2.startswith('Segment'):
                            break
                        else:
                            data = line2.split()
                            fluorescence.append(float(data[1]))
                            temperature.append(float(data[2]))

                    df = pd.DataFrame({'temperature':temperature,'signal'+str(well_num):fluorescence})
                    df.sort_values('temperature',inplace=True)
                    dfs.append(df)

        # Combine dataframes so we can obtain a vector of temperatures and a matrix of fluorescence signal
        merged = pd.merge_asof(dfs[0].dropna(),dfs[1].dropna(),on="temperature", direction='nearest',allow_exact_matches=True)

        for df in dfs[2:]:
            
            merged = pd.merge_asof(merged,df.dropna(),on="temperature", direction='nearest')       

        fluo   = np.array(merged.iloc[:, 1:]).astype('float')
        temp   = np.array(merged.iloc[:, 0]).astype('float')  

        # Reduce data so we can plot and fit the data faster.
        while len(temp) > 700:
            fluo = fluo[::2]
            temp = temp[::2]

        self.init_dictionary_to_store_fluo_temp_data()
        signal = "Fluorescence"

        self.conditions_original            = well_num
        self.conditions                     = well_num
        self.signal_data_dictionary[signal] = fluo
        self.temp_data_dictionary[signal]   = temperature_to_kelvin(temp)

        self.signals = np.array([signal])

        self.full_spectrum = False

        return None

    def load_csv_file(self,file):

        """
        Input: A csv file where the first column has the temperature and all the next columns the
        fluorescence data, header is required

        """

        self.init_dictionary_to_store_fluo_temp_data()
        signals = []

        encoding = detect_encoding(file)

        # Try common delimiters
        for delimiter in [',',';','\t']:

            try:

                dat  = pd.read_csv(file,delimiter=delimiter,encoding=encoding)

                # Convert non-numeric columns to NaN
                dat = dat.apply(pd.to_numeric, errors='coerce')

                # Produce error if we don't have 2 or more columns
                if len(dat.columns) < 2:
                    raise ValueError('File does not have enough columns')

                # Set the conditions names and start index for the signal data
                if 'time' in dat.columns[0].lower():
                    conditions = [str(c) for c in dat.columns[2:]]
                    idx_start = 1
                else:
                    conditions = [str(c) for c in dat.columns[1:]]
                    idx_start = 0

                signal_data      = np.array(dat.iloc[:, (idx_start + 1):]).astype('float')
                temperature_data = np.array(dat.iloc[:, idx_start]).astype('float')

                # Produce error if signal data is empty
                if signal_data.size == 0:
                    raise ValueError('Signal data is empty')

                # Produce error if temperature data  is non-numeric
                if np.isnan(temperature_data).any():
                    raise ValueError('Temperature data is non-numeric')

                break

            except:

                pass

        idx_to_remove = find_indexes_of_non_signal_conditions(signal_data,conditions)

        # Remove the elements from the conditions array
        conditions = [cond for i,cond in enumerate(conditions) if i not in idx_to_remove]

        # Remove the columns from the signal data
        signal_data = np.delete(signal_data,idx_to_remove,axis=1)

        temperature_data = temperature_to_kelvin(temperature_data)

        # Divide the conditions into groups, according to the presence of the words '350', '330', 'ratio', 'scattering'
        conditions_350nm = [cond for cond in conditions if '350' in cond.lower() and 'ratio' not in cond.lower()]
        conditions_330nm = [cond for cond in conditions if '330' in cond.lower() and 'ratio' not in cond.lower()]
        conditions_ratio = [cond for cond in conditions if 'ratio' in cond.lower()]
        conditions_scatt = [cond for cond in conditions if 'scattering' in cond.lower()]

        # Check that the number of conditions is the same for all groups
        conditions_lst = [conditions_350nm,conditions_330nm,conditions_ratio,conditions_scatt]
        sel_cond_lst   = [cond for cond in conditions_lst if len(cond) > 0]
        n_conditions   = [len(cond) for cond in conditions_lst if len(cond) > 0]

        possible_rep_conditions = len(np.unique(n_conditions)) == 1

        if possible_rep_conditions:

            # For each group, store the signal data and the temperature data
            if len(conditions_350nm) > 0:
                signal_350 = signal_data[:, [i for i, cond in enumerate(conditions) if cond in conditions_350nm]]
                self.signal_data_dictionary["350nm"] = signal_350
                self.temp_data_dictionary["350nm"] = temperature_data
                signals.append("350nm")

            if len(conditions_330nm) > 0:
                signal_330 = signal_data[:, [i for i, cond in enumerate(conditions) if cond in conditions_330nm]]
                self.signal_data_dictionary["330nm"] = signal_330
                self.temp_data_dictionary["330nm"] = temperature_data
                signals.append("330nm")

            if len(conditions_ratio) > 0:
                signal_ratio = signal_data[:, [i for i, cond in enumerate(conditions) if cond in conditions_ratio]]
                self.signal_data_dictionary["Ratio 350 nm / 330 nm"] = signal_ratio
                self.temp_data_dictionary["Ratio 350 nm / 330 nm"] = temperature_data
                signals.append("Ratio 350 nm / 330 nm")

            if len(conditions_scatt) > 0:
                signal_scatt = signal_data[:, [i for i, cond in enumerate(conditions) if cond in conditions_scatt]]
                self.signal_data_dictionary["Scattering"] = signal_scatt
                self.temp_data_dictionary["Scattering"] = temperature_data
                signals.append("Scattering")

            cond_temp =  sel_cond_lst[0]

            if len(cond_temp) > 1:

                repeated_words = find_repeated_words(cond_temp)
                cond_temp      = [remove_words_in_string(cond,repeated_words) for cond in cond_temp]
                conditions     = cond_temp

            else:

                conditions = cond_temp

        else:

            # If we only have one condition, use the condition as signal name
            if len(conditions) == 1:
                signal_name = conditions[0]
            # Default signal name
            else:
                signal_name = "Fluorescence"

            self.signal_data_dictionary[signal_name] = signal_data
            self.temp_data_dictionary[signal_name]   = temperature_data
            signals.append(signal_name)

        self.conditions            = conditions
        self.conditions_original   = conditions

        self.signals = np.array(signals)

        return None

    def load_supr_dsf(self, json_file):

        """
        Load Supr_DSF JSON export format and interpolate signals to a fixed temperature grid.

        Args:
            json_data = file.read()
        """

        self.init_dictionary_to_store_fluo_temp_data()

        # Read JSON data from a file
        with open(json_file, "r") as file:
            json_data = file.read()

        # Parse JSON data into a dictionary
        data_dict     = json.loads(json_data)

        samples_name = [item["SampleName"] for item in data_dict['Samples']]
        samples_well = [item["WellLocations"] for item in data_dict['Samples']]

        samples_name_simple = []
        samples_well_simple = []

        for sn,sw in zip(samples_name,samples_well):

            if ',' in sw:

                sw = sw.split(',')
                sn = [sn for _ in sw]

            else:

                sw = [sw]
                sn = [sn]

            samples_name_simple += sn
            samples_well_simple += sw

        name_df = pd.DataFrame({
            'well': samples_well_simple,
            'name': samples_name_simple})

        scans   = [item["_scans"] for item in data_dict['Wells']]
        n_scans = len(scans)

        wavelengths   = data_dict['Wavelengths']
        wavelengths   = np.round(wavelengths, decimals=1)
        n_wavelengths = len(wavelengths)

        temperatures = []
        signals      = []

        temperature_fixed = np.arange(5,110,0.5)

        well = [item["PhysicalLocation"] for item in data_dict['Wells']]

        # Create a categorical data type with the custom order
        cat_type = pd.CategoricalDtype(categories=well, ordered=True)

        # Convert the column to the categorical data type
        name_df['well'] = name_df['well'].astype(cat_type)

        # Sort the DataFrame based on the custom order
        name_df = name_df.sort_values(by='well')

        conditions = name_df['name'].values.astype(str)

        self.conditions_original =  conditions
        self.conditions          =  conditions

        for i in range(n_scans):

            temperatures.append([item['Temperature'] for item in scans[i]])
            signals.append([item['Signal']           for item in scans[i]])

        temperatures = np.array(temperatures).T
        temperatures = np.round(temperatures, decimals=1)

        signals = np.array(signals)

        self.wavelengths = wavelengths
        named_wls = [str(wl) + ' nm' for wl in wavelengths]

        for i in range(n_wavelengths):

            signals_temp  = signals[:,:,i].T
            signal_interp = []

            # Iterate over the columns of the arrays
            for ii in range(n_scans):

                x = temperatures[:, ii]
                y = signals_temp[:, ii]

                sorted_indices = np.argsort(x)
                x = x[sorted_indices]
                y = y[sorted_indices]

                y_interpolated = np.interp(temperature_fixed, x, y,left=np.nan,right=np.nan)

                # Linear interpolation every 0.5 degrees
                signal_interp.append(y_interpolated)

            fluo = np.array(signal_interp).T
            
            non_nas  = np.logical_not(np.isnan(fluo).any(axis=1))

            wl = named_wls[i]
            self.signal_data_dictionary[wl]   = fluo[non_nas,:]
            self.temp_data_dictionary[wl]     = temperature_to_kelvin(temperature_fixed[non_nas])

        barycenters = []

        for i, condition in enumerate(conditions):

            # We need to extract the signal data for this condition
            y = signals[i,:,:]
            x = temperatures[:, i]

            barycenter = np.apply_along_axis(get_barycenter, 1, y,wavelengths = wavelengths)

            # Interpolate every 0.5 degrees
            y_interpolated = np.interp(
                temperature_fixed,
                x, barycenter,
                left=np.nan, right=np.nan)

            barycenters.append(y_interpolated)

        # Create a matrix from the barycenters
        barycenters = np.array(barycenters).T

        non_nas = np.logical_not(np.isnan(barycenters).any(axis=1))

        self.signal_data_dictionary['BCM'] = barycenters[non_nas, :]
        self.temp_data_dictionary['BCM']   = temperature_to_kelvin(temperature_fixed[non_nas])

        named_wls = ['BCM'] + named_wls

        self.signals = np.array(named_wls)

        # Assign minimum and maximum wavelengths to self
        self.min_wavelength = np.floor(np.min(wavelengths))
        self.max_wavelength = np.ceil(np.max(wavelengths))

        self.full_spectrum = True

        # Now we try to add the ratio signal
        idx350 = np.argmin(np.abs(wavelengths - 350))
        idx330 = np.argmin(np.abs(wavelengths - 330))

        num_wl = named_wls[idx350]
        den_wl = named_wls[idx330]

        self.create_ratio_signal(num_wl,den_wl)

        return None

    def load_uncle_multi_channel(self,uncle_file):

        """
        Load UNCLE multichannel Excel file, extract per-wavelength matrices and interpolate to fixed temps.

        Args:
            uncle_file (str): Path to the UNCLE .xlsx file.

        Function to load the data from the UNCLE instrument
        Returns:
        Below is an example of the data format:

            Toms run.uni			Sample Name	1 mg/ml Protein
	        Temp :25, Time:134.9	Temp :25.49, Time:185.4
            Wavelength	Intensity	Intensity

            249.18	22.384	22.417
            249.66	13.144	14.261

        The line should be discarded
        The second line contains the time and temperature data. We're interested only in the temperature data
        The third line can be discarded
        The fourth line is empty
        The fifth line contains the first row of the signal data, with the wavelength data in the first column

        The file is a xlsx file with as many sheets as channels
        """

        # Get the names of the sheets
        sheet_names = get_sheet_names_of_xlsx(uncle_file)

        # Create empty dictionaries to store the data
        self.init_dictionary_to_store_fluo_temp_data()

        temperature_fixed = np.arange(5, 110, 0.5)

        wavelengths       = None
        temperatures      = []
        conditions        = []
        signals           = []

        # Loop through each sheet and read the data
        for sheet_name in sheet_names:

            try:

                # Read the data from the sheet
                data = pd.read_excel(
                    uncle_file,
                    sheet_name=sheet_name,
                    header=None,
                    skiprows=0
                )

                # Extract the sample name, from the first row, fifth column
                sample_name = data.iloc[0, 4]

                # Remove the first row
                data = data.iloc[1:, :]

                # Extract the time/temperature data
                temperature_data = data.iloc[0, 1:].values
                # Select the temperature data
                temperature_data = [x.split(',')[0] for x in temperature_data]
                temperature_data = [x.split(':')[1] for x in temperature_data]
                temperature_data = np.array(temperature_data,dtype=float)

                # Check that we have temperature data
                if len(temperature_data) < 10:
                    continue

                # Extract the signal data
                # It contains one column per temperature and one row per wavelength
                signal_data = np.array(data.iloc[3:, 1:].values,dtype=float)

                # Check that we have non nan signal data
                non_nas = np.logical_not(np.isnan(signal_data).any(axis=1))

                if np.sum(non_nas) < 100:
                    continue

                # Assign the wavelength data if wavelengths is None
                if wavelengths is None:
                    wavelengths = np.round(
                        np.array(data.iloc[3:, 0].values,dtype=float),
                        decimals=1)

                signals.append(signal_data)
                temperatures.append(temperature_data)
                conditions.append(sample_name)

            except:

                pass

        # Now we interpolate the signal data to the given fixed temperature vector
        # Iterate over the wavelengths
        self.wavelengths = wavelengths
        named_wls = [str(wl) + ' nm' for wl in wavelengths]

        # We require one signal matrix per wavelength
        # with one column per condition
        for i in range(len(wavelengths)):

            # Extract the signal data for the current wavelength
            # signal_temp has one row per condition and one column per wavelength
            signal_temp = np.array([signal[i,:] for signal in signals])

            signal_interp = []

            # Iterate over the conditions
            for row in range(signal_temp.shape[0]):

                # Extract the signal data for the current condition
                y = signal_temp[row, :]

                # Interpolate the signal data to the fixed temperature vector
                y_interpolated = np.interp(
                    temperature_fixed,
                    temperatures[row], y,
                    left=np.nan, right=np.nan)

                # Store the interpolated signal data
                signal_interp.append(y_interpolated)

            # Convert signal_interp to a 2D array
            # It should have one column per condition
            fluo = np.array(signal_interp).T

            non_nas = np.logical_not(np.isnan(fluo).any(axis=1))

            wl                              = named_wls[i]
            self.signal_data_dictionary[wl] = fluo[non_nas, :]
            self.temp_data_dictionary[wl]   = temperature_to_kelvin(temperature_fixed[non_nas])

        barycenters = []

        for i, condition in enumerate(conditions):

            # We need to extract the signal data for this condition
            signal_temp      = signals[i]
            temperature_temp = temperatures[i]

            barycenter = np.apply_along_axis(get_barycenter, 0, signal_temp,wavelengths = wavelengths)

            # Interpolate every 0.5 degrees
            y_interpolated = np.interp(
                temperature_fixed,
                temperature_temp, barycenter,
                left=np.nan, right=np.nan)

            barycenters.append(y_interpolated)

        # Create a matrix from the barycenters
        barycenters = np.array(barycenters).T

        non_nas = np.logical_not(np.isnan(barycenters).any(axis=1))

        self.signal_data_dictionary['BCM'] = barycenters[non_nas, :]
        self.temp_data_dictionary['BCM']   = temperature_to_kelvin(temperature_fixed[non_nas])

        named_wls = ['BCM'] + named_wls

        self.signals = np.array(named_wls)

        self.conditions          = np.array(conditions)
        self.conditions_original = np.array(conditions)

        # Assign minimum and maximum wavelengths to self
        self.min_wavelength = np.floor(np.min(wavelengths))
        self.max_wavelength = np.ceil(np.max(wavelengths))
        self.full_spectrum  = True

        # Now we try to add the ratio signal
        idx350 = np.argmin(np.abs(wavelengths - 350))
        idx330 = np.argmin(np.abs(wavelengths - 330))

        num_wl = named_wls[idx350]
        den_wl = named_wls[idx330]

        self.create_ratio_signal(num_wl,den_wl)

        return None

    def load_aunty_xlsx_file(self,file_path):

        """
        Load AUNTY-format multi-sheet Excel file where each sheet is a condition.

        Args:
            file_path (str): Path to the AUNTY .xls/.xlsx file.
        Import the AUNTY xlsx file which has the following format
        Many sheets where the name is the condition name
        The data is stored as follows:
            The first column has the temperature data
            Subsequent columns have the fluorescence data
            The first row indicates the wavelength
        Returns:
            EMPTY_CELL	wavelength
            temperature	250	    250.490005493164	250.979995727539	251.470001220703
            15.00	    98.22	58.0299987792969	49.9000015258789	93.1100006103516
        # Create empty dictionaries to store the data
        """
        # Get the names of the sheets
        sheet_names = get_sheet_names_of_xlsx(file_path)

        self.init_dictionary_to_store_fluo_temp_data()

        signals      = []
        conditions   = []
        wavelengths  = None
        temperature  = None

        for sheet_name in sheet_names:

            # Load the data
            data = pd.read_excel(
                file_path, sheet_name=sheet_name,
                header=None,skiprows=0)

            try:

                wavelength_cell = data.iloc[0, 1]
                temperature_cell = data.iloc[1, 0]
                corner_cell = data.iloc[0, 0]

                # Verify that the word wavelength is in the first row, second column
                condition1 = isinstance(wavelength_cell, str) and 'wavelength' in wavelength_cell.lower()

                # Verify that the word temperature is in the second row, first column
                condition2 = isinstance(temperature_cell, str) and 'temperature' in temperature_cell.lower()

                # Verify that the corner cell is empty (NaN)
                condition3 = np.isnan(corner_cell)

                if not (condition1 and condition2 and condition3):
                    continue

                fluorescence  = np.array(data.iloc[2:, 1:]).astype(float)
                signals.append(fluorescence)

                conditions.append(sheet_name)

                if wavelengths is None:

                    wavelengths    = np.round(np.array(data.iloc[1, 1:]).astype(float),2)
                    min_wavelength = np.min(wavelengths)
                    max_wavelength = np.max(wavelengths)

                    self.min_wavelength = np.floor(min_wavelength)
                    self.max_wavelength = np.ceil(max_wavelength)

                if temperature is None:

                    temperature   = np.round(np.array(data.iloc[2:, 0]).astype(float), 2)

            except:

                continue

        # Now we have one signal matrix per condition, but we need one signal matrix per wavelength, with
        # the conditions as columns
        # Therefore, for each wavelength, we create one dataframe per condition and then merge them
        self.wavelengths = wavelengths
        named_wls = [str(wl) + ' nm' for wl in wavelengths]

        for i, named_wl in enumerate(named_wls):

            wl_signal = [x[:,i] for x in signals]

            wl_signal_matrix = np.array(wl_signal).transpose().astype(float)

            self.signal_data_dictionary[named_wl] = wl_signal_matrix
            self.temp_data_dictionary[named_wl]   = temperature_to_kelvin(temperature)


        barycenters = []

        for i, condition in enumerate(conditions):

            # We need to extract the signal data for this condition
            signal_temp      = signals[i]

            barycenter = np.apply_along_axis(get_barycenter, 1, signal_temp,wavelengths = wavelengths)

            barycenters.append(barycenter)

        # Create a matrix from the barycenters
        barycenters = np.array(barycenters).T

        self.signal_data_dictionary['BCM'] = barycenters
        self.temp_data_dictionary['BCM']   = temperature_to_kelvin(temperature)

        named_wls = ['BCM'] + named_wls

        self.conditions          = np.array(conditions)
        self.conditions_original = np.array(conditions)

        self.signals = np.array(named_wls)

        self.full_spectrum  = True

        # Now we try to add the ratio signal
        idx350 = np.argmin(np.abs(wavelengths - 350))
        idx330 = np.argmin(np.abs(wavelengths - 330))

        num_wl = named_wls[idx350]
        den_wl = named_wls[idx330]

        self.create_ratio_signal(num_wl,den_wl)

        return None

    def create_ratio_signal(self,signal_num,signal_den):

        """
        Custom function to create a ratio signal from two signals.

        Args:
            signal_num (str): Numerator signal key in self.signal_data_dictionary.
            signal_den (str): Denominator signal key in self.signal_data_dictionary.

        """

        ratio_signal_name = 'Ratio ' + signal_num + ' / ' + signal_den

        fluo_num = self.signal_data_dictionary[signal_num]
        fluo_den = self.signal_data_dictionary[signal_den]

        fluo_ratio = fluo_num / fluo_den

        self.signal_data_dictionary[ratio_signal_name] = fluo_ratio
        self.temp_data_dictionary[ratio_signal_name] = self.temp_data_dictionary[signal_num]

        # To prevent error when inserting a longer string if the dtype is str
        self.signals = np.asarray(self.signals,dtype=object)

        self.signals = np.insert(self.signals, 0, ratio_signal_name)

        return None

    def set_signal(self,which):

        """
        Set the active signal for analysis by assigning self.fluo and self.temps.

        Args:
            which (str): Signal key present in self.signal_data_dictionary (e.g., '350nm').
        Assign self.fluo and self.temps according to the desired signal
        Returns:
        self.fluo is an n*m matrix where n is the number of measurements and m is the number of different samples
        self.temps is a vector of length n (same n as before)

        Assign self.dt (average value of the temperature step in self.temps),
        and set the signal we are analyzing, i.e., "350nm"
        """

        self.fluo   = self.signal_data_dictionary[which]
        self.temps  = self.temp_data_dictionary[which]

        self.set_signal_type(which)

        return None

    def set_min_max_temp(self):

        """
        Notes:
            Set self.max_temp and self.min_temp as the maximum and minimum of self.temps.
        """

        self.min_temp = np.min(self.temps)
        self.max_temp = np.max(self.temps)

        return None

    def set_dt(self):

        """
        Notes:
            Set self.dt as the average temperature step in self.temps.
        """

        if self.min_temp is None:

            self.set_min_max_temp()

        self.dt = ( self.max_temp - self.min_temp ) / (len(self.temps) - 1)

        return None

    def set_signal_type(self,which):
        
        """
        Set the signal we want to analyze, i.e., "DSF_RFU", "350nm", ...
        """

        self.signal_type = which

        return None

    def filter_by_temperature(self,min_temp,max_temp):

        """
        Filter fluorescence and temperature data to a specified temperature range.

        Args:
            min_temp (float): Minimum temperature threshold.
            max_temp (float): Maximum temperature threshold.
        Filter self.fluo and self.temps according to min_temp and max_temp
        """

        min_temp = temperature_to_kelvin(min_temp)
        max_temp = temperature_to_kelvin(max_temp)

        self.fluo = filter_fluo_by_temp(self.fluo,self.temps,min_temp,max_temp)
        self.temps = filter_temp_by_temp(self.temps,min_temp,max_temp)

        self.set_min_max_temp()

        self.set_dt()

        return None

    def sort_by_conditions_name(self,sort=False):

        """
        Optionally reorder samples according to their names using natural sorting.

        Args:
            sort (bool): If True, sort conditions; if False, revert to original order.
        Reorder the samples according to the sample names - check index_natsorted function from the natsort package
        """

        if sort:

            ind                      = index_natsorted(self.conditions_original)
            self.conditions_original = self.conditions_original[ind]
            self.conditions          = self.conditions_original
            self.signal_data_dictionary = {k: v[:,ind] for k, v in self.signal_data_dictionary.items()}

            # Get idx so we can reverse the sorting.
            self.idx_no_sorted   = [ind.index(i) for i in range(len(ind))]

        else:

            ind                      = self.idx_no_sorted
            self.conditions_original = self.conditions_original[ind]
            self.conditions          = self.conditions_original
            self.signal_data_dictionary = {k: v[:,ind] for k, v in self.signal_data_dictionary.items()}

        return None

    def set_colors(self,color_list):

        """
        Set custom colors for samples.

        Args:
            color_list (list): List of color specifications for each sample.
        Returns:
            None
        Notes:
            Sets self.all_colors
        """

        if len(color_list) == 0:
            return None

        # convert to list if it is not one
        if not isinstance(color_list, (list, np.ndarray)):
            color_list = [color_list]

        if len(color_list) != len(self.conditions_original):
            raise ValueError("Length of color_list must match number of conditions.")

        self.all_colors = color_list

        return None

    def set_conditions(self,condition_list):

        """
        Set custom condition names for samples.

        Args:
            condition_list (list): List of condition names for each sample.
        Returns:
            None
        Notes:
            Sets self.conditions
        """

        # convert to list if it is not one
        if not isinstance(condition_list, (list, np.ndarray)):
            condition_list = [condition_list]

        if len(condition_list) != len(self.conditions_original):
            raise ValueError("Length of condition_list must match number of conditions.")

        self.conditions = condition_list

        return None

    def select_conditions(self,boolean_mask):

        """
        Select a subset of sample columns for analysis.

        Args:
            boolean_mask (array-like): Boolean mask or index list selecting columns.
        Select a subset of samples that we want to analyze

        """

        # Convert to list if it is not one
        if not isinstance(boolean_mask, (list, np.ndarray)):
            boolean_mask = [boolean_mask]

        self.conditions = [x for i,x in enumerate(self.conditions) if boolean_mask[i]]

        if self.all_colors is not None:

            self.colors = [x for i,x in enumerate(self.all_colors) if boolean_mask[i]]

        if self.fluo is not None:

            self.fluo = self.fluo[:,boolean_mask]

        return None

    def median_filter(self,n_degree_window):

        """

        Use this function if the fluorescence curves present spikes. It applies a
        rolling median window filter

        Args:
            n_degree_window (int): Window size in temperature units for median filtering.
        Returns:
            None
        Notes:
            Sets self.fluo
        """

        self.fluo = np.apply_along_axis(
            median_filter_from_fluo_and_temp_vectors,
            0,
            self.fluo,
            self.temps,
            n_degree_window
        )

        return None

    def estimate_fluo_derivates(self,temp_window_length=8):

        """
        Estimate first and second derivatives of fluorescence using Savitzky-Golay.

        Args:
            temp_window_length (int): Approximate window length in temperature units used to compute the filter.
        Returns:
            None
        Notes:
            Sets self.derivative, self.derivative2, and self.tms_from_deriv
        """

        if self.dt is None:
            self.set_dt()

        odd_n_data_points_window_len         = np.ceil(temp_window_length / self.dt) // 2 * 2 + 1
        odd_n_data_points_window_len_2nd_der = np.ceil((temp_window_length+5) / self.dt) // 2 * 2 + 1

        self.derivative = savgol_filter(
            self.fluo,
            axis=0,
            window_length=odd_n_data_points_window_len,
            polyorder=4,
            deriv=1,
            mode="nearest"
        )

        self.derivative2 = savgol_filter(
            self.fluo,
            axis=0,
            window_length=odd_n_data_points_window_len_2nd_der,
            polyorder=4,
            deriv=2,
            mode="nearest"
        )

        der_temp_init = filter_fluo_by_temp(
            self.derivative,
            self.temps,
            self.min_temp+6,
            self.min_temp+11
        )

        der_temp_end  = filter_fluo_by_temp(
            self.derivative,
            self.temps,
            self.max_temp-11,
            self.max_temp-6
        )

        med_init  = np.median(der_temp_init,axis=0)
        med_end   = np.median(der_temp_end ,axis=0)

        mid_value = np.array([(x+y)/2 for x,y in zip(med_init,med_end)])

        der_temp = filter_fluo_by_temp(
            self.derivative,
            self.temps,
            self.min_temp+6,
            self.max_temp-6
        )

        mid_value = mid_value * np.where(mid_value>0,1,-1)

        der_temp = np.add(der_temp,mid_value) 

        temp_temp = filter_temp_by_temp(self.temps,self.min_temp+6,self.max_temp-6)

        max_der = np.amax(der_temp,axis=0)
        min_der = np.amin(der_temp,axis=0)

        der_direction_temp = [abs(maxd) > abs(mind) for maxd,mind in zip(max_der,min_der)]

        if  sum(der_direction_temp) > (len(der_direction_temp) / 2):

            self.tms_from_deriv = get_temp_at_maximum_of_derivative(temp_temp,der_temp)

        else:

            self.tms_from_deriv = get_temp_at_minimum_of_derivative(temp_temp,der_temp)

        return None

    def decompose_spectra(self,min_wl=0,max_wl=800):

        # Raise an error if full_spectrum is False
        if not self.full_spectrum:
            raise ValueError("Cannot decompose spectra if full_spectrum is False.")

        selected_rows = []

        keys = list(self.signal_data_dictionary.keys())

        for key in keys:

            condition1 = 'bcm' in key.lower()
            condition2 = 'ratio' in key.lower()
            condition3 = 'svd' in key.lower()

            if condition1 or condition2 or condition3:
                selected_rows.append(False)
            else:
                wl = float(key.split(' nm')[0])

                selected_rows.append(wl >= min_wl and wl <= max_wl)

        first_coeff_all = []
        second_coeff_all = []

        for id,_ in enumerate(self.conditions_original):

            # We need to have temperature as columns and wavelengths as rows
            signal_2D = generate_2D_signal_matrix(
                id,
                self.signal_data_dictionary,
                selected_rows
            )

            explained_variance, basis_spectra, coefficients = apply_svd(signal_2D)

            coeff_first_basis_spectrum = coefficients[0,:]
            coeff_second_basis_spectrum = coefficients[1,:]

            first_coeff_all.append(coeff_first_basis_spectrum)
            second_coeff_all.append(coeff_second_basis_spectrum)

        first_coeff_matrix = np.array(first_coeff_all).T
        second_coeff_matrix = np.array(second_coeff_all).T

        coeff_matrices = [second_coeff_matrix,first_coeff_matrix]
        coeff_names = ['SVD Second coefficient','SVD First coefficient']

        for svd_signal_name, coeff_matrix in zip(coeff_names, coeff_matrices):

            self.signal_data_dictionary[svd_signal_name] = coeff_matrix
            self.temp_data_dictionary[svd_signal_name]   = self.temp_data_dictionary[keys[0]]

            self.signals = np.asarray(self.signals,dtype=object)

            self.signals = np.insert(self.signals, 0, svd_signal_name)

        return None

    def set_baseline_types(self,poly_order_native=1,poly_order_unfolded=1):

        """
        Args:
            poly_order_native (int): Polynomial order for native baseline fitting (default 1).
            poly_order_unfolded (int): Polynomial order for unfolded baseline fitting (default 1
        Returns:
            None
        Notes:
            Sets self.poly_order_native, self.poly_order_unfolded,
            self.fit_kN, self.fit_kU, self.fit_qN, self.fit_q
        """

        self.poly_order_native   = poly_order_native
        self.poly_order_unfolded = poly_order_unfolded

        self.fit_kN = poly_order_native > 0
        self.fit_kU = poly_order_unfolded > 0

        self.fit_qN = poly_order_native > 1
        self.fit_qU = poly_order_unfolded > 1

        return None

    def estimate_baselines_parameters(self,baseline_degree_window=12):

        """
        Fit linear baselines to native and unfolded portions to estimate slopes & intercepts.

        Args:
            baseline_degree_window (float): Temperature window (in degrees) used to select baseline regions.

        Returns:
            None

        Notes:
            Sets self.bNs, self.bUs, self.kNs, self.kUs and self.tms_initial.
        """

        max_native_temp   = self.min_temp + baseline_degree_window
        min_unfolded_temp = self.max_temp - baseline_degree_window

        native_fluo       = filter_fluo_by_temp(
            self.fluo,
            self.temps,
            self.min_temp,
            max_native_temp
        )

        native_temps      = filter_temp_by_temp(
            self.temps,
            self.min_temp,
            max_native_temp
        )

        unfolded_fluo     = filter_fluo_by_temp(
            self.fluo,
            self.temps,
            min_unfolded_temp,
            self.max_temp
        )

        unfolded_temps    = filter_temp_by_temp(
            self.temps,
            min_unfolded_temp,
            self.max_temp
        )

        bNs = []
        bUs = []

        kNs = []
        kUs = []

        qNs = []
        qUs = []

        bIs = []

        native_temps_shifted = temperature_to_kelvin(native_temps) - temp_standard
        unfolded_temps_shifted = temperature_to_kelvin(unfolded_temps) - temp_standard

        for i in range(native_fluo.shape[1]):

            y_native = native_fluo[:,i]
            y_unfolded = unfolded_fluo[:,i]

            bN, kN, qN = None, None, None
            bU, kU, qU = None, None, None

            if self.poly_order_native == 0:

                bN = np.average(y_native)

                b_start = bN

            if self.poly_order_native == 1:

                kN, bN = fit_line_robust(native_temps_shifted,y_native)

                b_start = bN + kN * np.min(native_temps_shifted)

            if self.poly_order_native == 2:

                qN, kN, bN = fit_quadratic_robust(native_temps_shifted,y_native)

                min_t = np.min(native_temps_shifted)

                b_start = bN + kN * min_t + qN * min_t**2

            if self.poly_order_unfolded == 0:

                bU = np.average(y_unfolded)

                b_end = bU

            if self.poly_order_unfolded  == 1:

                kU, bU = fit_line_robust(unfolded_temps_shifted,y_unfolded)

                b_end = bU + kU * np.max(unfolded_temps_shifted)

            if self.poly_order_unfolded  == 2:

                qU, kU, bU = fit_quadratic_robust(unfolded_temps_shifted,y_unfolded)

                max_t = np.max(unfolded_temps_shifted)

                b_end = bU + kU * max_t + qU * max_t**2

            bNs.append(bN)
            bUs.append(bU)

            kNs.append(kN)
            kUs.append(kU)

            qNs.append(qN)
            qUs.append(qU)

            bI = (b_start + b_end) / 2

            bIs.append(bI)

        bNs = np.array(bNs)
        bUs = np.array(bUs)

        kNs = np.array(kNs)
        kUs = np.array(kUs)

        qNs = np.array(qNs)
        qUs = np.array(qUs)

        bIs = np.array(bIs)

        self.bNs = bNs
        self.bUs = bUs

        self.kNs = kNs
        self.kUs = kUs

        self.qNs = qNs
        self.qUs = qUs

        self.bIs = bIs

        return None

    def initialize_model(self,model_name,params_name):

        """
        Initialize internal storage before performing a fit with a model.

        Args:
            model_name (str): Name of the model (used for bookkeeping).
            params_name (list): List of parameter names corresponding to model parameters.

        Returns:
            None: Initializes lists used to store fit results and errors.
        """

        self.model_name  = model_name   # EquilibriumTwoState, EmpiricalTwoState, ...
        self.params_name = params_name  # i.e. ['kN', 'bN', 'kU', 'bU', 'dHm', 'Tm']

        #  fitted_conditions_indexes is a list of integers: indexes of the conditions that could be 'successfully' fitted
        self.fitted_conditions_indexes  = []
        
        self.params_all                 = [] # values of the fitted parameters
        self.errors_abs_all             = [] # std error of the fitted parameters
        self.errors_percentage_all      = [] # relative error of the fitted parameters
        self.fluo_predictions_all       = [] # Predicted signal
        self.std_error_estimate_all     = [] # std error of the fitting
        self.parameters_far_from_bounds = [] # Boolean, check that the parameters are far from the fitting bounds

        return None

    def fit_model_and_fill_params_and_errors(
            self,
            model_function,
            initial_estimates,
            low_bounds,
            high_bounds,
            fit_algorithm="trf"):

        """
        Fit model_function independently to each column (condition) in self.fluo and store results.

        Args:
            model_function (callable): Function f(T, *params) returning predicted fluorescence.
            initial_estimates (list): List of initial parameter tuples per condition.
            low_bounds (list): List of lower-bound lists per condition.
            high_bounds (list): List of upper-bound lists per condition.
            fit_algorithm (str): Optimization method passed to scipy.curve_fit (default 'trf').

        Returns:
            None: Populates self.params_all, self.errors_abs_all, self.errors_percentage_all,
                  self.fluo_predictions_all, self.std_error_estimate_all and self.fitted_conditions_indexes.
        """

        n = self.fluo.shape[1]

        for index in range(n):

            low_bound     =  low_bounds[index]
            high_bound    =  high_bounds[index]
            p0            =  initial_estimates[index]

            fluo_vec      =  self.fluo[:,index].flatten()

            try:

                params, cov = curve_fit(
                    model_function,
                    self.temps,
                    fluo_vec,
                    p0=p0,
                    bounds=(tuple(low_bound), tuple(high_bound)),
                    method=fit_algorithm
                )

                errors = np.sqrt(np.diag(cov))

                self.params_all.append(params)
                self.errors_abs_all.append(errors)
                self.errors_percentage_all.append(abs(errors / params) * 100)

                fluo_predictions      = model_function(self.temps,*params)

                sum_square_residuals = np.sum( ( fluo_predictions - fluo_vec ) ** 2)
                std_error_estimate   = np.sqrt( sum_square_residuals / len(fluo_vec))

                self.fluo_predictions_all.append(fluo_predictions)
                self.std_error_estimate_all.append(std_error_estimate)
                self.fitted_conditions_indexes.append(index)

            except:

                pass

        if self.fitted_conditions_indexes:

            self.fitted_conditions     = [self.conditions[index] for index in self.fitted_conditions_indexes]
            self.fitted_fluo           =  self.fluo[:,np.array(self.fitted_conditions_indexes)]

        return None


    def equilibrium_two_state(self,fit_algorithm="trf"):

        """
        Fit the thermodynamic equilibrium two-state model to the data.

        Args:
            fit_algorithm (str): Curve-fit algorithm passed to scipy (default 'trf').

        Returns:
            None

        Notes:
            Sets attributes describing the fit (params, errors, dG_std, dCp_component, T_onset).
        """

        model, params_name = get_fit_fx_two_state_signal_fit(
            fit_kN=self.fit_kN,
            fit_kU=self.fit_kU,
            fit_qN=self.fit_qN,
            fit_qU=self.fit_qU,
            type = "equilibrium"
        )

        self.initialize_model("EquilibriumTwoState",params_name)

        init_dH = 80 # kcal/mol
        
        low_bounds        = []
        high_bounds       = []
        initial_estimates = []

        for index in range(self.fluo.shape[1]):

            tm_init = self.tms_from_deriv[index]

            if tm_init > self.max_temp - 3:
                tm_init = self.max_temp - 10

            if tm_init < self.min_temp + 3:
                tm_init = self.min_temp + 10

            p0 = [init_dH, tm_init, self.bNs[index], self.bUs[index]]

            if self.fit_kN:
                p0.append(self.kNs[index])

            if self.fit_kU:
                p0.append(self.kUs[index])

            if self.fit_qN:
                p0.append(self.qNs[index])

            if self.fit_qU:
                p0.append(self.qUs[index])

            low_bound  = [0,  self.min_temp+3] +  [-np.inf for _ in p0[:-2]]
            high_bound = [600,self.max_temp-3] +  [ np.inf for _ in p0[:-2]]

            initial_estimates.append(p0)
            low_bounds.append(low_bound)
            high_bounds.append(high_bound)

        self.fit_model_and_fill_params_and_errors(
            model,initial_estimates,low_bounds,
            high_bounds,fit_algorithm
        )

        dHm_all,      Tm_all     = [p[0] for p in self.params_all], [p[1] for p in self.params_all]

        self.dG_std , self.dCp_component = get_two_state_deltaG(dHm_all, Tm_all, self.cp)
        self.T_onset = get_two_state_tonset(dHm_all, Tm_all)

        return None

    def equilibrium_three_state(self,t1min=25,t1max=90,t2min=25,t2max=90,fit_algorithm="trf"):

        """

        Fit an equilibrium three-state model (N ⇆ I ⇆ U) to the data.

        Args:
            t1min (float) : Lower bound for the T1 parameter
            t1max (float) : Upper bound for the T1 parameter
            t2min (float) : Lower bound for the T2 parameter
            t2max (float) : Upper bound for the T2 parameter
            fit_algorithm (str): Curve fit algorithm passed to scipy (default 'trf').

        Returns:
            None

        Notes:
            Sets fit result attributes including dG_comb_std and T_onset.
        """

        model, params_name = get_fit_fx_three_state_signal_fit(
            fit_kN=self.fit_kN,
            fit_kU=self.fit_kU,
            fit_qN=self.fit_qN,
            fit_qU=self.fit_qU,
            type = "equilibrium"
        )

        self.initialize_model("EquilibriumThreeState",params_name)

        t1min = temperature_to_kelvin(t1min)
        t1max = temperature_to_kelvin(t1max)

        t2min = temperature_to_kelvin(t2min)
        t2max = temperature_to_kelvin(t2max)

        init_dH = 80 # kcal/mol

        low_bounds        = []
        high_bounds       = []
        initial_estimates = []

        temp1_init = (t1min + t1max) / 2
        temp2_init = (t2min + t2max) / 2

        for index in range(self.fluo.shape[1]):

            p0 = [
                init_dH,
                temp1_init,
                init_dH,
                temp2_init,
                self.bNs[index],
                self.bUs[index],
                self.bIs[index]
            ]

            if self.fit_kN:
                p0.append(self.kNs[index])

            if self.fit_kU:
                p0.append(self.kUs[index])

            if self.fit_qN:
                p0.append(self.qNs[index])

            if self.fit_qU:
                p0.append(self.qUs[index])

            low_bound  = [5, t1min,5,t2min]    +  [-np.inf for _ in p0[:-4]]
            high_bound = [600,t1max,600,t2max] +  [ np.inf for _ in p0[:-4]]

            initial_estimates.append(p0)
            low_bounds.append(low_bound)
            high_bounds.append(high_bound)

        self.fit_model_and_fill_params_and_errors(model,initial_estimates,low_bounds,high_bounds,fit_algorithm)

        dHm1_all,      Tm1_all     = [p[0] for p in self.params_all], [p[1] for p in self.params_all]
        dHm2_all,      Tm2_all     = [p[2] for p in self.params_all], [p[3] for p in self.params_all]

        self.T_onset      = get_two_state_tonset(dHm1_all, Tm1_all)
        self.dG_comb_std  = get_eq_three_state_combined_deltaG(dHm1_all, Tm1_all, dHm2_all, Tm2_all)
       
        return None

    def empirical_two_state(self,fit_algorithm="trf",onset_threshold=0.01):

        """

        Fit an empirical two-state model using Tonset to describe transition steepness.

        Args:
            fit_algorithm (str): Optimization method for curve_fit (default 'trf').
            onset_threshold (float): Fraction unfolded that defines Tonset (default 0.01).

        Returns:
            None

        Notes:
            Sets fit attributes and computes an empirical score.
        """

        model, params_name = get_fit_fx_two_state_signal_fit(
            fit_kN=self.fit_kN,
            fit_kU=self.fit_kU,
            fit_qN=self.fit_qN,
            fit_qU=self.fit_qU,
            type = "empirical"
        )

        self.initialize_model("EmpiricalTwoState",params_name)

        low_bounds        = []
        high_bounds       = []
        initial_estimates = []

        for index in range(self.fluo.shape[1]):

            tm_init = self.tms_from_deriv[index]

            if tm_init > self.max_temp - 3:
                tm_init = self.max_temp - 10

            if tm_init < self.min_temp + 6:
                tm_init = self.min_temp + 10

            p0 = [tm_init-5, tm_init, self.bNs[index], self.bUs[index]]

            if self.fit_kN:
                p0.append(self.kNs[index])

            if self.fit_kU:
                p0.append(self.kUs[index])

            if self.fit_qN:
                p0.append(self.qNs[index])

            if self.fit_qU:
                p0.append(self.qUs[index])

            low_bound  = [self.min_temp,  self.min_temp+3] + [-np.inf for _ in p0[:-2]]
            high_bound = [self.max_temp-3,self.max_temp-3] + [ np.inf for _ in p0[:-2]]

            initial_estimates.append(p0)
            low_bounds.append(low_bound)
            high_bounds.append(high_bound)

        self.fit_model_and_fill_params_and_errors(model,initial_estimates,low_bounds,high_bounds,fit_algorithm)

        T_onset_all, Tm_all = [p[0] for p in self.params_all], [p[1] for p in self.params_all]

        self.score = get_empirical_two_state_score(Tm_all, T_onset_all)

        return None

    def empirical_three_state(self,t1min=20,t1max=70,t2min=40,t2max=90,fit_algorithm="trf",onset_threshold=0.01):

        """

        Fit an empirical three-state model using Tonset1/Tonset2 to describe transitions.

        Args:
            t1min (float): lower bound for the T1 parameter
            t1max (float): upper bound for the T1 parameter
            t2min (float): lower bound for the T2 parameter
            t2max (float): upper bound for the T2 parameter
            fit_algorithm (str): Optimization method for curve_fit (default 'trf').
            onset_threshold (float): Fraction unfolded that defines Tonset (default 0.01).

        Returns:
            None

        Notes:
            Sets fit attributes and computes an empirical three-state score.
        """

        model, params_name = get_fit_fx_three_state_signal_fit(
            fit_kN=self.fit_kN,
            fit_kU=self.fit_kU,
            fit_qN=self.fit_qN,
            fit_qU=self.fit_qU,
            type = "empirical"
        )

        self.initialize_model("EmpiricalThreeState",params_name)

        t1min = temperature_to_kelvin(t1min)
        t1max = temperature_to_kelvin(t1max)

        t2min = temperature_to_kelvin(t2min)
        t2max = temperature_to_kelvin(t2max)

        low_bounds        = []
        high_bounds       = []
        initial_estimates = []

        temp1_init = (t1min + t1max) / 2
        temp2_init = (t2min + t2max) / 2

        for index in range(self.fluo.shape[1]):

            p0 = [
                temp1_init-10,
                temp1_init,
                temp2_init-10,
                temp2_init,
                self.bNs[index],
                self.bUs[index],
                self.bIs[index]
            ]

            if self.fit_kN:
                p0.append(self.kNs[index])

            if self.fit_kU:
                p0.append(self.kUs[index])

            if self.fit_qN:
                p0.append(self.qNs[index])

            if self.fit_qU:
                p0.append(self.qUs[index])

            low_bound  = [t1min-8,t1min,t2min-8,t2min] +  [-np.inf for _ in p0[:-4]]
            high_bound = [t1max-8,t1max,t2max-8,t2max] +  [ np.inf for _ in p0[:-4]]

            initial_estimates.append(p0)
            low_bounds.append(low_bound)
            high_bounds.append(high_bound)

        self.fit_model_and_fill_params_and_errors(model,initial_estimates,low_bounds,high_bounds,fit_algorithm)

        T_onset1_all,      T1_all     = [p[0] for p in self.params_all], [p[1] for p in self.params_all]
        T_onset2_all,      T2_all     = [p[2] for p in self.params_all], [p[3] for p in self.params_all]

        self.T_eucl_comb = get_empirical_three_state_score(T_onset1_all, T1_all, T_onset2_all, T2_all)

        return None

    def irreversible_two_state(self,fit_algorithm="trf"):

        """
        Fit an irreversible two-state model using an ODE for fraction native versus temperature.

        Args:
            fit_algorithm (str): Curve-fit algorithm used by scipy (default 'trf').

        Returns:
            None

        Notes:
            Sets fit attributes and computes a pkd-like quantity (self.pkd).
        """

        model, params_name = get_fit_fx_two_state_signal_fit(
            fit_kN=self.fit_kN,
            fit_kU=self.fit_kU,
            fit_qN=self.fit_qN,
            fit_qU=self.fit_qU,
            type = "irreversible"
        )

        self.initialize_model("IrreversibleTwoState",params_name)

        low_bounds        = []
        high_bounds       = []
        initial_estimates = []

        init_Tf = self.min_temp + (self.max_temp - self.min_temp) / 2
        init_Ea = 20 # in kcal/mol

        for index in range(self.fluo.shape[1]):

            p0 = [init_Tf, init_Ea, self.bNs[index], self.bUs[index]]

            if self.fit_kN:
                p0.append(self.kNs[index])

            if self.fit_kU:
                p0.append(self.kUs[index])

            if self.fit_qN:
                p0.append(self.qNs[index])

            if self.fit_qU:
                p0.append(self.qUs[index])

            low_bound  = [self.min_temp-100, 1E-2] + [-np.inf for _ in p0[:-2]]
            high_bound = [self.max_temp+100, 1E5]  + [ np.inf for _ in p0[:-2]]

            initial_estimates.append(p0)
            low_bounds.append(low_bound)
            high_bounds.append(high_bound)

        self.fit_model_and_fill_params_and_errors(model,initial_estimates,low_bounds,high_bounds,fit_algorithm)

        Tf_all, Ea_all = [p[0] for p in self.params_all], [p[1] for p in self.params_all]

        self.pkd = get_irrev_two_state_pkd(Tf_all, Ea_all)

        return None

    def filter_by_relative_error(self,threshold_percentage):

        """
        Filter fitted conditions based on relative error threshold.

        Args:
            threshold_percentage (float): Maximum allowed relative error percentage for all parameters.

        Returns:

            boolean_mask (list): List of booleans indicating which conditions meet the error criteria.

        """

        boolean_mask = []

        for errors_percentage in self.errors_percentage_all:

            condition = all(err <= threshold_percentage for err in errors_percentage)

            boolean_mask.append(condition)

        return boolean_mask

    def filter_by_fitting_std_error(self,threshold_std_error):

        """
        Filter fitted conditions based on fitting standard error threshold.

        Args:
            threshold_std_error (float): Maximum allowed standard error of the fit.

        Returns:

            boolean_mask (list): List of booleans indicating which conditions meet the error criteria.

        """

        boolean_mask = []

        for std_error in self.std_error_estimate_all:

            boolean_mask.append(std_error <= threshold_std_error)

        return boolean_mask

    def filter_by_param_values(self,param_name,low_value,high_value):

        """
        Filter fitted conditions based on parameter value range.

        Args:
            param_name (str): Name of the parameter to filter on.
            low_value (float): Minimum acceptable value for the parameter.
            high_value (float): Maximum acceptable value for the parameter.

        Returns:

            boolean_mask (list): List of booleans indicating which conditions meet the parameter criteria.

        """

        if param_name not in self.params_name:
            raise ValueError(f"Parameter '{param_name}' not found in fitted parameters.")

        param_index = self.params_name.index(param_name)

        boolean_mask = []

        for params in self.params_all:

            param_value = params[param_index]

            condition = (param_value >= low_value) and (param_value <= high_value)

            boolean_mask.append(condition)

        return boolean_mask

    def update_conditions_after_filtering(self,filtered_indexes):

        """
        Update fitted conditions and associated data after filtering.

        Args:
            filtered_indexes (list): List of indexes of conditions to retain.
        Returns:
            None
        Notes:
            Updates self.fitted_conditions, self.fitted_fluo, and all fit result attributes.
        """

        self.fitted_conditions_indexes = [self.fitted_conditions_indexes[i] for i in filtered_indexes]
        self.fitted_conditions        = [self.fitted_conditions[i] for i in filtered_indexes]
        self.fitted_fluo              = self.fitted_fluo[:,filtered_indexes]

        self.params_all               = [self.params_all[i] for i in filtered_indexes]
        self.errors_abs_all           = [self.errors_abs_all[i] for i in filtered_indexes]
        self.errors_percentage_all    = [self.errors_percentage_all[i] for i in filtered_indexes]
        self.fluo_predictions_all     = [self.fluo_predictions_all[i] for i in filtered_indexes]
        self.std_error_estimate_all   = [self.std_error_estimate_all[i] for i in filtered_indexes]

        return None


test = False



if test:

    import matplotlib.pyplot as plt

    mp = DsfFitter()
    mp.load_nano_dsf_xlsx("/home/os/Downloads/demo.xlsx")

    mp.set_signal("350nm")

    n = 3

    bool_mask = [True for _ in range(n)] + [False for _ in range(mp.fluo.shape[1]-n)]

    mp.select_conditions(bool_mask)

    mp.estimate_fluo_derivates()

    mp.cp = 0
    mp.scan_rate = 1

    for combi in [(0,0),(1,1),(2,2)]:

        mp.set_baseline_types(combi[0],combi[1])
        mp.estimate_baselines_parameters()
        mp.empirical_two_state()

        mp.filter_by_relative_error(50)

        for index in range(mp.fitted_fluo.shape[1]):

            plt.plot(mp.temps,mp.fitted_fluo[:,index],label=mp.conditions[index],color='black',alpha=0.3)
            plt.plot(mp.temps,mp.fluo_predictions_all[index],label="fit",color='red')

        plt.legend()
        plt.show()

    quit()
