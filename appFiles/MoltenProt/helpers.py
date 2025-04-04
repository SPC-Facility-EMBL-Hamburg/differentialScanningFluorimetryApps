import numpy  as np
import codecs
import pandas as pd
import codecs
from openpyxl import load_workbook
from xlrd     import open_workbook

from collections import Counter

def temperature_to_kelvin(temp):
    """
    Convert temperature to Kelvin
    """

    if (np.max(temp) < 273.15):
        # If the temperature is in Celsius
        return temp + 273.15
    else:
        # If the temperature is already in Kelvin
        return temp

def get_sheet_names_of_xlsx(filepath):

    """

    Get the sheet names of a xls or xlsx file without loading it. 
    The open_workbook function is used so we can handle the error "openpyxl does not support the old .xls file format" 

    """

    try:

        wb = load_workbook(filepath, read_only=True, keep_links=False)
        sheetnames = wb.sheetnames
    except:

        xls = open_workbook(filepath, on_demand=True)
        sheetnames =  xls.sheet_names()
 
    return sheetnames

def file_is_of_type_uncle(xlsx_file):

    """
    Check if the file is an uncle file
    """

    try:

        # Read the first sheet of the xlsx file
        data = pd.read_excel(xlsx_file,skiprows=1,nrows=5,header=None)

        # Extrac the first row
        row = data.iloc[0,1:]

        # concatenate the row into a string
        row_str = ' '.join([str(x) for x in row])

        # count the number of times we have word 'Time' and 'Temp'
        count_time = row_str.count('Time')
        count_temp = row_str.count('Temp')

        if count_time > 20 and count_temp > 20 and count_temp == count_time:
            return True
        else:
            return False
    except:
        return False

def detect_encoding(file_path):
    encodings = ["utf-8", "latin1", "iso-8859-1", "cp1252"]
    for enc in encodings:
        try:
            with codecs.open(file_path, encoding=enc, errors="strict") as f:
                f.read()
            return enc
        except UnicodeDecodeError:
            continue
    return "Unknown encoding"

def isBlank (myString):

    return (not (myString and myString.strip())) or myString == "nan" 

def find_indexes_of_non_signal_conditions(signal_data,conditions):

    # Find the indexes of conditions with derivative in it, or 'S.D.' in it
    idx_to_remove1 = [i for i, cond in enumerate(conditions) if 'derivative' in cond.lower() or 's.d.' in cond.lower()]

    # Find indexes empty signal_data columns
    idx_to_remove2 = [i for i in range(signal_data.shape[1]) if np.isnan(signal_data[:, i]).all()]

    # Find columns where all values are NaN
    nan_columns = np.all(np.isnan(signal_data), axis=0)

    # Get the column indices
    idx_to_remove3 = np.where(nan_columns)[0].tolist()

    # Find unique indexes to remove
    idx_to_remove = list(set(idx_to_remove1 + idx_to_remove2 + idx_to_remove3))

    return idx_to_remove

def find_repeated_words(string_lst):

    """
    Given a list of strings, find the repeated words in the list

    Args:

        string_lst (lst): list of strings

    Returns:

        repeated_words (lst): list of repeated words

    """

    # Repeated words
    words = " ".join(string_lst).split()

    # Count the occurrences of each word
    word_counts = Counter(words)

    # Find words that are repeated (appear, at least, as many times as elements in string_lst)
    repeated_words = [word for word, count in word_counts.items() if count >= len(string_lst)]

    return repeated_words

def remove_words_in_string(input_string,word_list):

    # Split the string into words
    words = input_string.split()

    # Replace words from the list with an empty string
    filtered_words = [word for word in words if word not in word_list]

    # Join the words back into a single string
    output_string = ' '.join(filtered_words)

    return output_string

def subset_panta_data(data,columns_positions):

    dfs    = [pd.DataFrame(np.array(data.iloc[:,[pos-1,pos]]),
        columns = ["temperature","signal"+str(i)]) for i,pos in enumerate(columns_positions)]

    # Case N = 1 (1 capillary)
    if len(dfs) == 1:

        mat   = dfs[0].to_numpy(dtype='float')

        fluo = mat[:,1:]         # You must use '1:' and not '1' to return a 2D array!
        temp = mat[:,0] + 273.15 # To kelvin 

        return fluo, temp

    # Combine dataframes so we can obtain a vector of temperatures and a matrix of fluorescence signal
    merged = pd.merge_asof(dfs[0].dropna(),dfs[1].dropna(),on="temperature", direction='nearest',allow_exact_matches=True)

    for df in dfs[2:]:
        
        merged = pd.merge_asof(merged,df.dropna(),on="temperature", direction='nearest',allow_exact_matches=True)       
    
    fluo   = np.array(merged.iloc[:, 1:]).astype('float')
    temp   = np.array(merged.iloc[:, 0]).astype('float') + 273.15 # To kelvin 

    # Reduce data so we can plot and fit the data faster.

    while len(temp) > 700:
        fluo = fluo[::2]
        temp = temp[::2]

    return fluo, temp

def getStartLineQuantStudioFile(file_name):

    """
    Returns the number of the first line not starting with "*" + 1
    """

    with codecs.open(file_name, 'r', encoding='utf-8',errors='ignore') as rf:
        ls       = rf.read().splitlines()
        splitted = [l.split() for l in ls]
        
        for i,s in enumerate(splitted):
            if len(s) > 5 and "*" not in s[0]:
                start_row = i+1
                break

    return(start_row)

def generateQS3DataFrame(df):

    df = pd.DataFrame(np.array(df.iloc[:,[3,4]]),
        columns = ["temperature","signal"],dtype=str)

    df.replace(to_replace=',', value='',inplace=True,regex=True)

    df[["temperature", "signal"]] = df[["temperature", "signal"]].apply(pd.to_numeric)

    df.sort_values('temperature',inplace=True)

    return df

def generateMergedQuantStudioDataFrame(data):

    dfs = [generateQS3DataFrame(x) for _, x in data.groupby(1)]
    dfs = [df.rename(columns={"signal": 'signal'+str(i)}) for i,df in enumerate(dfs)]

    # Combine dataframes so we can obtain a vector of temperatures and a matrix of fluorescence signal
    merged = pd.merge_asof(dfs[0].dropna(),dfs[1].dropna(),on="temperature", direction='nearest',allow_exact_matches=True)

    for df in dfs[2:]:
        
        merged = pd.merge_asof(merged,df.dropna(),on="temperature", direction='nearest',allow_exact_matches=True)       

    fluo   = np.array(merged.iloc[:, 1:]).astype('float')
    temp   = np.array(merged.iloc[:, 0]).astype('float') + 273.15 # To kelvin 

    # Reduce data so we can plot and fit the data faster.

    nConditions = fluo.shape[1]

    if nConditions > 192:
        maxPoints = 150
    elif nConditions > 96:
        maxPoints = 300
    else:
        maxPoints = 700

    while len(temp) > maxPoints:
        fluo = fluo[::2]
        temp = temp[::2]

    return fluo, temp

def detect_txt_file_type(file):

    '''

    Decide if the file was generated by an Agilent's MX3005P qPCR instrument or a 
    QuantStudio qPCR instrument

    '''

    with codecs.open(file, 'r', encoding='utf-8',errors='ignore') as rf:
        ls       = rf.read().splitlines()

        for line in ls:
            if line.startswith('Segment') and 'Well' in line:
                return 'MX3005P'

            splittedLine = line.split()
            if 'Well' in splittedLine and 'Target' in splittedLine and 'Reading' in splittedLine:
                return 'QuantStudio'

    return None

def load_layout(xls_file):

    xls = pd.ExcelFile(xls_file)
    dat = pd.read_excel(xls)

    conditions =  []

    for row in range(dat.shape[0]):

        condition = dat.iloc[row,1:].values
        condition = ' '.join([str(x) for x in condition])
        condition = condition.replace('nan','')
        condition = " ".join(condition.split())

        conditions.append(condition)

    return np.array(conditions)

def filter_fluo_by_temp(np_fluo,np_temp,min_temp,max_temp):

	"""
	Filter the fluorescence matrix using a temperature range

	Requires: 
		
		1) The fluorescence matrix 'np_fluo' of dimensions
	n*m where n is the number of temperatures at which the fluorescence 
	was measured and m the number of conditions  
		2) The temperature vector 'np_temp' of length n.
		3) The lower bound 'min_temp'
		4) The upper bound 'max_temp'


	Returns the filtered fluorescence matrix 
	"""

	shape  = np_fluo.shape
	# Convert 1D array to a 2D numpy array of n rows and 1 column
	if len(shape) == 1:
		np_fluo = np.reshape(np_fluo, (-1, 1))

	tot = np_fluo.shape[1]

	np_tog = np.column_stack((np_fluo, np_temp))
	np_tog = np_tog[np_tog[:,tot ] >= min_temp]
	np_tog = np_tog[np_tog[:,tot] <= max_temp]
	np_tog = np.delete(np_tog, tot, 1)

	return(np_tog)

def filter_temp_by_temp(np_temp,min_temp,max_temp):

	"""
	Filter the temperature vector using a temperature range

	Requires: 
		
		1) The temperature vector 'np_temp' of length n.
		2) The lower bound 'min_temp'
		3) The upper bound 'max_temp'


	Returns the filtered temperature vector 
	"""

	np_temp_filter = np_temp[np_temp >= min_temp]
	np_temp_filter = np_temp_filter[np_temp_filter <= max_temp]

	return(np_temp_filter)


def median_filter_from_fluo_and_temp_vectors(fluo_vec,temp_vec,rolling_window):

    """

    Compute the median filter of the fluorescence vector using a temperature rolling window

    First, we convert the temperature into an integer vector and then 
        into time variable to take advantage of pandas function 
            rolling().median() 


	Requires: 

		1) The fluorescence 1D vector 'fluo_vec'
		2) The temperature  1D vector 'temp_vec'
		3) The size of the rolling window in celsius 'rolling_window'

	Returns the fluorescence vector passed through the median filter

    """

    scaling_factor = 10000

    temp_vec     =  np.multiply(temp_vec,scaling_factor).astype(int) 
    series       =  pd.Series(fluo_vec,index=temp_vec,dtype=float)
    series.index =  pd.to_datetime(series.index,unit='s')

    roll_window  = str(int(rolling_window*scaling_factor))+"s"

    fluo_df_median_filt = series.rolling(roll_window).median().to_numpy()

    return fluo_df_median_filt

def get_temp_at_maximum_of_derivative(temps,fluo_derivative):

    tms_derivative = np.take(temps, np.argmax(fluo_derivative,axis=0)) 
    return tms_derivative

def get_temp_at_minimum_of_derivative(temps,fluo_derivative):

    tms_derivative = np.take(temps, np.argmin(fluo_derivative,axis=0)) 
    return tms_derivative

def estimate_initial_tm_from_baseline_fitting(bUs,bNs,kNs,kUs,temps,fit_length,fluo_derivative):


    """
    Use an heuristic algorithm to get an initial estimate of the Tms

    Requires: 

        The intercept of the baseline fitting for the native and unfolded states 
            'bNs' 'bUs' (1D numpy arrays of length K where K is the number of conditions)

        The slope of the baseline fitting for the native and unfolded states 
            'kNs' 'kUs' (1D numpy arrays of length K where K is the number of conditions)
        
        The 1D temperature vector 
            'temps'

        A scalar 'fit_length' (The temperature range used to estimate 'bNs' 'bUs' 'kNs' 'kUs')
    
        The first derivative of the fluorescence 
            'fluo_derivative' (numpy array of dimension m*n where m is the same dimension
            as the temperature vector and n is the number of conditions)
                
    Returns:

        A 1D vector of Tms (numpy array of length K where K is the number of conditions)

    """

    dintersect = - (bUs - bNs) / (kNs - kUs)
    
    tmin = min(temps)
    tmax = max(temps)
    tmid = (tmin + tmax) / 2

    c1 = dintersect > (tmin+fit_length)
    c2 = dintersect < (tmax+fit_length)
    c3 = np.logical_and(c1,c2)

    tms_case1 = np.where(c3,tmid,0)

    c4 = np.logical_not(c3)

    b_diff = (kUs - kNs ) * tmid + (bUs - bNs)

    c5 = b_diff > 0
    c6 = np.logical_and(c4,c5)

    tms_case2 = get_temp_at_maximum_of_derivative(temps,fluo_derivative) * c6
    
    c7 = np.logical_not(c5)
    c8 = np.logical_and(c7,c4)

    tms_case3 = np.take(temps, np.argmin(fluo_derivative,axis=0)) * c8

    tms = tms_case1 + tms_case2 + tms_case3

    return tms

def check_good_parameters_estimation(params,low_bound,high_bound,params_name):

    """
        
    Check that the estimated parameters are far from the bounds, 
    if not there is a problem with the fitting.

    For the first  4 params 'kN', 'bN', 'kU', 'bU' we will normalize low_bound - high_bound range to 0-1
    and then verify that the parameters lie in the interval 0.02-0.98.
        
    """

    low_bound  = np.array(low_bound)
    high_bound = np.array(high_bound)
    params     = np.array(params)

    params_normalized = (params[:4] - low_bound[:4] ) / (high_bound[:4] - low_bound[:4])

    lie_in_correct_interval   = np.logical_and(params_normalized < 0.98,
        params_normalized > 0.02)

    """

    For Tm or T1 and T2 and T_onset or T_onset1 and T_onset2  we will check that they are 1 degree from the boundaries

    """

    for temp_param_name in ["Tm","T1","T2","T_onset","T_onset1","T_onset2"]:

        if temp_param_name in params_name:

            position_of_tm = params_name.index(temp_param_name) 
            tm             = params[position_of_tm]
            tm_bound_low   = low_bound[position_of_tm]
            tm_bound_up    = high_bound[position_of_tm]

            lie_in_correct_interval = np.append(lie_in_correct_interval, np.array([tm > (tm_bound_low+1) and tm < (tm_bound_up-1)]))

    """

    For Tm or dHm, dHm1 and dHm2 we will check that they are between 2.5 and 750 kcal/mol. 10466 and 3098230 in Joules

    """

    for dh_param_name in ["dHm","dHm1","dHm2"]:

        if dh_param_name in params_name:

            position_of_dh = params_name.index(dh_param_name) 
            dh             = params[position_of_dh]
            lie_in_correct_interval = np.append(lie_in_correct_interval, np.array([dh > 10466 and dh < 3098230]))

    """

    For Ki we will check that it lies between 1e-3 and 1e3

    """

    if "Ki" in params_name:
        position_of_ki = params_name.index("Ki") 
        ki             = params[position_of_ki]
        lie_in_correct_interval = np.append(lie_in_correct_interval, np.array([ki >= 0.001 and ki < 1000]))


    return all(lie_in_correct_interval)

def estimate_baseline_factor(kN_fit, bN_fit, kU_fit, bU_fit, Tm_fit,std_error_estimate):

    """
    Use the fluorescence dependence of the folded and unfolded states to estimate
    the baseline separation factor
    """

    baseline_factor = 1 - 6 * std_error_estimate / abs(kU_fit*Tm_fit + \
        bU_fit - (kN_fit * Tm_fit + bN_fit))

    baseline_factor = np.round(baseline_factor,2)

    return(baseline_factor)

def get_EquilibriumTwoState_model_dGstd(dHms,tms,cp,temp_standard=298.15):


    """

    Obtain 'dG_std' - Gibbs free energy of unfolding at standard temperature (298 K), 
        extrapolated using the values of dCp

    dG_std is extrapolated to standard temperature using the model described in Becktel and Schellman, 1987
    Tm is chosen as the reference temperature for equation (4), which also means that dS(Tm) = dH(Tm)/Tm

    """
    dHms = np.array(dHms)
    tms  = np.array(tms)

    dG_std = dHms * (1 - temp_standard / tms) - cp * (tms - temp_standard + temp_standard * np.log(temp_standard / tms))

    # Calculate contribution of dCp to the real dG 
        
    dCp_component = temp_standard - tms - temp_standard * np.log(temp_standard / tms)

    return dG_std, dCp_component

def get_EquilibriumTwoState_model_Tons(dHms,tms,dHms_sd,tms_sd,onset_threshold=0.01,R_gas_constant=8.314):

    """

    Compute onset temperature (and stdev using error propagation) 
        based on fitted dHm and Tm
    
    onset_threshold is the fraction unfolded that corresponds to onset of unfolded.
        The default value means that (0.01 - 1% must be unfolded)

    """

    dHms,    tms        = np.array(dHms),    np.array(tms) 
    dHms_sd, tms_sd     = np.array(dHms_sd), np.array(tms_sd)

    T_onset = 1 / (1 / tms - R_gas_constant / dHms * np.log(onset_threshold / (1 - onset_threshold))) 

    return T_onset

def get_EmpiricalTwoState_model_ranking(Tm,T_onset):

    """

    Imagine a two dimensional space with coordinates 
        x1 = Tmelting and 
        x2 = Tonset (where 1% of the protein is unfolded)

    Return the distance from the origin (0,0) to the point (Tm,Tonset)

    """

    return np.sqrt(np.array(Tm)**2 + np.array(T_onset)**2)

def get_EquilibriumThreeState_model_dG_comb_std(dHm1_fit,T1_fit,dHm2_fit,T2_fit,temp_standard=298.15):

    dHm1_fit,   dHm2_fit    = np.array(dHm1_fit), np.array(dHm2_fit)
    T1_fit,     T2_fit      = np.array(T1_fit),   np.array(T2_fit)

    dG_comb_std = dHm1_fit * (1 - temp_standard / T1_fit) + dHm2_fit * (1 - temp_standard / T2_fit)

    return dG_comb_std

def get_EmpiricalThreeState_model_T_eucl_comb(T_onset1_fit,T1_fit,T_onset2_fit,T2_fit):

    """

    Return the distance from the origin (0,0) to the point (T1_fit,Tonset1) +
           the distance from the origin (0,0) to the point (T2_fit,Tonset2) 

    """

    T_onset1_fit,   T_onset2_fit    = np.array(T_onset1_fit),   np.array(T_onset2_fit)
    T1_fit,         T2_fit          = np.array(T1_fit),         np.array(T2_fit)

    T_eucl_comb = np.sqrt(T1_fit**2 + T_onset1_fit**2) + np.sqrt(T2_fit**2 + T_onset2_fit**2)

    return T_eucl_comb

def arrhenius(T, Tf, Ea,R_gas_constant=8.314):
    """
    Arrhenius equiation: defines dependence of reaction rate constant k on temperature
    In this version of the equation we use Tf (a temperature of k=1)
    to get rid of instead of pre-exponential constant A
    """
    return np.exp(-Ea / R_gas_constant * (1 / T - 1 / Tf))

def get_IrrevTwoState_pkd(Tf, Ea,temp_standard=298.15):

	Tf = np.array(Tf)
	Ea = np.array(Ea)

	return -np.log10(arrhenius(temp_standard,Tf,Ea))

