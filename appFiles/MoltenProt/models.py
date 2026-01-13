# Collection of models  to fit thermal unfolding data
import numpy as np

from constants import temp_standard, R_gas
from helpers   import temperature_to_kelvin

from scipy.integrate  import solve_ivp


def baseline_signal(T,b,k=0,q=0):

    """
    For the dependence of the pre- and post-transition states.
    Args:
        T (array-like): temperature in Celsius or Kelvin
        b (float): intercept
        k (float): slope
        q (float): quadratic term

    Returns:
        array-like: signal at temperature T

    """

    dT = temperature_to_kelvin(T) - temp_standard

    return b + k*dT + q*(dT**2)

def eq_constant(T,DH1,T1,Cp=0):

    """
    Args:
    T (array-like): temperature in Celsius or Kelvin
    DH1 (float): enthalpy change at T1 in kcal/mol
    T1 (float): temperature at which ΔG(T) = 0 in Celsius or Kelvin
    Cp (float): change in heat capacity between the two states in kcal/(mol*K)
    """

    T  = temperature_to_kelvin(T)
    T1 = temperature_to_kelvin(T1)

    DG = DH1*(1 - T/T1) - Cp*(T1 - T + T*np.log(T/T1)) 
    K  = np.exp(-DG / (R_gas * T))

    return K

def eq_two_state_signal(T,dh,Tm,bN,bU,kN=0,kU=0,qN=0,qU=0,Cp=0):

    """
    Signal of a two-state unfolding model  N ⇔ U

    Args:
        T (array-like): temperature in Celsius or Kelvin
        DH (float): enthalpy change at Tm in kcal/mol
        Tm (float): melting temperature in Celsius or Kelvin
        bN (float): intercept of the native baseline
        bU (float): intercept of the unfolded baseline
        kN (float): slope of the native baseline
        kU (float): slope of the unfolded baseline
        qN (float): quadratic term of the native baseline
        qU (float): quadratic term of the unfolded baseline
        Cp (float): change in heat capacity between the two states in kcal/(mol*K)

    Returns:
        array-like: signal at temperature T
    """

    K   = eq_constant(T,dh,Tm,Cp)
    fn  = (1/(1 + K))
    fu  = 1 - fn

    y = fn*baseline_signal(T,bN,kN,qN) + fu*baseline_signal(T,bU,kU,qU)

    return y

def empirical_two_state_signal(T,Tonset,Tm,bN,bU,kN=0,kU=0,qN=0,qU=0):

    """
    Signal of a two-state unfolding model  N ⇔ U

    Args:
        T (array-like): temperature in Celsius or Kelvin
        Tonset (float): temperature at which unfolding starts in Celsius or Kelvin
        Tm (float): melting temperature in Celsius or Kelvin
        bN (float): intercept of the native baseline
        bU (float): intercept of the unfolded baseline
        kN (float): slope of the native baseline
        kU (float): slope of the unfolded baseline
        qN (float): quadratic term of the native baseline
        qU (float): quadratic term of the unfolded baseline
        Cp (float): change in heat capacity between the two states in kcal/(mol*K)

    Returns:
        array-like: signal at temperature T
    """

    Tm = temperature_to_kelvin(Tm)
    Tonset = temperature_to_kelvin(Tonset)
    T = temperature_to_kelvin(T)

    """
    Previous implementation:
        delta_g_2 = -R_gas * Tonset * np.log(0.01/0.99)

        delta_g = (Tm - T) * -delta_g_2 / (Tonset - Tm)

        K = np.exp(-delta_g / (R_gas * T))
    """

    # New simplified implementation
    exp_term = (Tm - T) * (Tonset * 4.595119) / ((Tonset - Tm) * T)

    K = np.exp(exp_term)

    fn = (1/(1 + K))
    fu = 1 - fn

    y = fn*baseline_signal(T,bN,kN,qN) + fu*baseline_signal(T,bU,kU,qU)

    return y

def arrhenius(T, Tf, Ea):

    """
    Arrhenius equation: defines dependence of reaction rate constant k on temperature
    In this version of the equation we use Tf (a temperature of k=1)
    to get rid of instead of pre-exponential constant A

    T, Tf, must be given in Kelvin, Ea in kcal units

    Args:
        T (array-like): temperature in Celsius or Kelvin
        Tf (float): reference temperature in Celsius or Kelvin
        Ea (float): activation energy in kcal/mol

    Returns:
        array-like: Arrhenius factor at temperature T

    """

    T  = temperature_to_kelvin(T)
    Tf = temperature_to_kelvin(Tf)

    return np.exp(-Ea / R_gas * (1 / T - 1 / Tf))


def ode_irrev(T, xn, Tf, Ea, scan_rate_v):

    """
    ordinary differential equation for the native fraction native versus temperature
    dxn/dT = -1/v*k(T)*xn

    Tf is the temperature at which k = 1/min

    Args:
        T (array-like): temperature in Celsius or Kelvin
        xn (array-like): fraction of native protein
        Tf (float): reference temperature in Celsius or Kelvin
        Ea (float): activation energy in kcal/mol
        scan_rate_v (float): scan rate in K/min
    Returns:
        array-like: derivative of native fraction with respect to temperature
    """

    return -1 / scan_rate_v * arrhenius(T, Tf, Ea) * xn

def irrev_signal(T, Tf,Ea,bN,bF,kN,kF,qN,qF,scan_rate_v):

    """
    Nomeclature from Jose M. Sanchez-Ruiz 1992
    T,t1, Tf should be in Kelvin!
    The variable
        'xn' refers to the folded  state
    """

    ivp_result = solve_ivp(
        ode_irrev,
        t_span=[np.min(T), np.max(T)],
        t_eval=T,
        y0=[1], # We assume that at the beginning all protein is in the native state
        args=(Tf,Ea,scan_rate_v),
        method="LSODA"
    )

    xn = ivp_result.y[0, :]

    return xn*baseline_signal(T,bN,kN,qN) + (1 - xn)*baseline_signal(T,bF,kF,qF)

def get_irrev_signal_fixed_scan_rate(scan_rate_v):

    """
    Returns a function to calculate the irreversible unfolding signal with a fixed scan rate.
    Args:
        scan_rate_v (float): scan rate in K/min
    Returns:
        function: function to calculate the irreversible unfolding signal
    """

    def fx(T, Tf,Ea,bF,bN,kN,kF,qN,qF):

        return irrev_signal(T, Tf,Ea,bF,bN,kN,kF,qN,qF,scan_rate_v)

    return fx

def get_fit_fx_two_state_signal_fit(
        fit_kN,fit_kU,fit_qN,fit_qU,
        type="equilibrium",
        scan_rate_v=1):

    """
    Returns a function to fit a two-state unfolding model with fixed or variable baselines.

    Args:
        fit_kN (bool): whether to fit the slope of the native baseline
        fit_kU (bool): whether to fit the slope of the unfolded baseline
        fit_qN (bool): whether to fit the quadratic term of the native baseline
        fit_qU (bool): whether to fit the quadratic term of the unfolded baseline
        type (str): type of the model to fit ("equilibrium", "empirical", "irreversible")
        scan_rate_v (float): scan rate in K/min (only used for irreversible model)

    Returns:
        function: function to fit a two-state unfolding model
        params_names (list): list of parameter names to fit
    """

    if type == "equilibrium":

        params_names = ["DH","Tm","bN","bU"]

    elif type == "empirical":

        params_names = ["T_onset","Tm","bN","bU"]

    elif type == "irreversible":

        params_names = ["Tf","Ea","bN","bU"]

    if fit_kN:
        params_names.append("kN")

    if fit_kU:
        params_names.append("kU")

    if fit_qN:
        params_names.append("qN")

    if fit_qU:
        params_names.append("qU")

    def fx(T,first_param,second_param,bN,bU, *args):

        idx = 0

        if fit_kN:
            kN = args[idx]
            idx += 1
        else:
            kN = 0

        if fit_kU:
            kU = args[idx]
            idx += 1
        else:
            kU = 0

        if fit_qN:
            qN = args[idx]
            idx += 1
        else:
            qN = 0

        if fit_qU:
            qU = args[idx]
            idx += 1
        else:
            qU = 0

        if type == "equilibrium":

            return eq_two_state_signal(T,first_param,second_param,bN,bU,kN,kU,qN,qU)

        elif type == "empirical":

            return empirical_two_state_signal(T,first_param,second_param,bN,bU,kN,kU,qN,qU)

        elif type == "irreversible":

            irrev_two_state_signal = get_irrev_signal_fixed_scan_rate(scan_rate_v)  # assuming scan rate of 1 K/min for fitting

            return irrev_two_state_signal(T,first_param,second_param,bN,bU,kN,kU,qN,qU)

    return fx, params_names

def eq_three_state_signal(T,DH1,T1,DH2,T2,bN,bU,bI,kN=0,kU=0,qN=0,qU=0,Cp1=0,CpTh=0):

    """
    Signal of a two-state unfolding model  N ⇔ I ⇔ U

    Args:
        T (array-like): temperature in Celsius or Kelvin
        DH1 (float): enthalpy change at T1 in kcal/mol
        T1 (float): melting temperature of first transition in Celsius or Kelvin
        DH2 (float): enthalpy change at T2 in kcal/mol
        T2 (float): melting temperature of second transition in Celsius or Kelvin
        bN (float): intercept of the native baseline
        bU (float): intercept of the unfolded baseline
        bI (float): intercept of the intermediate baseline
        kN (float): slope of the native baseline
        kU (float): slope of the unfolded baseline
        qN (float): quadratic term of the native baseline
        qU (float): quadratic term of the unfolded baseline
        Cp (float): change in heat capacity between the two states in kcal/(mol*K)

    Returns:
        array-like: signal at temperature T
    """

    A = eq_constant(T,DH1,T1,Cp1)
    B = eq_constant(T,DH2,T2,CpTh-Cp1)

    den = (1+A+A*B)

    fn, fi, fu = 1/den, A/den, A*B/den

    y = fn*baseline_signal(T,bN,kN,qN) + fu*baseline_signal(T,bU,kU,qU) + fi*bI

    return y

def empirical_three_state_signal(T,Tonset_1,T1,Tonset_2,T2,bN,bU,bI,kN=0,kU=0,qN=0,qU=0):

    """
    Signal of a two-state unfolding model  N ⇔ I ⇔ U

    Args:
        T (array-like): temperature in Celsius or Kelvin
        Tonset_1 (float): temperature at which unfolding starts in Celsius or Kelvin
        T1 (float): melting temperature of first transition in Celsius or Kelvin
        Tonset_2 (float): temperature at which second unfolding starts in Celsius or Kelvin
        T2 (float): melting temperature of second transition in Celsius or Kelvin
        bN (float): intercept of the native baseline
        bU (float): intercept of the unfolded baseline
        bI (float): intercept of the intermediate baseline
        kN (float): slope of the native baseline
        kU (float): slope of the unfolded baseline
        qN (float): quadratic term of the native baseline
        qU (float): quadratic term of the unfolded baseline
        Cp (float): change in heat capacity between the two states in kcal/(mol*K)

    Returns:
        array-like: signal at temperature T
    """

    T1 = temperature_to_kelvin(T1)
    T2 = temperature_to_kelvin(T2)
    Tonset_1 = temperature_to_kelvin(Tonset_1)
    Tonset_2 = temperature_to_kelvin(Tonset_2)

    T  = temperature_to_kelvin(T)

    """
    Previous implementation:

    delta_g_N_I_at_onset = -R_gas * Tonset_1 * np.log(0.01/0.99)
    delta_g_I_U_at_onset = -R_gas * Tonset_2 * np.log(0.01/0.99)

    delta_g_N_I = (T1 - T) * -delta_g_N_I_at_onset / (Tonset_1 - T1)
    delta_g_U_I = (T2 - T) * -delta_g_I_U_at_onset / (Tonset_2 - T2)

    K1 = np.exp(-delta_g_N_I / (R_gas * T))
    K2 = np.exp(-delta_g_U_I / (R_gas * T))
    """

    # New simplified implementation
    exp_term_1 = (T1 - T) * (Tonset_1 * 4.595119) / ((Tonset_1 - T1)*T)
    exp_term_2 = (T2 - T) * (Tonset_2 * 4.595119) / ((Tonset_2 - T2)*T)

    K1 = np.exp(exp_term_1)
    K2 = np.exp(exp_term_2)

    den = (1+K1+K1*K2)

    fn, fi, fu = 1/den, K1/den, K1*K2/den

    y = fn*baseline_signal(T,bN,kN,qN) + fu*baseline_signal(T,bU,kU,qU) + fi*bI

    return y

def get_fit_fx_three_state_signal_fit(fit_kN,fit_kU,fit_qN,fit_qU,type="equilibrium"):

    """
    Returns a function to fit a two-state unfolding model with fixed or variable baselines.

    Args:
        fit_kN (bool): whether to fit the slope of the native baseline
        fit_kU (bool): whether to fit the slope of the unfolded baseline
        fit_qN (bool): whether to fit the quadratic term of the native baseline
        fit_qU (bool): whether to fit the quadratic term of the unfolded baseline

    Returns:
        function: function to fit a two-state unfolding model
        params_names (list): list of parameter names to fit
    """

    if type == "equilibrium":

        params_names = ["DH1","T1","DH2","T2","bN","bU","bI"]

    elif type == "empirical":

        params_names = ["T_onset1","T1","T_onset2","T2","bN","bU","bI"]

    if fit_kN:
        params_names.append("kN")

    if fit_kU:
        params_names.append("kU")

    if fit_qN:
        params_names.append("qN")

    if fit_qU:
        params_names.append("qU")

    def fx(T,param1,param2,param3,param4,bN,bU,bI,*args):

        idx = 0

        if fit_kN:
            kN = args[idx]
            idx += 1
        else:
            kN = 0

        if fit_kU:
            kU = args[idx]
            idx += 1
        else:
            kU = 0

        if fit_qN:
            qN = args[idx]
            idx += 1
        else:
            qN = 0

        if fit_qU:
            qU = args[idx]
            idx += 1
        else:
            qU = 0

        if type == "equilibrium":

            return eq_three_state_signal(T,param1,param2,param3,param4,bN,bU,bI,kN,kU,qN,qU)

        elif type == "empirical":

            return empirical_three_state_signal(T,param1,param2,param3,param4,bN,bU,bI,kN,kU,qN,qU)

    return fx, params_names