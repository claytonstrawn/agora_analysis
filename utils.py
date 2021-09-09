import numpy as np

def convert_T_over_mu_to_T_func():
    def calc_mu_table_local(temperature):
        tt = np.array([1.0e+01, 1.0e+02, 1.0e+03, 1.0e+04, \
                       1.3e+04, 2.1e+04, 3.4e+04, 6.3e+04, \
                       1.0e+05, 1.0e+09])
        mt = np.array([1.18701555, 1.15484424, 1.09603514, 0.9981496, \
                       0.96346395, 0.65175895, 0.6142901, 0.6056833, \
                       0.5897776, 0.58822635])
        logttt= np.log(temperature)
        # linear interpolation in log-log space
        logmu = np.interp(logttt,np.log(tt),np.log(mt)) 
        return np.exp(logmu)

    temperature_values = []
    mu_values = []
    T_over_mu_values = []
    current_temperature = 1e1
    final_temperature = 1e9
    dlogT = 1.
    while current_temperature < final_temperature:
        temperature_values.append(current_temperature)
        current_mu = calc_mu_table_local(current_temperature)
        mu_values.append(current_mu)
        T_over_mu_values.append(current_temperature/current_mu)
        current_temperature = np.exp(np.log(current_temperature)+dlogT)
    def convert_T_over_mu_to_T(T_over_mu):
            logT_over_mu = np.log(T_over_mu)
             # linear interpolation in log-log space
            logT = np.interp(logT_over_mu, np.log(T_over_mu_values), \
                             np.log(temperature_values))
            return np.exp(logT)
    return convert_T_over_mu_to_T