import numpy as np
import os

class NoMetadataError(Exception):
    def __init__(self,message):
        self.message = message
        
class NotCloseEnoughError(Exception):
    def __init__(self,message):
        self.message = message

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

def read_metadata(code,simnum,redshift,metadata_type,athreshold = .1):
    a = 1/redshift+1
    zthreshold = athreshold/a**2
    _ROOT = os.path.abspath(os.path.dirname(__file__))
    path_to_metadata = os.path.join(_ROOT,"agora_metadata")
    foldernames = {'C1':'Cal1','C2':'Cal2','C3':'Cal3','CR':'CosmoRun'}
    foldername = foldernames[simnum]
    metadata_location_folder = os.path.join(path_to_metadata,foldername)
    metadata_file = "AGORA_%s_%s_%s.txt"%(code,simnum,metadata_type)
    metadata_location = os.path.join(metadata_location_folder,metadata_file)
    with open(metadata_location) as f:
        lines = [line for line in f.readlines() if line.strip()]
    if metadata_type == 'base':
        zs,snapnums,center_xs,center_ys,center_zs,Rvirs = [],[],[],[],[],[]
        all_quantities = ('z,snapnum,center_x,center_y,center_z,Rvir').split(',')
    elif metadata_type == 'derived':
        zs,snapnums = [],[]
        Lxs,Lys,Lzs = [],[],[]
        bvxs,bvys,bvzs = [],[],[]
        all_quantities = ('z,snapnum,Lx,Ly,Lz,bvx,bvy,bvz').split(',')
    if len(lines) == 1:
        raise NoMetadataError('No metadata found for simulation %s!'%metadata_file)
    for i,line in enumerate(lines):
        if i == 0:
            known_quantities = line[:-1].split(', ')
            remaining_quantites = list(set(all_quantities)-set(known_quantities))
        else:
            sections = line[:-1].split(', ')
            for j,q in enumerate(known_quantities):
                floateval = '' if (q in ['snapnum']) else 'float'
                exec('%ss.append(%s(sections[%s]))'%(q,floateval,j))
            for q in remaining_quantites:
                exec('%ss.append(np.nan)'%(q,))
    closest_z_index = np.argmin(np.abs(np.array(zs)-redshift))
    if np.abs(zs[closest_z_index]-redshift)>zthreshold:
        message = "Nothing in %s is within %.3f of z=%.3f"%(metadata_location,zthreshold,redshift)
        raise NotCloseEnoughError(message)
    dict_to_return = {}
    for q in all_quantities:
        dict_to_return[q] = eval('%ss[%s]'%(q,closest_z_index))
    return dict_to_return