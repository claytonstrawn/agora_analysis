import os
import numpy as np
from unyt import unyt_array,unyt_quantity
from numpy import array

radii = np.array([0.1,0.15,0.2,0.3,0.5,0.7,1.0,1.5,2,4,8])

def read_analysis_data(code,redshift,which_data,\
                            printfields = False):
    if not os.path.exists('agora_analysis_quantities/%s_data_%s_%.3f.txt'%\
                          (which_data,code,redshift)):
        return (None,) if which_data in ['mass','flux'] else (None,None)
    with open('agora_analysis_quantities/%s_data_%s_%.3f.txt'%\
                          (which_data,code,redshift),'r') as f:
        lines = f.read().replace('\n','')
        d = eval(lines)
    radii = sorted(d.keys())
    for i,r in enumerate(radii):
        rd = d[r]
        if i==0:
            if which_data == 'resolution':
                all_res_data = np.zeros((len(radii),len(rd['resolution_info'])))
                all_res_data2 = np.zeros((len(radii),len(rd['resolution_info'])))
                to_return = (all_res_data,all_res_data2)
            elif which_data in ['mass','flux']:
                fields = list(rd.keys())
                fields = sorted(fields)
                all_nonres_data = np.zeros((len(radii),len(fields)))
                if printfields:
                    print(fields)
                to_return = (all_nonres_data,)
        if which_data == 'resolution':
            all_res_data[i] = rd['resolution_info']
            all_res_data2[i] = rd['resolution_info2']
        elif which_data in ['mass','flux']:
            for j,f in enumerate(fields):
                all_nonres_data[i,j] = rd[f].v
    return to_return