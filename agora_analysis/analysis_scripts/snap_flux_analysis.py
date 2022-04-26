from agora_analysis.field_setup.metal_functions import calculate_inflow_outflow

from unyt import unyt_array,unyt_quantity
import numpy as np
import yt
yt.set_log_level(50)
import os

radii = np.array([0.1,0.15,0.2,0.3,0.5,0.7,1.0,1.5,2,4,8])

    
def snap_flux_analysis(snap,overwrite = False):
    code = snap.code
    redshift = snap.approx_redshift
    if os.path.exists('agora_analysis_quantities/flux_data_%s_%.3f.txt'%(code,redshift)):
        if not overwrite:
            print("file"+'agora_analysis_quantities/flux_data_%s_%.3f.txt'%(code,redshift)+\
                  " already exists! skip!")
            return False
    #want to scan concentric spheres
    all_data = {}
    
    if code == 'gadget':
        dt = 1e7
    else:
        dt = 1e6
    
    for i,r in enumerate(radii):
        data = {}
        #print('\r'+'flux progress bar: '+'.'*(i+1)+' '*(len(radii)-i-1)+'!',end = '')
        #record gas inflow/outflow through surface
        i,o = calculate_inflow_outflow(snap,sphere_size = r, dt = dt,metals = False,stars = False)
        data['gas_mass_inflow'] = i
        data['gas_mass_outflow'] = o
        #record gas metal inflow/outflow through surface
        i,o = calculate_inflow_outflow(snap,sphere_size = r, dt = dt,metals = True,stars = False)
        data['gas_metals_inflow'] = i
        data['gas_metals_outflow'] = o
        #record stellar inflow/outflow through surface
        i,o = calculate_inflow_outflow(snap,sphere_size = r, dt = dt,metals = False,stars = True)
        data['star_mass_inflow'] = i
        data['star_mass_outflow'] = o
        #record gas metal inflow/outflow through surface
        i,o = calculate_inflow_outflow(snap,sphere_size = r, dt = dt,metals = True,stars = True)
        data['star_metals_inflow'] = i
        data['star_metals_outflow'] = o
        all_data[r] = data
    with open('agora_analysis_quantities/flux_data_%s_%.3f.txt'%(code,redshift),'w') as f:
        f.write('%s\n'%(all_data))
        print()
    return True

