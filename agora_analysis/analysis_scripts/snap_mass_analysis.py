from agora_analysis.field_setup.main import load_necessary_fields
from agora_analysis.field_setup.metal_functions import calculate_inflow_outflow,\
                                                    RegionTooSmallError

import yt
import numpy as np
import gc
yt.set_log_level(50)
import agora_analysis
import os

resolution_bins = np.array([0.0,10**-3.0,10**-2.5,10**-2.0,\
                            10**-1.5,10.**-1.0,10.**-0.5,\
                            10.**0.0,10.**0.5,10.**1.0,\
                            10.**1.5,10.**2.0,np.inf])
radii = np.array([0.1,0.15,0.2,0.3,0.5,0.7,1.0,1.5,2,4,8])

    
def snap_mass_analysis(name,redshift,overwrite = False):
    if os.path.exists('agora_mass_data_files/mass_data_%s_%.3f.txt'%(name,redshift)):
        if not overwrite:
            print("file"+'agora_mass_data_files/mass_data_%s_%.3f.txt'%(name,redshift)+\
                  " already exists! skip!")
            return False
    try:
        snap = agora_analysis.AgoraSnapshot(name,redshift = redshift)
        snap.load_snapshot()
        list_of_used_fields = [('gas','agora_mass'),('agora_stars','agora_mass'),\
                               ('gas','metal_mass'),('agora_stars','metal_mass'),\
                               ('gas','agora_cell_volume'),('gas','agora_particle_volume'),\
                               ('agora_stars','metal_movement'),('gas','metal_movement'),\
                               ('agora_stars','star_mass_movement'),('gas','gas_movement')]
        load_necessary_fields(snap,list_of_used_fields)
        code = snap.code
    except agora_analysis.utils.NotCloseEnoughError:
        print('skipping %s z=%1.3f!'%(name,redshift))
        return False
    except agora_analysis.utils.NoMetadataError:
        print('skipping %s z=%1.3f!'%(name,redshift))
        return False
    except Exception as e:
        raise e
    #want to scan concentric spheres
    all_data = {}
    for i,r in enumerate(radii):
        print('\r'+'.'*(i+1)+' '*(len(radii)-i)+'.',end = '')
        sp = snap.ds.sphere(snap.center,r*snap.Rvir)
        data = {}
        #record total gas mass
        tot_gas_mass = np.sum(sp['gas','agora_mass']).in_units('Msun')
        data['tot_gas_mass'] = tot_gas_mass
        #record total star mass
        tot_star_mass = np.sum(sp['agora_stars','agora_mass']).in_units('Msun')
        data['tot_star_mass'] = tot_star_mass
        #record total gas metals
        tot_gas_metals = np.sum(sp['gas','metal_mass']).in_units('Msun')
        data['tot_gas_metals'] = tot_gas_metals
        #record total star metals
        tot_star_metals = np.sum(sp['agora_stars','metal_mass']).in_units('Msun')
        data['tot_star_metals'] = tot_star_metals
        #record gas mass in resolution bins
        bin_masses = np.zeros(len(resolution_bins)-1)
        for i in range(len(bin_masses)):
            low = resolution_bins[i]
            high = resolution_bins[i+1]
            if code in agora_analysis.grid_codes:
                f = ('gas','agora_cell_volume')
            else:
                f = ('gas','agora_particle_volume')
            in_bin = np.logical_and((sp[f]**(1./3.)).in_units('kpc').v>low,
                                    (sp[f]**(1./3.)).in_units('kpc').v<high)
            bin_masses[i] = np.sum(sp['gas','agora_mass'][in_bin])
        data['resolution_info'] = bin_masses
        bin_masses2 = np.zeros(len(resolution_bins)-1)
        for i in range(len(bin_masses2)):
            low = resolution_bins[i]
            high = resolution_bins[i+1]
            if code in agora_analysis.grid_codes:
                continue
            else:
                f = ('gas','smoothing_length')
            in_bin = np.logical_and((sp[f]**(1./3.)).in_units('kpc').v>low,
                                    (sp[f]**(1./3.)).in_units('kpc').v<high)
            bin_masses[i] = np.sum(sp['gas','agora_mass'][in_bin])
        data['resolution_info2'] = bin_masses2
        #record gas inflow/outflow through surface
        try:
            i,o = calculate_inflow_outflow(snap,sphere_size = r, d = 0.01,metals = False,stars = False)
        except RegionTooSmallError:
            i,o = np.nan,np.nan
        data['gas_mass_inflow'] = i
        data['gas_mass_outflow'] = o
        #record gas metal inflow/outflow through surface
        try:
            i,o = calculate_inflow_outflow(snap,sphere_size = r, d = 0.01,metals = True,stars = False)
        except RegionTooSmallError:
            i,o = np.nan,np.nan
        data['gas_metals_inflow'] = i
        data['gas_metals_outflow'] = o
        #record stellar inflow/outflow through surface
        try:
            i,o = calculate_inflow_outflow(snap,sphere_size = r, d = 0.05,metals = False,stars = True)
        except RegionTooSmallError:
            i,o = np.nan,np.nan
        data['star_mass_inflow'] = i
        data['star_mass_outflow'] = o
        #record gas metal inflow/outflow through surface
        try:
            i,o = calculate_inflow_outflow(snap,sphere_size = r, d = 0.05,metals = True,stars = True)
        except RegionTooSmallError:
            i,o = np.nan,np.nan
        data['star_metals_inflow'] = i
        data['star_metals_outflow'] = o
        all_data[r] = data
    with open('agora_mass_data_files/mass_data_%s_%.3f.txt'%(name,redshift),'w') as f:
        f.write('%s\n'%(all_data))
    return True

from numpy import array
from unyt import unyt_quantity
def read_mass_analysis_data(name,redshift,printfields = True,include_zeros = False):
    if not os.path.exists('agora_mass_data_files/mass_data_%s_%.3f.txt'%(name,redshift)):
        return None,None,None
    with open('agora_mass_data_files/mass_data_%s_%.3f.txt'%(name,redshift),'r') as f:
        lines = f.read().replace('\n','')
        d = eval(lines)
    import numpy as np
    if include_zeros:
        radii = [0.0]+sorted(d.keys())
    else:
        radii = sorted(d.keys())

    for i,r in enumerate(radii):
        if i==0:
            rd = d[radii[1]]
            fields = list(rd.keys())
            fields.remove('resolution_info')
            fields.remove('resolution_info2')
            fields = sorted(fields)
            if printfields:
                print(fields)
            all_nonres_data = np.zeros((len(radii),len(fields)))
            all_res_data = np.zeros((len(radii),len(rd['resolution_info'])))
            all_res_data2 = np.zeros((len(radii),len(rd['resolution_info2'])))
            if include_zeros:
                continue
        rd = d[r]
        for j,f in enumerate(fields):
            all_nonres_data[i,j] = rd[f].v
        all_res_data[i] = rd['resolution_info']
        all_res_data2[i] = rd['resolution_info2']
    print('all_nonres_data,all_res_data,all_res_data2')
    return all_nonres_data,all_res_data,all_res_data2