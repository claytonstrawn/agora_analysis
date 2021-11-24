from agora_analysis.field_setup.common_fields import add_metallicity_fields,\
                                        add_gas_fields, \
                                        add_star_fields,\
                                        add_resolution_fields,\
                                        grid_codes, \
                                        particle_codes
from agora_analysis.field_setup.metal_functions import add_flux_fields,\
                                                    calculate_inflow_outflow,\
                                                    add_star_metals,\
                                                    RegionTooSmallError

import yt
import numpy as np
import gc
yt.set_log_level(50)
import agora_analysis
import os

resolution_bins = np.array([0.0,10**-3.0,10**-2.5,10**-2.0,10**-1.5,10.**-1.0,10.**-0.5,10.**0.0,10.**0.5,10.**1.0,\
                                       10.**1.5,10.**2.0,np.inf])
radii = np.array([0.1,0.2,0.3,0.5,0.7,1.0,1.5,2,4,8])

    
def snap_mass_analysis(name,redshift,overwrite = False):
    if os.path.exists('agora_mass_data_files/mass_data_%s_%.3f.txt'%(name,redshift)):
        if not overwrite:
            print("file"+'agora_mass_data_files/mass_data_%s_%.3f.txt'%(name,redshift)+\
                  " already exists! skip!")
            return False
    try:
        snap = agora_analysis.AgoraSnapshot(name,redshift = redshift)
        snap.load_snapshot()
        code = snap.code
        add_gas_fields(snap)
        add_star_fields(snap)
        add_metallicity_fields(snap)
        add_resolution_fields(snap)
        add_flux_fields(snap)
        add_star_metals(snap)
    except agora_analysis.utils.NotCloseEnoughError:
        print('skipping %s z=%1.3f!'%(name,redshift))
        return False
    except agora_analysis.utils.NoMetadataError:
        print('skipping %s z=%1.3f!'%(name,redshift))
        return False
    #want to scan concentric spheres
    all_data = {}
    for r in radii:
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
            if code in grid_codes:
                f = ('index','agora_cell_volume')
            else:
                f = ('gas','agora_particle_volume')
            in_bin = np.logical_and((sp[f]**(1./3.)).in_units('kpc').v>low,
                                    (sp[f]**(1./3.)).in_units('kpc').v<high)
            bin_masses[i] = np.sum(sp['gas','agora_mass'][in_bin])
        data['resolution_info'] = bin_masses
        #record gas inflow/outflow through surface
        try:
            i,o = calculate_inflow_outflow(snap,sphere_size = r, d = 0.01,metals = False)
        except RegionTooSmallError:
            i,o = np.nan,np.nan
        data['gas_mass_inflow'] = i
        data['gas_mass_outflow'] = o
        #record gas metal inflow/outflow through surface
        try:
            i,o = calculate_inflow_outflow(snap,sphere_size = r, d = 0.01,metals = True)
        except RegionTooSmallError:
            i,o = np.nan,np.nan
        data['gas_metals_inflow'] = i
        data['gas_metals_outflow'] = o
            
        all_data[r] = data
    with open('agora_mass_data_files/mass_data_%s_%.3f.txt'%(name,redshift),'w') as f:
        f.write('%s\n'%(all_data))
    return True

from numpy import array
from unyt import unyt_quantity
def read_mass_analysis_data(name,redshift,printfields = False,include_zeros = False):
    if not os.path.exists('agora_mass_data_files/mass_data_%s_%.3f.txt'%(name,redshift)):
        return None,None
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
            fields = sorted(fields)
            if printfields:
                print(fields)
            all_nonres_data = np.zeros((len(radii),len(fields)))
            all_res_data = np.zeros((len(radii),len(rd['resolution_info'])))
            if include_zeros:
                continue
        rd = d[r]
        for j,f in enumerate(fields):
            all_nonres_data[i,j] = rd[f].v
        all_res_data[i] = rd['resolution_info']
    return all_nonres_data,all_res_data