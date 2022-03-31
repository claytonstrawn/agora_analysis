import agora_analysis
from unyt import unyt_array,unyt_quantity
import numpy as np
import yt
yt.set_log_level(50)
import os

resolution_bins = np.array([0.0,10**-3.0,10**-2.5,10**-2.0,\
                            10**-1.5,10.**-1.0,10.**-0.5,\
                            10.**0.0,10.**0.5,10.**1.0,\
                            10.**1.5,10.**2.0,np.inf])
radii = np.array([0.1,0.15,0.2,0.3,0.5,0.7,1.0,1.5,2,4,8])
    
def snap_resolution_analysis(snap,overwrite = False):
    code = snap.code
    redshift = snap.approx_redshift
    if os.path.exists('agora_analysis_quantities/resolution_data_%s_%.3f.txt'%(code,redshift)):
        if not overwrite:
            print("file"+'agora_analysis_quantities/resolution_data_%s_%.3f.txt'%(code,redshift)+\
                  " already exists! skip!")
            return False
    #want to scan concentric spheres
    all_data = {}
    oldsp = None
    bin_masses = unyt_array(np.zeros(len(resolution_bins)-1),'Msun')
    bin_masses2 = unyt_array(np.zeros(len(resolution_bins)-1),'Msun')

    for i,r in enumerate(radii):
        data = {}
        print('\r'+'res progress bar:  '+'.'*(i+1)+' '*(len(radii)-i-1)+'!',end = '')
        if oldsp is not None:
            sh = snap.ds.sphere(snap.center,r*snap.Rvir)-oldsp
        else:
            sh = snap.ds.sphere(snap.center,r*snap.Rvir)
        data = {}
        #record gas mass in resolution bins
        for j in range(len(resolution_bins)-1):
            low = resolution_bins[j]
            high = resolution_bins[j+1]
            if code in agora_analysis.grid_codes:
                f = ('gas','agora_cell_volume')
            else:
                f = ('gas','agora_particle_volume')
            in_bin = np.logical_and((sh[f]**(1./3.)).in_units('kpc').v>low,
                                    (sh[f]**(1./3.)).in_units('kpc').v<high)
            bin_masses[j] += np.sum(sh['gas','agora_mass'][in_bin])
            if code in agora_analysis.grid_codes:
                continue
            else:
                f = ('gas','smoothing_length')
            in_bin = np.logical_and((sh[f]).in_units('kpc').v>low,
                                    (sh[f]).in_units('kpc').v<high)
            bin_masses2[j] += np.sum(sh['gas','agora_mass'][in_bin])
        data['resolution_info'] = unyt_array(np.copy(bin_masses.v),bin_masses.units)
        data['resolution_info2'] = unyt_array(np.copy(bin_masses2.v),bin_masses2.units)
        oldsp = snap.ds.sphere(snap.center,r*snap.Rvir)
        all_data[r] = data
    with open('agora_analysis_quantities/resolution_data_%s_%.3f.txt'%(code,redshift),'w') as f:
        f.write('%s\n'%(all_data))
        print()
    return True

