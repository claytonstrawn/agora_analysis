from unyt import unyt_array,unyt_quantity
import numpy as np
import yt
yt.set_log_level(50)
import os

def load_for_massloadingfactor(code,redshift):
    snap = AgoraSnapshot(code,redshift)
    snap.load_snapshot()
    snap.add_correct_fields()
    return snap
    
def snap_mass_analysis(snap,overwrite = False):
    code = snap.code
    redshift = snap.approx_redshift
    if os.path.exists('agora_analysis_quantities/mass_data_%s_%.3f.txt'%(code,redshift)):
        if not overwrite:
            print("file"+'agora_analysis_quantities/massloading_data_%s_%.3f.txt'%(code,redshift)+\
                  " already exists! skip!")
            return False
    #want to scan concentric spheres
    all_data = {}
    oldsp = None
    
    tot_gas_mass = unyt_quantity(0,'Msun')
    tot_star_mass = unyt_quantity(0,'Msun')
    tot_gas_metals = unyt_quantity(0,'Msun')
    tot_star_metals = unyt_quantity(0,'Msun')
    for i,r in enumerate(radii):
        print('\r'+'mass progress bar: '+'.'*(i+1)+' '*(len(radii)-i-1)+'!',end = '')
        if oldsp is not None:
            sh = snap.ds.sphere(snap.center,r*snap.Rvir)-oldsp
        else:
            sh = snap.ds.sphere(snap.center,r*snap.Rvir)
        data = {}
        #record total gas mass
        tot_gas_mass += np.sum(sh['gas','agora_mass']).in_units('Msun')
        data['tot_gas_mass'] = unyt_quantity(np.copy(tot_gas_mass.v),tot_gas_mass.units)
        #record total star mass
        tot_star_mass += np.sum(sh['agora_stars','agora_mass']).in_units('Msun')
        data['tot_star_mass'] = unyt_quantity(np.copy(tot_star_mass.v),tot_star_mass.units)
        #record total gas metals
        tot_gas_metals += np.sum(sh['gas','metal_mass']).in_units('Msun')
        data['tot_gas_metals'] = unyt_quantity(np.copy(tot_gas_metals.v),tot_gas_metals.units)
        #record total star metals
        tot_star_metals += np.sum(sh['agora_stars','metal_mass']).in_units('Msun')
        data['tot_star_metals'] = unyt_quantity(np.copy(tot_star_mass.v),tot_star_mass.units)
        all_data[r] = data
    with open('agora_analysis_quantities/mass_data_%s_%.3f.txt'%(code,redshift),'w') as f:
        f.write('%s\n'%(all_data))
        print()
    return True

