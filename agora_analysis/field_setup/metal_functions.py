import numpy as np
    
def add_metal_mass_fields(snap):
    ds = snap.ds
    sampling_type = snap.sampling_type
    def metal_mass(field,data):
        return data['gas','agora_mass']*data['gas','agora_metallicity']
    ds.add_field(('gas','metal_mass'),function = metal_mass,
                    units = 'g',sampling_type = sampling_type,force_override = True)
    def star_metal_mass(field,data):
        return data['agora_stars','agora_mass']*data['agora_stars','agora_metallicity']
    ds.add_field(('agora_stars','metal_mass'),function = star_metal_mass,
                    units = 'g',sampling_type = 'particle',force_override = True)

def add_flux_fields(snap):
    ds = snap.ds
    sampling_type = snap.sampling_type
    def gas_movement(field,data):
        radial_velocity = data['gas','radial_velocity']
        gas_mass = data['gas','agora_mass']
        return radial_velocity*gas_mass
    ds.add_field(('gas','gas_movement'),function = gas_movement,
                 units = 'g*m*s**-1',sampling_type = sampling_type,force_override = True)
    def metal_movement(field,data):
        radial_velocity = data['gas','radial_velocity']
        metal_mass = data['gas','metal_mass']
        return radial_velocity*metal_mass
    ds.add_field(('gas','metal_movement'),function = metal_movement,
                 units = 'g*m*s**-1',sampling_type = sampling_type,force_override = True)
    def star_mass_movement(field,data):
        radial_velocity = data['agora_stars','radial_velocity']
        star_mass = data['agora_stars','agora_mass']
        return radial_velocity*star_mass
    ds.add_field(('agora_stars','star_mass_movement'),function = star_mass_movement,
                 units = 'g*m*s**-1',sampling_type = 'particle',force_override = True)
    def star_metal_movement(field,data):
        radial_velocity = data['agora_stars','radial_velocity']
        metal_mass = data['agora_stars','metal_mass']
        return radial_velocity*metal_mass
    ds.add_field(('agora_stars','metal_movement'),function = star_metal_movement,
                 units = 'g*m*s**-1',sampling_type = 'particle',force_override = True) 

def calculate_inflow_outflow(snap,sphere_size=1,d=0.01,metals = True,stars = False):
    ds =snap.ds
    dr = d*snap.Rvir*sphere_size

    sp = ds.sphere(snap.center, snap.Rvir*sphere_size)
    sp_surf = sp - ds.sphere(snap.center, snap.Rvir*sphere_size-dr)
    if 1:
    #try:
        if metals and stars:
            mass_over_sp = sp_surf['agora_stars','metal_movement'].in_units('Msun*yr**-1*kpc')
        elif metals and (not stars):
            mass_over_sp = sp_surf['gas','metal_movement'].in_units('Msun*yr**-1*kpc')
        elif (not metals) and stars:
            mass_over_sp = sp_surf['agora_stars','star_mass_movement'].in_units('Msun*yr**-1*kpc')
        elif (not metals) and (not stars):
            mass_over_sp = sp_surf['gas','gas_movement'].in_units('Msun*yr**-1*kpc')
    else:
    #except ValueError:
        print('d=%.2f gives surface with 0 cells!'%d)
        return np.nan,np.nan
    positive_flux = mass_over_sp>0
    negative_flux = mass_over_sp<0
    outflowing_mass = np.sum(mass_over_sp[positive_flux]/dr)
    inflowing_mass = np.sum(mass_over_sp[negative_flux]/dr)
    return inflowing_mass,outflowing_mass