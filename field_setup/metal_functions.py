import numpy as np

def add_flux_fields(snap):
    ds = snap.ds
    code = snap.code

    grid_codes = ['art','enzo','ramses']
    particle_codes = ['gadget','gear','gizmo','changa']

    sampling_type = 'grid' if code in grid_codes else 'particle'
    def radial_distance(field,data):
        xdist = data['gas','x']-snap.center_x
        ydist = data['gas','y']-snap.center_y
        zdist = data['gas','z']-snap.center_z
        return np.sqrt(xdist**2+ydist**2+zdist**2)
    ds.add_field(('gas','radial_distance'),function = radial_distance,
                 units = 'kpc',sampling_type = sampling_type,force_override = True)

    def metal_mass(field,data):
        return data['gas','agora_mass']*data['gas','agora_metallicity']
    ds.add_field(('gas','metal_mass'),function = metal_mass,
                    units = 'g',sampling_type = sampling_type,force_override = True)

    def radial_velocity(field,data):
        xdist = data['gas','x']-snap.center_x
        ydist = data['gas','y']-snap.center_y
        zdist = data['gas','z']-snap.center_z
        xvel = data['gas','velocity_x']
        yvel = data['gas','velocity_y']
        zvel = data['gas','velocity_z']
        dot = xdist*xvel+ydist*yvel+zdist*zvel
        return dot/data['gas','radial_distance']
    ds.add_field(('gas','radial_velocity'),function = radial_velocity,
                 units = 'km/s',sampling_type = sampling_type,force_override = True)
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
    
def calculate_inflow_outflow(snap,sphere_size=1,d=0.01,metals = True):
    ds =snap.ds
    dr = d*snap.Rvir*sphere_size

    sp = ds.sphere(snap.center, snap.Rvir*sphere_size)
    sp_surf = sp - ds.sphere(snap.center, snap.Rvir*sphere_size-dr)

    if metals:
        mass_over_sp = sp_surf['gas','metal_movement'].in_units('Msun*yr**-1*kpc')
    else:
        mass_over_sp = sp_surf['gas','gas_movement'].in_units('Msun*yr**-1*kpc')
    positive_flux = mass_over_sp>0
    negative_flux = mass_over_sp<0
    outflowing_mass = np.sum(mass_over_sp[positive_flux]/dr)
    inflowing_mass = np.sum(mass_over_sp[negative_flux]/dr)
    return inflowing_mass,outflowing_mass

def add_star_metals(snap):
    ds = snap.ds
    def metal_mass(field,data):
        return data['agora_stars','agora_mass']*data['agora_stars','agora_metallicity']
    ds.add_field(('agora_stars','metal_mass'),function = metal_mass,
                    units = 'g',sampling_type = 'particle',force_override = True)