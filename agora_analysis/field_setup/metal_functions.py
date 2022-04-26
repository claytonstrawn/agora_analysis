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

def calculate_inflow_outflow(snap,sphere_size=1,outer_sphere = 10,dt=1e6,metals = True,stars = False):
    print('stars=%s,metals=%s'%(stars,metals))
    ds = snap.ds
    if isinstance(dt,float):
        dt = ds.quan(dt,'yr')
    if stars:
        part_type = 'agora_stars' 
        sampling_type_str = 'particles'
    else:
        part_type = 'gas'
        if snap.sampling_type == 'particle':
            sampling_type_str = 'particles'
        else:
            sampling_type_str = 'cells'
    position_type = 'radial_distance'
    velocity_type = 'radial_velocity'
    
    r = ds.quan(snap.Rvir*sphere_size)
    R = ds.quan(snap.Rvir*outer_sphere)
    small_sp = ds.sphere(snap.center, r)
    large_sp = ds.sphere(snap.center, R)
    radial_velocities = large_sp[part_type,velocity_type]
    minimum = np.amin(radial_velocities)
    maximum = np.amax(radial_velocities)
    print(minimum,maximum)
    if r-maximum*dt>0.0*R:
        cannot_outflow = ds.sphere(snap.center,r-maximum*dt)
        inner = small_sp - cannot_outflow
    else:
        inner = small_sp
    print("outflows can be between %f and %f Rvir"%((r-maximum*dt)/snap.Rvir,r/snap.Rvir))
    if r-minimum*dt<1.0*R:
        cannot_inflow = large_sp - ds.sphere(snap.center,min(r-minimum*dt,R))
        outer = large_sp - small_sp - cannot_inflow
    else:
        outer = large_sp - small_sp
    print("inflows can be between %f and %f Rvir"%(r/snap.Rvir,(r-minimum*dt)/snap.Rvir))
    will_outflow = (inner[part_type,velocity_type]*dt + inner[part_type,position_type]) >= r
    will_inflow = (outer[part_type,velocity_type]*dt + outer[part_type,position_type]) <= r
    if metals:
        mass_type = 'metal_mass'
    else:
        mass_type = 'agora_mass'
    outflowing_mass = np.sum(inner[part_type,mass_type][will_outflow]).in_units('Msun')
    inflowing_mass = np.sum(outer[part_type,mass_type][will_inflow]).in_units('Msun')
    print('total number of %s outflowing = %d out of %d (%.2f%%), with total mass %.3e'%\
          (sampling_type_str,np.sum(will_outflow),len(will_outflow),\
           np.sum(will_outflow)/len(will_outflow)*100, outflowing_mass))
    print('total number of %s inflowing = %d out of %d (%.2f%%), with total mass %.3e'%\
          (sampling_type_str,np.sum(will_inflow),len(will_inflow),\
           np.sum(will_inflow)/len(will_inflow)*100, inflowing_mass))
    outflow = outflowing_mass/dt
    inflow = inflowing_mass/dt
    print(inflow.in_units('Msun/yr'),outflow.in_units('Msun/yr'))
    return inflow.in_units('Msun/yr'),outflow.in_units('Msun/yr')