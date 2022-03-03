from agora_analysis.field_setup.common_fields import add_star_fields,\
                                                        add_metallicity_fields,\
                                                        add_gas_fields,\
                                                        add_resolution_fields,\
                                                        add_temperature_fields,\
                                                        add_pressure_fields,\
                                                        add_new_star_ages,\
                                                        add_radial_distance_fields
from agora_analysis.field_setup.relative_velocity_fields import add_cylindrical_velocity_fields,\
                                                        add_radial_velocity_fields
from agora_analysis.field_setup.metal_functions import add_metal_mass_fields,\
                                                        add_flux_fields
from agora_analysis.field_setup.ionization_fields import add_PI_CI_fields,\
                                                            is_ion_type,\
                                                            add_number_density_fields
from functools import partial

def load_all_fields(snap):
    add_star_fields(snap)
    add_gas_fields(snap)
    add_metallicity_fields(snap)
    add_resolution_fields(snap)
    add_temperature_fields(snap)
    add_pressure_fields(snap)
    add_new_star_ages(snap)
    add_radial_distance_fields(snap)
    add_cylindrical_velocity_fields(snap)
    add_radial_velocity_fields(snap)
    add_metal_mass_fields(snap)
    add_flux_fields(snap)
    
def add_to_required(toadd, required):
    if toadd not in required:
        required.append(toadd)
    return required        
    
def detect_necessary_fields(field,required = None):
    if required is None:
        required = []
    if isinstance(field,str):
        ftype,fname = 'gas',field
    else:
        ftype, fname = field
    if ftype == 'gas':
        add_to_required(add_gas_fields,required)
        add_to_required(add_metallicity_fields,required)
        add_to_required(add_metallicity_fields,required)
        add_to_required(add_temperature_fields,required)
    if ftype == 'agora_stars':
        add_to_required(add_star_fields,required)
    if fname in ['agora_cell_volume','agora_cell_volume_inv2',\
                 'agora_particle_volume','agora_particle_volume_inv2']:
        add_to_required(add_resolution_fields,required)
    if fname in ['radial_distance','radial_velocity',\
                 'cylindrical_velocity_z','gas_movement',\
                 'metal_movement','gas_movement',\
                 'star_mass_movement']:
        add_to_required(add_radial_distance_fields,required)
    if fname in ['radial_velocity','cylindrical_velocity_z','gas_movement',
                 'metal_movement','gas_movement','star_mass_movement']:
        add_to_required(add_radial_velocity_fields,required)
    if fname in ['metal_mass','gas_movement','metal_movement','gas_movement','star_mass_movement']:
        add_to_required(add_metal_mass_fields,required)
    if fname in ['gas_movement','metal_movement','gas_movement','star_mass_movement']:
        add_to_required(add_flux_fields,required)
    if is_ion_type(fname):
        ion = is_ion_type(fname)
        partial_add_PI_CI_fields = partial(add_PI_CI_fields,ion=ion)
        partial_add_number_density_fields = partial(add_number_density_fields,ion=ion)
        add_to_required(partial_add_PI_CI_fields,required)
        add_to_required(partial_add_number_density_fields,required)
    if 'cylindrical_velocity' in fname:
        add_to_required(add_cylindrical_velocity_fields,required)
    return required

def load_necessary_fields(snap,fields,override = False,print_requirements = True):
    if isinstance(fields,list):
        fields = fields
    else:
        fields = [fields]
    required = None
    for field in fields:
        required = detect_necessary_fields(field,required)
    if print_requirements:
        required_to_print = [f.__name__ for f in required]
        print('required field additions are %s to instantiate %s'%(required_to_print,fields))
    for func in required:
        func(snap)
    if override is True:
        override_default_fields(snap)
    elif isinstance(override,str):
        override_default_fields(snap,which = override)

def override_default_fields(snap,which = 'all'):
    if which == 'all':
        which = ['temperature','metallicity','density']
    elif isinstance(which,str):
        which = [which]
    print('replacing default %s fields with agora defined versions...'%which)
    for i,fname in enumerate(which):
        units = {'temperature':'K','metallicity':'','density':'g/cm**3'}
        def recreate_field(field,data):
            return data['gas','agora_%s'%fname]
        snap.ds.add_field(('gas', fname), 
                         function=recreate_field, 
                         force_override=True,
                         sampling_type = snap.sampling_type,
                         display_name=fname.capitalize(), 
                         take_log=True, 
                         units="")