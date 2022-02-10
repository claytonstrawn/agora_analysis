from agora_analysis.field_setup.common_fields import add_gas_fields,\
                                                        add_star_fields,\
                                                        add_metallicity_fields,\
                                                        add_resolution_fields,\
                                                        add_pressure_field
from agora_analysis.field_setup.metal_functions import add_flux_fields,\
                                                        add_star_metals
from agora_analysis.field_setup.ionization_fields import add_PI_CI_fields,is_ion_type,\
                                                        add_number_density_fields
from functools import partial

def load_all_fields(snap):
    add_gas_fields(snap)
    add_star_fields(snap)
    add_metallicity_fields(snap)
    add_resolution_fields(snap)
    add_temperature_fields(snap)
    add_pressure_field(snap)
    add_flux_fields(snap)
    add_star_metals(snap)
    
def detect_necessary_fields(field,required = None):
    if required is None:
        required = []
    if isinstance(field,str):
        ftype,fname = 'gas',field
    else:
        ftype, fname = field
    if ftype == 'gas':
        required.append(add_gas_fields)
        required.append(add_metallicity_fields)
        required.append(add_pressure_field)
    if fname in ['agora_cell_volume','agora_cell_volume_inv2',\
                 'agora_particle_volume','agora_particle_volume_inv2']:
        required.append(add_resolution_fields)
    if fname in ['radial_distance','metal_mass','radial_velocity',\
                 'gas_movement','metal_movement']:
        required.append(add_flux_fields)
    if ftype == 'agora_stars':
        required.append(add_star_fields)
        required.append(add_star_metals)
    if is_ion_type(fname):
        ion = is_ion_type(fname)
        partial_add_PI_CI_fields = partial(add_PI_CI_fields,ion=ion)
        partial_add_number_density_fields = partial(add_number_density_fields,ion=ion)
        required.append(partial_add_PI_CI_fields)
        required.append(partial_add_number_density_fields)
    return required

def load_necessary_fields(snap,fields):
    if isinstance(fields,list):
        fields = fields
    else:
        fields = [fields]
    required = None
    for field in fields:
        required = detect_necessary_fields(field,required)
    for func in required:
        func(snap)

def reload_fields_for_trident(snap,which = 'all'):
    if which == 'all':
        which = ['metallicity']
    elif isinstance(which,str):
        which = [which]
    for i,fname in enumerate(which):
        units = {'temperature':'K','metallicity':'','density':'g/cm**3'}
        def recreate_field(field,data):
            return data['gas','agora_'+fname]
        snap.ds.add_field(('gas', fname), 
                         function=recreate_field, 
                         force_override=True,
                         sampling_type = snap.sampling_type,
                         display_name=fname.capitalize(), 
                         take_log=True, 
                         units="")