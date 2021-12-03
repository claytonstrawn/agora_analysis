from agora_analysis.field_setup.common_fields import *
from agora_analysis.field_setup.metal_functions import *

def load_all_fields(snap):
    add_gas_fields(snap)
    add_star_fields(snap)
    add_metallicity_fields(snap)
    add_resolution_fields(snap)
    add_temperature_fields(snap)
    add_flux_fields(snap)
    add_star_metals(snap)