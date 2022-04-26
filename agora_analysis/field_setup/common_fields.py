from yt.data_objects.particle_filters import add_particle_filter
from yt.utilities.physical_constants import mp, kb
from yt.utilities.cosmology import Cosmology
from yt import YTArray
from agora_analysis.utils import convert_T_over_mu_to_T_func
import numpy as np

PartType_Gas_to_use_dict = {'art':'gas','ramses':'gas','enzo':'gas','gadget':'PartType0',
                       'gear':'PartType0','changa':'Gas','gizmo':'PartType0'}
PartType_Star_to_use_dict = {'art':'stars','ramses':'star','enzo':'Stars','gadget':'PartType4',
                       'gear':'PartType1','changa':'Stars','gizmo':'PartType4'}
PartType_StarBeforeFiltered_to_use_dict = {'art':'stars','ramses':'all','enzo':'all',
                                          'gadget':'PartType4','gear':'PartType1',
                                           'changa':'Stars',
                                           'gizmo':'PartType4'}
GasMassType_to_use_dict = {'art':'cell_mass','ramses':'cell_mass',\
                           'enzo':'cell_mass','gadget':'Masses',\
                       'gear':'Masses','changa':'Mass','gizmo':'Masses'}
GasDensityType_to_use_dict = {'art':'density','ramses':'density',\
                           'enzo':'density','gadget':'Density',\
                       'gear':'Density','changa':'Density','gizmo':'Density'}
StarMassType_to_use_dict = {'art':'particle_mass','ramses':'particle_mass',
                            'enzo':'particle_mass','gadget':'Masses',
                       'gear':'Masses','changa':'Mass','gizmo':'Masses'}
FormationTimeType_to_use_dict = {'art':'particle_creation_time',\
                                 'ramses':'particle_birth_time','enzo':'creation_time',\
                                 'gadget':'StellarFormationTime','gear':'StarFormationTime',\
                                 'changa':'FormationTime','gizmo':'StellarFormationTime'}
StarMetallicityType_to_use_dict = {'art':'particle_metallicity','ramses':'particle_metallicity',\
                              'enzo':'metallicity_fraction','gadget':'metallicity',\
                              'gear':'StarMetals','changa':'metallicity','gizmo':'metallicity'}
            
grid_codes = ['art','enzo','ramses']
particle_codes = ['gadget','gear','gizmo','changa']

def add_gas_fields(snap):
    pf = snap.ds
    code = snap.code
    PartType_Gas_to_use = PartType_Gas_to_use_dict[code]
    GasMassType_to_use = GasMassType_to_use_dict[code]
    GasDensityType_to_use = GasDensityType_to_use_dict[code]
    sampling_type = 'cell' if code in grid_codes else 'particle'
    def gas_mass(field, data):
        return data[PartType_Gas_to_use,GasMassType_to_use]
    pf.add_field(('gas','agora_mass'),
                    sampling_type = sampling_type,
                    function = gas_mass,
                    force_override = True,
                    units = 'Msun')
    def gas_density(field, data):
        return data[PartType_Gas_to_use,GasDensityType_to_use]
    pf.add_field(('gas','agora_density'),
                    sampling_type = sampling_type,
                    function = gas_density,
                    force_override = True,
                    units = 'g/cm**3')
    def density_squared(field, data):
        return data[(PartType_Gas_to_use, GasDensityType_to_use)]**2
    pf.add_field(("gas", "agora_density_squared"),
                 function=density_squared, 
                 sampling_type = sampling_type,
                 units="g**2/cm**6")


def add_star_fields(snap):
    pf = snap.ds
    code = snap.code
    PartType_Star_to_use = PartType_Star_to_use_dict[code]
    PartType_StarBeforeFiltered_to_use = PartType_StarBeforeFiltered_to_use_dict[code]
    FormationTimeType_to_use = FormationTimeType_to_use_dict[code]
    StarMetallicityType_to_use = StarMetallicityType_to_use_dict[code]
    StarMassType_to_use = StarMassType_to_use_dict[code]

    if code == 'art':
        def agora_stars_art(pfilter, data):
            return (data[(pfilter.filtered_type, StarMassType_to_use)] > 0)
        add_particle_filter('agora_stars', function=agora_stars_art, \
                        filtered_type=PartType_StarBeforeFiltered_to_use, \
                        requires=[FormationTimeType_to_use])
    else:
        def agora_stars(pfilter, data):
            return (data[(pfilter.filtered_type, FormationTimeType_to_use)] > 0)
        add_particle_filter('agora_stars', function=agora_stars, \
                filtered_type=PartType_StarBeforeFiltered_to_use, \
                requires=[FormationTimeType_to_use])
    pf.add_particle_filter('agora_stars')
    if code == 'art':
        pf.add_particle_filter('agora_stars')
        def star_metallicity_art(field,data):
            return data['agora_stars','particle_metallicity1']+\
                    data['agora_stars','particle_metallicity2']
        pf.add_field(('agora_stars','agora_metallicity'),
                    sampling_type = 'particle',
                    function = star_metallicity_art,
                    force_override = True,
                    units = '')
    elif code == 'gear':
        def gear_star_metallicity(field, data):
            if len(data['agora_stars', StarMetallicityType_to_use].shape) == 1:
                return data['agora_stars',StarMetallicityType_to_use].in_units("")
            else:
                return data['agora_stars', StarMetallicityType_to_use][:,9].in_units("") 
            # in_units("") turned out to be crucial!; otherwise code_metallicity 
            #will be used and it will mess things up
            # We are creating ("Gas", "Metallicity") here, 
            #different from ("Gas", "metallicity") which is 
            #auto-generated by yt but doesn't work properly
        pf.add_field(('agora_stars', 'agora_metallicity'), 
                     function=gear_star_metallicity,
                     force_override=True,
                     sampling_type = 'particle',
                     display_name="Metallicity",
                     take_log=True, 
                     units="")
    else:
        def star_metallicity(field, data):
            return data['agora_stars',StarMetallicityType_to_use]
        pf.add_field(('agora_stars','agora_metallicity'),
                        sampling_type = 'particle',
                        function = star_metallicity,
                        force_override = True,
                        units = '')

    def star_mass(field, data):
        return data['agora_stars',StarMassType_to_use]
    pf.add_field(('agora_stars','agora_mass'),
                sampling_type = 'particle',
                function = star_mass,
                force_override = True,
                units = 'Msun')

def add_metallicity_fields(snap):
    pf = snap.ds
    code = snap.code
    sampling_type = snap.sampling_type
    PartType_Gas_to_use = PartType_Gas_to_use_dict[code]
    if code == 'art': 
        def art_gas_metallicity(field, data):
            return (data["gas", "metal_ii_density"] + data["gas", "metal_ia_density"]) / \
                    data["gas", "agora_density"]
        pf.add_field(("gas", "agora_metallicity"), 
                     function=art_gas_metallicity, 
                     force_override=True, 
                     display_name="Metallicity", 
                     sampling_type = 'cell',
                     take_log=True, 
                     units="")
    elif code == 'gadget':
        def gadget_metallicity(field, data):
            return data[(PartType_Gas_to_use, "Metallicity")]+YTArray(1e-4*0.02041,'')
        pf.add_field(('gas', 'agora_metallicity'), 
                     function=gadget_metallicity, 
                     force_override=True,
                     take_log=True, 
                     sampling_type = 'particle',
                     display_name="Metallicity", 
                     units="")
    elif code == 'gear':
        def gear_metallicity(field, data):
            if len(data[PartType_Gas_to_use, "metallicity"].shape) == 1:
                return data[PartType_Gas_to_use,"metallicity"].in_units("")
            else:
                return data[PartType_Gas_to_use, "metallicity"][:,9].in_units("") 
            # in_units("") turned out to be crucial!; otherwise code_metallicity 
            #will be used and it will mess things up
            # We are creating ("Gas", "Metallicity") here, 
            #different from ("Gas", "metallicity") which is 
            #auto-generated by yt but doesn't work properly
        pf.add_field(('gas', 'agora_metallicity'), 
                     function=gear_metallicity,
                     force_override=True,
                     sampling_type = 'particle',
                     display_name="Metallicity",
                     take_log=True, 
                     units="")
    elif code == 'enzo':
        def enzo_metallicity(field, data):
            return data["gas", "metal_density"] / data["gas", "density"]
        pf.add_field(("gas", "agora_metallicity"), 
                     function=enzo_metallicity, 
                     force_override=True, 
                     display_name="Metallicity", 
                     sampling_type = 'cell',
                     take_log=True, 
                     units="")
    else:
        def gas_metallicity(field, data):
            return data['gas',"metallicity"]
        pf.add_field(('gas','agora_metallicity'),
                    sampling_type = sampling_type,
                    function = gas_metallicity,
                    force_override = True,
                    units = '')


def add_resolution_fields(snap):
    pf = snap.ds
    code = snap.code
    if code in grid_codes:
        def cell_volume(field,data):
            return data[("index", "cell_volume")]
        pf.add_field(("gas", "agora_cell_volume"),function=cell_volume,sampling_type = 'cell',
                     units='kpc**3', display_name="Resolution $\Delta$ x", take_log=True )
        def inverse_cell_volume_squared(field,data):
            return data[("index", "cell_volume")]**-2
        pf.add_field(("gas", "agora_cell_volume_inv2"), function=inverse_cell_volume_squared, \
                     units='pc**(-6)',sampling_type = 'cell',
                     display_name="Inv2CellVolumeCode", take_log=True)
    elif code in particle_codes:
        def particle_volume(field, data):
            return data[('gas', 'agora_mass')]/data[('gas', "agora_density")]
        pf.add_field(('gas', "agora_particle_volume"), 
                     function=particle_volume, 
                     units="kpc**3", 
                     display_name="Resolution $\Delta$ x", 
                     sampling_type = 'particle',
                     take_log=True)
        def inverse_particle_size_squared(field, data):
            return (data[('gas', 'agora_mass')]/data[('gas', "agora_density")])**(-2.)
        pf.add_field(('gas', "agora_particle_volume_inv2"), 
                     function=inverse_particle_size_squared, 
                     units="pc**(-6)", 
                     display_name="Inv2ParticleVolume", 
                     sampling_type = 'particle',
                     take_log=True)
            
def add_temperature_fields(snap):
    pf = snap.ds
    code = snap.code
    convert_T_over_mu_to_T = convert_T_over_mu_to_T_func()
    PartType_Gas_to_use = PartType_Gas_to_use_dict[code]
    if code == 'ramses':
        def ramses_temp(field, data):
            T_J = YTArray(0.0,'K')  # in K
            n_H = YTArray(8.0,'1/cm**3') \
            # Should be density threshold for star formation
            del_star=200.
            rhoc=YTArray(1.88e-29,'g/cm**3') # Critical density
            Omega_b=0.042
            hubble=pf.hubble_constant
            aexp=1/(pf.current_redshift+1)
            gamma_0 = 2.0
            x_H = 0.76
            mH = YTArray(1.660539e-24,'g')
            # from pymses/utils/constants/__init__.py  (vs. in yt, mass_hydrogen_cgs = 
            # 1.007947*amu_cgs = 1.007947*1.660538921e-24 = 1.6737352e-24)
            kB = YTArray(1.3806494e-16,'cm**2*g/s**2/K')
            # from pymses/utils/constants/__init__.py  (vs. in yt, 
            #boltzmann_constant_cgs = 1.3806488e-16)
            n_J = max(YTArray(rhoc*del_star*Omega_b*(hubble)**2/aexp**3*x_H/mH,'1/cm**3'),n_H) 
            #max value n_ISM vs. n_SF
            T_over_mu = data[PartType_Gas_to_use, "pressure"].in_units('g/s**2/cm')/\
                        (data[PartType_Gas_to_use, "density"].in_units('g/cm**3')) *\
                        mH / kB - \
                        T_J * (data[PartType_Gas_to_use, "density"].in_units('g/cm**3') *
                               x_H / mH / n_J)**(gamma_0 - 1.0) # T/mu = T2 in Ramses
            return YTArray(convert_T_over_mu_to_T(T_over_mu), 'K') # now T
        pf.add_field(("gas", "agora_temperature"), 
                     function = ramses_temp, 
                     sampling_type = 'cell',
                     force_override=True, 
                     display_name="Temperature",
                     take_log=True,
                     units="K")
    elif code == 'gear':
        def gear_temp(field, data):
            # Assume cosmic abundances
            gamma = 5.0/3.0
            T_over_mu_GEAR= data[PartType_Gas_to_use, "InternalEnergy"].in_units('J/g')*\
                            (gamma-1)*mp.in_units('g')/kb.in_units('J/K')
            return YTArray(convert_T_over_mu_to_T(T_over_mu_GEAR), 'K') # now T
        pf.add_field(('gas', 'agora_temperature'),
                     function=gear_temp,
                     force_override=True,
                     sampling_type = 'particle',
                     display_name="Temperature",
                     take_log=True, 
                     units="K")
    elif code == 'gadget':
        HYDROGEN_MASSFRAC = 0.76
        gamma=5.0/3.0
        GAMMA_MINUS1=gamma-1.
        PROTONMASS=mp.in_units('g')
        BOLTZMANN=kb.in_units('J/K')
        u_to_temp_fac=(4.0 / (8.0 - 5.0 * (1.0 - HYDROGEN_MASSFRAC))) *\
                        PROTONMASS / BOLTZMANN * GAMMA_MINUS1
        def gadget_temp(field, data):
            # Assume cosmic abundances
            gamma = 5.0/3.0
            T_over_mu_GADGET= data[PartType_Gas_to_use, "InternalEnergy"].in_units('J/g')*\
                            (gamma-1)*mp.in_units('g')/kb.in_units('J/K')
            return YTArray(convert_T_over_mu_to_T(T_over_mu_GADGET), 'K') # now T
        pf.add_field(('gas', 'agora_temperature'), 
                     function=gadget_temp,
                     force_override=True,
                     sampling_type = 'particle',
                     display_name="Temperature",
                     take_log=True, 
                     units="K")
    elif code == 'gizmo':
        def gizmo_temp(field, data):
            # Assume cosmic abundances
            x_H = 0.76
            gamma = 5.0/3.0
            if data.has_field_parameter(('gas',"mean_molecular_weight")):
                mu = data.get_field_parameter(('gas',"mean_molecular_weight"))
            else:
                # Assume zero ionization
                mu = 4.0 / (3.0 * x_H + 1.0)
            ret = data[PartType_Gas_to_use, "InternalEnergy"]*(gamma-1)*mu*mp/kb
            return ret.in_units('K')
        pf.add_field(('gas', 'agora_temperature'), 
                     function=gizmo_temp,
                     force_override=True,
                     sampling_type = 'particle',
                     display_name="Temperature",
                     take_log=True,
                     units="K")
    else:
        sampling_type = 'cell' if code in grid_codes else 'particle'
        def gas_temperature(field,data):
            return data[PartType_Gas_to_use,'temperature']
        pf.add_field(('gas','agora_temperature'),
                        sampling_type = sampling_type,
                        function = gas_temperature,
                        force_override = True,
                        display_name="Temperature",
                        take_log=True,
                        units = 'K')

def add_pressure_fields(snap):
    pf = snap.ds
    code = snap.code
    def gas_pressure(field,data):
        return data['gas','temperature']*data['gas','number_density']
    pf.add_field(('gas','agora_pressure'),
                    sampling_type = snap.sampling_type,
                    function = gas_pressure,
                    force_override = True,
                    display_name=r"$\rm{Pressure}/k_B$",
                    take_log=True,
                    units = 'K/cm**3')

def add_new_star_ages(snap):
    pf = snap.ds
    code = snap.code
    FormationTimeType_to_use=FormationTimeType_to_use_dict[code]
    PartType_Star_to_use=PartType_Star_to_use_dict[code]
    def new_age(field,data):
        formation_z=1./pf.arr(data[(PartType_Star_to_use,FormationTimeType_to_use)],'')-1.
        yt_cosmo=Cosmology()
        return pf.current_time.in_units('Gyr')-yt_cosmo.t_from_z(formation_z).in_units('Gyr')
    pf.add_field(('agora_stars','particle_ages'), 
                 function=new_age,  
                 sampling_type = 'particle',
                 units="Gyr")
        
def add_radial_distance_fields(snap):
    ds = snap.ds
    def radial_distance(field,data):
        xdist = data['gas','x']-snap.center_x
        ydist = data['gas','y']-snap.center_y
        zdist = data['gas','z']-snap.center_z
        return np.sqrt(xdist**2+ydist**2+zdist**2)
    ds.add_field(('gas','radial_distance'),function = radial_distance,
                 units = 'kpc',sampling_type = snap.sampling_type,force_override = True)
    def star_radial_distance(field,data):
        xdist = data['agora_stars','particle_position_x']-snap.center_x
        ydist = data['agora_stars','particle_position_y']-snap.center_y
        zdist = data['agora_stars','particle_position_z']-snap.center_z
        return np.sqrt(xdist**2+ydist**2+zdist**2)
    #ds.add_field(('agora_stars','radial_distance'),function = star_radial_distance,
    #             units = 'kpc',sampling_type = 'particle',force_override = True)