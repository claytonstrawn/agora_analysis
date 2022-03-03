default_proj_cmaps = {
                     ('gas','density'): None,
                     ('gas','temperature'): 'magma',
                     ('gas','number_density'): None,
                    }  
default_proj_zlims = {
                     ('gas','PI_OVI'):(0,1),
                     #('gas','OVI_PI_dens'):(1e13,1e16),
                     ('gas','OVI_CI_dens'):(1e13,1e16),
                     ('gas','density'): (1e-4,1e2),
                     ('gas','temperature'): (1e3,1e7),
                     ('gas','number_density'): (1e20,1e26),
                    }  
default_proj_weights = {
                        ('gas','density'): None,
                        ('gas','temperature'): ('gas','density'),
                        ('gas','number_density'): None,
                       }  



default_slc_cmaps = {('gas','PI_OVI_number_density'):('Spectral', 'Diverging', 7),
                     ('gas','CI_OVI_number_density'):('Spectral', 'Diverging', 7),
                     ('gas','number_density'): ('Spectral', 'Diverging', 8),
                     ('gas','temperature'): 'RED TEMPERATURE',
                     ('gas','agora_metallicity'): 'fake_PRGn',
                     ('gas','agora_pressure'): 'RdYlGn',
                     ('gas','density'): None,
                     ('gas','O_p5_ion_fraction'):('PuOr', 'Diverging', 8),
                     ('gas','radial_velocity'):'seismic',
                     ('gas','cylindrical_velocity_z'):'seismic',
                    }  
default_slc_zlims = {
                     ('gas','PI_OVI_number_density'):(1e-15,1e-8),
                     ('gas','CI_OVI_number_density'):(1e-15,1e-8),
                     ('gas','O_p5_ion_fraction'):(1e-8,1e0),
                     ('gas','density'): None,
                     ('gas','temperature'): (1e3,1e7),
                     ('gas','number_density'): (1e-5,1e-1),
                     ('gas','agora_metallicity'): (1e-3,1e-1),
                     ('gas','agora_pressure'): (1e1,1e6),
                     ('gas','radial_velocity'): (-1000,1000),
                     ('gas','cylindrical_velocity_z'):(-1000,1000),
                    }  

default_thin_weights = {('gas','temperature'): ('gas','density'),
                        ('gas','agora_metallicity'): ('gas','density'),
                        ('gas','agora_pressure'): ('gas','density'),
                        ('gas','O_p5_ion_fraction'): ('gas','density'),
                        ('gas','CI_OVI'):('gas','O_p5_number_density'),
                        ('gas','radial_velocity'): ('gas','density'),
                        ('gas','density'): None,
                        ('gas','temperature'): ('gas','density'),
                        ('gas','number_density'): None,
                       }

default_thin_cmaps = {('gas','PI_OVI_number_density'):('Spectral', 'Diverging', 7),
                     ('gas','CI_OVI_number_density'):('Spectral', 'Diverging', 7),
                     ('gas','number_density'): ('Spectral', 'Diverging', 8),
                     ('gas','temperature'): 'RED TEMPERATURE',
                     ('gas','agora_metallicity'): 'fake_PRGn',
                     ('gas','agora_pressure'): 'RdYlGn',
                     ('gas','density'): None,
                     ('gas','O_p5_ion_fraction'):('PuOr', 'Diverging', 8),
                     ('gas','CI_OVI'):'coolwarm',
                     ('gas','radial_velocity'):'PiYG',
                    }  
default_thin_zlims = {('gas','number_density'):(1e18,1e22),
                     ('gas','temperature'): (1e3,1e7),
                     ('gas','agora_metallicity'): (1e-3,1e-1),
                     ('gas','agora_pressure'): (1e1,1e6),                     
                    ('gas','PI_OVI_number_density'):(1e8,1e15),
                    ('gas','CI_OVI_number_density'):(1e8,1e15),
                    ('gas','O_p5_ion_fraction'):(1e-8,1e0),
                    ('gas','CI_OVI'):(0,1),
                    ('gas','radial_velocity'): (-400,400),
                    }  

default_phase_cmaps = {('gas','mass'):None,
                    }  
default_phase_weights = {('gas','mass'):None,
                    }  
default_phase_fractional = {('gas','mass'):True,
                    }  
default_phase_xlims = {('gas','number_density'):(4e-5,1e3),
                    }  
default_phase_ylims = {('gas','temperature'):(5,8e7),
                    }  

def choose_default(field,dict_type,v,printing = True):
    if v!= 'default':
        return v
    dictionary = eval(dict_type)
    try:
        to_ret = dictionary[field]
    except KeyError:
        if printing:
            print('field %s not found in default dictionary %s! Consider adding...'%\
                  (field,dict_type))
        to_ret = None
    if to_ret == 'fake_PRGn':
        make_fake_PRGn()
    return to_ret


def make_fake_PRGn():
    print('adding "fake_PRGn" to colormap list...')
    import yt
    palettablestrlist = "#762a83,#9970ab,#c2a5cf,#e7d4e8,#f7f7f7,#d9f0d3,#a6dba0,#5aae61,#1b7837".split(',')

    def hex2rgb(h):
        h = h.lstrip('#')
        return tuple(int(h[i:i+2], 16)/256.0 for i in (0, 2, 4))

    colormaplist=[None]*9
    for i in range(9):
        colormaplist[i] = (hex2rgb(palettablestrlist[i]),1)

    yt.make_colormap(colormaplist,
                     name='fake_PRGn', interpolate=False)