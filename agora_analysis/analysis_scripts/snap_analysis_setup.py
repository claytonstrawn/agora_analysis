import agora_analysis
from agora_analysis.field_setup.main import load_necessary_fields
import os 

def load_snapshot_for_analysis(name,redshift,load_for = 'all',overwrite = False):
    code = name.split('_')[1]
    if load_for == 'all':
        check = ['mass','flux','resolution']
    elif isinstance(load_for,str):
        check = [load_for]
    else:
        check = load_for
    skip = True
    for item in check:
        if not os.path.exists('agora_analysis_quantities/%s_data_%s_%.3f.txt'%(item,code,redshift)):
            skip = False
    if overwrite or skip:
        return None
    try:
        snap = agora_analysis.AgoraSnapshot(name,redshift = redshift)
        snap.load_snapshot()
        list_of_used_fields = [('gas','agora_mass'),('agora_stars','agora_mass'),\
                               ('gas','metal_mass'),('agora_stars','metal_mass'),\
                               ('gas','agora_cell_volume'),('gas','agora_particle_volume'),\
                               ('agora_stars','metal_movement'),('gas','metal_movement'),\
                               ('agora_stars','star_mass_movement'),('gas','gas_movement')]
        load_necessary_fields(snap,list_of_used_fields)
    except agora_analysis.utils.NotCloseEnoughError:
        print('skipping %s z=%1.3f!'%(name,redshift))
        return False
    except agora_analysis.utils.NoMetadataError:
        print('skipping %s z=%1.3f!'%(name,redshift))
        return False
    except Exception as e:
        raise e
    return snap
        