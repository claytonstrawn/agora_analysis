from agora_analysis import AgoraSnapshot
from agora_analysis.field_setup.common_fields import add_gas_fields,\
                                                        add_temperature_fields
import numpy as np
import yt

class BadCalculationError(Exception):
    def __init__(self,message):
        self.message = message

def cold_ang_mom(snap):
    ds = snap.ds
    Rvir = snap.Rvir
    center = snap.center
    sp = ds.sphere(center,0.15*Rvir)
    cold = ds.cut_region(sp, ["obj[('gas', 'agora_temperature')] < 1e3"])
    L = cold.quantities.angular_momentum_vector(use_gas = True, use_particles= False)
    if np.isnan(L).all():
        L = sp.quantities.angular_momentum_vector(use_gas = True, use_particles= False)
    print(L)
    return L
    
def calculate_L(code,redshift,interested_redshifts = 'all',close_enough =0.01,\
                snap=None,throw_errors = True):
    if interested_redshifts != 'all' and \
            not min(np.abs(redshift-np.array(interested_redshifts))<close_enough):
        return np.array([np.nan]*3)
    if snap is None:
        snap = AgoraSnapshot('AGORA_%s_CR'%code,redshift)
    if ('L' in snap.__dict__) and (not np.isnan(snap.L).all()):
        return snap.L.v*snap.Lmag.v
    if 'ds' not in snap.__dict__:
        snap.load_snapshot()
        add_gas_fields(snap)
        add_temperature_fields(snap)
    try:
        L = cold_ang_mom(snap).v
    except:
        L = np.array([np.nan]*3)
        if throw_errors:
            raise BadCalculationError('failed to create L for %s at z=%.2f'%(code,redshift))
        else:
            print('failed to create L for %s at z=%.2f'%(code,redshift))
    return L

def bulk_velocity(snap):
    ds = snap.ds
    Rvir = snap.Rvir
    center = snap.center
    sp = ds.sphere(center,Rvir)
    if snap.sampling_type == 'cell':
        bv = sp.quantities.bulk_velocity(use_gas = True, use_particles= True)
    elif snap.sampling_type == 'particle':
        bv = sp.quantities.bulk_velocity(use_gas = True, use_particles= True)
    print(bv.in_units('km/s'))
    return bv

def calculate_bulk_vel(code,redshift,interested_redshifts = 'all',close_enough =0.01,snap = None):
    if interested_redshifts != 'all' and \
            not min(np.abs(redshift-np.array(interested_redshifts))<close_enough):
        return np.array([np.nan]*3)    
    if snap is None:
        snap = AgoraSnapshot('AGORA_%s_CR'%code,redshift)
    if ('bv' in snap.__dict__) and (not np.isnan(snap.bv).any()):
        return snap.bv.v
    if 'ds' not in snap.__dict__:
        snap.load_snapshot()
    try:
        bv = bulk_velocity(snap).in_units('km/s').v
    except:
        print('failed to create bv for %s at z=%.2f'%(code,redshift))
        bv = np.array([np.nan]*3)
    return bv