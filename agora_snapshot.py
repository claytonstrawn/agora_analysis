import os
import numpy as np
from unyt import unyt_array,unyt_quantity
from agora_analysis.agora_metadata.all_file_locations import locations
from agora_analysis.utils import read_metadata,NotCloseEnoughError,NoMetadataError
import yt



class NotImplementedError(Exception):
    def __init__(self,message):
        self.message = message
        
class NoSimulationError(Exception):
    def __init__(self,message):
        self.message = message

grid_codes = ['art','enzo','ramses']
particle_codes = ['gadget','gear','gizmo','changa']
codes = grid_codes+particle_codes

class AgoraSnapshot(object):
    def __init__(self,fullname,redshift,Rvir_method = 'average',yt_log_level = 'default'):
        if yt_log_level != 'default':
            yt.set_log_level(yt_log_level)
        self.fullname = fullname
        self.simulation = fullname.split('_')[0]
        assert self.simulation == 'AGORA','Simulation name must start with "AGORA"!'
        self.code = fullname.split('_')[1]
        self.simnum = fullname.split('_')[2]
        if self.simnum in ['C3','CR']:
            self.stars_in = True
        elif self.simnum in ['C1','C2']:
            self.stars_in = False
        path_to_metadata = os.path.expanduser("~/agora_analysis/agora_metadata")
        foldernames = {'C1':'Cal1','C2':'Cal2','C3':'Cal3','CR':'CosmoRun'}
        foldername = foldernames[self.simnum]
        self.metadata_location_folder = os.path.join(path_to_metadata,foldername)
        self.approx_redshift = redshift
        self.lookup_metadata()
        self.lookup_Rvir(Rvir_method = Rvir_method)
        #self.Rvir = self.lookup_Rvir(fullname,self.redshift,Rvir_method)
        #can add other calculated variables (Mstar, sfr, Mgas, bulk_velocity, angular_momentum, etc.)
        #can make loading the actual dataset optional -- you can plot stuff like Rvir, Mstar by code 
        #& redshift without loading the full snapshots.
        
    def __repr__(self):
        toret = self.fullname+' at redshift '+str(self.redshift)
        if 'ds' not in self.__dict__.keys():
            toret+= ' (not loaded)'
        return toret
        
    def load_snapshot(self,path_to_snap = None):
        if path_to_snap is None:
            path_to_snap = self.lookup_snap_path()
        self.ds = self.ytload(path_to_snap)
        self.edit_center_units()
        self.set_up_fields(self.code,self.ds)

    def lookup_metadata(self):
        z,snapnum,center_x,center_y,center_z,Rvir = read_metadata(self.code,self.simnum,self.approx_redshift)
        self.redshift = unyt_quantity(z,'')
        self.snapnum = float(snapnum)
        self.center_x = unyt_quantity(center_x,'kpc')
        self.center_y = unyt_quantity(center_y,'kpc')
        self.center_z = unyt_quantity(center_z,'kpc')
        self.center = unyt_array([self.center_x,self.center_y,self.center_z])
        self.specific_Rvir = unyt_quantity(Rvir,'kpc')
        
    def lookup_Rvir(self,Rvir_method = 'average'):
        Rvirs = []
        for code in codes:
            try:
                _,_,_,_,_,Rvir = read_metadata(code,self.simnum,self.approx_redshift)
                Rvirs.append(Rvir)
            except NotCloseEnoughError:
                continue
            except NoMetadataError:
                continue
        if Rvir_method == 'average':
            Rvir = unyt_quantity(np.average(Rvirs),'kpc')
            self.Rvir = Rvir
        else:
            self.Rvir = self.specific_Rvir
        
    def lookup_snap_path(self):
        root_loc = "/project/projectdirs/agora/paper_CGM"
        #name of the folder where the final version of the code is stored
        location = locations["%s_%s"%(self.code,self.simnum)]
        if location is None:
            raise NoSimulationError('Location for "%s" files is unknown!'%\
                                   self.fullname)
        #this translates the different code file naming conventions
        #it does not account for auxiliary files (see `ytload` below)
        self.snapnum = float(self.snapnum)
        if self.code == 'art':
            #there are two methods for loading art snaps, by 
            #expansion factor a (<1) or by snap number (>=1)
            if self.snapnum < 1:
                filename = '10MpcBox_csf512_a%s.d'%self.snapnum
            else:
                filename = '10MpcBox_csf512_%05d.d'%self.snapnum
        elif self.code == 'changa':
            name = {'CR':'11-l-cal-IV','C1':'cal0'}[self.simnum]
            filename = '%s.%06d'%(name,self.snapnum)
        elif self.code == 'enzo':
            filename = 'RD%04d/RD%04d'%(self.snapnum,self.snapnum)
        elif self.code == 'gadget':
            filename = 'snapshot_%03d/snapshot_%03d.0.hdf5'%(self.snapnum,self.snapnum)
        elif self.code == 'gear':
            filename = 'snapshot_%04d.hdf5'%self.snapnum
        elif self.code == 'gizmo':
            filename = 'snapshot_%03d.hdf5'%self.snapnum
        elif self.code == 'ramses':
            filename = 'output_%05d/info_%05d.txt'%(self.snapnum,self.snapnum)
        path_to_snap = os.path.join(root_loc,location,filename)
        return path_to_snap
    
    def ytload(self,path_to_snap):
        print("Loading snapshot %s, redshift %1.3f stored at %s"%(self.fullname,self.redshift,path_to_snap))
        if self.code == 'art':
            projectdir = path_to_snap.split('10M')[0]
            timestep = path_to_snap.split("_")[-1][:-2]
            h = projectdir+"PMcrd_%s.DAT"%timestep
            d = projectdir+"PMcrs0_%s.DAT"%timestep
            s = projectdir+"stars_%s.dat"%timestep
            if 'a' in timestep:
                h,d,s = h.replace('_a','a'),d.replace('_a','a'),s.replace('_a','a')
            if not self.stars_in:
                s = None
            ds = yt.load(path_to_snap,file_particle_header=h,\
                                file_particle_data=d,\
                                file_particle_stars=s)
        elif self.code == 'changa':
            ds = yt.load(path_to_snap)
        elif self.code == 'enzo':
            ds = yt.load(path_to_snap)
        elif self.code == 'gadget':
            unit_base = {'length':(1.0, "Mpccm/h")}
            ds = yt.load(path_to_snap,unit_base = unit_base)
        elif self.code == 'gear':
            ds = yt.load(path_to_snap)
        elif self.code == 'gizmo':
            ds = yt.load(path_to_snap)
        elif self.code == 'ramses':
            ds = yt.load(path_to_snap)
        return ds
    
    def recenter(self,approx_center,recentering_distance,use_gas = True,use_particles = False):
        sp = self.ds.sphere(approx_center, (recentering_distance,'kpc'))
        approx_center = sp.quantities.center_of_mass(use_gas=use_gas,use_particles=use_particles)
        sp = self.ds.sphere(approx_center, (recentering_distance*0.5,'kpc'))
        approx_center = sp.quantities.center_of_mass(use_gas=use_gas,use_particles=use_particles)
        sp = self.ds.sphere(approx_center, (recentering_distance*0.25,'kpc'))
        center = sp.quantities.center_of_mass(use_gas=use_gas,use_particles=use_particles)
        return center.in_units('kpc')
    
    def edit_center_units(self):
        #using ds.arr instead of unyt_array
        #means the units can access "unitary" units
        #instead of just physical units
        cx,cy,cz = self.center_x.in_units('kpc'),self.center_y.in_units('kpc'),self.center_z.in_units('kpc')
        c = self.center.in_units('kpc')
        rv = self.Rvir.in_units('kpc')
        self.center_x = self.ds.quan(cx.v,'kpc')
        self.center_y = self.ds.quan(cy.v,'kpc')
        self.center_z = self.ds.quan(cz.v,'kpc')
        self.center = self.ds.arr(c.v,'kpc')
        self.Rvir = self.ds.quan(rv.v,'kpc')
    
    def set_up_fields(self,code,ds):
        #codes refer to things by different names
        #we need standard names, add aliases for each
        #e.g. 'stars','particle_mass' instead of 'star','particle_mass' in RAMSES
        #now all plotting scripts can refer to 'agora_stars','agora_mass'
        pass  