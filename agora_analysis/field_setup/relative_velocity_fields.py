import numpy as np
from unyt import unyt_array
def add_cylindrical_velocity_fields(snap):
    def cylindrical_velocity_z(field,data):
        axis = snap.L
        x_comp_vel = (data['gas','velocity_x']-snap.bvx)*axis[0]
        y_comp_vel = (data['gas','velocity_y']-snap.bvy)*axis[1]
        z_comp_vel = (data['gas','velocity_z']-snap.bvz)*axis[2]
        total_vel = x_comp_vel+y_comp_vel+z_comp_vel
        x_comp_pos = (data['gas','x']-snap.center_x)*axis[0]
        y_comp_pos = (data['gas','y']-snap.center_y)*axis[1]
        z_comp_pos = (data['gas','z']-snap.center_z)*axis[2]
        above_disk = ((x_comp_pos+y_comp_pos+z_comp_pos>0.0)*2)-1
        return total_vel*above_disk
    snap.ds.add_field(('gas','cylindrical_velocity_z'),function = cylindrical_velocity_z,
                 units = 'km/s',sampling_type = snap.sampling_type,force_override = True)
    def star_cylindrical_velocity_z(field,data):
        axis = snap.L
        x_comp_vel = (data['agora_stars','particle_velocity_x']-snap.bvx)*axis[0]
        y_comp_vel = (data['agora_stars','particle_velocity_y']-snap.bvy)*axis[1]
        z_comp_vel = (data['agora_stars','particle_velocity_z']-snap.bvz)*axis[2]
        total_vel = x_comp_vel+y_comp_vel+z_comp_vel
        x_comp_pos = (data['agora_stars','particle_position_x']-snap.center_x)*axis[0]
        y_comp_pos = (data['agora_stars','particle_position_y']-snap.center_y)*axis[1]
        z_comp_pos = (data['agora_stars','particle_position_z']-snap.center_z)*axis[2]
        above_disk = ((x_comp_pos+y_comp_pos+z_comp_pos>0.0)*2)-1
        return total_vel*above_disk
    snap.ds.add_field(('agora_stars','cylindrical_velocity_z'),function = cylindrical_velocity_z,
                 units = 'km/s',sampling_type = 'particle',force_override = True)
    
def add_radial_velocity_fields(snap):
    ds = snap.ds
    def radial_velocity(field,data):
        xdist = data['gas','x']-snap.center_x
        ydist = data['gas','y']-snap.center_y
        zdist = data['gas','z']-snap.center_z
        xvel = data['gas','velocity_x']-snap.bvx
        yvel = data['gas','velocity_y']-snap.bvy
        zvel = data['gas','velocity_z']-snap.bvz
        dot = xdist*xvel+ydist*yvel+zdist*zvel
        return dot/data['gas','radial_distance']
    ds.add_field(('gas','radial_velocity'),function = radial_velocity,
                 units = 'km/s',sampling_type = snap.sampling_type,force_override = True)
    def star_radial_velocity(field,data):
        xdist = data['agora_stars','particle_position_x']-snap.center_x
        ydist = data['agora_stars','particle_position_y']-snap.center_y
        zdist = data['agora_stars','particle_position_z']-snap.center_z
        xvel = data['agora_stars','particle_velocity_x']-snap.bvx
        yvel = data['agora_stars','particle_velocity_y']-snap.bvy
        zvel = data['agora_stars','particle_velocity_z']-snap.bvz
        dot = xdist*xvel+ydist*yvel+zdist*zvel
        return dot/data['agora_stars','radial_distance']
    ds.add_field(('agora_stars','radial_velocity'),function = star_radial_velocity,
                 units = 'km/s',sampling_type = 'particle',force_override = True)