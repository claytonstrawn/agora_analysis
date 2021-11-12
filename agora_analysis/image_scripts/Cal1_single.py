import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from matplotlib import cm
import copy
from yt import ProjectionPlot,YTArray
from yt.utilities.exceptions import YTFieldNotFound
from agora_analysis.agora_snapshot import AgoraSnapshot,NoSimulationError
from agora_analysis.common_fields import add_gas_fields,add_temperature_fields,add_resolution_fields
from mpl_toolkits.axes_grid1 import AxesGrid

add_nametag                      = 1           # 0/1       = OFF/ON

figure_width                     = 200          # in kpc

orig_codes_labels = {'art':'ART-I','ramses':'RAMSES','enzo':'ENZO',
                     'gadget':'GADGET-3','gear':'GEAR','changa':'CHANGA',
                     'gizmo':'GIZMO'}
orig_codes = ['art','ramses','enzo','gadget','gear','changa','gizmo']

plots = ['density','temperature']

grid_codes = ['art','enzo','ramses']
particle_codes = ['gadget','gear','gizmo','changa']

fig_density_maps = [None]*3
grid_density_maps = [None]*3
fig_temp_maps = [None]*3
grid_temp_maps = [None]*3
fig_res_maps = [None]*3
grid_res_maps = [None]*3

def Cal1_single(save = True,show = True,codes = orig_codes,n_axes = 3):
    for ax_plot in range(0,n_axes):
        """fig_density_maps[ax_plot] = plt.figure(figsize=(100.*len(codes)/len(orig_codes),20))
        grid_density_maps[ax_plot] = AxesGrid(fig_density_maps[ax_plot], (0.01,0.01,0.99,0.99), 
                                        nrows_ncols = (1, len(codes)), axes_pad = 0.02, 
                                        add_all = True, share_all = True,
                                        label_mode = "1", cbar_mode = "single", 
                                        cbar_location = "right", cbar_size = "2%", 
                                        cbar_pad = 0.02)
        fig_temp_maps[ax_plot] = plt.figure(figsize=(100.*len(codes)/len(orig_codes),20))
        grid_temp_maps[ax_plot] = AxesGrid(fig_temp_maps[ax_plot], (0.01,0.01,0.99,0.99), 
                                        nrows_ncols = (1, len(codes)), axes_pad = 0.02, 
                                        add_all = True, share_all = True,
                                        label_mode = "1", cbar_mode = "single", 
                                        cbar_location = "right", cbar_size = "2%", 
                                        cbar_pad = 0.02)"""
        fig_res_maps[ax_plot] = plt.figure(figsize=(100.*len(codes)/len(orig_codes),20))
        grid_res_maps[ax_plot] = AxesGrid(fig_res_maps[ax_plot], (0.01,0.01,0.99,0.99), 
                                        nrows_ncols = (1, len(codes)), axes_pad = 0.02, 
                                        add_all = True, share_all = True,
                                        label_mode = "1", cbar_mode = "single", 
                                        cbar_location = "right", cbar_size = "2%", 
                                        cbar_pad = 0.02)

    for code_num,code in enumerate(codes):
        fullname = 'AGORA_%s_C1'%code
        redshift = 7.0
        ####################################
        #        PRE-ANALYSIS STEPS        #
        ####################################

        # LOAD DATASETS
        try:
            snap = AgoraSnapshot(fullname,redshift,Rvir_method = 'average')
            snap.load_snapshot()
            pf = snap.ds
            add_gas_fields(snap)
            add_temperature_fields(snap)
            add_resolution_fields(snap)

        except PermissionError:
            print('unable to load snap %s (permission denied), skipping...'%fullname)
            continue
        except NoSimulationError:
            print('unable to load snap %s (simulation file missing), skipping...'%fullname)
            continue
        except:
            return snap
        
        pf.coordinates.x_axis[1] = 0
        pf.coordinates.y_axis[1] = 2
        pf.coordinates.x_axis['y'] = 0
        pf.coordinates.y_axis['y'] = 2
        
        if code in grid_codes:
            fw = figure_width
            proj_region = pf.box(snap.center - YTArray([fw, fw, fw], 'kpc'),
                          snap.center + YTArray([fw, fw, fw], 'kpc')) 
            # projected images made using a (2*figure_width)^3 box for AMR codes
        else:
            proj_region = pf.all_data()
        
        if draw_density_DF >= 1:
            sp = pf.sphere(snap.center, (0.5*figure_width, "kpc"))
            p6 = ProfilePlot(sp, ("gas", "agora_density"),  ("gas", "agora_mass"), 
                             weight_field=None, n_bins=50, x_log=True, accumulation=False)
            p6.set_log("agora_mass", True)
            p6.set_xlim(1e-29, 1e-20)
            density_DF_xs[time].append(p6.profiles[0].x.in_units('g/cm**3').d)
            density_DF_profiles[time].append(p6.profiles[0]["agora_mass"].in_units('Msun').d)
            """else:
                # Because ParticleProfilePlot doesn't exist, I will do the following trick.  
                p6 = ProfilePlot(sp, (PartType_Gas_to_use, "Density_2"),  
                                 (PartType_Gas_to_use, "Mass_2"), 
                                 weight_field=None, n_bins=50, 
                                 x_log=True, accumulation=False)
                p6.set_log((PartType_Gas_to_use,"Mass_2"), True)
                p6.set_xlim(1e-29, 1e-20)
                density_DF_xs[time].append(p6.profiles[0].x.in_units('g/cm**3').d)
                density_DF_profiles[time].append(p6.profiles[0]["Mass_2"].in_units('Msun').d)
            """
            ####################################
            #        POST-ANALYSIS STEPS       #
            ####################################

    # SAVE FIGURES
    if save:
        for ax_plot in range(0,n_axes):
            fig_density_map = fig_density_maps[ax_plot]
            fig_density_map.savefig("Cal1_single_density_z%d_ax%d"%(redshift,ax_plot), 
                                  bbox_inches='tight',
                                  pad_inches=0.03, dpi=300)
            fig_temp_map = fig_temp_maps[ax_plot]
            fig_temp_map.savefig("Cal1_single_temp_z%d_ax%d"%(redshift,ax_plot), 
                                  bbox_inches='tight',
                                  pad_inches=0.03, dpi=300)
            fig_res_map = fig_res_maps[ax_plot]
            fig_res_map.savefig("Cal1_single_resolution_z%d_ax%d"%(redshift,ax_plot), 
                                  bbox_inches='tight',
                                  pad_inches=0.03, dpi=300)
    # SHOW FIGURES (if interactive)
    if show:
        plt.show()

if __name__ == '__main__':
    Cal1_single(save = True,show = False)