######################################################################
#
#  UNIFIED ANALYSIS SCRIPT FOR COSMORUN SIMULATION FOR THE AGORA PROJECT
#  (Based on "disk script" by Ji-hoon Kim)
#
######################################################################


import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from matplotlib import cm
import copy
from yt import ProjectionPlot,YTArray
from yt.utilities.exceptions import YTFieldNotFound
from agora_analysis.agora_snapshot import AgoraSnapshot,NoSimulationError
from agora_analysis.common_fields import add_gas_fields,add_temperature_fields
from mpl_toolkits.axes_grid1 import AxesGrid

add_nametag                      = 1           # 0/1       = OFF/ON

figure_width                     = 200          # in kpc

codes_labels = ['ART-I','RAMSES','ENZO','GADGET-3','GEAR','CHANGA','GIZMO']
codes = ['art','ramses','enzo','gadget','gear','changa','gizmo']

grid_codes = ['art','enzo','ramses']
particle_codes = ['gadget','gear','gizmo','changa']

fig_short_maps = [None]*3
grid_short_maps = [None]*3

def Cal1_multipanel(save = True,show = True):
    for ax_plot in range(0,3):
        fig_short_maps[ax_plot]  = plt.figure(figsize=(100,20))
        grid_short_maps[ax_plot] = AxesGrid(fig_short_maps[ax_plot],(0.01,0.01,0.99,0.99), \
                                      nrows_ncols = (2, len(codes)), axes_pad = 0.02,\
                                      add_all = True, share_all = True,\
                                      label_mode = "1", cbar_mode = "edge", \
                                      cbar_location = "right", cbar_size = "6%", \
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
        except PermissionError:
            print('unable to load snap %s (permission denied), skipping...'%fullname)
            continue
        except NoSimulationError:
            print('unable to load snap %s (simulation file missing), skipping...'%fullname)
            continue
        
        pf.coordinates.x_axis[1] = 0
        pf.coordinates.y_axis[1] = 2
        pf.coordinates.x_axis['y'] = 0
        pf.coordinates.y_axis['y'] = 2
        
        for ax_plot in range(0,3):
            fig_short_map = fig_short_maps[ax_plot]
            grid_short_map = grid_short_maps[ax_plot]


        ####################################
        #      MAIN ANALYSIS ROUTINES      #
        ####################################

            if code in grid_codes:
                fw = figure_width
                proj_region = pf.box(snap.center - YTArray([fw, fw, fw], 'kpc'),
                              snap.center + YTArray([fw, fw, fw], 'kpc')) 
                # projected images made using a (2*figure_width)^3 box for AMR codes
            else:
                proj_region = pf.all_data()

            rvir_circle = snap.Rvir

            #first plot
            my_cmap = copy.copy(cm.get_cmap('viridis'))
            my_cmap.set_bad(my_cmap(0)) # cmap range [0, 1)
            my_cmap.set_under(my_cmap(0))
            p11 = ProjectionPlot(pf, ax_plot, 
                                     ("gas", "agora_density"), 
                                     center = snap.center, 
                                     data_source=proj_region, 
                                     width = (figure_width, 'kpc'), 
                                     weight_field = None, 
                                     fontsize=12)
            p11.set_zlim(("gas", "agora_density"), 1.5e-4, 1e-1)
            p11.set_cmap(("gas", "agora_density"), my_cmap)
            p11.set_colorbar_label(("gas", "agora_density"), 
                                   "$\Sigma_{\mathrm{gas}}$"+
                                   "$\left(\\frac{\mathrm{g}}{\mathrm{cm}^2}\\right)$")
            plot = p11.plots[("gas", "agora_density")]
            plot.figure = fig_short_map
            plot.axes = grid_short_map[code_num].axes
            if code_num == 0: 
                plot.cax = grid_short_map.cbar_axes[0]
            p11._setup_plots()
            circle2 = plt.Circle((0, 0), rvir_circle, color='k', 
                                 linewidth=1,fill=False, ls='--')
            grid_short_map[code_num].axes.add_artist(circle2)
            if add_nametag:
                at = AnchoredText("%s" % (codes_labels[code_num]), 
                                  loc=2, prop=dict(size=9), frameon=True)
                grid_short_map[code_num].axes.add_artist(at)

            #second plot
            my_cmap = copy.copy(cm.get_cmap('magma'))
            my_cmap.set_bad(my_cmap(0.6))
            # (log(1e4) - log(1e1)) / (log(1e6) - log(1e1)) = 0.6
            my_cmap.set_under(my_cmap(0))
            p12 = ProjectionPlot(pf, ax_plot, 
                                 ("gas", "agora_temperature"), 
                                 center = snap.center, 
                                 data_source=proj_region, 
                                 width = (figure_width, 'kpc'), 
                                 weight_field = ("gas", "agora_density"), fontsize=12)
            p12.set_zlim(("gas", "agora_temperature"), 5e3, 1e6)
            p12.set_cmap(("gas", "agora_temperature"), my_cmap)
            p11.set_colorbar_label(("gas", "agora_temperature"), 
                       "$\mathrm{Temperature}\left(\mathrm{K}\right)")
            plot2 = p12.plots[("gas", "agora_temperature")]
            plot2.figure = fig_short_map
            plot2.axes = grid_short_map[len(codes)+code_num].axes
            circle2 = plt.Circle((0, 0), rvir_circle, color='k', 
                                 linewidth=1,fill=False, ls='--')
            grid_short_map[len(codes)+code_num].axes.add_artist(circle2)
            if code_num == 0: 
                plot2.cax = grid_short_map.cbar_axes[1]
            p12._setup_plots()
            circle2 = plt.Circle((0, 0), rvir_circle, color='k', 
                                 linewidth=1,fill=False, ls='--')
            grid_short_map[len(codes)+code_num].axes.add_artist(circle2)
        ####################################
        #        POST-ANALYSIS STEPS       #
        ####################################

    # SAVE FIGURES
    if save:
        for ax_plot in range(0,3):
            fig_short_map = fig_short_maps[ax_plot]
            fig_short_map.savefig("Cal1_multipanel_z%d_ax%d"%(redshift,ax_plot), 
                                  bbox_inches='tight',
                                  pad_inches=0.03, dpi=300)
    # SHOW FIGURES (if interactive)
    if show:
        plt.show()

if __name__ == '__main__':
    Cal2_multipanel(save = True,show = False)