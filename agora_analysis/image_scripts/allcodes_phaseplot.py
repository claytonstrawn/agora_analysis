from agora_analysis import AgoraSnapshot,NotCloseEnoughError
from agora_analysis.field_setup.main import load_necessary_fields
from unyt import unyt_array
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import yt,trident
from agora_analysis.image_scripts.default_cmaps_zlims import choose_default,\
                                                            default_proj_cmaps,\
                                                            default_proj_zlims,\
                                                            default_proj_weights,\
                                                            default_thin_cmaps,\
                                                            default_thin_zlims,\
                                                            default_thin_weights,\
                                                            default_slc_cmaps,\
                                                            default_slc_zlims
from agora_analysis import codes,official_names



def _plot_phaseplot_allcodes(region,redshift,xfield = 'number_density',yfield = 'temperature',\
                             zfield = 'mass',test_one_code = False,throw_errors = 'warn',\
                             simnum = 'CR'):
    if isinstance(test_one_code,bool):
        if test_one_code:
            code_list = [None,None,None,None,None,None,'art']
        else:
            code_list = codes
    elif isinstance(test_one_code,int):
        code_list = [None]*(8-test_one_code)+codes[:test_one_code]
    temp_fields = [xfield,yfield,zfield]
    fields = []
    for field in temp_fields:
        if not isinstance(field,tuple):
            field = ('gas',field)
            fields.append(field)
    cmap = choose_default(fields[2],'default_phase_cmaps','default')
    weight = choose_default(fields[2],'default_phase_weights','default')
    frac = choose_default(fields[2],'default_phase_fractional','default')
    xlim = choose_default(fields[0],'default_phase_xlims','default')
    ylim = choose_default(fields[1],'default_phase_ylims','default')
    default_name = 'phaseplot_%s_%s_%.1f.png'%(zfield,region,redshift)
    fig  = plt.figure(figsize=(1,4))
    grid = AxesGrid(fig,(0.01,0.5,0.99,3.), (len(code_list),1), \
                     axes_pad= 0.05,add_all = True, share_all = True,\
                      label_mode= "L", cbar_mode = "single",\
                      cbar_location = "right", cbar_size = "0.5%", \
                      cbar_pad = 0.005,aspect=False)
    for i,code in enumerate(code_list):
        if code is None:
            continue
        else:
            print('plotting fields %s for code "%s" in region %s'%(zfield,code,region))
        try:
            snap = AgoraSnapshot('AGORA_%s_%s'%(code,simnum),redshift)
        except NotCloseEnoughError:
            continue
        snap.load_snapshot()
        load_necessary_fields(snap,field)
        ds = snap.ds
        if region == 'CGM':
            data_source = ds.sphere(snap.center,snap.Rvir)-ds.sphere(snap.center,0.15*snap.Rvir)
        elif region == 'gal':
            data_source = ds.sphere(snap.center,0.15*snap.Rvir)
        elif region == 'zoom-in':
            data_source = ds.sphere(snap.center,5*snap.Rvir)
        elif isinstance(region,(float,int)):
            data_source = ds.sphere(snap.center,region*snap.Rvir)
        elif isinstance(region,(tuple)):
            data_source = ds.sphere(snap.center,region[1]*snap.Rvir)-\
                            ds.sphere(snap.center,region[0]*snap.Rvir)
        try:
            p = yt.PhasePlot(data_source, xfield, yfield, zfield,\
                             weight_field = weight,fractional = frac)
        except Exception as e:
            if throw_errors == 'warn':
                print('unable to plot %s field %s for code "%s"'%(proj_or_slc,field,code))
            elif throw_errors == True:  
                raise e
            continue
        p.set_cmap(zfield, cmap)
        p.set_xlim(xlim[0],xlim[1])
        p.set_ylim(ylim[0],ylim[1])
        plot = p.plots[zfield]
        plot.figure = fig
        plot.axes = grid[i].axes
        if i == len(code_list)-1:
            plot.cax = grid.cbar_axes[0]
        p._setup_plots()
        grid[i].axes.set_ylabel(official_names[code])
            
    return fig,default_name

def phaseplot_allcodes(region,redshift,**kwargs):
    return _plot_phaseplot_allcodes(region,redshift,**kwargs)