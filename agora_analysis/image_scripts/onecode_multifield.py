from agora_analysis import AgoraSnapshot,NotCloseEnoughError
from agora_analysis.field_setup.main import load_necessary_fields
from unyt import unyt_array
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from quasarscan.preprocessing.code_specific_setup import load_and_setup
from quasarscan.utils import ion_lists
import yt,trident
from agora_analysis.image_scripts.default_cmaps_zlims import default_proj_cmaps,\
                                                            default_proj_zlims,\
                                                            default_proj_weights,\
                                                            default_slc_cmaps,\
                                                            default_slc_zlims,\
                                                            default_thin_cmaps,\
                                                            default_thin_zlims,\
                                                            default_thin_weights,\
                                                            choose_default

def _plot_1code_4fields(proj_or_slc,code,fields,redshift,circles=[1.0],axis = 0,
                      cmaps = 'default',zlims = 'default',weights = 'default',
                      width = 'default',throw_errors = 'warn',
                      textsize = 30, textcolor = 'white',circlecolor = 'white',
                     simnum = 'CR',offset_center = None,load_method = 'AGORA',
                     depth_frac = 0.1,test_one_code = False):    
    for i,field in enumerate(fields):
        if isinstance(field,str):
            fields[i] = ('gas',field)
        elif field is None:
            fields[i] = ('None','None')
    if cmaps=='default':
        cmaps=['default']*4
    if zlims=='default':
        zlims=['default']*4
    if weights=='default':
        weights=['default']*4

    default_name = '%s_%d_%.1f_%s_%s%s%s%s.png'%(proj_or_slc,axis,redshift,code,fields[0][1],\
                                                fields[1][1],fields[2][1],\
                                                fields[3][1])
    assert len(fields) == 4,'_plot_1code_4fields only accepts 4 fields'
    
    snap = AgoraSnapshot('AGORA_%s_%s'%(code,simnum),redshift)
    snap.load_snapshot()
    if load_method == 'AGORA':
        load_necessary_fields(snap,fields)
        ds = snap.ds
    elif load_method == 'quasarscan':
        path = snap.lookup_snap_path()
        ds,_ = load_and_setup(path,code,ions = ion_lists.agoraions)
        
    fig  = plt.figure(figsize=(100,20))
    grid = AxesGrid(fig,(0.01,0.01,0.99,0.99), \
                      nrows_ncols = (2,2), axes_pad = 0.02,\
                      add_all = True, share_all = True,\
                      label_mode = "1", cbar_mode = "each", \
                      cbar_location = "right", cbar_size = "6%", \
                      cbar_pad = 0.02)    
    
    for i,field in enumerate(fields):
        if field == ('None','None'):
            continue
        print('plotting %s field %s for code "%s"'%(proj_or_slc,field,code))
        cmap = choose_default(field,'default_%s_cmaps'%proj_or_slc,cmaps[i])
        zlim = choose_default(field,'default_%s_zlims'%proj_or_slc,zlims[i])
        if proj_or_slc in ['proj','thin']:
            weight = choose_default(field,'default_%s_weights'%proj_or_slc,weights[i])
        if not isinstance(width,unyt_array):
            width = snap.Rvir if width == 'default' else width * snap.Rvir
        if offset_center is None:
            offset_center = unyt_array([0,0,0],'kpc')
        elif not isinstance(offset_center,unyt_array):
            offset_center = unyt_array(offset_center,'kpc')
        try:
            if proj_or_slc == 'proj':
                proj_region = ds.box(snap.center+offset_center - width,snap.center+offset_center + width)
                p = yt.ProjectionPlot(ds, axis, field, center = snap.center+offset_center, 
                                   weight_field = weight,width = width, fontsize=9,
                                   data_source = proj_region)
            elif proj_or_slc == 'slc':
                p = yt.SlicePlot(ds, axis, field, center = snap.center+offset_center, 
                                  width = width, fontsize=9)
            elif proj_or_slc == 'thin':
                depth = width*depth_frac
                to_add = unyt_array([width.v,width.v,width.v],width.units)
                to_add[axis] = depth.v
                proj_region = ds.box(snap.center+offset_center - to_add,\
                                     snap.center+offset_center + to_add)
                p = yt.ProjectionPlot(ds, axis, field, center = snap.center+offset_center, 
                                   weight_field = weight,width = width, fontsize=9,
                                   data_source = proj_region)
            for c in circles:
                circle_radius = ds.arr(c*snap.Rvir.v,'kpc')
                p.annotate_sphere(snap.center,circle_radius,circle_args={'color':circlecolor})
            p.annotate_text((0.2, 0.8), field[1], coord_system="axis",
                                text_args = {'size':textsize,'color':textcolor})
            if i==0:
                p.annotate_text((0.2, 0.2), code, coord_system="axis",
                                text_args = {'size':textsize,'color':textcolor})
        except Exception as e:
            if throw_errors == 'warn':
                print('unable to plot %s field %s for code "%s"'%(proj_or_slc,field,code))
            elif throw_errors == True:  
                raise e
            continue
        if zlim == (0,1) or field[1] in ['radial_velocity']:
            p.set_log(field,False)
        if zlim is not None:
            p.set_zlim(field,zlim[0],zlim[1])
        p.set_cmap(field, cmap)
        plot = p.plots[field]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
        p._setup_plots()
    return fig,default_name

def slc_1code_4fields(code,fields,redshift,**kwargs):
    return _plot_1code_4fields('slc',code,fields,redshift,**kwargs)
    
def proj_1code_4fields(code,fields,redshift,**kwargs):
    return _plot_1code_4fields('proj',code,fields,redshift,**kwargs)
    
def thin_1code_4fields(code,fields,redshift,**kwargs):
    return _plot_1code_4fields('thin',code,fields,redshift,**kwargs)