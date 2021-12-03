from agora_analysis import AgoraSnapshot,NotCloseEnoughError
from agora_analysis.field_setup.main import load_all_fields
from unyt import unyt_array
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import yt

default_proj_cmaps = {
                     ('gas','density'): None,
                     ('gas','temperature'): 'magma',
                     ('gas','number_density'): None,
                    }  
default_proj_zlims = {
                     ('gas','density'): (1e-4,1e2),
                     ('gas','temperature'): (1e3,1e7),
                     ('gas','number_density'): (1e20,1e26),
                    }  
default_proj_weights = {
                        ('gas','density'): None,
                        ('gas','temperature'): ('gas','density'),
                        ('gas','number_density'): None,
                       }  
default_slc_cmaps = {
                     ('gas','density'): None,
                     ('gas','temperature'): 'magma',
                     ('gas','number_density'): None,
                    }  
default_slc_zlims = {
                     ('gas','density'): None,
                     ('gas','temperature'): (1e3,1e7),
                     ('gas','number_density'): (1e-4,1e3),
                    }  


def choose_default(field,dict_type,v,printing = True):
    if v!= 'default':
        return v
    dictionary = eval(dict_type)
    try:
        to_ret = dictionary[field]
    except KeyError:
        if printing:
            print('field %s not found in default dictionary %s! Consider adding...'%(field,dict_type))
        to_ret = None
    return to_ret

def _plot_1field_2rows(proj_or_slc,field,redshift,width,circles=[1],axis = 0,
                      cmap = 'default',zlims = 'default',weight = 'default',
                      textsize = 30, textcolor = 'white',circlecolor = 'white',
                     simnum = 'CR',offset_center = None):
    code_list = ['art','enzo','ramses',None,'gadget','gear','gizmo','changa']
    if not isinstance(field,tuple):
        field = ('gas',field)
    default_name = 'proj_%d_%.1f_%s.png'%(axis,redshift,field[1])
    
    cmap = choose_default(field,'default_%s_cmaps'%proj_or_slc,cmap)
    zlims = choose_default(field,'default_%s_zlims'%proj_or_slc,zlims)
    if proj_or_slc == 'proj':
        weight = choose_default(field,'default_proj_weights',weight)
    fig  = plt.figure(figsize=(100,20))
    grid = AxesGrid(fig,(0.01,0.01,0.99,0.99), \
                      nrows_ncols = (2,4), axes_pad = 0.02,\
                      add_all = True, share_all = True,\
                      label_mode = "1", cbar_mode = "single", \
                      cbar_location = "right", cbar_size = "6%", \
                      cbar_pad = 0.02)
    for i,code in enumerate(code_list):
        if code is None:
            continue
        try:
            snap = AgoraSnapshot('AGORA_%s_%s'%(code,simnum),redshift)
        except NotCloseEnoughError:
            continue
        snap.load_snapshot()
        load_all_fields(snap)
        if not isinstance(width,unyt_array):
            width = width * snap.Rvir
        if offset_center is None:
            offset_center = unyt_array([0,0,0],'kpc')
        elif not isinstance(offset_center,unyt_array):
            offset_center = unyt_array(offset_center,'kpc')
        if proj_or_slc == 'proj':
            proj_region = snap.ds.box(snap.center - width,snap.center + width)
            p = yt.ProjectionPlot(snap.ds, axis, field, center = snap.center+offset_center, 
                               weight_field = weight,width = width, fontsize=9,
                               data_source = proj_region)
        elif proj_or_slc == 'slc':
            p = yt.SlicePlot(snap.ds, axis, field, center = snap.center+offset_center, 
                              width = width, fontsize=9)
        for c in circles:
            circle_radius = c*snap.Rvir
            p.annotate_sphere(snap.center,circle_radius,circle_args={'color':circlecolor})
            p.annotate_text((0.2, 0.8), code, coord_system="axis",
                            text_args = {'size':textsize,'color':textcolor})
        if zlims is not None:
            p.set_zlim(field,zlims[0],zlims[1])
        p.set_cmap(field, cmap)
        plot = p.plots[field]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[0]
        p._setup_plots()
    return fig,default_name

def slc_1field_2rows(field,redshift,width,**kwargs):
    _plot_1field_2rows('slc',field,redshift,width,**kwargs)
    
def proj_1field_2rows(field,redshift,width,**kwargs):
    _plot_1field_2rows('proj',field,redshift,width,**kwargs