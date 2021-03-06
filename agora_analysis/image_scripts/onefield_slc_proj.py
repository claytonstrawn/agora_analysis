from agora_analysis import AgoraSnapshot,NotCloseEnoughError,official_names
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
                                                            default_slc_zlims


def add_PI_CI_dens_field(code,ds):
    def _PI_dens(field,data):
        return data['gas','O_p5_number_density'] * data['gas','PI_OVI'] 
    def _CI_dens(field,data):
        return data['gas','O_p5_number_density'] * data['gas','CI_OVI']
    if code in ['art','enzo','ramses']:
        sampling_type = 'cell'
    else:
        sampling_type = 'particle'
    ds.add_field(('gas','OVI_PI_dens'),
                           sampling_type=sampling_type,
                           function=_PI_dens,
                           units='cm**-3')
    ds.add_field(('gas','OVI_CI_dens'),
                           sampling_type=sampling_type,
                           function=_CI_dens,
                           units='cm**-3')

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
                      textsize = 30, textcolor = 'white',circlecolor = 'white',throw_errors = True,
                     simnum = 'CR',offset_center = None,load_method = 'AGORA',test_one_code = False):
    code_list = [None,'enzo','art','ramses','gadget','gear','gizmo','changa']
    if test_one_code:
        code_list = [None,None,None,None,None,None,None,'art']

    
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
        else:
            print('plotting %s field %s for code "%s"'%(proj_or_slc,field,code))
        try:
            snap = AgoraSnapshot('AGORA_%s_%s'%(code,simnum),redshift)
        except NotCloseEnoughError:
            continue
        snap.load_snapshot()
        if load_method == 'AGORA':
            load_necessary_fields(snap,field)
            ds = snap.ds
        elif load_method == 'quasarscan':
            path = snap.lookup_snap_path()
            ds,_ = load_and_setup(path,code,ions = ion_lists.agoraions)
        if not isinstance(width,unyt_array):
            width = width * snap.Rvir
        if offset_center is None:
            offset_center = unyt_array([0,0,0],'kpc')
        elif not isinstance(offset_center,unyt_array):
            offset_center = unyt_array(offset_center,'kpc')
        try:
            if proj_or_slc == 'proj':
                proj_region = ds.box(snap.center - width,snap.center + width)
                p = yt.ProjectionPlot(ds, axis, field, center = snap.center+offset_center, 
                                   weight_field = weight,width = width, fontsize=9,
                                   data_source = proj_region)
            elif proj_or_slc == 'slc':
                p = yt.SlicePlot(ds, axis, field, center = snap.center+offset_center, 
                                  width = width, fontsize=9)
            for c in circles:
                circle_radius = ds.arr(c*snap.Rvir.v,'kpc')
                p.annotate_sphere(snap.center,circle_radius,circle_args={'color':circlecolor})
            p.annotate_text((0.2, 0.8), code, coord_system="axis",
                                text_args = {'size':textsize,'color':textcolor})
        except Exception as e:
            if throw_errors == 'warn':
                print('unable to plot %s field %s for code "%s" because of %s'%(proj_or_slc,field,code,e))
                continue
            elif throw_errors == True:
                raise e
        if zlims == (0,1):
            p.set_log(field,False)
        if zlims is not None:
            p.set_zlim(field,zlims[0],zlims[1])
        p.set_cmap(field, cmap)
        plot = p.plots[field]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[0]
        p._setup_plots()
        grid[i].axes.text(0,0,official_names[code])
    return fig,default_name

def slc_1field_2rows(field,redshift,width,**kwargs):
    _plot_1field_2rows('slc',field,redshift,width,**kwargs)
    
def proj_1field_2rows(field,redshift,width,**kwargs):
    _plot_1field_2rows('proj',field,redshift,width,**kwargs)
    
