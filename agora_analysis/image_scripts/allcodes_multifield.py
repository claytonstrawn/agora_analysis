from agora_analysis import AgoraSnapshot,NotCloseEnoughError
from agora_analysis.field_setup.main import load_necessary_fields
from unyt import unyt_array
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from quasarscan.preprocessing.code_specific_setup import load_and_setup
from quasarscan.utils import ion_lists
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

def _plot_nfields_allcodes(proj_or_slc,fields,redshift,width,circles=[1],axis = 0,
                      textsize = 30, textcolor = 'white',circlecolor = 'white',
                      simnum = 'CR',offset_center = None,load_method = 'AGORA',
                      depth_frac = 0.1,test_one_code = False,throw_errors = 'warn',**kwargs):
    if isinstance(test_one_code,int):
        code_list = [None]*(8-test_one_code)+codes[:test_one_code]
    elif test_one_code:
        code_list = [None,None,None,None,None,None,None,'art']
    else:
        code_list = codes
    temp_fields = list(fields)
    fields,fields_names,cmaps,zlims,weights = [],'',[],[],[]
    for field in temp_fields:
        if not isinstance(field,tuple):
            field = ('gas',field)
            fields.append(field)
        fields_names += field[1]
        cmaps.append(choose_default(field,'default_%s_cmaps'%proj_or_slc,'default'))
        zlims.append(choose_default(field,'default_%s_zlims'%proj_or_slc,'default'))
        if proj_or_slc in ['proj','thin']:
            weights.append(choose_default(field,'default_%s_weights'%proj_or_slc,'default'))
        else:
            weights.append(None)
    default_name = '%s_%d_%.1f_%s.png'%(proj_or_slc,axis,redshift,fields_names)
    fig  = plt.figure(figsize=(100,20))
    grid = AxesGrid(fig,(0.01,0.01,0.99,0.99), (len(fields),len(code_list)), \
                      **kwargs)
    for i,code in enumerate(code_list):
        if code is None:
            continue
        else:
            print('plotting %s fields %s for code "%s"'%(proj_or_slc,fields,code))
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
        for j,field in enumerate(fields):
            cmap = cmaps[j]
            zlim = zlims[j]
            weight = weights[j]
            try:
                if proj_or_slc == 'proj':
                    proj_region = ds.box(snap.center - width,snap.center + width)
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
            plot.axes = grid[i+j*len(codes)].axes
            if i == (len(code_list) - 1):
                plot.cax = grid.cbar_axes[i+j*len(codes)]
            p._setup_plots()
            if j==0:
                grid[i+j*len(codes)].axes.set_title(official_names[code])
            
    return fig,default_name

def slc_nfields_allcodes(fields,redshift,width,**kwargs):
    _plot_nfields_allcodes('slc',fields,redshift,width,**kwargs)
    
def proj_nfields_allcodes(fields,redshift,width,**kwargs):
    _plot_nfields_allcodes('proj',fields,redshift,width,**kwargs)
    
def thin_nfields_allcodes(fields,redshift,width,**kwargs):
    _plot_nfields_allcodes('thin',fields,redshift,width,**kwargs)