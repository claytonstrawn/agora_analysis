from agora_analysis import AgoraSnapshot,NotCloseEnoughError
from agora_analysis.field_setup.main import load_necessary_fields
from unyt import unyt_array
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
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

def _run_plot(proj_or_slc,snap,width,axis,weight,field,depth_frac):
    ds = snap.ds
    if proj_or_slc == 'thin':
        depth = width*depth_frac
    else:
        depth = width
    if isinstance(axis,(int,float)):
        func = {'proj':yt.ProjectionPlot,'slc':yt.SlicePlot,'thin':yt.ProjectionPlot}[proj_or_slc]
        north = None
        proj_ax = np.zeros(3)
        proj_ax[axis] = 1.0
    elif axis in ['face','edge']:
        func = {'proj':yt.OffAxisProjectionPlot,
                'slc':yt.OffAxisSlicePlot,
                'thin':yt.OffAxisProjectionPlot
               }[proj_or_slc]
        if snap.sampling_type == 'particle' and proj_or_slc == 'slc': 
            func = yt.OffAxisProjectionPlot
            depth = 0.005*width
        if axis == 'face':
            axis = snap.L
            north = None
        elif axis == 'edge':
            axis = np.cross(np.array([0,0,1]),snap.L)
            north = snap.L
        proj_ax = axis
    proj_region = ds.disk(snap.center,proj_ax,2*width,depth)
    args = {'center' : snap.center, 
            'width': (float(width.v),width.units), 
            'fontsize':9,
            'data_source': proj_region}
    if proj_or_slc!= 'slc':
        args['weight_field'] = weight
    elif (proj_or_slc == 'slc') and (func is yt.OffAxisProjectionPlot):
        args['weight_field'] = ('gas','density')
    if north is not None:
        args['north_vector'] = north
    return func(ds, axis, field, **args)

def _plot_nfields_allcodes(proj_or_slc,fields,redshift,width,circles=[1],axis = 0,
                      textsize = 30, textcolor = 'white',circlecolor = 'white',
                      simnum = 'CR',depth_frac = 0.1,test_one_code = False,throw_errors = 'warn'):
    if isinstance(test_one_code,bool):
        if test_one_code:
            code_list = [None,None,None,None,None,None,'art']
        else:
            code_list = codes
    elif isinstance(test_one_code,int):
        code_list = [None]*(7-test_one_code)+codes[:test_one_code]
    elif isinstance(test_one_code,str):
        code_list = [None,None,None,None,None]+[test_one_code]
    elif isinstance(test_one_code,list):
        code_list = [None]*(7-len(test_one_code))+test_one_code
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
    default_name = '%s_%s_%.1f_%s.png'%(proj_or_slc,axis,redshift,fields_names)
    fig  = plt.figure(figsize=(100,20))    
    grid = AxesGrid(fig,(0.01,0.01,0.99,0.99), (len(fields),len(code_list)), \
                     axes_pad= 0.02,add_all = True, share_all = True,\
                      label_mode= "L", cbar_mode = "edge",\
                      cbar_location = "right", cbar_size = "6%", \
                      direction = 'row',\
                      cbar_pad = 0.02)
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
        load_necessary_fields(snap,fields)
        ds = snap.ds
        if not isinstance(width,unyt_array):
            width = width * snap.Rvir
        for j,field in enumerate(fields):
            cmap = cmaps[j]
            zlim = zlims[j]
            weight = weights[j]
            try:
                p = _run_plot(proj_or_slc,snap,width,axis,weight,field,depth_frac)
                for c in circles:
                    circle_radius = ds.arr(c*snap.Rvir.v,'kpc')
                    p.annotate_sphere(snap.center,circle_radius,circle_args={'color':circlecolor})
            except Exception as e:
                if throw_errors == 'warn':
                    print('unable to plot %s field %s for code "%s"'%(proj_or_slc,field,code))
                elif throw_errors == True:  
                    raise e
                continue
            if zlim == (0,1) or field[1] in ['radial_velocity','cylindrical_velocity_z']:
                p.set_log(field,False)
            if zlim is not None:
                p.set_zlim(field,zlim[0],zlim[1])
            p.set_cmap(field, cmap)
            plot = p.plots[field]
            plot.figure = fig
            plot.axes = grid[i+j*len(code_list)].axes
            if i == len(code_list)-1:
                plot.cax = grid.cbar_axes[j]
            p._setup_plots()
            if j==0:
                grid[i+j*len(code_list)].axes.set_title(official_names[code])
    return fig,default_name

def slc_nfields_allcodes(fields,redshift,width,**kwargs):
    _plot_nfields_allcodes('slc',fields,redshift,width,**kwargs)
    
def proj_nfields_allcodes(fields,redshift,width,**kwargs):
    _plot_nfields_allcodes('proj',fields,redshift,width,**kwargs)
    
def thin_nfields_allcodes(fields,redshift,width,**kwargs):
    _plot_nfields_allcodes('thin',fields,redshift,width,**kwargs)