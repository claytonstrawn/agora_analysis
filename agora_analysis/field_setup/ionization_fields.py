import numpy as np
import yt,trident
from yt.utilities.physical_constants import mh
from trident import add_ion_fields
import os
import roman
from copy import copy

def read_table(filename = "default"):
    if filename == 'default':
        _ROOT = os.path.abspath(os.path.dirname(__file__))
        path = os.path.join(_ROOT,'CI_PI_cutoff_tables.txt')
        f=open(path,'r')
    else:
        f=open(filename,'r')
    lines = f.readlines()[6:]
    f.close()
    ions = []
    redshifts = []
    rhos = {}
    ts = {}
    for line in lines:
        if len(line.split())==2:
            current_ion = line.replace('\n','')
            ions.append(current_ion)
        elif len(line.split())>2:
            redshift = line.split(',')[0]
            if not redshift in redshifts:
                redshifts.append(redshift)
            rho_or_t = line.split()[1].replace(':','')
            start_of_list = line.index('[')
            if rho_or_t == 'rho':
                rhos[(current_ion,redshift)] = np.fromstring(line[start_of_list+1:-2],sep = ' ')
            elif rho_or_t == 't':
                ts[(current_ion,redshift)] = np.fromstring(line[start_of_list+1:-2],sep = ' ')
    return ions,redshifts,rhos, ts
table = read_table()

def cutoffs_for_ion_at_redshift(ion,redshift):
    if isinstance(redshift,str):
        redshift = float(redshift)
    ions,redshifts,rhos,ts = table
    if ion in ions and str(redshift) in redshifts:
        return rhos[(ion,str(redshift))],ts[(ion,str(redshift))]
    else:
        for i,z in enumerate(redshifts):
            if float(z)>redshift:
                break
            if i == len(redshifts)-1:
                return rhos[(ion,str(z))],ts[(ion,str(z))]
        z_below = redshifts[i-1]
        z_above = redshifts[i]
        fraction_z_below = (float(z_above)-redshift)/(float(z_above)-float(z_below))
        fraction_z_above = 1-fraction_z_below
        if len(rhos[(ion,z_below)]) > len(rhos[(ion,z_above)]):
            adjusted_rho_below = np.append(rhos[(ion,z_below)][:len(rhos[(ion,z_above)])-1],rhos[(ion,z_below)][-1])
            adjusted_rho_above = rhos[(ion,z_above)]
            adjusted_t_below = np.append(ts[(ion,z_below)][:len(ts[(ion,z_above)])-1],ts[(ion,z_below)][-1])
            adjusted_t_above = ts[(ion,z_above)]
        if len(rhos[(ion,z_below)]) < len(rhos[(ion,z_above)]):
            adjusted_rho_below = rhos[(ion,z_below)]
            adjusted_rho_above = np.append(rhos[(ion,z_above)][:len(rhos[(ion,z_below)])-1],rhos[(ion,z_above)][-1])
            adjusted_t_below = ts[(ion,z_below)]
            adjusted_t_above = np.append(ts[(ion,z_above)][:len(ts[(ion,z_below)])-1],ts[(ion,z_above)][-1])
        if len(rhos[(ion,z_below)]) == len(rhos[(ion,z_above)]):
            adjusted_rho_below = rhos[(ion,z_below)]
            adjusted_rho_above = rhos[(ion,z_above)]
            adjusted_t_below = ts[(ion,z_below)]
            adjusted_t_above = ts[(ion,z_above)]
        rho_final = adjusted_rho_below*fraction_z_below+adjusted_rho_above*fraction_z_above
        t_final = adjusted_t_below*fraction_z_below+adjusted_t_above*fraction_z_above
        return rho_final,t_final

def make_PI_CI_funcs(ion,redshift):
    rhos,ts = cutoffs_for_ion_at_redshift(ion,redshift)
    def PI_ion(field, data):
        #0 if CI, 1 if PI
        temps = data['gas','temperature']
        comp = np.zeros(temps.shape)
        for i in range(1,len(ts)):
            if i==1:
                mask = temps<ts[i]
                comp[mask] = np.inf
            elif i < len(ts)-1:
                mask = np.logical_and(temps<ts[i],temps>ts[i-1])
                comp[mask] = np.sqrt(rhos[i]*rhos[i-1])
            else: 
                mask = temps>ts[i]
                comp[mask] = -np.inf
        tr = yt.YTArray((data['gas','density']/mh)<comp).astype(float)
        return tr
    def CI_ion(field,data):
        PI_field_name = ('gas','PI_%s'%ion.replace(' ',''))
        return 1-(data[PI_field_name])
    return PI_ion,CI_ion

def binary_field_PI_CI(ion,field):
    if isinstance(field,str):
        field = ('gas',field)
    def PI_ion_field(f,data):
        return data[field]*data['gas','PI_%s'%ion.replace(' ','')]
    def CI_ion_field(f,data):
        return data[field]*data['gas','CI_%s'%ion.replace(' ','')]
    return PI_ion_field,CI_ion_field

atoms = ['He', 'Li', 'Be','Ne', 'Na', 'Mg', 'Al', 'Si',\
          'Cl', 'Ar', 'Ca', 'Sc', 'Ti', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',\
         'Cu', 'Zn', 'H', 'B', 'C', 'N', 'O', 'F', 'P', 'S','K','V']
def is_ion_type(fname):
    toret = []
    words = fname.split('_')
    if words[0] in ['PI','CI']:
        word = words[1]
    else:
        word = words[0]
    for atom in atoms:
        if word.startswith(atom):
            element = atom
            try:
                ionization = roman.fromRoman(word.split(element)[1])
                return element+' '+word.split(element)[1]
            except roman.InvalidRomanNumeralError:
                return False
    return False

def ion_str_to_tup(s):
    s = s.replace(' ','')
    for atom in atoms:
        if s.startswith(atom):
            element = atom
            ionization = roman.fromRoman(s.split(element)[1])
            return element,ionization
        
def add_PI_CI_fields(snap,ion):
    PI_field_name = 'PI_%s'%(ion.replace(' ',''))
    CI_field_name = 'CI_%s'%(ion.replace(' ',''))
    PI_ion,CI_ion = make_PI_CI_funcs(ion,snap.redshift)
    snap.ds.add_field(('gas',PI_field_name),
                   sampling_type=snap.sampling_type,
                   function=PI_ion,
                   units='')
    snap.ds.add_field(('gas',CI_field_name),
                   sampling_type=snap.sampling_type,
                   function=CI_ion,
                   units='')

def add_number_density_fields(snap,ion,trident_method = 'default'):
    atom,ion_num = ion_str_to_tup(ion)
    if trident_method == 'uniform':
        high_priority_fields = ['H_nuclei_density',"%s_p%d_number_density" % (atom, ion_num-1),
                        "%s_nuclei_density" % atom,"%s_nuclei_mass" % atom,
                        "%s_metallicity" % atom,'H_nuclei_mass']
        old_field_list = copy(snap.ds.derived_field_list)
        for field in high_priority_fields:
            if field in old_field_list:
                snap.ds.derived_field_list.remove(field)
    trident.add_ion_fields(snap.ds,atom)
    
    if trident_method == 'uniform':
        snap.ds.derived_field_list = old_field_list
        
    def number_density_rename(field,data):
        return data['gas','%s_p%d_number_density'%(atom,ion_num-1)]
    snap.ds.add_field(('gas','%s_number_density'%ion.replace(' ','')),
                   sampling_type=snap.sampling_type,
                   function=number_density_rename,
                   units='1/cm**3')
    f_PI,f_CI = binary_field_PI_CI(ion,'%s_number_density'%ion.replace(' ',''))
    snap.ds.add_field(('gas','PI_%s_number_density'%ion.replace(' ','')),
                   sampling_type=snap.sampling_type,
                   function=f_PI,
                   units='1/cm**3')
    snap.ds.add_field(('gas','CI_%s_number_density'%ion.replace(' ','')),
                   sampling_type=snap.sampling_type,
                   function=f_CI,
                   units='1/cm**3')