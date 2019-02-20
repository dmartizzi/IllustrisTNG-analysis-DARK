import sys
import os.path
import numpy as np
import scipy
import matplotlib as mpl
import requests
import h5py

#def get(path, params=None):
    # make HTTP GET request to path
#    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
#    r.raise_for_status()

#    if r.headers['content-type'] == 'application/json':
#        return r.json() # parse json responses automatically
#    return r

def get(path, params=None):
    # make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically

    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r

def set_simulation(base_url,sim_name):
    r = get(base_url)
    names = [sim['name'] for sim in r['simulations']]
    if sim_name in names:
        i = names.index(sim_name)
        sim = get( r['simulations'][i]['url'] )
        return sim
    else:
        print "The selected simulation name doesn't exist!!!"
        sys.exit()

def set_redshift(snap_url,sim_z):
    # Find the appropriate snapshot for the selected redshift
    isnap = len(snaps)-1
    while sim_z >= snaps[isnap]["redshift"]:
	isnap -= 1
    return isnap, snaps[isnap]["number"], snaps[isnap]["redshift"]

def get_subhalos(snap,mass_min,mass_max,lim):
    subs = get(snap["subhalos"], {"limit" : lim})
    # Set a stellar mass range
    search_query = "?mass__gt=" + str(mass_min) + "&mass__lt=" + str(mass_max)
    search_url = snap["url"]+"/subhalos/"+search_query
    subhalos = get(search_url, {"limit" : lim})
    return subhalos
    
if __name__ == "__main__":
    # Input parameters 
    # Select simulation and redshift 
    sim_name = sys.argv[1] ## Name of the simulation
    sim_z = float(sys.argv[2]) ## Redshift of the snapshot
    mass_min = 10**float(sys.argv[3]) / 1e10 / 0.6774 ## Minimum stellar mass in 1e10 Msun/h
    mass_max = 10**float(sys.argv[4]) / 1e10 / 0.6774 ## Maximum stellar mass in 1e10 Msun/h
    lim = int(sys.argv[5]) ## Maximum number of subhalos in list
    
    # read API key
    if os.path.exists('.api_key'):
        api_key = np.loadtxt('.api_key',dtype="S50")
        api_key = str(api_key)
    else:
        print "Please, create file .api_key with your API key."
        sys.exit()
        
    # URL to the API
    base_url = "http://www.tng-project.org/api/"
    headers = {"api-key": api_key}

    # Select the simulation
    print("Reading %s simulation."%sim_name)
    sim = set_simulation(base_url,sim_name)

    # List of simulation snapshots
    snaps = get( sim["snapshots"] )
    
    # Find the appropriate snapshot for the selected redshift
    snap_url = sim["snapshots"]
    isnap, inumber, snap_z = set_redshift(snap_url,sim_z)
    snap = get(snaps[isnap]["url"])
    print("Selecting snapshot number %d at redshift z = %e."%(inumber,snap_z))

    # Get subhalos in a given stellar mass range 
    subhalos = get_subhalos(snap,mass_min,mass_max,lim)
    nselected = subhalos["count"]
    print("# of selected subhalos with %f < log10(Mstar) < %f = %d."%(np.log10(1e10*mass_min*0.6774),np.log10(1e10*mass_max*0.6774),nselected))
    # Subhalo ids
    ids = [ subhalos['results'][i]['id'] for i in range(0,len(subhalos['results'][:])) ]
    # Loop on all selected subhalo ids
    for idsub in ids:
        sub_url = snap["url"]+"/subhalos/"+str(idsub)
        subhalo = get(sub_url)
        # If you want the subhalo cutout
        vdisp = subhalo["veldisp"]
        print idsub, vdisp
        cutout_url = subhalo["cutouts"]["subhalo"]
        # If you want the parent halo cutout
        #cutout = subhalo["cutouts"]["parent_halo"]
        print cutout_url

        cutout_request = {'stars':'Coordinates,Masses'}
        cutout_hdf5 = get(cutout_url, cutout_request)

        dic = {"gas" : "0", "stars" : "4", "dm" : "1", "bhs" : "5"}
        for ptype in ["gas", "stars", "dm", "bhs"]:
            with h5py.File(cutout_hdf5,'r') as f:
                try :
                    x = f['PartType'+dic[ptype]]['Coordinates'][:,0] 
                    y = f['PartType'+dic[ptype]]['Coordinates'][:,1] 
                    dens = np.log10(f['PartType'+dic[ptype]]['Masses'][:])
                    print x
                    print y 
                except :                
                    print("Particle type %s does not exist for halo %d"%(ptype,idsub))
        del cutout_url
        del cutout_hdf5
