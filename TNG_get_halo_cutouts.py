import sys
import os.path
import numpy as np
import scipy
import matplotlib as mpl
import matplotlib.pyplot as plt
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
        print("Box size = %f kpc/h"%sim["boxsize"])
        print("DM particle mass = %f x 1e10 Msun/h"%sim["mass_dm"])
        print("Gas particle mass = %f x 1e10 Msun/h"%sim["mass_gas"])
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

def get_subhalos(snap,mass_min,mass_max,vmax,lim):
    subs = get(snap["subhalos"], {"limit" : lim})
    # Set a stellar mass range
    search_query = "?mass__gt=" + str(mass_min) + "&mass__lt=" + str(mass_max) + "&vmax__gt=" + str(vmax)
    search_url = snap["url"]+"/subhalos/"+search_query
    subhalos = get(search_url, {"limit" : lim})
    return subhalos

def make_projected_prof(cutout_hdf5,ptype_dic,rc,proj_axis,rmin,rmax,nbins,mdmpart):
    # Initialize profile
    nvars = int(1+3*len(ptype_dic.keys()))
    shape = (nbins,nvars)
    projected_prof = np.zeros(shape)
    iprof_dic = {"dm":1,"gas":2,"stars":3,"bhs":4}
    
    # Radial Bins
    rmin = np.log10(rmin)
    rmax = np.log10(rmax)
    dlogr = (rmax-rmin)/float(nbins)
    for i in range(0,nbins):
        projected_prof[i,0] = rmin + 0.5*dlogr + float(i)*dlogr
    
    # Loop on all particle types 
    for ptype in ptype_dic.keys():
        with h5py.File(cutout_hdf5,'r') as f:
	    rnorm = np.array([(f['PartType'+ptype_dic[ptype]]['Coordinates'][:,0]-rc[0])*proj_axis[0], \
                              (f['PartType'+ptype_dic[ptype]]['Coordinates'][:,1]-rc[1])*proj_axis[1], \
	                      (f['PartType'+ptype_dic[ptype]]['Coordinates'][:,2]-rc[2])*proj_axis[2]])
            rnorm = rnorm.T
            rproj = np.sqrt((f['PartType'+ptype_dic[ptype]]['Coordinates'][:,0]-rc[0]-rnorm[:,0])**2+ \
                            (f['PartType'+ptype_dic[ptype]]['Coordinates'][:,1]-rc[1]-rnorm[:,1])**2+ \
                            (f['PartType'+ptype_dic[ptype]]['Coordinates'][:,2]-rc[2]-rnorm[:,2])**2)
            rproj = np.log10(rproj)
            del rnorm
            try :
                mass = f['PartType'+ptype_dic[ptype]]['Masses'][:]
            except :
                #print("Particle type %s does not exist for halo %d"%(ptype,idsub))
                mass = np.array(len(rproj[:])*[mdmpart])
                
            # Surface densities in bins
            for i in range(0,nbins):
                #print 'bin ',i, iprof_dic[ptype]
                ishell = np.logical_and(projected_prof[i,0]-0.5*dlogr < rproj, rproj <= projected_prof[i,0]+0.5*dlogr)
                icum = rproj <= projected_prof[i,0]
                shell = np.pi*(10**(2*(projected_prof[i,0]+0.5*dlogr))-10**(2*(projected_prof[i,0]-0.5*dlogr)))*1e6 # (pc/h)**2
                surf = np.pi*(10**(2*projected_prof[i,0]))*1e6 # (pc/h)**2
                ivar = iprof_dic[ptype]
		projected_prof[i,ivar] = np.sum(mass[icum])*1e10/surf # Msun/h/(pc/h)**2
                projected_prof[i,ivar+len(ptype_dic.keys())] = np.sum(mass[ishell])*1e10/shell # Msun/h/(pc/h)**2
                projected_prof[i,ivar+2*len(ptype_dic.keys())] = projected_prof[i,ivar] - projected_prof[i,ivar+len(ptype_dic.keys())]

            del rproj,mass,ishell,icum
        
    return projected_prof 



if __name__ == "__main__":
    # Input parameters 
    # Select simulation and redshift 
    sim_name = sys.argv[1] ## Name of the simulation
    sim_z = float(sys.argv[2]) ## Redshift of the snapshot
    mass_min = 10**float(sys.argv[3]) / 1e10 / 0.6774 ## Minimum stellar mass in 1e10 Msun/h
    mass_max = 10**float(sys.argv[4]) / 1e10 / 0.6774 ## Maximum stellar mass in 1e10 Msun/h
    vmax = float(sys.argv[5]) # Vmax cutoff in km/s
    lim = int(sys.argv[5]) ## Maximum number of subhalos in list
    
    # read API key
    if os.path.exists('.api_key'):
        api_key = np.loadtxt('.api_key',dtype="S50")
        api_key = str(api_key)
    elif os.path.exists('~/.api_key'):
        api_key = np.loadtxt('~/.api_key',dtype="S50")
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

    # DM particle mass
    mdmpart = sim["mass_dm"]
    
    # List of simulation snapshots
    snaps = get( sim["snapshots"] )
    
    # Find the appropriate snapshot for the selected redshift
    snap_url = sim["snapshots"]
    isnap, inumber, snap_z = set_redshift(snap_url,sim_z)
    snap = get(snaps[isnap]["url"])
    print("Selecting snapshot number %d at redshift z = %e."%(inumber,snap_z))

    # Get subhalos in a given stellar mass range 
    subhalos = get_subhalos(snap,mass_min,mass_max,vmax,lim)
    nselected = subhalos["count"]
    print("# of selected subhalos with %f < log10(Mstar) < %f = %d."%(np.log10(1e10*mass_min*0.6774),np.log10(1e10*mass_max*0.6774),nselected))
    # Subhalo ids
    ids = [ subhalos['results'][i]['id'] for i in range(0,len(subhalos['results'][:])) ]
    flag = 0
    cutout_dic = {}
    # Loop on all selected subhalo ids
    for idsub in ids:
        sub_url = snap["url"]+"/subhalos/"+str(idsub)
        subhalo = get(sub_url)
        # If you want the subhalo cutout
        #vdisp = subhalo["veldisp"]
        #print idsub, vdisp
        #cutout_url = subhalo["cutouts"]["subhalo"]
        # If you want the parent halo cutout
        cutout_url = subhalo["cutouts"]["parent_halo"]

        # Subhalo center
        rc = np.array([subhalo["pos_x"], subhalo["pos_y"], subhalo["pos_z"]])

        cutout_exists = os.path.isfile('cutout_'+str(subhalo["grnr"])+'.hdf5')
        if cutout_exists:
            cutout_dic[subhalo["grnr"]] = 'cutout_'+str(subhalo["grnr"])+'.hdf5'
        
        if subhalo["grnr"] not in cutout_dic.keys():
	    cutout_request = {'gas':'Coordinates,Masses','stars':'Coordinates,Masses','dm':'Coordinates','bhs':'Coordinates,Masses'}        
            cutout_hdf5 = get(cutout_url, cutout_request)
            cutout_dic[subhalo["grnr"]] = cutout_hdf5
        else :
            cutout_hdf5 = cutout_dic[subhalo["grnr"]]

        print cutout_hdf5
        
        # Compute surface density profile for this subhalo
        ptype_dic = {"gas" : "0", "stars" : "4", "dm" : "1", "bhs" : "5"}
        nbins = 40 # 10 bins per dex
        rmin = 1.0e1 #kpc/h
        rmax = 1.0e5 #kpc/h

        # Project along z axis
        proj_axis = [0.0,0.0,1.0]
        projected_profz = make_projected_prof(cutout_hdf5,ptype_dic,rc,proj_axis,rmin,rmax,nbins,mdmpart)
        if flag == 0:
            prof_tot = np.copy(projected_profz)
            count = np.copy(projected_profz)
            izero = count <= 0
            ione = count > 0
            count[izero] = 1.0e-20
            count[ione] = 1
            del ione, izero
            flag = 1
        else :
            prof_tot = prof_tot + projected_profz
            ione = projected_profz > 0
            count[ione] = count[ione] + 1
            del ione
            
        del cutout_url
        del cutout_hdf5

    prof_tot = prof_tot/count

    params={'backend': 'ps','text.fontsize': 28,'axes.labelsize': 28,'text.fontsize': 28,'legend.fontsize': 18,
            'xtick.labelsize': 28,'ytick.labelsize': 28,'lines.linewidth': 4,'axes.linewidth': 3,
            'xtick.major.width': 3,'ytick.major.width': 3,'xtick.minor.width': 3,'ytick.minor.width': 3,
            'xtick.major.size': 10,'ytick.major.size': 10,'xtick.minor.size': 5,'ytick.minor.size': 5,
            'lines.markeredgewidth' : 3, 'lines.markersize': 6,
            'text.usetex': True,'text.dvipnghack': True,'figure.figsize': [11,10]}
    mpl.rcParams.update(params)
    
    plt.clf()
    plt.xscale('log')
    plt.yscale('linear')
    plt.xlabel(r'${\rm R \, [Mpc/h]}$')
    plt.ylabel(r'${\rm R \Delta \Sigma \, [Mpc \,M_{\odot}/pc^2]}$')
    plt.xlim([0.02,20])
    plt.ylim([0,10])
    plt.plot(10**prof_tot[:,0]/1e3, 10**prof_tot[:,0]/1e3*prof_tot[:,9],'k-',label=sim_name+' - z='+str(sim_z)+' - DM')
    plt.plot(10**prof_tot[:,0]/1e3, 10**prof_tot[:,0]/1e3*prof_tot[:,10],'g-',label=sim_name+' - z='+str(sim_z)+' - Gas')
    plt.plot(10**prof_tot[:,0]/1e3, 10**prof_tot[:,0]/1e3*prof_tot[:,11],'r-',label=sim_name+' - z='+str(sim_z)+' - Stars')
    plt.plot(10**prof_tot[:,0]/1e3, 10**prof_tot[:,0]/1e3*prof_tot[:,12],'b-',label=sim_name+' - z='+str(sim_z)+' - BHs')
    plt.plot(10**prof_tot[:,0]/1e3, 10**prof_tot[:,0]/1e3*(prof_tot[:,9]+prof_tot[:,10]+prof_tot[:,11]+prof_tot[:,12]),'k--',label=sim_name+' - z='+str(sim_z)+' - Total')
    plt.legend(loc=1,frameon=False)
    plt.savefig('DeltaSigma.png')
    plt.show()
