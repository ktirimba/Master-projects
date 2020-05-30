import matplotlib.pyplot as plt
from astropy.io.votable import parse_single_table
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy import units as u

from astroquery.gaia import Gaia
import numpy as np
import pandas as pd


maser_list= pd.read_csv('masers.dat',sep = ' ', header =0)
for i in range(len(maser_list)):
    #if maser_list['Refs'][i] =='42,66,65':
    if maser_list['arm'][i]=='Per':
        #print (maser_list['name'][i])
        coords = str(maser_list['ra'][i]) + str(maser_list['dec'][i])  ##getting maser coordinates in a string
        print (coords)
        physical_size = 125.0  ##physical size set
        parallax = maser_list['parallax'][i]

        parallax_error = maser_list['parallax_error'][i]
        pmra = maser_list['pmra'][i]
        pmra_err = maser_list['pmra_error'][i]
        # print ('error par',parallax)
        radial_distance = (1/parallax) * 1000
        print ('radial distance:',radial_distance)
        rad_dist_max = radial_distance/1000 + physical_size/1000 + (parallax_error/parallax**2)
        rad_dist_min = radial_distance/1000 - physical_size/1000 - (parallax_error/parallax**2)
        # print (radial_distance)
        print (rad_dist_max, 'kpc')
        print (rad_dist_min, 'kpc')
        # print  (maser_list['parallax_error'][i])
        radius = (206265.0* (physical_size) /radial_distance) *(1/3600.)
        #print ('radius',radius,'degrees')
        # print ((physical_size* u.parsec).to(u.parallax()))
        #print (parallax)
        min_par = (1/rad_dist_max)# + 3*maser_list['parallax_error'][i]# + maser_list['parallax'][i] / maser_list['parallax_error'][i] **2
        max_par = (1/rad_dist_min)#3*maser_list['parallax_error'][i] #parallax - maser_list['parallax'][i] / maser_list['parallax_error'][i] **2

        print (min_par,max_par)
        # print (max_par)
        # print ('radius', radius)
        c_deg= SkyCoord(coords, unit=(u.hourangle, u.degree)).to_string().split()## conversion of maser coordinates to degree
        c_deg = [float(i) for i in c_deg]
'''
        #print (c_deg)
        # k="SELECT * FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',{},{},{}))=1 AND (gaiadr2.gaia_source.parallax>=({}-gaiadr2.gaia_source.parallax_error)) " \
        #   "AND (gaiadr2.gaia_source.parallax<=({}+gaiadr2.gaia_source.parallax_error)) AND (gaiadr2.gaia_source.parallax_over_error >= 5)".format(c_deg[0],c_deg[1],radius,min_par,max_par)
        k="SELECT * FROM gaiadr2.gaia_source AS gaia INNER JOIN gaiadr1.tmass_best_neighbour AS xmatch ON gaia.source_id = xmatch.source_id INNER JOIN gaiadr1.tmass_original_valid AS tmass ON tmass.tmass_oid = xmatch.tmass_oid" \
          " WHERE CONTAINS(POINT('ICRS',gaia.ra,gaia.dec),CIRCLE('ICRS',{},{},{}))=1 AND (gaia.parallax>=({}-gaia.parallax_error)) AND (gaia.parallax<=({}+gaia.parallax_error)) AND (gaia.parallax_over_error >= 5)".format(c_deg[0],c_deg[1],radius,min_par,max_par)
        # #
        print (k)
        # with open("gaia_queries.txt", "a") as myfile:
        #     myfile.write(k + '\n')
        #     myfile.close()
        try:
            job = Gaia.launch_job_async(k,name=maser_list['name'][i])


            # r=job.get_results()
            #print (r['parallax'])
            # print (r['parallax_error'])
            # propx = np.array(r['pmra'])
            # print ('number of objects',len(r['pmra']))
            # propy = np.array(r['pmdec'])
            # abs_mag =  np.array(r['phot_g_mean_mag']) +5 - 5* np.log10(1/(np.array(r['parallax'])/1e3))

            #print('parallax over error:', np.array(r['parallax_over_error']))
            # plt.plot(r['bp_rp'],abs_mag,'.',c='r')
            # plt.xlabel('G$_{BP}$ - G$_{RP}$')
            # plt.ylabel('M$_G$')
            # plt.gca().invert_yaxis()
            # plt.show()
            # # print ('pmra list:',propx)
            # # print ('pmdec list:',propy)
            # plt.plot(propx,propy,'.')
            #
            # plt.plot(maser_list['prop_motion_x'][i],maser_list['prop_motion_y'][i],'.',c='r')
            #
            # plt.show()

        except ValueError:
            pass




    # j = Gaia.cone_search_async(c_deg,radius)
    # r = j.get_results()
    # vals = np.array([r['ra'],r['dec'],r['parallax'],r['parallax_error'],r['pmra'],r['pmra_error'],r['pmdec'],r['pmdec_error']])
    #print (np.shape(vals))'''
