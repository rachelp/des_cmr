import string
from pylab import *
import slr_zeropoint_shiftmap
import sys, getopt
import ezgal
import os
import io
from scipy import stats
import pdb
from PIL import Image
import numpy as np
from astropy import wcs
import pyfits as pyfits
import astropy.coordinates as coord
import astropy.units as u
import astropy.coordinates.angles as angles
import pyspherematch as psm
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
from matplotlib.patches import RegularPolygon
import matplotlib.pyplot as plt
from mpl_multiaxes import MultiPane_Subplot 
import time

class des_image():
      def __init__(self, catname,fitsname,jpegname,tracername,membername, boxname, bcgname, bcgname2, notesname, nedspec, photspec,choose_gal, readbcgs, xcs_z):

        im = imread(jpegname)
        scaled_image_data = im/255.
        self.im = scaled_image_data
        imfit = pyfits.open(fitsname)
	w = wcs.WCS(imfit[0].header)
        self.wcs = w
#        dt = [('ra',float), ('dec',float), ('mag_auto_g',float), ('magerr_auto_g',float), ('flags_g',float), ('class_star_g',float), ('spread_model_g',float), \
#                 ('mag_auto_r',float), ('magerr_auto_r',float), ('flags_r',float), ('class_star_r',float), ('spread_model_r',float),\
#                 ('mag_auto_i',float), ('magerr_auto_i',float), ('flags_i',float), ('class_star_i',float), ('spread_model_i',float),\
#                 ('mag_auto_z',float), ('magerr_auto_z',float), ('flags_z',float), ('class_star_z',float), ('spread_model_z',float)]
        dt = [('ra',float), ('dec',float), ('mag_auto_g',float), ('magerr_auto_g',float), ('mag_auto_r',float), ('magerr_auto_r',float), ('mag_auto_i',float), ('magerr_auto_i',float), \
                 ('mag_auto_z',float), ('magerr_auto_z',float), \
                 ('mag_detmodel_g',float), ('magerr_detmodel_g',float), ('mag_detmodel_r',float), ('magerr_detmodel_r',float), ('mag_detmodel_i',float), ('magerr_detmodel_i',float), \
                 ('mag_detmodel_z',float), ('magerr_detmodel_z',float)]

        cnr1 = w.wcs_pix2world(imfit[0].header['NAXIS1'],imfit[0].header['NAXIS2'],1)
        cnr2 = w.wcs_pix2world(imfit[0].header['NAXIS1'],0,1)
        cnr3 = w.wcs_pix2world(0,imfit[0].header['NAXIS2'],1)
        cnr4 =  w.wcs_pix2world(0,0,1)
        cnrra = [cnr1[0],cnr2[0],cnr3[0],cnr4[0]]
        cnrdec = [cnr1[1],cnr2[1],cnr3[1],cnr4[1]]
        maxra = max(cnrra)
        minra = min(cnrra)
        maxdec = max(cnrdec)
        mindec = min(cnrdec)
       
        # cut specz ad photoz galaxies
        try:
           if (nedspec != 0):
              wkeep = where((nedspec['ra'] >= minra) & (nedspec['ra'] <= maxra) & (nedspec['dec'] >= mindec) & (nedspec['dec'] <= maxdec))
              self.zra = nedspec[wkeep]['ra']
              self.zdec = nedspec[wkeep]['dec']
              self.zz = nedspec[wkeep]['z']
        except:
           self.zra = np.array([])
           self.zdec = np.array([])
           self.zz = np.array([])

        try:
           if (photspec != 0):
              wkeep = where((photspec.field('ra') >= minra) & (photspec.field('ra') <= maxra) & (photspec.field('dec') >= mindec) & (photspec.field('dec') <= maxdec))
              self.zra = np.append(self.zra, photspec[wkeep].field('ra'))
              self.zdec = np.append(self.zdec, photspec[wkeep].field('dec'))
              self.zz = np.append(self.zz, photspec[wkeep].field('z'))
        except:
           self.zra = np.array([])
           self.zdec = np.array([])
           self.zz = np.array([])

        # read in all galaxies
        galaxy = np.loadtxt(catname, dtype=dt, skiprows=(1), delimiter=',')
        inside = ((galaxy['ra'] < maxra) & (galaxy['ra'] > minra) & (galaxy['dec'] < maxdec) & (galaxy['dec'] > mindec))
        galaxy = galaxy[inside]
        self.gr_mag_all = galaxy['mag_auto_r']
        self.ri_mag_all = galaxy['mag_auto_r']
        self.iz_mag_all = galaxy['mag_auto_i']
        self.rz_mag_all = galaxy['mag_auto_i']
        self.gi_mag_all = galaxy['mag_auto_i']
        self.gr_color_all = galaxy['mag_detmodel_g'] - galaxy['mag_detmodel_r']
        self.gi_color_all = galaxy['mag_detmodel_g'] - galaxy['mag_detmodel_i']
        self.ri_color_all = galaxy['mag_detmodel_r'] - galaxy['mag_detmodel_i']
        self.iz_color_all = galaxy['mag_detmodel_i'] - galaxy['mag_detmodel_z']
        self.rz_color_all = galaxy['mag_detmodel_r'] - galaxy['mag_detmodel_z']
        self.ra_all =  galaxy['ra']
        self.dec_all =  galaxy['dec']

       # open tracer files. 
        dt = [('ra',float), ('dec',float), ('g_mag',float)]
        self.ra = np.array([])
        self.dec = np.array([])
        self.gr_color = np.array([])
        self.ri_color = np.array([])
        self.iz_color = np.array([])
        self.rz_color = np.array([])
        self.gi_color = np.array([])
        self.g_mag = np.array([])
        self.r_mag = np.array([])
        self.i_mag = np.array([])
        self.z_mag = np.array([])
	try:
 	  with open(membername):
                cmrgalaxy = np.loadtxt(membername, dtype=dt, skiprows=(0), delimiter=' ')
                m1 = np.zeros_like(cmrgalaxy['ra'], dtype=bool)
                m2 = np.zeros_like(cmrgalaxy['dec'], dtype=bool)
                m1[np.unique(cmrgalaxy['ra'], return_index=True)[1]] = True
                m2[np.unique(cmrgalaxy['dec'], return_index=True)[1]] = True
                keep = where((m1 == True) & (m2 == True))
                if size(keep) > 0:
                 for kk in range (0,size(keep)):
                   try:
                     result = psm.spherematch(cmrgalaxy['ra'][keep][kk],cmrgalaxy['dec'][keep][kk], self.ra_all, self.dec_all, nnearest=1)
                     if (result[0] < 0.1/3600.):
                       self.ra = np.append(self.ra,self.ra_all[result[1].astype(int)])
                       self.dec = np.append(self.dec,self.dec_all[result[1].astype(int)])
	               self.gr_color = np.append(self.gr_color,self.gr_color_all[result[1].astype(int)])
       		       self.gi_color = np.append(self.gi_color,self.gi_color_all[result[1].astype(int)])
        	       self.ri_color = np.append(self.ri_color,self.ri_color_all[result[1].astype(int)])
        	       self.iz_color = np.append(self.iz_color,self.iz_color_all[result[1].astype(int)])
        	       self.rz_color = np.append(self.rz_color,self.rz_color_all[result[1].astype(int)])
        	       self.g_mag = np.append(self.g_mag,galaxy['mag_auto_g'][result[1].astype(int)])
        	       self.r_mag = np.append(self.r_mag,galaxy['mag_auto_r'][result[1].astype(int)])
        	       self.i_mag = np.append(self.i_mag,galaxy['mag_auto_i'][result[1].astype(int)])
        	       self.z_mag = np.append(self.z_mag,galaxy['mag_auto_z'][result[1].astype(int)])
                   except:
                     pass
	except:
           pass
        dt = [('ra',float), ('dec',float), ('g_mag',float), ('r_mag',float),('i_mag',float),('z_mag',float),('gr_color',float), ('gi_color',float),('ri_color',float),('iz_color',float),('rz_color',float)]
        try:
          with open(tracername):
                cmrgalaxy = np.loadtxt(tracername, dtype=dt, skiprows=(0), delimiter=' ')
                m1 = np.zeros_like(cmrgalaxy['ra'], dtype=bool)
                m2 = np.zeros_like(cmrgalaxy['dec'], dtype=bool)
                m1[np.unique(cmrgalaxy['ra'], return_index=True)[1]] = True
                m2[np.unique(cmrgalaxy['dec'], return_index=True)[1]] = True
                keep = where((m1 == True) & (m2 == True))
                if size(keep) > 0:
                 for kk in range (0,size(keep)):
                   try:
                     result = psm.spherematch(cmrgalaxy['ra'][keep][kk],cmrgalaxy['dec'][keep][kk], self.ra_all, self.dec_all, nnearest=1)
                     if (result[0] < 0.1/3600.):
                       self.ra = np.append(self.ra,self.ra_all[result[1].astype(int)])
                       self.dec = np.append(self.dec,self.dec_all[result[1].astype(int)])
                       self.gr_color = np.append(self.gr_color,self.gr_color_all[result[1].astype(int)])
                       self.gi_color = np.append(self.gi_color,self.gi_color_all[result[1].astype(int)])
                       self.ri_color = np.append(self.ri_color,self.ri_color_all[result[1].astype(int)])
                       self.iz_color = np.append(self.iz_color,self.iz_color_all[result[1].astype(int)])
                       self.rz_color = np.append(self.rz_color,self.rz_color_all[result[1].astype(int)])
                       self.g_mag = np.append(self.g_mag,galaxy['mag_auto_g'][result[1].astype(int)])
                       self.r_mag = np.append(self.r_mag,galaxy['mag_auto_r'][result[1].astype(int)])
                       self.i_mag = np.append(self.i_mag,galaxy['mag_auto_i'][result[1].astype(int)])
                       self.z_mag = np.append(self.z_mag,galaxy['mag_auto_z'][result[1].astype(int)])
                   except:
                     pass
        except:
           pass
        if len(self.ra)>1:
           m1 = np.zeros_like(np.array(self.ra), dtype=bool)
           m2 = np.zeros_like(np.array(self.ra), dtype=bool)
           m1[np.unique(self.ra, return_index=True)[1]] = True
           m2[np.unique((np.array(self.ra))*(np.array(self.dec)), return_index=True)[1]] = True
           keep = where((m1 == True) & (m2 == True))
           self.ra=self.ra[keep]
           self.dec=self.dec[keep]
           self.gr_color=self.gr_color[keep]
           self.gi_color=self.gi_color[keep]
           self.ri_color=self.ri_color[keep]
           self.iz_color=self.iz_color[keep]
           self.rz_color=self.rz_color[keep]
           self.g_mag=self.g_mag[keep]
           self.r_mag=self.r_mag[keep]
           self.i_mag=self.i_mag[keep]
           self.z_mag=self.z_mag[keep]

        #initialize BCG list
        dt = [('ra',float), ('dec',float), ('g_mag',float), ('r_mag',float),('i_mag',float),('z_mag',float),('gr_color',float), ('gi_color',float),('ri_color',float),('iz_color',float),('rz_color',float)]
        self.bcg_xpix = np.array([])
        self.bcg_ypix = np.array([])
        self.bcg_ra  = np.array([])
        self.bcg_dec = np.array([])
        self.bcg_gr_color = np.array([])
        self.bcg_gi_color = np.array([])
        self.bcg_ri_color = np.array([])
        self.bcg_iz_color = np.array([])
        self.bcg_rz_color = np.array([])
        self.bcg_g_mag = np.array([])
        self.bcg_r_mag = np.array([])
        self.bcg_i_mag = np.array([])
        self.bcg_z_mag = np.array([])
        if readbcgs:
           try:
             with open(bcgname):
               bcggals=np.loadtxt(bcgname, dtype=dt, skiprows=(0), delimiter=' ')
               xpixtemp, ypixtemp= w.wcs_world2pix(bcggals['ra'], bcggals['dec'],1)
               self.bcg_xpix=np.append(self.bcg_xpix, xpixtemp)
               self.bcg_ypix=np.append(self.bcg_ypix, ypixtemp)
               self.bcg_ra=np.append(self.bcg_ra, bcggals['ra'])
               self.bcg_dec=np.append(self.bcg_dec, bcggals['dec'])
               self.bcg_gr_color=np.append(self.bcg_gr_color, bcggals['gr_color'])
               self.bcg_gi_color=np.append(self.bcg_gi_color, bcggals['gi_color'])
               self.bcg_ri_color=np.append(self.bcg_ri_color, bcggals['ri_color'])
               self.bcg_iz_color=np.append(self.bcg_iz_color, bcggals['iz_color'])
               self.bcg_rz_color=np.append(self.bcg_rz_color, bcggals['rz_color'])
               self.bcg_g_mag=np.append(self.bcg_g_mag, bcggals['g_mag'])
               self.bcg_r_mag=np.append(self.bcg_r_mag, bcggals['r_mag'])
               self.bcg_i_mag=np.append(self.bcg_i_mag, bcggals['i_mag'])
               self.bcg_z_mag=np.append(self.bcg_z_mag, bcggals['z_mag'])
           except:
             pass
           try:
             with open(bcgname2):
               bcggals=np.loadtxt(bcgname2, dtype=dt, skiprows=(0), delimiter=' ')
               xpixtemp, ypixtemp= w.wcs_world2pix(bcggals['ra'], bcggals['dec'],1)
               self.bcg_xpix=np.append(self.bcg_xpix, xpixtemp)
               self.bcg_ypix=np.append(self.bcg_ypix, ypixtemp)
               self.bcg_ra=np.append(self.bcg_ra, bcggals['ra'])
               self.bcg_dec=np.append(self.bcg_dec, bcggals['dec'])
               self.bcg_gr_color=np.append(self.bcg_gr_color, bcggals['gr_color'])
               self.bcg_gi_color=np.append(self.bcg_gi_color, bcggals['gi_color'])
               self.bcg_ri_color=np.append(self.bcg_ri_color, bcggals['ri_color'])
               self.bcg_iz_color=np.append(self.bcg_iz_color, bcggals['iz_color'])
               self.bcg_rz_color=np.append(self.bcg_rz_color, bcggals['rz_color'])
               self.bcg_g_mag=np.append(self.bcg_g_mag, bcggals['g_mag'])
               self.bcg_r_mag=np.append(self.bcg_r_mag, bcggals['r_mag'])
               self.bcg_i_mag=np.append(self.bcg_i_mag, bcggals['i_mag'])
               self.bcg_z_mag=np.append(self.bcg_z_mag, bcggals['z_mag'])
           except:
             pass
        if len(self.bcg_ra)>1:
           m1 = np.zeros_like(np.array(self.bcg_ra), dtype=bool)
           m2 = np.zeros_like(np.array(self.bcg_ra), dtype=bool)
           m1[np.unique(self.bcg_ra, return_index=True)[1]] = True
           m2[np.unique((np.array(self.bcg_ra))*(np.array(self.bcg_dec)), return_index=True)[1]] = True
           keep = where((m1 == True) & (m2 == True))
           self.bcg_xpix=self.bcg_xpix[keep]
           self.bcg_ypix=self.bcg_ypix[keep]
           self.bcg_ra=self.bcg_ra[keep]
           self.bcg_dec=self.bcg_dec[keep]
           self.bcg_gr_color=self.bcg_gr_color[keep]
           self.bcg_gi_color=self.bcg_gi_color[keep]
           self.bcg_ri_color=self.bcg_ri_color[keep]
           self.bcg_iz_color=self.bcg_iz_color[keep]
           self.bcg_rz_color=self.bcg_rz_color[keep]
           self.bcg_g_mag=self.bcg_g_mag[keep]
           self.bcg_r_mag=self.bcg_r_mag[keep]
           self.bcg_i_mag=self.bcg_i_mag[keep]
           self.bcg_z_mag=self.bcg_z_mag[keep]
           

        #open box files
        dt = [('xbot',float), ('xtop',float), ('ybot',float), ('ytop',float)]
        self.gr_xbot = np.array([])
        self.gr_xtop = np.array([])
        self.ri_xbot = np.array([])
        self.ri_xtop = np.array([])
        self.iz_xbot = np.array([])
        self.iz_xtop = np.array([])
        self.gi_xbot = np.array([])
        self.gi_xtop = np.array([])
        self.rz_xbot = np.array([])
        self.rz_xtop = np.array([])
        self.gr_ybot = np.array([])
        self.gr_ytop = np.array([])
        self.ri_ybot = np.array([])
        self.ri_ytop = np.array([])
        self.iz_ybot = np.array([])
        self.iz_ytop = np.array([])
        self.gi_ybot = np.array([])
        self.gi_ytop = np.array([])
        self.rz_ybot = np.array([])
        self.rz_ytop = np.array([])
        try:
          with open(boxname):
                 boxit = np.loadtxt(boxname, dtype=dt, skiprows=(0))
                 self.gr_xbot = boxit[0]['xbot']
                 self.gr_xtop = boxit[0]['xtop']
                 self.gr_ybot = boxit[0]['ybot']
                 self.gr_ytop = boxit[0]['ytop']
                 self.ri_xbot = boxit[1]['xbot']
                 self.ri_xtop = boxit[1]['xtop']
                 self.ri_ybot = boxit[1]['ybot']
                 self.ri_ytop = boxit[1]['ytop']
                 self.iz_xbot = boxit[2]['xbot']
                 self.iz_xtop = boxit[2]['xtop']
                 self.iz_ybot = boxit[2]['ybot']
                 self.iz_ytop = boxit[2]['ytop']
                 self.gi_xbot = boxit[3]['xbot']
                 self.gi_xtop = boxit[3]['xtop']
                 self.gi_ybot = boxit[3]['ybot']
                 self.gi_ytop = boxit[3]['ytop']
                 self.rz_xbot = boxit[4]['xbot']
                 self.rz_xtop = boxit[4]['xtop']
                 self.rz_ybot = boxit[4]['ybot']
                 self.rz_ytop = boxit[4]['ytop']
        except IOError:
                 self.gr_xbot = 0
                 self.gr_xtop = 0
                 self.gr_ybot = 0
                 self.gr_ytop = 0
                 self.ri_xbot = 0
                 self.ri_xtop = 0
                 self.ri_ybot = 0
                 self.ri_ytop = 0
                 self.iz_xbot = 0
                 self.iz_xtop = 0
                 self.iz_ybot = 0
                 self.iz_ytop = 0
                 self.gi_xbot = 0
                 self.gi_xtop = 0
                 self.gi_ybot = 0
                 self.gi_ytop = 0
                 self.rz_xbot = 0
                 self.rz_xtop = 0
                 self.rz_ybot = 0
                 self.rz_ytop = 0
        
        # trim member, specz and photoz galaxies
        if (len(self.ra) > 1):
            a = len(self.ra)
            xpix = np.zeros([len(self.ra)])
            ypix = np.zeros([len(self.ra)])
            for i in range(0,a):
                xpix[i],ypix[i] = w.wcs_world2pix(self.ra[i],self.dec[i],1)
        else:
            xpix = np.array([0])
            ypix = np.array([0])
        m1 = np.zeros_like(xpix, dtype=bool)
        m2 = np.zeros_like(ypix, dtype=bool)
        m1[np.unique(xpix, return_index=True)[1]] = True
        m2[np.unique(ypix, return_index=True)[1]] = True
        keep = where((m1 == True) & (m2 == True) & (xpix != 0))
        self.xpix = xpix[keep]
        self.ypix = ypix[keep]
        self.ra = self.ra[keep]
        self.dec = self.dec[keep]
        self.gr_color = self.gr_color[keep]
        self.gi_color = self.gi_color[keep]
        self.ri_color = self.ri_color[keep]
        self.iz_color = self.iz_color[keep]
        self.rz_color = self.rz_color[keep]
        self.g_mag = self.g_mag[keep]
        self.r_mag = self.r_mag[keep]
        self.i_mag = self.i_mag[keep]
        self.z_mag = self.z_mag[keep]


        self.z_gr_color = np.array([])
        self.z_gi_color = np.array([])
        self.z_ri_color = np.array([])
        self.z_iz_color = np.array([])
        self.z_rz_color = np.array([])
        self.z_gr_mag = np.array([])
        self.z_gi_mag = np.array([])
        self.z_ri_mag = np.array([])
        self.z_iz_mag = np.array([])
        self.z_rz_mag = np.array([])
        self.z_member = np.array([])
        self.z_member_z = np.array([])
        if (len(self.zra) > 1 and len(self.ra_all) > 1):
            a = len(self.zra)
            for i in range(0,a):
#                xpix,ypix = w.wcs_world2pix(self.zra[i],self.zdec[i],1)
#                world = w.wcs_pix2world(xpix,667-ypix,1)
#                result = psm.spherematch(world[0],world[1], self.ra_all, self.dec_all, nnearest=1)
                result = psm.spherematch(self.zra[i],self.zdec[i], self.ra_all, self.dec_all, nnearest=1)
                if (result[0] < 1./3600.):
                  self.z_gr_color = append(self.z_gr_color,self.gr_color_all[result[1].astype(int)])
                  self.z_gi_color = append(self.z_gi_color,self.gi_color_all[result[1].astype(int)])
                  self.z_ri_color = append(self.z_ri_color,self.ri_color_all[result[1].astype(int)])
                  self.z_iz_color = append(self.z_iz_color,self.iz_color_all[result[1].astype(int)])
                  self.z_rz_color = append(self.z_rz_color,self.rz_color_all[result[1].astype(int)])
                  self.z_gr_mag = append(self.z_gr_mag,self.gr_mag_all[result[1].astype(int)])
                  self.z_gi_mag = append(self.z_gi_mag,self.gi_mag_all[result[1].astype(int)])
                  self.z_ri_mag = append(self.z_ri_mag,self.ri_mag_all[result[1].astype(int)])
                  self.z_iz_mag = append(self.z_iz_mag,self.iz_mag_all[result[1].astype(int)])
                  self.z_rz_mag = append(self.z_rz_mag,self.rz_mag_all[result[1].astype(int)])
                  self.z_member = append(self.z_member,result[1].astype(int))
                  self.z_member_z = append(self.z_member_z, self.zz[i])


        #start making figures        
        ion()
        fig = plt.figure(num = 1, figsize=(12,12))
        fig.clf()
        ax = MultiPane_Subplot(fig, subplot_pos=(1, 1, 1), nrows_ncols = (2, 2), n_pane=4, pane_direction="row", axes_pad_inch=0.0)
        fig.add_subplot(ax)
        ax[0].imshow(scaled_image_data)
        ax[1].imshow(scaled_image_data[:,:,0], cmap=plt.cm.Reds_r)
        ax[2].imshow(scaled_image_data[:,:,1], cmap=plt.cm.Greens_r)
        ax[3].imshow(scaled_image_data[:,:,2], cmap=plt.cm.Blues_r)
        fig2 = plt.figure(num = 2, figsize=(18,9))
        fig2.clf()
        bx = MultiPane_Subplot(fig2, subplot_pos=(1, 1, 1), nrows_ncols = (2, 3), n_pane=6, pane_direction="row", axes_pad_inch=0.0)
        fig2.add_subplot(bx)
        self.temp_xpix=-100
        self.temp_ypix=-100
        self.temp_gr_color=-100
        self.temp_gi_color=-100
        self.temp_ri_color=-100
        self.temp_iz_color=-100
        self.temp_rz_color=-100
        self.temp_r_mag=-100
        self.temp_i_mag=-100
        def display_ax0():
            ax[0].cla()
            ax[0].imshow(scaled_image_data)
            for i in range(0,len(self.xpix)):
               ellipse = Ellipse(xy=(self.xpix[i], 667-self.ypix[i]), width=30, height=30, edgecolor='g', fc='None', lw=2)
               ax[0].add_patch(ellipse)
            if (len(self.zra) > 0):
               a = len(self.zra)
               for i in range(a):
                  xpix,ypix = w.wcs_world2pix(self.zra[i],self.zdec[i],1)
                  left = xpix - 10
                  bottom = 667 - (ypix + 10)
                  rectangle = Rectangle((left,bottom), 20,20,edgecolor='c', fc='None', lw=1)
                  ax[0].add_patch(rectangle)
            if (len(self.bcg_ra)>0):
               a = len(self.bcg_ra)
               for i in range(a):
                  xpix,ypix = w.wcs_world2pix(self.bcg_ra[i],self.bcg_dec[i],1)
                  left = xpix - 20
                  bottom = 667 - (ypix + 20)
                  rectangle = Rectangle((left,bottom), 40,40,edgecolor='r', fc='None', lw=1)
                  ax[0].add_patch(rectangle)
            if self.temp_gr_color !=-100:
               left = self.temp_xpix - 15
               bottom = 667 - (self.temp_ypix + 15)
               rectangle = Rectangle((left,bottom), 30, 30,edgecolor='g', fc='none', lw=1)
               ax[0].add_patch(rectangle)
            fig.suptitle('Measured Redshift = '+'{: .3f}'.format(xcs_z), fontsize=20)
            fig.canvas.draw()
        def display_bx():
            bx[0].cla()
            bx[1].cla()
            bx[3].cla()
            bx[4].cla()
            bx[5].cla()
            bx[0].text(20, 0.0, r'g-r', fontsize=15)
            bx[1].text(20, 0.0, r'g-i', fontsize=15)
            bx[3].text(20, 2, r'r-i', fontsize=15)
            bx[4].text(20, 2, r'i-z', fontsize=15)
            bx[5].text(20, 2, r'r-z', fontsize=15)
            bx[0].scatter(self.gr_mag_all, self.gr_color_all,color="black", marker='.',s=1)
            bx[1].scatter(self.gi_mag_all, self.gi_color_all,color="black", marker='.',s=1)
            bx[3].scatter(self.ri_mag_all, self.ri_color_all,color="black", marker='.',s=1)
       	    bx[4].scatter(self.iz_mag_all, self.iz_color_all,color="black", marker='.',s=1)
            bx[5].scatter(self.rz_mag_all, self.rz_color_all,color="black", marker='.',s=1)
       	    if len(self.ra) > 1:
               bx[0].scatter(self.r_mag,self.gr_color,marker='o',s=40,color="blue")
               bx[1].scatter(self.i_mag,self.gi_color,marker='o',s=40,color="green")
               bx[3].scatter(self.r_mag,self.ri_color,marker='o',s=40,color="orange")
               bx[4].scatter(self.i_mag,self.iz_color,marker='o',s=40,color="red")
               bx[5].scatter(self.i_mag,self.rz_color,marker='o',s=40,color="magenta")    
            if len(self.zra) > 1:
               bx[0].scatter(self.z_gr_mag,self.z_gr_color,color="cyan",marker='s',s=5)
               bx[1].scatter(self.z_gi_mag,self.z_gi_color,color="cyan",marker='s',s=5)
               bx[3].scatter(self.z_ri_mag,self.z_ri_color,color="cyan",marker='s',s=5)
               bx[4].scatter(self.z_iz_mag,self.z_iz_color,color="cyan",marker='s',s=5)
               bx[5].scatter(self.z_rz_mag,self.z_rz_color,color="cyan",marker='s',s=5)
            if len(self.bcg_ra>1):
               bx[0].scatter(self.bcg_r_mag,self.bcg_gr_color,color="red",marker=(5, 0),s=250, alpha=0.2)
               bx[1].scatter(self.bcg_i_mag,self.bcg_gi_color,color="red",marker=(5, 0),s=250, alpha=0.2)
               bx[3].scatter(self.bcg_r_mag,self.bcg_ri_color,color="red",marker=(5, 0),s=250, alpha=0.2)
               bx[4].scatter(self.bcg_i_mag,self.bcg_iz_color,color="red",marker=(5, 0),s=250, alpha=0.2)
               bx[5].scatter(self.bcg_i_mag,self.bcg_rz_color,color="red",marker=(5, 0),s=250, alpha=0.2)
            if self.temp_gr_color !=-100:
               bx[0].scatter(self.temp_r_mag, self.temp_gr_color, color="black", marker=(5, 2), s=250, alpha=0.4)
               bx[1].scatter(self.temp_i_mag, self.temp_gi_color, color="black", marker=(5, 2), s=250, alpha=0.4)
               bx[3].scatter(self.temp_r_mag, self.temp_ri_color, color="black", marker=(5, 2), s=250, alpha=0.4)
               bx[4].scatter(self.temp_i_mag, self.temp_iz_color, color="black", marker=(5, 2), s=250, alpha=0.4)
               bx[5].scatter(self.temp_i_mag, self.temp_rz_color, color="black", marker=(5, 2), s=250, alpha=0.4)
            if (self.gr_xbot != 0):
               rect = Rectangle( ( 0,0 ), 1, 1, alpha = 0.2, ec = "gray", fc = "CornflowerBlue", visible = True, axes=bx[0])
               rect.set_width(self.gr_xtop - self.gr_xbot)
               rect.set_height(self.gr_ytop - self.gr_ybot)
               rect.set_xy((self.gr_xbot, self.gr_ybot))
               bx[0].add_patch(rect)
            if (self.gi_xbot != 0):
               rect = Rectangle( ( 0,0 ), 1, 1, alpha = 0.2, ec = "gray", fc = "CornflowerBlue", visible = True, axes=bx[1])
               rect.set_width(self.gi_xtop - self.gi_xbot)
               rect.set_height(self.gi_ytop - self.gi_ybot)
               rect.set_xy((self.gi_xbot, self.gi_ybot))
               bx[1].add_patch(rect)
            if (self.ri_xbot != 0):
               rect = Rectangle( ( 0,0 ), 1, 1, alpha = 0.2, ec = "gray", fc = "CornflowerBlue", visible = True, axes=bx[3])
               rect.set_width(self.ri_xtop - self.ri_xbot)
               rect.set_height(self.ri_ytop - self.ri_ybot)
               rect.set_xy((self.ri_xbot, self.ri_ybot))
               bx[3].add_patch(rect)
            if (self.iz_xbot != 0):
               rect = Rectangle( ( 0,0 ), 1, 1, alpha = 0.2, ec = "gray", fc = "CornflowerBlue", visible = True, axes=bx[4])
               rect.set_width(self.iz_xtop - self.iz_xbot)
               rect.set_height(self.iz_ytop - self.iz_ybot)
               rect.set_xy((self.iz_xbot, self.iz_ybot))
               bx[4].add_patch(rect)
            if (self.rz_xbot != 0):
               rect = Rectangle( ( 0,0 ), 1, 1, alpha = 0.2, ec = "gray", fc = "CornflowerBlue", visible = True, axes=bx[5])
               rect.set_width(self.rz_xtop - self.rz_xbot)
               rect.set_height(self.rz_ytop - self.rz_ybot)
               rect.set_xy((self.rz_xbot, self.rz_ybot))
               bx[5].add_patch(rect)  
            bx.set_xlim(16, 24)
            bx.set_ylim(-0.2, 3.5)
            bx.axes_llc.set_xlabel("magnitude")
            bx.axes_llc.set_ylabel("color")
            fig2.suptitle('Measured Redshift = '+'{: .3f}'.format(xcs_z), fontsize=20)
            fig2.canvas.draw()  
        display_ax0()    
        display_bx()

        #start procedures to enable mouse actions
        def on_press(event):
            self.x0 = event.xdata
            self.y0 = event.ydata
            self.x0a = event.x
            self.y0a = event.y
        def on_release(event):
            x0=self.x0
            y0=self.y0
            x0a=self.x0a
            y0a=self.y0a
            x1 = event.xdata
            y1 = event.ydata
            x1a = event.x
            y1a = event.y
            self.x1=x1
            self.y1=y1
            self.x1a=x1a
            self.y1a=y1a
            if ( y1>y0 ):
               ytop = y1
               ybot = y0
            else:
               ybot = y1
               ytop = y0
            if ( x1>x0 ):
               xtop = x1
               xbot = x0
            else:
               xbot = x1
               xtop = x0            
            bbox = fig2.get_window_extent().transformed(fig2.dpi_scale_trans.inverted())
            fullwidth = bbox.width*fig2.dpi
            fullheight = bbox.height*fig2.dpi
            bbox = bx.get_window_extent()
            bxwidth = bbox.width
            bxheight = bbox.height
            padx = (fullwidth - bxwidth*3)/2.0
            pady = (fullheight - bxheight*2)/2.0
            panel0xmin = padx
            panel0xmax = padx + bxwidth
            panel0ymin = pady + bxheight
            panel0ymax = pady + 2*bxheight
            panel1xmin = padx + bxwidth
            panel1xmax = padx + 2*bxwidth
            panel1ymin = pady + bxheight
            panel1ymax = pady + 2*bxheight
            panel3xmin = padx 
            panel3xmax = padx + bxwidth
            panel3ymin = pady
            panel3ymax = pady + bxheight
            panel4xmin = padx + bxwidth
            panel4xmax = padx + 2*bxwidth
            panel4ymin = pady
            panel4ymax = pady + bxheight
            panel5xmin = padx + 2*bxwidth
            panel5xmax = padx + 3*bxwidth
            panel5ymin = pady
            panel5ymax = pady + bxheight
            if ((x0a < panel0xmax) & (x0a > panel0xmin) & (y0a > panel0ymin) & (y0a < panel0ymax)):
               self.gr_xbot = xbot
               self.gr_xtop = xtop
               self.gr_ybot = ybot
               self.gr_ytop = ytop
            elif ((x0a < panel3xmax) & (x0a > panel3xmin) & (y0a > panel3ymin) & (y0a < panel3ymax)):
               self.ri_xbot = xbot
               self.ri_xtop = xtop
               self.ri_ybot = ybot
               self.ri_ytop = ytop
            elif ((x0a < panel1xmax) & (x0a > panel1xmin) & (y0a > panel1ymin) & (y0a < panel1ymax)):
               self.gi_xbot = xbot
               self.gi_xtop = xtop
               self.gi_ybot = ybot
               self.gi_ytop = ytop
            elif ((x0a < panel4xmax) & (x0a > panel4xmin) & (y0a > panel4ymin) & (y0a < panel4ymax)):
               self.iz_xbot = xbot
               self.iz_xtop = xtop
               self.iz_ybot = ybot
               self.iz_ytop = ytop
            elif ((x0a < panel5xmax) & (x0a > panel5xmin) & (y0a > panel5ymin) & (y0a < panel5ymax)):
               self.rz_xbot = xbot
               self.rz_xtop = xtop
               self.rz_ybot = ybot
               self.rz_ytop = ytop
            display_bx()

        def onclick(event):
            if event.button !=3:
               coords = np.vstack((event.xdata,667-event.ydata)).T
               world = w.wcs_pix2world(coords,1)
               result = psm.spherematch(world[:,0],world[:,1], self.ra_all, self.dec_all, nnearest=1)
               if (result[0] < 1./3600.):
                  xpix,ypix = w.wcs_world2pix(self.ra_all[result[1].astype(int)], self.dec_all[result[1].astype(int)],1)
                  if (((xpix == self.xpix) & (ypix == self.ypix)).any() == False):
                     self.xpix = np.append(self.xpix,xpix)
                     self.ypix = np.append(self.ypix,ypix)
                     self.ra = np.append(self.ra, self.ra_all[result[1].astype(int)])
                     self.dec = np.append(self.dec, self.dec_all[result[1].astype(int)])
                     self.gr_color = np.append(self.gr_color, self.gr_color_all[result[1].astype(int)])
                     self.gi_color = np.append(self.gi_color, self.gi_color_all[result[1].astype(int)])
                     self.ri_color = np.append(self.ri_color, self.ri_color_all[result[1].astype(int)])
                     self.iz_color = np.append(self.iz_color, self.iz_color_all[result[1].astype(int)])
                     self.rz_color = np.append(self.rz_color, self.rz_color_all[result[1].astype(int)])
                     self.g_mag = np.append(self.g_mag, galaxy['mag_auto_g'][result[1].astype(int)])
                     self.r_mag = np.append(self.r_mag, galaxy['mag_auto_r'][result[1].astype(int)])
                     self.i_mag = np.append(self.i_mag, galaxy['mag_auto_i'][result[1].astype(int)])
                     self.z_mag = np.append(self.z_mag, galaxy['mag_auto_z'][result[1].astype(int)])
                  else:
                     keep = where(((xpix == self.xpix) & (ypix == self.ypix))==False)
                     self.xpix = self.xpix[keep]
                     self.ypix = self.ypix[keep]
                     self.ra = self.ra[keep]
                     self.dec = self.dec[keep]
              	     self.gr_color = self.gr_color[keep]
                     self.gi_color = self.gi_color[keep]
                     self.ri_color = self.ri_color[keep]
                     self.iz_color = self.iz_color[keep]
                     self.rz_color = self.rz_color[keep]
                     self.g_mag = self.g_mag[keep]
                     self.r_mag = self.r_mag[keep]
                     self.i_mag = self.i_mag[keep]
                     self.z_mag = self.z_mag[keep]
                  self.temp_xpix=xpix
                  self.temp_ypix=ypix
                  self.temp_gr_color=self.gr_color_all[result[1].astype(int)]
                  self.temp_gi_color=self.gi_color_all[result[1].astype(int)]
                  self.temp_ri_color=self.ri_color_all[result[1].astype(int)]
                  self.temp_iz_color=self.iz_color_all[result[1].astype(int)]
                  self.temp_rz_color=self.rz_color_all[result[1].astype(int)]
                  self.temp_r_mag=galaxy['mag_auto_r'][result[1].astype(int)]
                  self.temp_i_mag=galaxy['mag_auto_i'][result[1].astype(int)]
                  display_ax0()
                  display_bx()
                
        def onrightclick(event):
            if event.button==3:
               coords = np.vstack((event.xdata,667-event.ydata)).T
               world = w.wcs_pix2world(coords,1)
               result = psm.spherematch(world[:,0],world[:,1], self.ra_all, self.dec_all, nnearest=1)
               if (result[0] < 1./3600.):
                  xpix,ypix = w.wcs_world2pix(self.ra_all[result[1].astype(int)], self.dec_all[result[1].astype(int)],1)
                  if (((xpix == self.bcg_xpix) & (ypix == self.bcg_ypix)).any() == False):
                     self.bcg_xpix = np.append(self.bcg_xpix,xpix)
                     self.bcg_ypix = np.append(self.bcg_ypix,ypix)
                     self.bcg_ra =  np.append(self.bcg_ra, self.ra_all[result[1].astype(int)])
                     self.bcg_dec = np.append(self.bcg_dec, self.dec_all[result[1].astype(int)])
                     self.bcg_gr_color = np.append(self.bcg_gr_color, self.gr_color_all[result[1].astype(int)])
                     self.bcg_gi_color = np.append(self.bcg_gi_color, self.gi_color_all[result[1].astype(int)])
                     self.bcg_ri_color = np.append(self.bcg_ri_color, self.ri_color_all[result[1].astype(int)])
                     self.bcg_iz_color = np.append(self.bcg_iz_color, self.iz_color_all[result[1].astype(int)])
                     self.bcg_rz_color = np.append(self.bcg_rz_color, self.rz_color_all[result[1].astype(int)])
                     self.bcg_g_mag = np.append(self.bcg_g_mag, galaxy['mag_auto_g'][result[1].astype(int)])
                     self.bcg_r_mag = np.append(self.bcg_r_mag, galaxy['mag_auto_r'][result[1].astype(int)])
                     self.bcg_i_mag = np.append(self.bcg_i_mag, galaxy['mag_auto_i'][result[1].astype(int)])
                     self.bcg_z_mag = np.append(self.bcg_z_mag, galaxy['mag_auto_z'][result[1].astype(int)])
                  else:
                     keep = where(((xpix == self.bcg_xpix) & (ypix == self.bcg_ypix))==False)
                     self.bcg_xpix = self.bcg_xpix[keep]
                     self.bcg_ypix = self.bcg_ypix[keep]
                     self.bcg_ra = self.bcg_ra[keep]
               	     self.bcg_dec = self.bcg_dec[keep]
                     self.bcg_gr_color = self.bcg_gr_color[keep]
               	     self.bcg_gi_color = self.bcg_gi_color[keep]
               	     self.bcg_ri_color = self.bcg_ri_color[keep]
                     self.bcg_iz_color = self.bcg_iz_color[keep]
                     self.bcg_rz_color = self.bcg_rz_color[keep]
                     self.bcg_g_mag = self.bcg_g_mag[keep]
	             self.bcg_r_mag = self.bcg_r_mag[keep]
                     self.bcg_i_mag = self.bcg_i_mag[keep]
         	     self.bcg_z_mag = self.bcg_z_mag[keep]
               display_bx()
               display_ax0()

        #link mouse action procedures with figures
	cid = fig.canvas.mpl_connect('button_press_event', onclick)
        cid = fig.canvas.mpl_connect('button_press_event', onrightclick)
        cid = fig2.canvas.mpl_connect('button_press_event', on_press)
        cid = fig2.canvas.mpl_connect('button_release_event', on_release)
        
        # Print existing note if there is one.
        if os.path.isfile(notesname) and os.path.getsize(notesname) > 0:
        	notesexist = True
        	notefile = open(notesname,'a+')
        	notefile.seek(0)
        	print "\nExisting notes:\n" + notefile.read()
        else:
        	notesexist = False
        	
        #decide if figures will be shown or not
        if (choose_gal == True):
        	# Add a note
        	valid = 0
        	while (valid ==0):
        		newnote = raw_input('Add a note? (y/n) -> ')
        		if newnote == 'y' or newnote == 'Y':
        			if notesexist: #if note exists already
        				notetext = raw_input('Type your note below then press enter. Your text will be appended to the existing note: \n')
        				notefile.write(time.strftime("%x") + ': ' + notetext + '\n')
        				valid = 1
        			else:
        				notetext = raw_input('Type your note below then press enter. \n')
        				if len(notetext) > 0:
        					notefile = open(notesname, 'a+')
        					notefile.write(time.strftime("%x") + ': ' + notetext + '\n')
        				valid = 1
        		elif newnote == 'n' or newnote == 'N':
        			valid = 1
        		else:
        			print 'Invalid option.\n'
        			
        	# Flag the galaxy
        	valid = 0
        	while (valid == 0):
        		raiseflag = raw_input('Flag this cluster? (y/n) -> ')
        		if raiseflag == 'n' or raiseflag == 'N':
        			self.flagged = False
        			valid = 1
        		elif raiseflag == 'y' or raiseflag == 'Y':
        			self.flagged = True
        			valid = 1
        		else:
        			print 'Invalid option.\n'
        	raw_input('Press ENTER when you are ready to move on to the next cluster.')
        	show()

        		

           
        else:
        	plt.close(fig)
        	plt.close(fig2)
        	
        if os.path.isfile(notesname):
        	notefile.close()
        


# Main starts here
total = len(sys.argv)
cmdargs = str(sys.argv)

inflagtype = 'flagged'
outflagtype = 'skewedxray'


woz=True
choose = True
owbcgs=False
readbcgs=True
output_dir = './'
input_dir = './'
file_list1='des_xcs_good_jan28_coords_sql.lis'
file_list2='des_xcs_good_jan28_spectro_photoz_redshifts.lis'


try:
   opts, args = getopt.getopt(sys.argv[1:],"i:o:f:h:", ['nochoose', 'overwrite_bcgs', 'noreadbcgs', 'wz'])
except getopt.GetoptError:
   print 'test.py -i <inputfile> -o <outputfile>'
   sys.exit(2)

for opt, arg in opts:
   if opt == '-h':
     print 'test.py -i <inputfile> -o <outputfile> -f <file>'
     sys.exit()
   elif opt in ("-i", "--icatdir"):
     input_dir = arg
   elif opt in ("-o", "--odir"):
     output_dir = arg
   elif opt in ("-f", "--file"):
     file_list1 = arg
     file_list2 = arg
for arg in args:
   if (arg ==  'nochoose'):
     choose = False
   elif (arg =='overwrite_bcgs'):
     owbcgs=True
   elif (arg =='noreadbcgs'):
     readbcgs=False
   elif (arg =='wz'):
     woz=False

photzspecname ='des_specz_master.fits'
nedspecname ='ned_gals_spec.dat'
#root_xcs, root_tile, clust_ra, clust_dec = np.loadtxt('des_xcs_good_sep18_coords_sql_modeltest.lis', dtype=str).T
#root_xcs, root_tile, clust_ra, clust_dec = np.loadtxt('des_xcs_good_sep18_coords_sql.lis', dtype=str).T
if woz:
   root_xcs, root_tile, clust_ra, clust_dec  = np.loadtxt(file_list1, dtype=str).T
   xcs_zs=[99.0]*len(root_xcs)
else:
   dt = [('XCSNAME',np.object), ('DESTILE',np.object), ('RA',np.object), ('DEC',np.object), ('z_lit',float), ('ref',np.object), ('type',float), ('z_myspec', float),('vdisp_my_spec', float) ,('num_my_spec', float), ('photoz',float)]
   xcsclusters = np.loadtxt(file_list2, dtype=dt, skiprows=(1)).T
   root_xcs=xcsclusters['XCSNAME']
   root_tile=xcsclusters['DESTILE']
   clust_ra=xcsclusters['RA']
   clust_dec=xcsclusters['DEC']
   xcs_zs=xcsclusters['photoz']

#read in input and ouput file names

fitsname = np.array(str,dtype='S50')
jpegname = np.array(str,dtype='S50')
incatname = np.array(str,dtype='S50')
inresname = np.array(str,dtype='S50')
inzname = np.array(str,dtype='S50')
inspeczname = np.array(str,dtype='S50')
intracername = np.array(str,dtype='S50')
inmembername = np.array(str,dtype='S50')
inboxname = np.array(str,dtype='S50')
inbcgname = np.array(str,dtype='S50')
innotesname = np.array(str, dtype = 'S50')

outbcgname = np.array(str,dtype='S50')

flagname = np.array(str, dtype = 'S50')

a = len(root_tile)
for i in range(0,a):
      name = input_dir + '/jpegs/desdm_xcs_xcontours_' + root_xcs[i] + '.jpg'
      jpegname = np.append(jpegname,name)
      
      name = input_dir + '/fits/' + root_tile[i] + '_z_' + root_xcs[i] + '.fits'
      fitsname = np.append(fitsname,name)
      
      name = input_dir + '/cats/' + root_tile[i] + '_' + root_xcs[i] + '.cat'
#      name = input_dir + '/cats_rev_stargal/' + root_tile[i] + '_' + root_xcs[i] + '.cat'
#      name = input_dir + '/cats_normal/' + root_tile[i] + '_' + root_xcs[i] + '.cat'
#      name = input_dir + '/cats_diskonly/' + root_tile[i] + '_' + root_xcs[i] + '_diskonly.cat'
#      name = input_dir + '/cats_diskbulge/' + root_tile[i] + '_' + root_xcs[i] + '_diskbulge.cat'
#      name = input_dir + '/cats_deblended/' + root_tile[i] + '_' + root_xcs[i] + '_undeblended.cat'
      incatname = np.append(incatname,name)
      
      name = input_dir + '/cmr/' + root_tile[i] + '_' + root_xcs[i] + '.res'
      inresname = np.append(inresname,name)
      
      name = input_dir + '/cmr/' + root_tile[i] + '_' + root_xcs[i] + '.cmr'
      intracername = np.append(intracername,name)
      
      name = input_dir + '/cmr/' + root_tile[i] + '_' + root_xcs[i] + '.bb'
      inboxname = np.append(inboxname,name)
      
      name = input_dir + '/cmr/' + root_tile[i] + '_' + root_xcs[i] + '.photoz'
      inzname = np.append(inzname,name)
      
      name = input_dir + '/cmr/' + root_tile[i] + '_' + root_xcs[i] + '.specz'
      inspeczname = np.append(inspeczname,name)
      
      name = '/malbec1/des-sv/sva1/XCS/cmr/' + root_tile[i] + '_' + root_xcs[i] + '.member'
      inmembername = np.append(inmembername,name)
      
      name = input_dir + '/bcg/' + root_tile[i] + '_' + root_xcs[i] + '.bcg'
      inbcgname = np.append(inbcgname,name)
     
      name = output_dir + '/bcg2/' + root_tile[i] + '_' + root_xcs[i] + '.bcg'
      outbcgname = np.append(outbcgname,name)
      
      name = input_dir + '/notes/' + root_tile[i] + '_' + root_xcs[i] + '.txt'
      innotesname = np.append(innotesname,name)      

# read in specZ
dt = [('ra',float), ('dec',float), ('z',float), ('ref',str)]
try:
  with open(nedspecname):
     nedspec = np.loadtxt(nedspecname, dtype=dt, skiprows=(0), delimiter=',')
except IOError:
  nedspec = 0
try:
  with open(photzspecname):
     photspec1 = pyfits.open(photzspecname)
     photspec = photspec1[1].data 
except:
  photspec = 0
  
# read in flags
inflagname = input_dir + '/flags/' + inflagtype + '.txt'
outflagname = input_dir + '/flags/' + outflagtype + '.txt'
flag_ind = np.loadtxt(inflagname).astype(int)
flags = np.array([])

a = len(jpegname)
#for i in range(1, a):
for i in flag_ind:
        print '\nWorking on: ', incatname[i],fitsname[i],jpegname[i], i, ' of ', a
        print 'Coordinates of current cluster:', coord.Longitude(clust_ra[i-1], unit=u.hour).degree,coord.Latitude(clust_dec[i-1], unit=u.degree).degree
        #go into the class des_image(), output as c
	c = des_image(incatname[i],fitsname[i],jpegname[i],intracername[i], inmembername[i], inboxname[i], inbcgname[i], outbcgname[i], innotesname[i], nedspec, photspec,choose, readbcgs, xcs_zs[i-1])


        if (os.path.isfile(outbcgname[i])  & owbcgs):
            while True:
                    overWrite = raw_input("Would you like to (over)write the BCG file? (y/n)->")
                    if overWrite.upper() == "Y":
                          np.savetxt(outbcgname[i], np.transpose([np.transpose(c.bcg_ra), np.transpose(c.bcg_dec), np.transpose(c.bcg_g_mag), np.transpose(c.bcg_r_mag), np.transpose(c.bcg_i_mag), np.transpose(c.bcg_z_mag),  \
                               np.transpose(c.bcg_gr_color),np.transpose(c.bcg_gi_color),np.transpose(c.bcg_ri_color),np.transpose(c.bcg_iz_color),np.transpose(c.bcg_rz_color)]),delimiter = ' ')
                          break
                    elif overWrite.upper() == "N":
                          break
                    else:
                          print "Invalid Option\n"
        elif not (os.path.isfile(outbcgname[i])):
             np.savetxt(outbcgname[i], np.transpose([np.transpose(c.bcg_ra), np.transpose(c.bcg_dec), np.transpose(c.bcg_g_mag), np.transpose(c.bcg_r_mag), np.transpose(c.bcg_i_mag), np.transpose(c.bcg_z_mag),  \
                        np.transpose(c.bcg_gr_color),np.transpose(c.bcg_gi_color),np.transpose(c.bcg_ri_color),np.transpose(c.bcg_iz_color),np.transpose(c.bcg_rz_color)]),delimiter = ' ')
        if c.flagged:
        	flags = np.append(flags, i)
open(outflagname, 'a+')
np.savetxt(outflagname, flags, delimiter= '\n', fmt = '%i') 


