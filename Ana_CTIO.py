from Sky_Brightness_With_Moonlight import *
from astropy.table import Table,vstack
from astropy.table.table import Column
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from Throughputs import Throughputs

ADD_OpSim=False
ANA_RESU=False
Compare_Throughputs=True

if ADD_OpSim:

    filename='SIMLIB_DUMP_DES-griz_mod.DAT'

    table_ctio=ascii.read(filename,fast_reader=False)

    print table_ctio


    table_opsim=Table(names=('SKYOPSIM','MOONPHASE','MOONZD_DEG','SKYOPSIM_NOFCORR'),dtype=('f8','f8','f8','f8'))

    for val in table_ctio:
        date=0.
        mysky=SkyBrightness(np.deg2rad(val['RA']),np.deg2rad(val['DEC']),val['MJD'],date,val['BAND'],0.)
        print mysky.new_skybrightness(),val['SKYMAG'],val['BAND']
    
        table_opsim.add_row((mysky.new_skybrightness(),mysky.moonPhase,mysky.moonZD_DEG,mysky.totBr))
    #break

    table_ctio.add_columns(table_opsim.columns.values())
    
    print table_ctio
    ascii.write(table_ctio,filename.replace('mod','OpSim'))

if ANA_RESU:

    def Plot(sel,varx='MOONPHASE',vary=['SKYOPSIM','SKYMAG']):
        theval=85.
        selb=sel[np.where(sel['MOONZD_DEG']<theval)]
        figb, axb = plt.subplots(ncols=2, nrows=2, figsize=(14,10))
    
        dict_pos={}
        dict_pos['u']=(1.,-5.5)
        dict_pos['g']=(90.,-5.)
        dict_pos['r']=(90.,-3.5)
        dict_pos['i']=(15.,-4)
        dict_pos['z']=(12.,-2.5)
        dict_pos['y']=(12.,-1.6)
        
        for j,band in enumerate(['g','r','i','z']):
            sela=sel[np.where(sel['BAND']==band)]
            selac=selb[np.where(selb['BAND']==band)]
            if j<2:
                k=0
            if j>= 2 and j < 4:
                k=1
            if j>=4:
                k=2
            tot_label=[]
            tot_label.append(axb[k][j%2].plot(sela[varx],sela[vary[0]]-sela[vary[1]],'k.',label='Moon ZD > '+str(theval)+'$^{o}$'))
            tot_label.append(axb[k][j%2].plot(selac[varx],selac[vary[0]]-selac[vary[1]],'r.',label='Moon ZD < '+str(theval)+'$^{o}$'))
            axb[k][j%2].set_xlabel(r'Moon Phase',{'fontsize': 14.})
            axb[k][j%2].set_ylabel(r'$m_{sky}^{Opsim}-m_{sky}^{CTIO}$',{'fontsize': 14.})
        #axb[k][j%2].legend(loc='center left',prop={'size':12})
            axb[k][j%2].legend(loc='best',prop={'size':10})
            axb[k][j%2].text(dict_pos[band][0], dict_pos[band][1], band, style='italic',
                             bbox={'facecolor':'yellow', 'alpha':0.5, 'pad':10})
            axb[k][j%2].set_title(band,loc='left')
    def Plot_Hist(sel):
        theval=85.

        selb=sel[np.where(sel['MOONZD_DEG']<theval)]
        sela=sel[np.where(sel['MOONZD_DEG']>=theval)]
        figb, axb = plt.subplots(ncols=2, nrows=2, figsize=(14,10))
    
        dict_pos={}
        dict_pos['u']=(1.,-5.5)
        dict_pos['g']=(90.,-5.)
        dict_pos['r']=(90.,-3.5)
        dict_pos['i']=(15.,-4)
        dict_pos['z']=(12.,-2.5)
        dict_pos['y']=(12.,-1.6)
        
        for j,band in enumerate(['g','r','i','z']):
            selab=sela[np.where(sela['BAND']==band)]
            selac=selb[np.where(selb['BAND']==band)]
            if j<2:
                k=0
            if j>= 2 and j < 4:
                k=1
            if j>=4:
                k=2

            selcc=sela[np.where(np.abs(sela['SKYOPSIM']-sela['SKYMAG'])<2.)]
            
            print band,np.mean(selcc['SKYOPSIM']-selcc['SKYMAG']),np.std(selcc['SKYOPSIM']-selcc['SKYMAG'])
            axb[k][j%2].hist(selab['SKYOPSIM']-selab['SKYMAG'],bins=40.,histtype='step',range=[-2.,2.],color='b',label='Moon ZD > '+str(theval)+'$^{o}$')
            axb[k][j%2].hist(selac['SKYOPSIM']-selac['SKYMAG'],bins=40.,histtype='step',range=[-2.,2.],color='r',label='Moon ZD <='+str(theval)+'$^{o}$')
            axb[k][j%2].set_ylabel(r'Number of Entries',{'fontsize': 14.})
            axb[k][j%2].set_xlabel(r'$m_{sky}^{Opsim}-m_{sky}^{CTIO}$',{'fontsize': 14.})
            axb[k][j%2].legend(loc='best',prop={'size':10})
            axb[k][j%2].set_title(band,loc='left')
    table_ctio=ascii.read('SIMLIB_DUMP_DES-griz_OpSim.DAT',fast_reader=False)
    sel=table_ctio


    def Plot_OpSim_Corr(phases,filterOffset):
        
        figb, axb = plt.subplots(ncols=2, nrows=3, figsize=(14,10))
        for j,band in enumerate(['u','g','r','i','z']):
            if j<2:
                k=0
            if j>= 2 and j < 4:
                k=1
            if j>=4:
                k=2
                
            vals=[]
            for val in phases:
                vals.append(filterOffset[band,val])
            axb[k][j%2].plot(phases,vals,'ko')
            axb[k][j%2].set_ylabel(r'OpSim corr [mag]',{'fontsize': 14.})
            axb[k][j%2].set_xlabel(r'Moon Phase',{'fontsize': 14.})

    Plot(sel)
    Plot(sel,varx='MOONPHASE',vary=['SKYOPSIM_NOFCORR','SKYMAG'])
    Plot_Hist(sel)


    #OpSim filter corrections:

    filterOffset = {}
        
        # Corrections for moonPhase = 0 percent (new moon)
    filterOffset['u', 0.] = 0.66
    filterOffset['g', 0.] = 0.41
    filterOffset['r', 0.] = -0.28
    filterOffset['i', 0.] = -1.36
    filterOffset['z', 0.] = -2.15
    
        # Corrections for moonPhase = 18 percent
    filterOffset['u', 18.] = 0.28
    filterOffset['g', 18.] = 0.30
    filterOffset['r', 18.] = -0.19
    filterOffset['i', 18.] = -1.17
    filterOffset['z', 18.] = -1.99
    
        # Corrections for moonPhase = 50 percent
    filterOffset['u', 50.] = -1.05
    filterOffset['g', 50.] = 0.03
    filterOffset['r', 50.] = 0.02
    filterOffset['i', 50.] = -0.96
    filterOffset['z', 50.] = -1.78
    
        # Corrections for moonPhase = 80 percent
    filterOffset['u', 80.] = -1.83
    filterOffset['g', 80.] = -0.08
    filterOffset['r', 80.] = 0.10
    filterOffset['i', 80.] = -0.78
    filterOffset['z', 80.] = -1.54
    
        # Corrections for moonPhase = 100 percent (full moon)
    filterOffset['u', 100.] = -2.50
    filterOffset['g', 100.] = -0.35
    filterOffset['r', 100.] = 0.31
    filterOffset['i', 100.] = -0.47
    filterOffset['z', 100.] = -1.16
    
    
    phases=[0.,18.,50.,80.,100.]
    
    Plot_OpSim_Corr(phases,filterOffset)



    plt.show()

if Compare_Throughputs:
    
    #this file was taken from http://www.ctio.noao.edu/noao/node/1033
    throughputs_ctio=ascii.read('DECam_filters.csv',fast_reader=False)

    print throughputs_ctio
    throughputs_lsst=Throughputs()

    style=[',',',',',',',']
    filtercolors = {'u':'b', 'g':'c', 'r':'g', 'i':'y', 'z':'r', 'y':'m'}
    band_corr={'u':'u','g':'g','r':'r','i':'i','z':'z','y':'Y'}
    for i,band in enumerate(['u','g','r','i','z','y']):
        plt.plot(throughputs_lsst.lsst_system[band].wavelen,throughputs_lsst.lsst_system[band].sb,linestyle='--',color=filtercolors[band], label='LSST')
        #sel=throughputs_ctio[np.where(throughputs_ctio['']==band)]
        plt.plot(throughputs_ctio['wavelength'],throughputs_ctio[band_corr[band]]/throughputs_ctio['atm'],linestyle='-',color=filtercolors[band], label='DES')
    plt.xlabel('Wavelength (nm)')    
    plt.ylabel('Sb (0-1)')
    plt.title('System throughput')
    plt.legend(loc='upper right', fontsize='smaller', fancybox=True, numpoints=1)
    plt.show()
