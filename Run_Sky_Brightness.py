from Sky_Brightness_With_Moonlight import * 
from Throughputs import *
from Parameters import *

parser = OptionParser()
parser.add_option("-f", "--fieldname", type="string", default="DD", help="filter [%default]")
parser.add_option("-n", "--fieldid", type="int", default=290, help="filter [%default]")

opts, args = parser.parse_args()

To_Process=False
Plot_Brightness=True
Plot=False
Test=False
Test_airmass=False
 
def Calc_Integ(bandpass):
        resu=0.
        dlam=0
        for i,wave in enumerate(bandpass.wavelen):
            if i < len(bandpass.wavelen)-1:
                dlam=bandpass.wavelen[i+1]-wave
                resu+=dlam*bandpass.sb[i]/wave
            #resu+=dlam*bandpass.sb[i]

        return resu  

def Plot_Filters(sela,dict_posd,whata=[],whatb=[],legx='',legy=''):
       

    fige, axe = plt.subplots(ncols=2, nrows=3, figsize=(14,10))
   
    for j,band in enumerate(['u','g','r','i','z','y']):
        sela=sel[np.where(sel['filter']==band)]
        selac=selb[np.where(selb['filter']==band)]
        if j<2:
            k=0
        if j>= 2 and j < 4:
            k=1
        if j>=4:
            k=2

        axe[k][j%2].plot(sela[whata[0]],sela[whatb[0]]-sela[whatb[1]],'k.')
        axe[k][j%2].set_xlabel(r''+legx,{'fontsize': 10.})
        axe[k][j%2].set_ylabel(r''+legy,{'fontsize': 10.})
        axe[k][j%2].text(dict_posd[band][0], dict_posd[band][1], band, style='italic',
                         bbox={'facecolor':'yellow', 'alpha':0.5, 'pad':10})

def Plot_3D(sel):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(sel['moonPhase'],sel['filtSkyBrightness']-sel['mbsky_through'],sel['moonZD_DEG'], c='r', marker='o')


if To_Process:
    filename='MoonLight.pkl'
    List=['Observations_'+opts.fieldname+'_'+str(opts.fieldid)+'.pkl']
    
    outdir='Obs_minion_1016'
    thedict={}
    fieldid={}
    for i,name in enumerate(List):
        pkl_file = open(outdir+'/'+name,'rb')
        thedict[i]=pkl.load(pkl_file)
        fieldid[i]=np.int(name.split('_')[2].split('.')[0])

    print thedict[0].keys(),thedict[0]['dataSlice'].dtype.names

#('obsHistID', 'filtSkyBrightness', 'airmass', 'moonPhase', 'fieldRA', 'fieldDec', 'visitExpTime', 'expDate', 'filter', 'fieldID', 'fiveSigmaDepth', 'ditheredDec', 'expMJD', 'ditheredRA', 'rawSeeing')

    data=thedict[0]['dataSlice']

    transmission=Throughputs()

    mytype=[('filter',np.dtype('a15')),('filtSkyBrightness', np.float),('mbsky_through', np.float),('mbsky_through_moon', np.float),('moonPhase', np.float),('moonDiam', np.float),('year',np.int),('month',np.int),('day',np.int),('hour',np.int),('min',np.int),('sec',np.int),('MJD', np.float),('year_sunset',np.int),('month_sunset',np.int),('day_sunset',np.int),('hour_sunset',np.int),('min_sunset',np.int),('sec_sunset',np.int),('year_sunrise',np.int),('month_sunrise',np.int),('day_sunrise',np.int),('hour_sunrise',np.int),('min_sunrise',np.int),('sec_sunrise',np.int),('moonZD_DEG', np.float),('moon_airmass', np.float),('targetZD_DEG', np.float),('target_airmass', np.float),('fiveSigmaDepth', np.float),('m5_with_moon', np.float),('m5_darksky', np.float),('m5_recalc', np.float),('m5_recalc_atmthrough', np.float),('katm',np.float),('katm_aero',np.float),('katm_opsim',np.float),('airmass',np.float),('moonBr',np.float),('totBr',np.float),('twilight',np.float),('distance2moon',np.float)]

    tab_resu=np.zeros((60,1),dtype=[type for type in mytype])

    num=-1
    io=0
    reftime=0
    for i in range(len(data)):
    #print data['fieldRA'][i],data['fieldDec'][i],data['expMJD'][i],data['moonPhase'][i],data['airmass'][i],data['expDate'][i],data['filtSkyBrightness'][i],data['filter'][i]
        mysky=SkyBrightness(data['fieldRA'][i],data['fieldDec'][i],data['expMJD'][i],data['expDate'][i],data['filter'][i],data['airmass'][i])
    #mysky=SkyBrightness(data['fieldRA'][i],data['fieldDec'][i],49355.199664,190051,data['filter'][i])
        filtre=data['filter'][i]
        transmission.Load_Atmosphere(data['airmass'][i])
        myup=transmission.darksky.calcInteg(transmission.lsst_system[filtre])
               
        Tb=Calc_Integ(transmission.lsst_atmos[filtre])
        Sigmab=Calc_Integ(transmission.lsst_system[filtre])
        katm=-2.5*np.log10(Tb/Sigmab)
        
        Tb_aero=Calc_Integ(transmission.lsst_atmos_aerosol[filtre])
        katm_aero=-2.5*np.log10(Tb_aero/Sigmab)
    
        mbsky_through=-2.5*np.log10(myup/(3631.*Sigmab))

    #print 'yes here',filtre,mysky.new_skybrightness(),mbsky_through,mysky.moonPhase
        seeing=data['rawSeeing'][i]
        photParams = PhotometricParameters()
        param=parameters()
        Filter_Wavelength_Correction = np.power(500.0 / param.filterWave[filtre], 0.3)
        Airmass_Correction = math.pow(data['airmass'][i],0.6)
        FWHM_Sys = param.FWHM_Sys_Zenith * Airmass_Correction
        FWHM_Atm = seeing * Filter_Wavelength_Correction * Airmass_Correction
        finSeeing = param.scaleToNeff * math.sqrt(np.power(FWHM_Sys,2) + param.atmNeffFactor * np.power(FWHM_Atm,2))
        FWHMeff = SignalToNoise.FWHMgeom2FWHMeff(finSeeing)
        
        wavelen_min, wavelen_max, wavelen_step=transmission.lsst_system[filtre].getWavelenLimits(None,None,None)
        flatSed = Sed()
        flatSed.setFlatSED(wavelen_min, wavelen_max, wavelen_step)
        #flux0=np.power(10.,-0.4*mysky.new_skybrightness())
        flux0=np.power(10.,-0.4*data['filtSkyBrightness'][i])
        flatSed.multiplyFluxNorm(flux0)
        m5_calc=SignalToNoise.calcM5(flatSed,transmission.lsst_atmos_aerosol[filtre],transmission.lsst_system[filtre],photParams=photParams,FWHMeff=FWHMeff)
        m5_calc_old=SignalToNoise.calcM5(transmission.darksky,transmission.lsst_atmos_aerosol[filtre],transmission.lsst_system[filtre],photParams=photParams,FWHMeff=FWHMeff)
        Tscale = data['visitExpTime'][i]/ 30.0 * np.power(10.0, -0.4*(data['filtSkyBrightness'][i] - param.msky[filtre]))
        dCm = param.dCm_infinity[filtre] - 1.25*np.log10(1 + np.power(10.,0.8*param.dCm_infinity[filtre]- 1.)/Tscale)

        m5_recalc=dCm+param.Cm[filtre]+0.5*(data['filtSkyBrightness'][i]-21.)+2.5*np.log10(0.7/finSeeing)-param.kAtm[filtre]*(data['airmass'][i]-1.)+1.25*np.log10(data['visitExpTime'][i]/30.)
        
        m5_recalc_atmthrough=dCm+param.Cm[filtre]+0.5*(data['filtSkyBrightness'][i]-21.)+2.5*np.log10(0.7/finSeeing)-katm_aero*(data['airmass'][i]-1.)+1.25*np.log10(data['visitExpTime'][i]/30.)
        if filtre=='g' and data['expMJD'][i]>=59758 and data['expMJD'][i]<=59770.:
            diff=np.abs(data['filtSkyBrightness'][i]-mysky.new_skybrightness())
            print 'hello',m5_calc,m5_recalc_atmthrough,m5_calc_old,data['fiveSigmaDepth'][i],m5_recalc,data['filtSkyBrightness'][i],mysky.new_skybrightness(),data['filtSkyBrightness'][i]-mysky.new_skybrightness(),data['fieldRA'][i],data['fieldDec'][i],mbsky_through,24.*3600.*(data['expMJD'][i]-reftime)-34.
            if io == 0:
                reftime=data['expMJD'][i]
            io+=1
        #print katm,katm_aero,param.kAtm[filtre]
            #if io >9:
                #break
        num+=1
        if len(tab_resu) <= num:
            tab_resu=np.resize(tab_resu,(len(tab_resu)+100,1))

        (yy, mm, dd,hh,mn,sec) = mysky.mjd2gre(data['expMJD'][i])[:6]
        (yy_sunset, mm_sunset, dd_sunset,hh_sunset,mn_sunset,sec_sunset) = mysky.mjd2gre(mysky.sunSetTwilMJD)[:6]
        (yy_sunrise, mm_sunrise, dd_sunrise,hh_sunrise,mn_sunrise,sec_sunrise) = mysky.mjd2gre(mysky.sunRiseTwilMJD)[:6]
        moonZD_DEG=mysky.moonZD_DEG
        moon_airmass=mysky.moon_airmass
        targetZD_DEG=mysky.targetZD_DEG
        target_airmass=mysky.target_airmass


        tab_resu['filter'][num]=filtre
        tab_resu['filtSkyBrightness'][num]=data['filtSkyBrightness'][i]
        tab_resu['mbsky_through'][num]=mbsky_through
        tab_resu['mbsky_through_moon'][num]=mysky.new_skybrightness()
        tab_resu['moonPhase'][num]=mysky.moonPhase
        tab_resu['moonDiam'][num]=mysky.moonDiam
        tab_resu['moonBr'][num]=mysky.moonBr
        tab_resu['totBr'][num]=mysky.totBr
        tab_resu['twilight'][num]=mysky.twilight
        tab_resu['distance2moon'][num]=mysky.distance2moon_DEG
        

        tab_resu['year'][num]=yy
        tab_resu['month'][num]=mm
        tab_resu['day'][num]=dd
        tab_resu['hour'][num]=hh
        tab_resu['min'][num]=mn
        tab_resu['sec'][num]=sec
        tab_resu['MJD'][num]=data['expMJD'][i]

        tab_resu['year_sunset'][num]=yy_sunset
        tab_resu['month_sunset'][num]=mm_sunset
        tab_resu['day_sunset'][num]=dd_sunset
        tab_resu['hour_sunset'][num]=hh_sunset
        tab_resu['min_sunset'][num]=mn_sunset
        tab_resu['sec_sunset'][num]=sec_sunset

        tab_resu['year_sunrise'][num]=yy_sunrise
        tab_resu['month_sunrise'][num]=mm_sunrise
        tab_resu['day_sunrise'][num]=dd_sunrise
        tab_resu['hour_sunrise'][num]=hh_sunrise
        tab_resu['min_sunrise'][num]=mn_sunrise
        tab_resu['sec_sunrise'][num]=sec_sunrise
         
        tab_resu['moonZD_DEG'][num]=moonZD_DEG
        tab_resu['moon_airmass'][num]=moon_airmass
        tab_resu['targetZD_DEG'][num]=targetZD_DEG
        tab_resu['target_airmass'][num]=target_airmass

        tab_resu['fiveSigmaDepth'][num]=data['fiveSigmaDepth'][i]
        tab_resu['m5_with_moon'][num]=m5_calc
        tab_resu['m5_darksky'][num]=m5_calc_old
        tab_resu['m5_recalc'][num]=m5_recalc
        tab_resu['m5_recalc_atmthrough'][num]=m5_recalc_atmthrough

        tab_resu['airmass'][num]=data['airmass'][i]
        tab_resu['katm'][num]=katm
        tab_resu['katm_aero'][num]=katm_aero
        tab_resu['katm_opsim'][num]=param.kAtm[filtre]

    #if i > 100:
    #   break

    tab_resu=np.resize(tab_resu,(num+1,1))
    
    
    pkl_file_res = open('MoonLight.pkl','wb')
    pkl.dump(tab_resu, pkl_file_res)
    pkl_file_res.close()

if Plot_Brightness:
    filename='MoonLight.pkl'
    pkl_file = open(filename,'rb')
    tab_load=pkl.load(pkl_file)

    sel=tab_load

    theval=85.
    selb=sel[np.where(sel['moonZD_DEG']<theval)]
    
    figb, axb = plt.subplots(ncols=2, nrows=3, figsize=(14,10))
    
    dict_pos={}
    dict_pos['u']=(1.,-5.5)
    dict_pos['g']=(90.,-5.)
    dict_pos['r']=(90.,-3.5)
    dict_pos['i']=(15.,-4)
    dict_pos['z']=(12.,-2.5)
    dict_pos['y']=(12.,-1.6)

    for j,band in enumerate(['u','g','r','i','z','y']):
        sela=sel[np.where(sel['filter']==band)]
        selac=selb[np.where(selb['filter']==band)]
        if j<2:
            k=0
        if j>= 2 and j < 4:
            k=1
        if j>=4:
            k=2
        tot_label=[]
        tot_label.append(axb[k][j%2].plot(sela['moonPhase'],sela['filtSkyBrightness']-sela['mbsky_through'],'k.',label='Moon ZD > '+str(theval)+'$^{o}$'))
        tot_label.append(axb[k][j%2].plot(selac['moonPhase'],selac['filtSkyBrightness']-selac['mbsky_through'],'r.',label='Moon ZD < '+str(theval)+'$^{o}$'))
        axb[k][j%2].set_xlabel(r'Moon Phase',{'fontsize': 14.})
        axb[k][j%2].set_ylabel(r'$m_{sky}^{Opsim}-m_{sky}^{darksky}$',{'fontsize': 14.})
        #axb[k][j%2].legend(loc='center left',prop={'size':12})
        axb[k][j%2].legend(loc='best',prop={'size':10})
        axb[k][j%2].text(dict_pos[band][0], dict_pos[band][1], band, style='italic',
                         bbox={'facecolor':'yellow', 'alpha':0.5, 'pad':10})

    figa, axa = plt.subplots(ncols=2, nrows=3, figsize=(14,10))
    for j,band in enumerate(['u','g','r','i','z','y']):
        sela=sel[np.where(np.logical_and(sel['filter']==band,sel['moonBr']>=0.))]
        selac=selb[np.where(np.logical_and(selb['filter']==band,selb['moonBr']>=0))]
        if j<2:
            k=0
        if j>= 2 and j < 4:
            k=1
        if j>=4:
            k=2
        tot_label=[]
        tot_label.append(axa[k][j%2].plot(sela['distance2moon'],sela['filtSkyBrightness']-sela['mbsky_through'],'k.',label='Moon ZD > '+str(theval)+'$^{o}$'))
        tot_label.append(axa[k][j%2].plot(selac['distance2moon'],selac['filtSkyBrightness']-selac['mbsky_through'],'r.',label='Moon ZD < '+str(theval)+'$^{o}$'))
        axa[k][j%2].set_xlabel(r'Distance to Moon [deg]',{'fontsize': 14.})
        axa[k][j%2].set_ylabel(r'$m_{sky}^{Opsim}-m_{sky}^{darksky}$',{'fontsize': 14.})
        #axb[k][j%2].legend(loc='center left',prop={'size':12})
        axa[k][j%2].legend(loc='best',prop={'size':10})
        axa[k][j%2].text(dict_pos[band][0], dict_pos[band][1], band, style='italic',
                         bbox={'facecolor':'yellow', 'alpha':0.5, 'pad':10})

    figc, axc = plt.subplots(ncols=2, nrows=3, figsize=(14,10))
    for j,band in enumerate(['u','g','r','i','z','y']):
        sela=sel[np.where(sel['filter']==band)]
        selac=selb[np.where(selb['filter']==band)]
        if j<2:
            k=0
        if j>= 2 and j < 4:
            k=1
        if j>=4:
            k=2
        tot_label=[]
        tot_label.append(axc[k][j%2].plot(sela['moonPhase'],sela['totBr']-sela['mbsky_through'],'k.',label='Moon ZD > '+str(theval)+'$^{o}$'))
        tot_label.append(axc[k][j%2].plot(selac['moonPhase'],selac['totBr']-selac['mbsky_through'],'r.',label='Moon ZD < '+str(theval)+'$^{o}$'))
        axc[k][j%2].set_xlabel(r'Moon Phase',{'fontsize': 14.})
        axc[k][j%2].set_ylabel(r'$m_{sky}^{Opsim}-m_{sky}^{darksky}$',{'fontsize': 14.})
        #axb[k][j%2].legend(loc='center left',prop={'size':12})
        axc[k][j%2].legend(loc='best',prop={'size':10})
        axc[k][j%2].text(dict_pos[band][0], dict_pos[band][1], band, style='italic',
                         bbox={'facecolor':'yellow', 'alpha':0.5, 'pad':10})


if Plot:
    filename='MoonLight.pkl'
    pkl_file = open(filename,'rb')
    tab_load=pkl.load(pkl_file)

    sel=tab_load
    #sel=tab_load[np.where(np.logical_and(tab_load['MJD']>59750.418928,tab_load['MJD']<59887.071327))]
    #sel=tab_load[np.where(np.logical_and(tab_load['MJD']>=59758,tab_load['MJD']<=59770.))]
#sel=sel[np.where(sel['hour']<1)]
#plt.plot(sel['hour_sunset']-sel['hour'],sel['filtSkyBrightness']-sel['mbsky_through'],'ko')
    figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    #axa.hist(sel['m5_with_moon']-sel['fiveSigmaDepth'],bins=20)
    axa.plot(sel['moonPhase'],sel['m5_with_moon']-sel['fiveSigmaDepth'],'ko')

    #Plot_3D(sel)

    theval=80.
    selb=sel[np.where(sel['moonZD_DEG']<theval)]
    
    figb, axb = plt.subplots(ncols=2, nrows=3, figsize=(14,10))
    
    dict_pos={}
    dict_pos['u']=(1.,-5.5)
    dict_pos['g']=(90.,-5.)
    dict_pos['r']=(90.,-3.5)
    dict_pos['i']=(15.,-4)
    dict_pos['z']=(12.,-2.5)
    dict_pos['y']=(12.,-1.6)

    for j,band in enumerate(['u','g','r','i','z','y']):
        sela=sel[np.where(sel['filter']==band)]
        selac=selb[np.where(selb['filter']==band)]
        if j<2:
            k=0
        if j>= 2 and j < 4:
            k=1
        if j>=4:
            k=2
        tot_label=[]
        tot_label.append(axb[k][j%2].plot(sela['moonPhase'],sela['filtSkyBrightness']-sela['mbsky_through'],'k.',label='Moon ZD > '+str(theval)+'$^{o}$'))
        tot_label.append(axb[k][j%2].plot(selac['moonPhase'],selac['filtSkyBrightness']-selac['mbsky_through'],'r.',label='Moon ZD < '+str(theval)+'$^{o}$'))
        axb[k][j%2].set_xlabel(r'Moon Phase',{'fontsize': 14.})
        axb[k][j%2].set_ylabel(r'$m_{sky}^{Opsim}-m_{sky}^{darksky}$',{'fontsize': 14.})
        #axb[k][j%2].legend(loc='center left',prop={'size':12})
        axb[k][j%2].legend(loc='best',prop={'size':10})
        axb[k][j%2].text(dict_pos[band][0], dict_pos[band][1], band, style='italic',
                         bbox={'facecolor':'yellow', 'alpha':0.5, 'pad':10})
    #axb.plot(selb['moonPhase'],selb['filtSkyBrightness']-selb['mbsky_through'],'ko')
 
    dict_posb={}
    dict_posb['u']=(17.5,-0.003)
    dict_posb['g']=(17.5,-0.003)
    dict_posb['r']=(17.5,-0.003)
    dict_posb['i']=(16.5,-0.003)
    dict_posb['z']=(17.2,-0.003)
    dict_posb['y']=(17.05,-0.03)
    
    figc, axc = plt.subplots(ncols=2, nrows=3, figsize=(14,10))
    selb=sel[np.where(sel['moonZD_DEG']<85.)]
    coeff=24.*3600.
    for j,band in enumerate(['u','g','r','i','z','y']):
        sela=sel[np.where(sel['filter']==band)]
        selac=selb[np.where(selb['filter']==band)]
        if j<2:
            k=0
        if j>= 2 and j < 4:
            k=1
        if j>=4:
            k=2
        #axc[k][j%2].hist(sela['filtSkyBrightness']-sela['mbsky_through_moon'],bins=40)
        print 'mean rms',band,np.mean(sela['filtSkyBrightness']-sela['mbsky_through_moon']),np.std(sela['filtSkyBrightness']-sela['mbsky_through_moon'])
        #axc[k][j%2].plot(sela['MJD'],sela['filtSkyBrightness']-sela['mbsky_through_moon'],'k.')
        if len(sela) > 0.:
            #axc[k][j%2].plot(coeff*(sela['MJD']-sela['MJD'].min()),sela['target_airmass'],'k.')
            #axc[k][j%2].plot(sela['airmass'],sela['target_airmass'],'k.')
            axc[k][j%2].plot(sela['filtSkyBrightness'],sela['filtSkyBrightness']-sela['mbsky_through_moon'],'k.')  
        axc[k][j%2].text(dict_posb[band][0], dict_posb[band][1], band, style='italic',
                         bbox={'facecolor':'yellow', 'alpha':0.5, 'pad':10})
        axc[k][j%2].set_xlabel(r'$m_{sky}^{Opsim}$',{'fontsize': 10.})
        axc[k][j%2].set_ylabel(r'$m_{sky}^{Opsim}-m_{sky}^{throughput+Moon}$',{'fontsize': 10.})

    Plot_Filters(sela,dict_posb,whata=['airmass'],whatb=['airmass','target_airmass'],legx='$airmass^{Opsim}$',legy='$airmass^{Opsim}-airmass^{throughput+Moon}$')

    Plot_Filters(sela,dict_posb,whata=['filtSkyBrightness'],whatb=['filtSkyBrightness','mbsky_through_moon'],legx='$m_{sky}^{Opsim}$',legy='$m_{sky}^{Opsim}-m_{sky}^{throughput+Moon}$')

    dict_posc={}
    dict_posc['u']=(21.5,0.02)
    dict_posc['g']=(22.5,0.09)
    dict_posc['r']=(22.2,0.11)
    dict_posc['i']=(21.7,0.10)
    dict_posc['z']=(21.6,0.11)
    dict_posc['y']=(21.05,0.09)
    Plot_Filters(sela,dict_posc,whata=['fiveSigmaDepth'],whatb=['fiveSigmaDepth','m5_with_moon'],legx='$m_{5}^{Opsim}$',legy='$m_{5}^{Opsim}-m_{5}^{Troughput+Moon}$')

    dict_posd={}
    dict_posd['u']=(0.10,0.25)
    dict_posd['g']=(0.05,0.16)
    dict_posd['r']=(0.03,0.10)
    dict_posd['i']=(0.02,0.08)
    dict_posd['z']=(0.015,0.07)
    dict_posd['y']=(0.04,0.16)


    fige, axe = plt.subplots(ncols=2, nrows=3, figsize=(14,10))
   
    for j,band in enumerate(['u','g','r','i','z','y']):
        sela=sel[np.where(sel['filter']==band)]
        selac=selb[np.where(selb['filter']==band)]
        if j<2:
            k=0
        if j>= 2 and j < 4:
            k=1
        if j>=4:
            k=2

        axe[k][j%2].plot(sela['katm_opsim']*(sela['airmass']-1),sela['katm_aero']*(sela['airmass']-1),'k.')
        axe[k][j%2].set_xlabel(r'$katm^{Opsim}*(airmass-1)$',{'fontsize': 10.})
        axe[k][j%2].set_ylabel(r'$katm^{Throughput}*(airmass-1)$',{'fontsize': 10.})
        axe[k][j%2].text(dict_posd[band][0], dict_posd[band][1], band, style='italic',
                         bbox={'facecolor':'yellow', 'alpha':0.5, 'pad':10})


if Test:
    filename='MoonLight.pkl'
    List=['Observations_'+opts.fieldname+'_'+str(opts.fieldid)+'.pkl']
    
    outdir='Obs_minion_1016'
    thedict={}
    fieldid={}
    for i,name in enumerate(List):
        pkl_file = open(outdir+'/'+name,'rb')
        thedict[i]=pkl.load(pkl_file)
        fieldid[i]=np.int(name.split('_')[2].split('.')[0])

    print thedict[0].keys(),thedict[0]['dataSlice'].dtype.names

#('obsHistID', 'filtSkyBrightness', 'airmass', 'moonPhase', 'fieldRA', 'fieldDec', 'visitExpTime', 'expDate', 'filter', 'fieldID', 'fiveSigmaDepth', 'ditheredDec', 'expMJD', 'ditheredRA', 'rawSeeing')

    data=thedict[0]['dataSlice']

    shift=100.
    mjd_ref=data['expMJD'][i]+shift
    for i in range(len(data)):
        if abs(data['expMJD'][i]-mjd_ref)<100.:
            print data['fieldRA'][i],data['fieldDec'][i],data['expMJD'][i],data['moonPhase'][i],data['airmass'][i],data['expDate'][i],data['filtSkyBrightness'][i],data['filter'][i],data['expMJD'][i]-mjd_ref
            (yy, mm, dd,hh,mn,sec) = mjd2gre_test(data['expMJD'][i])[:6]
            (yy_ref, mm_ref, dd_ref,hh_ref,mn_ref,sec_ref) = mjd2gre_test(mjd_ref)[:6]
            print yy, mm, dd,hh,mn,sec
            print yy_ref, mm_ref, dd_ref,hh_ref,mn_ref,sec_ref
            """
            ra=2.093429 
            dec=-1.082474
            """
            """
            ra=4.189756 
            dec=-1.082474
            mjd= 59874.32285+119.57636424
            """
            ra=data['fieldRA'][i]
            dec=data['fieldDec'][i]
            mjd=data['expMJD'][i]
            
            mysky=SkyBrightness(ra,dec,mjd,data['expDate'][i],data['filter'][i],data['airmass'][i])
            mysky=SkyBrightness(ra,dec,mjd_ref,data['expDate'][i],data['filter'][i],data['airmass'][i])
            break


if Test_airmass:

    ra=2.093429 
    dec=-1.082474

    mjd_orig=59000.

    mjd=[]
    airmass=[]
    airmass_pal=[]
    hours=[]
    while mjd_orig <= 59100.:
        mysky=SkyBrightness(ra,dec,mjd_orig,-1.,'g',-1.)
        mjd_orig +=0.01
        mjd.append(mjd_orig)
        airmass.append(mysky.target_airmass)
        airmass_pal.append(mysky.airmass_pal)
        (yy, mm, dd,hh,mn,sec) = mjd2gre_test(mjd_orig)[:6]
        hours.append(hh+mn/60.+sec/3600.)
        

    
    plt.plot(mjd,airmass,'b.')
    #plt.plot(mjd,airmass_pal,'r.')
    
    #plt.plot(hours,airmass,'b.')
    #plt.plot(hours,airmass_pal,'r.')
    


"""
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(sel['moonPhase'],sel['filtSkyBrightness']-sel['mbsky_through'],sel['moonZD_DEG'], c='r', marker='o')
"""
#plt.plot(sel['moon_airmass'],sel['filtSkyBrightness']-sel['mbsky_through'],'ko')
#plt.plot(sel['hour_sunset']+sel['min_sunset']/60.,sel['moonPhase'],'ko')

plt.show()
