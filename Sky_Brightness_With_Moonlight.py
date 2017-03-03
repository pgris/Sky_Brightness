import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u
import math
import healpy as hp
from astropy.table import Table
from optparse import OptionParser
import palpy as pal

from Throughputs import Throughputs

DEG2RAD = math.pi / 180.    # radians = degrees * DEG2RAD
RAD2DEG = 180. / math.pi    # degrees = radians * RAD2DEG
TWOPI = 2 * math.pi
DAY = 86400.
  
def Calc_Integ(bandpass):
    resu=0.
    dlam=0
    for i,wave in enumerate(bandpass.wavelen):
        if i < len(bandpass.wavelen)-1:
            dlam=bandpass.wavelen[i+1]-wave
            resu+=dlam*bandpass.sb[i]/wave
            #resu+=dlam*bandpass.sb[i]

    return resu  

class SkyBrightness:
    def __init__(self,ra,dec,mjd,date,filtre,extinction=0.172, skyBrightness=21.587):
        
        self.ra_RAD=ra
        self.dec_RAD=dec
        self.mjd=mjd
        self.k=extinction
        self.skyBrightness=skyBrightness
        self.date=date
        self.filter=filtre
        self.lon_RAD = -70.7494 *DEG2RAD #obs site long (from sims_operations/conf/system/SiteCP.conf) 
        self.lat_RAD = -30.2444* DEG2RAD #obs site lat (from sims_operations/conf/system/SiteCP.conf)
        self.lst_RAD  = self.Get_Sidereal_Time_at_Site() #local sidereal time at site (radians)
        # TWILIGHT LIMITS
        # Altitude of the Sun in degrees that define the twilight.
        # When the sun is above this limit and below the night limit, a special
        # twilight factor is included in the sky brightness model
        self.SunAltitudeTwilightLimit = -18.0
        self.SunAltitudeNightLimit = -12.0 
        self.n1dp92104 = (1. / 0.92104)
        self.log34p08 = math.log(34.08)
        self.simEpoch=59580
        #self.simEpoch=53371
        self.twilightBrightness = 17.3
        self.fluxTwilight = 10. ** ((-self.twilightBrightness) / 2.5)

        self.MoonProfile()
        self.Moon_Distance()
        alpha=self.Get_alpha()
        self.distance2moon_RAD = pal.dsep(self.moonRA_RAD, self.moonDec_RAD, self.ra_RAD, self.dec_RAD)
        self.distance2moon_DEG = self.distance2moon_RAD * RAD2DEG
        lha_RAD = self.lst_RAD - self.ra_RAD
        (az_RAD, d1, d2, alt_RAD, d4, d5, pa_RAD, d7, d8) = pal.altaz(lha_RAD, self.dec_RAD,self.lat_RAD)
        # Altitude -> Zenith distance (degrees)
        targetZD_DEG = 90. - (alt_RAD* RAD2DEG)
        # Altitude -> Zenith distance (radian)
        zd_RAD = 1.5707963 - alt_RAD
        # Airmass
        #am = slalib.sla_airmas (zd_RAD)
        am = pal.airmas(zd_RAD)
        
        # Compute Raileigh scattering
        rs = (10. ** 5.36) * (1.06 + (math.cos(self.distance2moon_RAD)) ** 2.)

        # Compute Mie scattering
        if self.distance2moon_DEG > 10.:
            ms = 10. ** (6.15 - (self.distance2moon_DEG / 40.))
        else:
            ms = 6.2E7 / (self.distance2moon_DEG ** 2.)

        # Compute illumninace of the Moon
        i = 10. ** (-0.4 * (3.84 + 0.026 * abs(alpha) + 4.E-9 * (alpha) ** 4.))

        # Compute optical pathlength (in units of airmass) as a
        # function of the Zenith Distance (ZD)
        x = lambda z: (1. - 0.96 * (math.sin(z * DEG2RAD)) ** 2.) ** -0.5

        #print 'airmass here',x(self.moonZD_DEG),x(targetZD_DEG),am

         # Put all together (nanoLamberts)!

        moonBr = (rs + ms) * i * 10. ** (-0.4 * self.k * x(self.moonZD_DEG))
        moonBr *= (1. - 10. ** (-0.4 * self.k * x(targetZD_DEG)))

        #  Taper brightness to 0, at moon altitude = self.sunAltitudeNightLimit
        if self.moonAlt_DEG < self.SunAltitudeNightLimit:
            moonBr = 0.
        elif self.moonAlt_DEG < 0. and self.moonAlt_DEG >= self.SunAltitudeNightLimit:
            moonBr *= (1 - self.moonAlt_DEG / self.SunAltitudeNightLimit)

        # Now, compute the dark sky brightness (in nanoLamberts) at
        # the Zenith from the same value in mag/arcsec^s
        skyBrightness = 34.08 * math.exp(20.7233 - 0.92104 * skyBrightness)

        # ... and at our position in the sky (targetZD_DEG)
        skyBr = skyBrightness * x(targetZD_DEG) * 10. ** (-0.4 * self.k * (x(targetZD_DEG) - 1.))

        # Add the brightness of the moon to that of the dark sky (nanoLamberts)
        totBr = moonBr + skyBr

        # Transform it into mag/arcsec^2
        # totBr = (1. / 0.92104) * math.log (34.08 * math.exp (20.7233) / totBr)
        totBr = self.n1dp92104 * (self.log34p08 + 20.7233 - math.log(totBr))

        #print 'totbr',totBr

        #Get the Twilight profile

        self.TwilightProfile()
      
        #if self.date < self.sunSetTwil or self.date > self.sunRiseTwil:
        (yy, mm, dd,hh,min,sec) = self.mjd2gre(mjd)[:6]
        if (hh < 20 and self.date > self.sunRiseTwil) or (hh>20  and self.date < self.sunSetTwil):
            totBr = -2.5 * math.log10(10. ** ((-totBr) / 2.5) + self.fluxTwilight)
            
        
        #print 'after twilight',totBr,self.date,self.sunSetTwil,self.sunRiseTwil

        #print 'filterskybrightness',self.FilterSkyBrightness(totBr)
        
       
        self.skybrightness_moon=self.FilterSkyBrightness(totBr)

    def new_skybrightness(self):

        return self.skybrightness_moon

    def Get_Sidereal_Time_at_Site(self):

        lst_RAD = pal.gmst(self.mjd) + self.lon_RAD
        if lst_RAD < 0:
            lst_RAD += TWOPI

        return lst_RAD

    def MoonProfile(self):
        (self.moonRA_RAD, self.moonDec_RAD, self.moonDiam) = pal.rdplan(self.mjd, 3, self.lon_RAD, self.lat_RAD)
        #print 'hello moon',self.moonRA_RAD, self.moonDec_RAD, self.moonDiam
        
        #getting the moon phase

        v6moon = pal.dmoon(self.mjd)
        
        Rmatrix = pal.prenut(2000.0, self.mjd)
        
        xyzMoon2000 = pal.dmxv(Rmatrix, v6moon)
       
        (moonra, moondec) = pal.dcc2s(xyzMoon2000)
        moonra = pal.dranrm(2.0 * math.pi + moonra)
        sun12 = pal.evp(self.mjd, 2000.0)
        sun3heliocentric = sun12[3]
        for i in range(3):
            sun3heliocentric[i] *= -1
        (sunra, sundec) = pal.dcc2s(sun3heliocentric)
        sunra = pal.dranrm(2.0 * math.pi + sunra)
        moonsunsep = pal.dsep(sunra, sundec, moonra, moondec)
        self.moonPhase = ((1.0 - math.cos(moonsunsep)) / 2.0) * 100

        #print 'and the phase',self.moonPhase

    def Moon_Distance(self):
        # Compute moon altitude in radians
        self.moonha_RAD = self.lst_RAD - self.moonRA_RAD

        (moonAz_RAD, d1, d2, self.moonAlt_RAD, d4, d5, d6, d7, d8) = \
            pal.altaz(self.moonha_RAD, self.moonDec_RAD, self.lat_RAD)
        self.moonAlt_DEG = self.moonAlt_RAD * RAD2DEG
        
        # compute moon's zenith distance
        if self.moonAlt_RAD < 0.:
            moonZD_RAD = 1.5707963
        else:
            moonZD_RAD = 1.5707963 - self.moonAlt_RAD
            
        self.moonZD_DEG = moonZD_RAD * RAD2DEG

    def Get_alpha(self):
        # Compute the phase angle given the illuminated fraction
        # formula from http://astro.ft.uam.es/TJM/tjm/webpaginas/practicas/ephemeris/moon.js
        alpha = math.acos(2. * self.moonPhase  / 100. - 1.) * RAD2DEG
        alpha = self.normalize(alpha, min=0., max=180., degrees=True)
        
        return alpha

    def normalize(self,angle, min=0., max=None, degrees=True):
        """
        Given an angle, make sure that its value is within the range
        [min, max].
        
        If degrees=True the ammount to add/subtract is 360.
        If degrees=False the ammount to add/subtract is 2*pi
        
        If min/max is not specified, it reverts to its default value
        (which depends on the value of degrees):
        min = 0
        max = 360./2*pi
        
        Return the normalized angle in its original units
        """
        angle = float(angle)
        
        if min is None:
            min = 0.
        if degrees:
            if max is None:
                max = 360.
            addit = 360.
        else:
            if max is None:
                max = TWOPI
            addit = TWOPI

        while angle < min:
            angle += addit
        while angle > max:
            angle -= addit
        return angle


    def TwilightProfile(self):
         (sunRise, sunSet, sunRiseMJD, sunSetMJD, self.sunRiseTwil,
            self.sunSetTwil) = self.getTwilightSunriseSunset(self.date)

         #print 'TwilightProfile','sunRiseMJD=',sunRiseMJD,'sunSetMJD=',sunSetMJD,'sunRiseTwil=',self.sunRiseTwil,'sunSetTwil=',self.sunSetTwil

    def getTwilightSunriseSunset(self, date):
        """
        Compute the time of the astronomical twilight for Sun raise
        and set times at the location on Earth specified by
        (self.latitude_RAD, self.longitude_RAD).

        Input
        date:   date in seconds from Jan 1 of the simulated year.

        Return
        A four-element array:
        (sunriseTime, sunsetTime, sunriseMJD, sunsetMJD)
        sunriseTime & sunsetTime are in seconds from Jan 1 of simulated year;
        sunriseMJD and sunsetMJD are MJD equivalents.
        """
        # Convert date in MJD
        mjd = (float(self.date) / float(DAY)) + self.simEpoch

        #mjd=49353.
        # MJD -> calendar date
        #mjd=self.mjd
        (yy, mm, dd,hh,min,sec) = self.mjd2gre(mjd)[:6]

        #print 'hello',yy,mm,dd,hh,min,sec
        # Compute sunset and sunrise at twilight. Return values are in
        # decimal hours.
        import Sun
        s = Sun.Sun()
        #(sunRise, sunSet) = s.sunRiseSet(yy, mm, dd,
        #                        self.longitude_RAD * RAD2DEG, self.latitude_RAD * RAD2DEG)
        # Following set to nauticalTwilight
        (sunRise, sunSet) = s.__sunriset__(yy, mm, dd, self.lon_RAD * RAD2DEG,
                                           self.lat_RAD * RAD2DEG, self.SunAltitudeNightLimit, 0)
        # -12.0, 0)
        (sunRiseTwil, sunSetTwil) = s.__sunriset__(yy, mm, dd, self.lon_RAD * RAD2DEG,
                                                   self.lat_RAD * RAD2DEG, self.SunAltitudeTwilightLimit,
                                                   0)
        
        # -18.0, 0)

        # Compute MJD values for sunrise and sunset
        sunSetMJD = int(mjd) + (sunSet / 24.)
        sunRiseMJD = int(mjd) + (sunRise / 24.)
        self.sunSetTwilMJD = int(mjd) + (sunSetTwil / 24.)
        self.sunRiseTwilMJD = int(mjd) + (sunRiseTwil / 24.)
        #print 'ici', sunSetTwilMJD,sunRiseTwilMJD
        # MJD -> simulated seconds
        sunsetDate = (sunSetMJD - self.simEpoch) * float(DAY)
        sunriseDate = (sunRiseMJD - self.simEpoch) * float(DAY)
        sunsetTwilDate = (self.sunSetTwilMJD - self.simEpoch) * float(DAY)
        sunriseTwilDate = (self.sunRiseTwilMJD - self.simEpoch) * float(DAY)
        


        #print 'twilight',self.mjd2gre(sunSetTwilMJD)[:6],self.mjd2gre(sunRiseTwilMJD)[:6]
        return (sunriseDate, sunsetDate, sunRiseMJD, sunSetMJD, sunriseTwilDate, sunsetTwilDate)

    def mjd2gre(self,mjd):
        # Use SLALIB to convert from MJD to (year, month, day)
    #(year, month, day, fraction, error) = slalib.sla_djcl (mjd)
        (year, month, day, fraction) = pal.djcl(mjd)
        
    #if (error):
    #    msg = 'Fatal error: MJD must correspond to a date later than 4701BC March 1'
    #    raise (SyntaxError, msg)
        
    # Now take care of the fraction of day
        hh = math.floor(24. * fraction)
        mm = math.floor(1440. * (fraction - hh / 24.))
        ss = 86400. * (fraction - hh / 24. - mm / 1440.)
        return (year, month, day, int(hh), int(mm), ss)

    def Corrections_Filters_for_moonPhase(self):
        
        self.skyBrightKeys = [0, 18, 50, 80, 100]
        self.filterOffset = {}
        
        # Corrections for moonPhase = 0 percent (new moon)
        self.filterOffset['u', 0.] = 0.66
        self.filterOffset['g', 0.] = 0.41
        self.filterOffset['r', 0.] = -0.28
        self.filterOffset['i', 0.] = -1.36
        self.filterOffset['z', 0.] = -2.15
        
        # Corrections for moonPhase = 18 percent
        self.filterOffset['u', 18.] = 0.28
        self.filterOffset['g', 18.] = 0.30
        self.filterOffset['r', 18.] = -0.19
        self.filterOffset['i', 18.] = -1.17
        self.filterOffset['z', 18.] = -1.99
        
        # Corrections for moonPhase = 50 percent
        self.filterOffset['u', 50.] = -1.05
        self.filterOffset['g', 50.] = 0.03
        self.filterOffset['r', 50.] = 0.02
        self.filterOffset['i', 50.] = -0.96
        self.filterOffset['z', 50.] = -1.78
        
        # Corrections for moonPhase = 80 percent
        self.filterOffset['u', 80.] = -1.83
        self.filterOffset['g', 80.] = -0.08
        self.filterOffset['r', 80.] = 0.10
        self.filterOffset['i', 80.] = -0.78
        self.filterOffset['z', 80.] = -1.54
        
        # Corrections for moonPhase = 100 percent (full moon)
        self.filterOffset['u', 100.] = -2.50
        self.filterOffset['g', 100.] = -0.35
        self.filterOffset['r', 100.] = 0.31
        self.filterOffset['i', 100.] = -0.47
        self.filterOffset['z', 100.] = -1.16


    def FilterSkyBrightness(self,skyBrightness):
        
        self.Corrections_Filters_for_moonPhase()
        moonPhase_PERCENT=self.moonPhase 

        if self.filter == 'y':
            filterSkyBright = 17.3
        else:      # g,r,i,z,u
            # If moon below horizon, use new moon offset for filter
            # brightness - MM
            if math.degrees(self.moonAlt_RAD) <= -6.0:
                adjustBright = self.filterOffset[self.filter, 0.]

            # Interpolate if needed. Note: moonPhase is a float not int
            elif moonPhase_PERCENT not in self.skyBrightKeys:
                i = 0
                while self.skyBrightKeys[i] < moonPhase_PERCENT:
                    i += 1

                # find upper and lower bound
                upperMoonPhase = self.skyBrightKeys[i]
                lowerMoonPhase = self.skyBrightKeys[i - 1]
                lowerAdjustBright = self.filterOffset[self.filter, lowerMoonPhase]
                upperAdjustBright = self.filterOffset[self.filter, upperMoonPhase]
                # linear interpolation
                diffMoonPhase = upperMoonPhase - lowerMoonPhase
                diffAdjustBright = upperAdjustBright - lowerAdjustBright
                moonPhaseAdj = moonPhase_PERCENT - lowerMoonPhase
                adjustBright = lowerAdjustBright + (((moonPhaseAdj) * (diffAdjustBright)) / (diffMoonPhase))

            else:          # moon not set and moon phase is key
                adjustBright = self.filterOffset[self.filter, moonPhase_PERCENT]
            #print 'adjust',adjustBright
            filterSkyBright = skyBrightness + adjustBright

            # z sky brightness should never be under 17.0
            if self.filter == 'z' and filterSkyBright < 17.0:
                filterSkyBright = 17.0

        # If twilight, set brightness for z and y
        (yy, mm, dd,hh,min,sec) = self.mjd2gre(self.mjd)[:6]
        if (hh < 20 and self.date > self.sunRiseTwil) or (hh>20  and self.date < self.sunSetTwil):
            if self.filter == 'z' or self.filter == 'y':
                filterSkyBright = 17.0

        return filterSkyBright


parser = OptionParser()
parser.add_option("-f", "--fieldname", type="string", default="DD", help="filter [%default]")
parser.add_option("-n", "--fieldid", type="int", default=290, help="filter [%default]")

opts, args = parser.parse_args()

To_Process=False
filename='MoonLight.pkl'

if To_Process:
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

    mytype=[('filter',np.dtype('a15')),('filtSkyBrightness', np.float),('mbsky_through', np.float),('mbsky_through_moon', np.float),('moonPhase', np.float),('moonDiam', np.float),('year',np.int),('month',np.int),('day',np.int),('hour',np.int),('min',np.int),('sec',np.int),('MJD', np.float),('year_sunset',np.int),('month_sunset',np.int),('day_sunset',np.int),('hour_sunset',np.int),('min_sunset',np.int),('sec_sunset',np.int),('year_sunrise',np.int),('month_sunrise',np.int),('day_sunrise',np.int),('hour_sunrise',np.int),('min_sunrise',np.int),('sec_sunrise',np.int)]

    tab_resu=np.zeros((60,1),dtype=[type for type in mytype])

    num=-1
    for i in range(len(data)):
    #print data['fieldRA'][i],data['fieldDec'][i],data['expMJD'][i],data['moonPhase'][i],data['airmass'][i],data['expDate'][i],data['filtSkyBrightness'][i],data['filter'][i]
        mysky=SkyBrightness(data['fieldRA'][i],data['fieldDec'][i],data['expMJD'][i],data['expDate'][i],data['filter'][i])
    #mysky=SkyBrightness(data['fieldRA'][i],data['fieldDec'][i],49355.199664,190051,data['filter'][i])
        filtre=data['filter'][i]
        transmission.Load_Atmosphere(data['airmass'][i])
        myup=transmission.darksky.calcInteg(transmission.lsst_system[filtre])
               
        Tb=Calc_Integ(transmission.lsst_atmos[filtre])
        Sigmab=Calc_Integ(transmission.lsst_system[filtre])
        katm=-2.5*np.log10(Tb/Sigmab)
    
        mbsky_through=-2.5*np.log10(myup/(3631.*Sigmab))

    #print 'yes here',filtre,mysky.new_skybrightness(),mbsky_through,mysky.moonPhase
       
        num+=1
        if len(tab_resu) <= num:
            tab_resu=np.resize(tab_resu,(len(tab_resu)+100,1))

        (yy, mm, dd,hh,mn,sec) = mysky.mjd2gre(data['expMJD'][i])[:6]
        (yy_sunset, mm_sunset, dd_sunset,hh_sunset,mn_sunset,sec_sunset) = mysky.mjd2gre(mysky.sunSetTwilMJD)[:6]
        (yy_sunrise, mm_sunrise, dd_sunrise,hh_sunrise,mn_sunrise,sec_sunrise) = mysky.mjd2gre(mysky.sunRiseTwilMJD)[:6]


        tab_resu['filter'][num]=filtre
        tab_resu['filtSkyBrightness'][num]=data['filtSkyBrightness'][i]
        tab_resu['mbsky_through'][num]=mbsky_through
        tab_resu['mbsky_through_moon'][num]=mysky.new_skybrightness()
        tab_resu['moonPhase'][num]=mysky.moonPhase
        tab_resu['moonDiam'][num]=mysky.moonDiam
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
        


    #if i > 100:
    #   break

    tab_resu=np.resize(tab_resu,(num+1,1))
    
    
    pkl_file_res = open('MoonLight.pkl','wb')
    pkl.dump(tab_resu, pkl_file_res)
    pkl_file_res.close()


pkl_file = open(filename,'rb')
tab_load=pkl.load(pkl_file)


sel=tab_load[np.where(tab_load['filter']=='g')]
#sel=sel[np.where(sel['hour']<1)]
#plt.plot(sel['hour_sunset']-sel['hour'],sel['filtSkyBrightness']-sel['mbsky_through'],'ko')
#plt.plot(sel['moonPhase'],sel['filtSkyBrightness']-sel['mbsky_through'],'ko')
plt.plot(sel['hour_sunset']+sel['min_sunset']/60.,sel['moonPhase'],'ko')

plt.show()

