import numpy as np
import matplotlib.pyplot as plt

import OQSrtk.SiteModel as sm

#--------------------------------------------------------------
# Manually creating a site model using dictionaries

site = sm.Site1D(Id=1, X=10., Y=20., Z=0.)

site.AddModel()
site.Mod[0].AddLayer({'Hl': 10.,'Vp': 300.,'Vs': 200.,'Dn': 1900.,'Qp': 50.,'Qs': 20.})
site.Mod[0].AddLayer({'Hl': 10.,'Vp': 500.,'Vs': 300.,'Dn': 1900.,'Qp': 50.,'Qs': 20.})
site.Mod[0].AddLayer({'Hl': 0.,'Vp': 1000.,'Vs': 800.,'Dn': 2100.,'Qp': 100.,'Qs': 50.})

#--------------------------------------------------------------
# Alternative method using lists (fixed format)

site = sm.Site1D()

site.AddModel()
site.Mod[0].AddLayer([10., 300., 200., 1900., 50., 20.])
site.Mod[0].AddLayer([10., 500., 300., 1900., 50., 20.])
site.Mod[0].AddLayer([0., 1000., 800., 2100., 100., 50.])

#--------------------------------------------------------------
# Parsing standard CSV data format

site = sm.Site1D()

site.ImportModel('Data/site01.csv')

#--------------------------------------------------------------
# Parsing Geopsy data format

site = sm.Site1D()

site.ImportModel('Data/site01.mod', FileType='Geopsy')

#--------------------------------------------------------------
# Parsing user-defined data format

site = sm.Site1D()

site.ImportModel('Data/site01.mod',
                  Header=['Hl','Vp','Vs','Dn','Qp','Qs'],
                  Delimiter=' ',
                  SkipLine=1)

#--------------------------------------------------------------
# Add sites to a regional container

zone = sm.SiteDb(Id=1, Info='Ferrara')

zone.AddSite(Site=site)
zone.AddSite(Site=site)
zone.AddSite(Site=site)

#--------------------------------------------------------------
# Computing Travel-Time average velocity

# Direct call (using default values)
site.ComputeTTAV()
print 'Vs30:', site.Mod[0].Eng['Vz'][30.], 'm/s'

# Arbitrary values
site.ComputeTTAV('Vp',50.0)
print 'Vp50:', site.Mod[0].Eng['Vz'][50.], 'm/s'

#--------------------------------------------------------------
# Compute Geotechnical soil class

site.ComputeGTClass('EC8')

print 'EC8 Site Class:', site.Mod[0].Eng['Gc']

#--------------------------------------------------------------
# Frequency axis

site.FrequencyAxis(1. ,100. ,1000, Log=True)
Freq = site.Freq

#--------------------------------------------------------------
# Computing Quarter-Wavelength parameters

site.ComputeQWL()

#--------------------------------------------------------------
# Computing Impedance Amplification (from Qwl)

site.ComputeImpAmp()
QA = site.Mod[0].Amp['Imp']

#--------------------------------------------------------------
# Computing SH transfer-function

# Anelastic Amplification
site.ComputeSHTF()
AA = site.Mod[0].Amp['Stf']

# Elastic Amplification
site.ComputeSHTF(Elastic=True)
EA = site.Mod[0].Amp['Stf']

# Elastic Amplification with variable angle of incidence
site.ComputeSHTF(Iang=45., Elastic=True)
IA = site.Mod[0].Amp['Stf']

#--------------------------------------------------------------
# Identify resonance frequencies and their amplitude

site.ComputeFnRes()

print 'F0 is:', site.Mod[0].Amp['Res']['Fn'][0], 'Hz'
print 'A0 is:', site.Mod[0].Amp['Res']['An'][0]

#--------------------------------------------------------------
# Compute Kappa0

# Using whole profile (down to last layer)
site.ComputeKappa()
print 'Kappa0 is:', site.Mod[0].Eng['K0'], 's'

# Using arbitrary depth
site.ComputeKappa(Z=100.)
print 'Kappa0 is:', site.Mod[0].Eng['K0'], 's'

# Spectral attenuation function
site.ComputeAttFun()
AT = site.Mod[0].Amp['Att']

#--------------------------------------------------------------
# Plot

plt.figure()

plt.semilogx(Freq,np.abs(AA),'k',label='SHTF Anelastic')
plt.semilogx(Freq,np.abs(EA),'r',label='SHTF Elastic')
plt.semilogx(Freq,np.abs(IA),'g',label='SHTF (45 Deg. Incidence)')
plt.semilogx(Freq,QA,'b',label='QWL Amplification')
plt.semilogx(Freq,AT,'r--',label='Kappa Attenuation')
plt.show(block=False)

plt.legend(loc='upper right')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplification')

plt.show(block=False)
plt.savefig('Pictures/AmpFunc.png',dpi=300)

#--------------------------------------------------------------
# Plot (2)

f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

ax1.semilogx(site.Freq,site.Mod[0].Eng['Qwl']['Hl'],'r-')
ax1.set_xlabel('Frequency (Hz)')
ax1.set_ylabel('Thickness (m)')
ax1.axis([1, 100, 0, 160])
ax1.grid('on')

ax2.semilogx(site.Freq,site.Mod[0].Eng['Qwl']['Vs'],'b-')
ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Average Velocity (m/s)')
ax2.axis([1, 100, 100, 650])
ax2.grid('on')

plt.show(block=False)
plt.savefig('Pictures/QwlPar.png',dpi=300)