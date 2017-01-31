import numpy as np
import matplotlib.pyplot as plt

import SiteModel as sm

#--------------------------------------------------------------
# Manually creating a site model using dictionaries

site = sm.Site1D(Id=1, X=10, Y=20)

site.AddLayer({'Hl': 10.,'Vp': 300.,'Vs': 200.,'Dn': 1900.,'Qp': 50.,'Qs': 20.})
site.AddLayer({'Hl': 10.,'Vp': 500.,'Vs': 300.,'Dn': 1900.,'Qp': 50.,'Qs': 20.})
site.AddLayer({'Hl': 0.,'Vp': 1000.,'Vs': 800.,'Dn': 2100.,'Qp': 100.,'Qs': 50.})

#--------------------------------------------------------------
# Alternative method using lists (fixed format)

site = sm.Site1D()

site.AddLayer([10., 300., 200., 1900., 50., 20.])
site.AddLayer([10., 500., 300., 1900., 50., 20.])
site.AddLayer([0., 1000., 800., 2100., 100., 50.])

#--------------------------------------------------------------
# Parsing standard CSV data format

site = sm.Site1D()

site.ImportModel('data/site01.csv')

#--------------------------------------------------------------
# Parsing Geopsy data format

site = sm.Site1D()

site.ImportModel('data/site01.mod',
                  Header=['Hl','Vp','Vs','Dn','Qp','Qs'],
                  Delimiter=' ',
                  Skipline=1)

#--------------------------------------------------------------
# Add sites to a regional container

zone = sm.SiteDb(Id=1, Name='Lodi')

zone.AddSite(site)
zone.AddSite(site)
zone.AddSite(site)

#--------------------------------------------------------------
# Computing Travel-Time average velocity

# Direct call (using default values)
site.ComputeTTAV()
print 'Vs30:', site.Eng['Vz']['30.0'], 'm/s'

# Arbitrary values
site.ComputeTTAV('Vp',50.0)
print 'Vp50:', site.Eng['Vz']['50.0'], 'm/s'

#--------------------------------------------------------------
# Compute Geotechnical soil class

site.ComputeGTClass('EC8')

print 'EC8 Site Class:', site.Eng['EC8']

#--------------------------------------------------------------
# Frequency axis

site.GenFreqAx(1 ,100. ,1000, Log=True)
Freq = site.Amp['Freq']

#--------------------------------------------------------------
# Computing Quarter-Wavelength parameters

site.ComputeQWL()
QA = site.Amp['Qwl']

#--------------------------------------------------------------
# Computing SH transfer-function

# Anelastic Amplification
site.ComputeSHTF()
AA = site.Amp['Shtf']

# Elastic Amplification
site.ComputeSHTF(Elastic=True)
EA = site.Amp['Shtf']

# Elastic Amplification with variable angle of incidence
site.ComputeSHTF(Iang=45., Elastic=True)
IA = site.Amp['Shtf']

#--------------------------------------------------------------
# Identify resonance frequencies and their amplitude

site.ComputeFnRes()

print 'F0 is:', site.Amp['Fn'][0], 'Hz'
print 'A0 is:', site.Amp['An'][0]

#--------------------------------------------------------------
# Compute Kappa0

# Using whole profile (down to last layer)
site.ComputeKappa0()
print 'Kappa0 is:', site.Eng['K0'], 's'

# Using arbitrary depth
site.ComputeKappa0(Z=100.)
print 'Kappa0 is:', site.Eng['K0'], 's'

# Spectral attenuation function
site.ComputeAttFun()
AT = site.Amp['Attf']

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
plt.savefig('AmpFunc.png',dpi=300)

#--------------------------------------------------------------
# Plot (2)

f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

ax1.semilogx(site.Amp['Freq'],site.Eng['Qwl']['Hl'],'r-o')
ax1.set_xlabel('Frequency (Hz)')
ax1.set_ylabel('Thickness (m)')
ax1.axis([1, 100, 0, 160])
ax1.grid('on')

ax2.semilogx(site.Amp['Freq'],site.Eng['Qwl']['Vs'],'b-o')
ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Average Velocity (m/s)')
ax2.axis([1, 100, 100, 650])
ax2.grid('on')

plt.show(block=False)
plt.savefig('Qwl.png',dpi=300)