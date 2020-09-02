#!/usr/bin/env python
from ncosmo3 import *

import matplotlib.pyplot as plt

#global cee
# global cee, h, h0, omega_l, omega_m, omega_nu, omega_r, omega_k
omega_nu = 0.

eps = 1.e-6


#model = 'Wright_2006_massless_neutrinos'
model = 'Planck2015'

h, h0, omega_k, omega_l, omega_m, omega_r = cosmological_model(model)
print(h)



print("h, h0, omega_l, omega_m, omega_r",h, h0, omega_l, omega_m, omega_r)

z  = 1.0
print(dc(z))

print(dl(z))


print(da(z))

solid_angle = np.pi*4.0
zmin   = eps
zmax   =  z
volume = zvolume(solid_angle, zmin, zmax)
print(volume/1e09)

print('lookback time ', tlookback(z), th_in_gyr(tlookback(z)))
print('age           ',    tform(z), th_in_gyr(tform(z)))

print('age of the Universe ', tform(zmin), th_in_gyr(tform(zmin)))

dist = distance_modulus(z)
print("distance modulus ", dist)

zz = zdistmodn(dist)
print ("recovered redshift ", zz)

radius = 1.00
angle_rad, angle_arcsec  = getangle(radius,z)
print("angle in radians", angle_rad," in arc sec ", angle_arcsec)

arcsec = np.pi/(3600.0*180.)
radius = getnproj(arcsec, z)
radius_kpc = radius * 1.e3
print("projected distance corresponding to 1 arc sec ",radius, radius_kpc)
angle_rad, angle_arcsec  = getangle(radius,z)
print("angle in radians", angle_rad," in arc sec ", angle_arcsec)

lsol = 1.0e9
wavelength = 2.0 ;
observed   = flux_observed(lsol, wavelength, z)
print("flux_observed", observed, "Jy")

emitted = flux_emitted(observed, wavelength, z)
print("flux_emitted", emitted," Lsol")


exit(0)

dh = cee/h0
npts = 200

zz   = np.zeros(npts)
ddm1 = np.zeros(npts)
ddm2 = np.zeros(npts)
ddm3 = np.zeros(npts)
ddm4 = np.zeros(npts)

model = 'Hogg_1999_high_lambda'
cosmological_model(model)
for ii in range(0,npts):
    zz[ii] = (ii+1) * 5.0/200.
    ddm1[ii] = dm(zz[ii])/dh
plt.xlim(0., 5.0)
plt.ylim(0., 3.0)
plt.plot(zz, ddm1,'g')


model = 'Hogg_1999_low_density'
cosmological_model(model)

for ii in range(0,npts):
    ddm2[ii] = dm(zz[ii])/dh

plt.plot(zz, ddm2,'r')


model = 'Hogg_1999_Einstein_deSitter'
cosmological_model(model)

for ii in range(0,npts):
    ddm3[ii] = dm(zz[ii])/dh

plt.plot(zz, ddm3,'b')


model = 'Planck2015'
cosmological_model(model)
for ii in range(0,npts):
    ddm4[ii] = dm(zz[ii])/dh

plt.plot(zz, ddm4,'y')
plt.xlabel('redshift')
plt.ylabel('proper motion distance DM(z)/DH')
plt.show()

#
#------------------------------------------------
#

model = 'Hogg_1999_high_lambda'
cosmological_model(model)
for ii in range(0,npts):
    ddm1[ii] = da(zz[ii])/dh
plt.xlim(0., 5.0)
plt.ylim(0., 0.5)
plt.xlabel('redshift')
plt.ylabel('angular diameter distance DA(z)/DH')


plt.plot(zz, ddm1,'g')


model = 'Hogg_1999_low_density'
cosmological_model(model)

for ii in range(0,npts):
    ddm2[ii] = da(zz[ii])/dh

plt.plot(zz, ddm2,'r')


model = 'Hogg_1999_Einstein_deSitter'
cosmological_model(model)

for ii in range(0,npts):
    ddm3[ii] = da(zz[ii])/dh

plt.plot(zz, ddm3,'b')


model = 'Planck2015'
cosmological_model(model)
for ii in range(0,npts):
    ddm4[ii] = da(zz[ii])/dh

plt.plot(zz, ddm4,'y')
plt.show()

#
#------------------------------------------------
#
model = 'Hogg_1999_high_lambda'
cosmological_model(model)

for ii in range(0,npts):
    ddm1[ii] = dl(zz[ii])/dh
plt.xlim(0., 5.0)
plt.ylim(0., 17.)
plt.xlabel('redshift')
plt.ylabel('Luminosity distance DL(z)/DH')


plt.plot(zz, ddm1,'g')


model = 'Hogg_1999_low_density'
cosmological_model(model)

for ii in range(0,npts):
    ddm2[ii] = dl(zz[ii])/dh

plt.plot(zz, ddm2,'r')


model = 'Hogg_1999_Einstein_deSitter'
cosmological_model(model)

for ii in range(0,npts):
    ddm3[ii] = dl(zz[ii])/dh

plt.plot(zz, ddm3,'b')


model = 'Planck2015'
cosmological_model(model)
for ii in range(0,npts):
    ddm4[ii] = dl(zz[ii])/dh

plt.plot(zz, ddm4,'y')
plt.show()
#
#------------------------------------------------
#
model = 'Hogg_1999_high_lambda'
cosmological_model(model)
solid_angle = 1.0

for ii in range(0,npts):
    ddm1[ii] = vol_elm(zz[ii], solid_angle)/dh**3
plt.xlim(0., 5.0)
plt.ylim(0., 1.2)
plt.xlabel('redshift')
plt.ylabel('volume element DL(z)/DH')


plt.plot(zz, ddm1,'g')


model = 'Hogg_1999_low_density'
cosmological_model(model)

for ii in range(0,npts):
    ddm2[ii] = vol_elm(zz[ii], solid_angle)/dh**3

plt.plot(zz, ddm2,'r')


model = 'Hogg_1999_Einstein_deSitter'
cosmological_model(model)

for ii in range(0,npts):
    ddm3[ii] = vol_elm(zz[ii], solid_angle)/dh**3

plt.plot(zz, ddm3,'b')


model = 'Planck2015'
cosmological_model(model)
for ii in range(0,npts):
    ddm4[ii] = vol_elm(zz[ii], solid_angle)/dh**3

plt.plot(zz, ddm4,'y')
plt.show()


#
#------------------------------------------------
#
model = 'Hogg_1999_high_lambda'
cosmological_model(model)

for ii in range(0,npts):
    ddm1[ii] = tlookback(zz[ii])
plt.xlim(0., 5.0)
plt.ylim(0., 1.2)
plt.xlabel('redshift')
plt.ylabel('Lookback time T/tH')


plt.plot(zz, ddm1,'g')


model = 'Hogg_1999_low_density'
cosmological_model(model)

for ii in range(0,npts):
    ddm2[ii] =tlookback(zz[ii])

plt.plot(zz, ddm2,'r')


model = 'Hogg_1999_Einstein_deSitter'
cosmological_model(model)

for ii in range(0,npts):
    ddm3[ii] =tlookback(zz[ii])

plt.plot(zz, ddm3,'b')


model = 'Planck2015'
cosmological_model(model)
for ii in range(0,npts):
    ddm4[ii] =tlookback(zz[ii])

plt.plot(zz, ddm4,'y')
plt.show()


#
#------------------------------------------------
#
model = 'Hogg_1999_high_lambda'
cosmological_model(model)

for ii in range(0,npts):
    ddm1[ii] = tform(zz[ii])
plt.xlim(0., 5.0)
plt.ylim(0., 1.2)
plt.xlabel('redshift')
plt.ylabel('Age T/tH')


plt.plot(zz, ddm1,'g')


model = 'Hogg_1999_low_density'
cosmological_model(model)

for ii in range(0,npts):
    ddm2[ii] =tform(zz[ii])

plt.plot(zz, ddm2,'r')


model = 'Hogg_1999_Einstein_deSitter'
cosmological_model(model)

for ii in range(0,npts):
    ddm3[ii] =tform(zz[ii])

plt.plot(zz, ddm3,'b')


model = 'Planck2015'
cosmological_model(model)
for ii in range(0,npts):
    ddm4[ii] =tform(zz[ii])

plt.plot(zz, ddm4,'y')
plt.show()


