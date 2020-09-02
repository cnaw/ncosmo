#!/usr/bin/env python
""
#     Verify if cosmology parameters are set and select mode to calculate:
#     i.e., if reproducing the plots (and values) obtained using 
#     Hogg, D. W. 1999, astro-ph/990511v4 from 2000-12-16 (mode 1)
#     or from Ned Wright's Cosmology calculator (modes 2 and 3), described
#     in  Wright, E. L. 2006, PASP, 118, 1711
#
#     Definitions following M. Tanabashi et al. Phys Rev D 98 030001 (2018),
#     Table 2.1 
#     also available at http://pdglive.lbl.gov/Particle.action?node=S066&init=0
#     h         = Hubble parameter/100 
#     omega_b   = baryon density of the Universe
#     omega_cm  = cold dark matter density of the Universe + baryons
#     omega_m   = pressureless matter density of the Universe
#     omega_r   = CMB radiation density of the Universe
#     omega_l   = dark energy density of the Universe
#     omega_nu  = neutrino density of the Universe
#     omega_k   = curvature
#     w         = dark energy equation of state parameter
#     wprime    = derivative of dark energy equation of state parameter
#
#     options
#
#     1. lambda  = 0, no neutrinos       (e.g., Hogg 1999)
#     2. lambda != 0, no neutrinos       (e.g., Hogg 1999)
#     3. massless neutrinos (original calculator)
#     4. neutrino masses + equation of state
#
#     In option 1 omega_r = 0,  omega_l = 0 omega_cm = omega_m 
#     In option 2 omega_r = 0               omega_cm = omega_m 
#     In option 3 omega_r = 4.165d-05/h**2, omega_cm = omega_m
#     in option 4 omega_r = 2.473d-05/h**2, omega_m = omega_nu + omega_cm
#     where omega_m = omega_cm + omega_nu
""

import numpy as np
from  scipy.integrate import romberg
from  scipy.integrate import quad
from  scipy.integrate import trapz
from  scipy.optimize import brentq
#
global cee, omega_nu
cee          = 299792.458 # km/sec

#
#----------------------------------------------------------------------
#
def cosmological_model(model):
    global h, h0, omega_k, omega_l, omega_m, omega_nu, omega_r
    omega_nu = 0.0
    print("cosmological model is ", model)

    if(model == 'CfA1'):
        h       = 1.0
        w       = -1.0
        wprime  = 0.0
        omega_r = 0.0
        omega_l = 0.0
        omega_m = 0.05
        omega_k  = 1.0 - omega_m - omega_l

    if(model == 'Hogg_1999_Einstein_deSitter'):
        h       = 0.7
        w       = -1.0
        wprime  = 0.0
        omega_r = 0.0
        omega_l = 0.0
        omega_m = 1.0
        omega_k  = 1.0 - omega_m - omega_l

    if(model == 'Hogg_1999_low_density'):
        h       = 0.7
        w       = -1.0
        wprime  = 0.0
        omega_r = 0.0
        omega_l = 0.0
        omega_m = 0.05
        omega_k  = 1.0 - omega_m - omega_l

    if(model == 'Hogg_1999_high_lambda'):
        h       = 0.7
        omega_r = 0.0
        omega_m = 0.2
        omega_l = 0.8
            
    if(model == 'standard'):
        h       = 0.7
        omega_r = 0.0
        omega_m = 0.3
        omega_l = 0.7
            
    if(model == 'Wright_2006_massless_neutrinos'):
        h       = 0.696
        omega_r = 4.165e-05/h**2
        omega_m = 0.286
        omega_l = 1.0 - omega_m - omega_r

    if(model == 'Wright_2006'):
        h       = 0.696
        omega_r = 2.473e-05/h**2
        omega_m = omega_nu + omega_cm
        omega_l = 1.0 - omega_cm


# using parameters from Tanabashi et al. Phys Rev D 98 030001 (2018)
    if(model == 'Planck2015'):
        h        = 0.678
        tnow     = 2.7255 # K  
        omega_r  = 2.473e-05/h**2  # at present day
        omega_b  = 0.02226/h**2
        omega_cm = 0.1186/h**2
        omega_nu = 0.68/(93.14*h*h) # upper limit of sum of neutrino masses
        omega_m  = omega_nu + omega_cm + omega_b
        omega_l  = 0.692
        
    omega_k = 1 - omega_l - omega_m  - omega_r
    eps = 1.e-6
    if(abs(omega_k) < eps):
        omega_k = 0.0

    h0 = h * 100.0
    return h, h0, omega_k, omega_l, omega_m, omega_r

#
#-----------------------------------------------------------------------
#
#     Calculate the neutrino density over rest mass density
#     Translation of N. Wright's AAC calculator and 
#     Wright, E. L. 2006, PASP, 118, 1711, eqs. 16 (and 18).
#
def nurho (mnu):
    beta     = 1.842
    inv_beta = 1.0/beta
    y        = 1.0 + (mnurel/mnu)**beta
    value    = y**inv_beta
    return value  
#
#----------------------------------------------------------------------
#
def omega_neutrino(h, tnow, t0, z) :
#
#     Calculate the contribution to the matter density due to the
#     Cosmic Neutrino Background. This is a translation of
#     E. L. Wright's Advanced Cosmology Calculator at
#     http://www.astro.ucla.edu/~wright/ACC.html
#     and described in  Wright, E. L. 2006, PASP, 118, 1711
#
#     Mangano et al. (2005) :
#     omega_nu = rho_nu/rho_critical = 3 * m_nu/(93.14*h**2 eV) (eq. 18)
#     contribution of neutrinos to present energy density is n_eff:
#     rho_rad = [ 1 + (7/8) * (4/11)**(4/3) * n_eff] * rho_photons (eq. 1)
#     where rho_photons is the present energy density of photons coming
#     from CMB temperature
    m_e         = 0.001     # Neutrino mass in eV/cee**2
    m_mu        = 0.009     # eV/cee**2
    m_tau       = 0.049    # ev/cee**2
    one_ev      = 11604.5050

    m_nu        = np.zeros(3)
    m_nu[0]     = m_e
    m_nu[1]     = m_mu
    m_nu[2]     = m_tau
#
#     The  present day neutrino temperature relative to the photon 
#     temperature given by the CMB (Wright 2006, Section 5):
#
    nu_temp = t0 * (4.0/11.0)**0.3333333333330
#
#     Mass of a neutrino "that is just now relativistic"
#
    momentum = 3.151 * nu_temp 
    m_rel    = (momentum/one_ev) * (t0/tnow)
#
#     The values printed below match those in paper (Section 5)
#      print(' nu_temp, momentum, m_rel',nu_temp, momentum, m_rel)
#
#     While the paper quotes a value of 93.14 eV for the normalisation (eq.18),
#     the ACC itself uses different values for the electron and mu and  tau
#     neutrinos (93.64, 93.90, 93.90); 
#     The value 0f 93.14h**2 eV is used by Lesgourgues & Verde in Chapter 25 of 
#     Tanabashi et al. Phys Rev D 98 030001 (2018) 
#
    a    = 1.0/(1.0+z)
    omega_nu = 0.0
    for i in range(0,len(m_nu)):
        result = nurho(m_rel, m_nu[i]*a)
        omega_nu = omega_neutrino + result * m_nu(i)/93.14
        
    omega_nu = omega_nu *(t0/tnow)**3
    omega_nu = omega_nu/(h*h)
    return omega_nu
#
#----------------------------------------------------------------------
def arcsinh(z):
    arc_sin_h = np.log(z + np.sqrt(1+z**2))
    return arc_sin_h
#
#
#----------------------------------------------------------------------
#
#    Angular diameter distance (Hogg eq. 18)
#
def da(z):
    angular_diameter_distance = dm(z)/(1+z)
    return angular_diameter_distance
    
#
#----------------------------------------------------------------------
#
#    Line of sight comoving distance (Hogg eqs. 15 and 17)
def dc(z):
    dh = cee/(h*100.0)
    if(z == 0):
        comoving_distance = 0.0
        return comoving_distance
        
    zmin = 0.0
    zmax = z
    ezinverse = lambda z : 1/ez(z)
#    comoving_distance = dh * romberg(ezinverse, zmin, zmax)
# quad is about 4 x faster than romberg
    integral  =  quad(ezinverse, zmin, zmax)
    comoving_distance = dh * integral[0]
    return comoving_distance
#
#----------------------------------------------------------------------
#
# Calculate the distance modulus in Mpc
#
def distance_modulus(z):
    if (z == 0.00) :
        dmod = 0.0
    else:
        dmod = 5.00 * np.log10(dl(z)) + 25.0
#    print("distance_modulus z, dmod", z, dmod)
    return dmod
#
#----------------------------------------------------------------------
#
#    Luminosity distance (Hogg eq. 21)
#
def dl(z):
    luminosity_distance = dm(z) * (1+z)
    return luminosity_distance
    
#
#----------------------------------------------------------------------
#
#    transverse comoving distance (Hogg eq. 16)
#
def dm(z):

    dh  = cee/h0

    if(omega_k == 0):
        transverse_distance = dc(z)

    if(omega_k < 0) :
        x = np.sqrt(np.abs(omega_k)) * dc(z)/dh
        transverse_distance = dh * np.sin(x)/np.sqrt(abs(omega_k))

    if(omega_k > 0) :
        x = np.sqrt(omega_k) * dc(z)/ dh
        transverse_distance = abs((np.exp(-x) - np.exp(x))/2.0)
        transverse_distance = dh * transverse_distance/np.sqrt(omega_k)
        
    return transverse_distance
#
#----------------------------------------------------------------------
#
# Eq 14 in Hogg 1999 for the function ez(z) which is 
# "proportional to the time derivative of the logarithm of the scale factor"
# and used to calculate the line of sight comoving distance (Hogg 1999, eq. 15)
#
def ez(z):
    zplus1 = 1.0 + z
#    if(omega_nu == 0) :
    ezterm = omega_m * zplus1**3 + omega_k * zplus1**2 + omega_l + omega_r * zplus1**4
    ezterm = np.sqrt(ezterm)
    return ezterm
#
#----------------------------------------------------------------------
#
def flux_emitted(f_nu, wavelength, z):
#
#     f_nu in Jy, lambda in microns, where lambda is the 
#     characteristic wavelength of the filter. Emitted flux
#     in solar luminosities
#     Constants from 
#     http://www.nist.gov/pml/div684/fcdc/codata.cfm
#     published in 2010, Reviews of Modern Physics, 84, 1527
#     Available at:
#     http://www.nist.gov/pml/div684/fcdc/upload/Wall-chart-3.pdf
#     and 
#     http://pdg.lbl.gov/2014/reviews/rpp2014-rev-astrophysical-constants.pdf
#     http://pdg.lbl.gov/2014/reviews/contents_sports.html
#     Code needs verification 2019-10-01    
    lsol      = 3.828e33        # solar luminosity in erg/s
    r_mpc     = 3.08567758149e24 # Mpc --> cm
    cee_cm    = 2.99792458e10   # cm/s
    four_pi   = 4.00 * np.pi
    jy_erg    = 1.0e23          # Jy per erg/s
    micron_cm = 1.0e-04
    f_nu_lsol = f_nu/(jy_erg*lsol)
    dist_cm   = dl(z) * r_mpc
    constant  = four_pi * dist_cm * dist_cm
    constant  = constant * cee_cm/(wavelength*micron_cm)

    result    = f_nu_lsol * constant
    return result
#
#
#----------------------------------------------------------------------
#
def flux_observed(luminosity, wavelength, z):
#
#     find flux in Jy for a luminosity  in solar units, and 
#     characteristic [rest] wavelength of the filter is in microns.
#     Constants from 
#     http://www.nist.gov/pml/div684/fcdc/codata.cfm
#     published in 2010, Reviews of Modern Physics, 84, 1527
#     Available at:
#     http://www.nist.gov/pml/div684/fcdc/upload/Wall-chart-3.pdf
#     and 
#     http://pdg.lbl.gov/2014/reviews/rpp2014-rev-astrophysical-constants.pdf
#     http://pdg.lbl.gov/2014/reviews/contents_sports.html
#     needs verification 2019-10-01
#    
    lsol      = 3.828e33         # solar luminosity in erg/s
    r_mpc     = 3.08567758149e24 # Mpc --> cm
    cee_cm    = 2.99792458e10    # cm/s
    four_pi   = 4.0 * np.pi
    jy_erg    = 1.0e23           # Jansky per erg
    micron_cm = 1.0e-04          # 1 micron in cm
#
    lum_jy    = luminosity * lsol * jy_erg
    dist_cm   = dl(z) * r_mpc
    constant  = four_pi * dist_cm * dist_cm 
    constant  = wavelength*micron_cm/(cee_cm*constant)
    result    = lum_jy*constant
    return result
#
#----------------------------------------------------------------------
#
def get_mabs(apm, z):
    result = apm - distance_modulus(z)
    return result
#
#----------------------------------------------------------------------
#
#    Find projected distance in Mpc for an angle in radians
#    Based on original routine written by David Reiss, UCSC, 1992
#
def getnproj(angle, z):
    dist = angle * da(z)
    return dist
#
#----------------------------------------------------------------------
#
#    Find angle on the sky given redshift and projected distance in MPC
#
def getangle(dist, z):
    radians = dist/da(z)
    arcsec  = radians*180.0*3600.0/np.pi
    return radians, arcsec
#
#----------------------------------------------------------------------
#
#    Find velocity (km/s) as a function of redshift
#
def getvel(z):
    if(zz == 0.0):
        vel = 0.00
    else :
        vel = cee*((2.0*z + z**2.0)/(2.0 + 2.0*z + z**2.0))
    return vel
#
#----------------------------------------------------------------------
#
#    Find z as a function of velocity (km/s)
#
def getz(v):
    if(v == 0.0):
        zz = 0.00
    else :
        zz = ((1.0 + v/cee)/(1.0-v/cee))**0.50 - 1.0

    return zz
#
#----------------------------------------------------------------------
#
# Calculate the contribution to Omega due to dark energy
# (Wright, E. 2006 PASP 118, 1711, eq. 19)
#
def rho_de(z,w, wprime, omega_l):

    a = 1.0/(1 + z)
    term1 = -(3 + 3*w + 6* wprime)
    term2 = -6 * wprime(1-a)
    rd    = omega_l * a**term1*np.exp(term2)
    return rd
#
#----------------------------------------------------------------------
#
#   Age since formation
#
def tform(z):
    thubble = 1.0/h0
    zmax = 3000.0
    
    if(z > zmax):
        tformation = 0
        return tformation
    
    zmin = z
    tformation =  quad(zplus1ezi, zmin, zmax) # thubble 
    return tformation[0]
#
#----------------------------------------------------------------------
#
# Hubble time in gigayears
#
#
def th_in_gyr(thubble):
    r_mpc     = 3.08567758149e19 # km
    year      =  31556925.2 # tropical year in seconds for 2011
    gyr = r_mpc/h0/year/1.e9 * thubble
    return gyr
#
#----------------------------------------------------------------------
#
#    Lookback time (Hogg eq. 30), referencing 
#    Peebles, 1993, Principles of Physical Cosmology, Princeton Univ. Press,
#    pages 313-315
#
def tlookback(z):
    thubble = 1.0/h0
    if(z == 0):
        tl = 0.0
        return
        
    zmin = 0.0
    zmax = z
#    tl = romberg(zplus1ezi, zmin, zmax)
    tl = quad(zplus1ezi, zmin, zmax)
    
    tl = tl[0]  # * thubble
    return tl
#
#----------------------------------------------------------------------
#
def vol_elm(z, solid_angle):
#    
#     Comoving Volume element Hogg (2000) eq. 28, referenced to
#     Weinberg 1972, Gravitation and Cosmolgy: Principles and Applications
#     of the General Theory of Relativity, John Wiley & Sons, New York, pg 486
#     and 
#     Peebles, 1993, Principles of Physical Cosmology, Princeton Univ. Press,
#     pages 331-333

    if(z <= 0):
        volume_element = 0.0
        return volume_element

    dh   = cee/h0
    daz  = da(z)
    
    volume_element = dh * (1 + z)**2 * daz * daz * solid_angle/ez(z)
    return volume_element
#
#----------------------------------------------------------------------
#
def zdistmodn(dmod):
#
#     Obtain redshift from the distance modulus in magnitudes
#     using the bisection method of
#     R. Brent, 1973 in "Algorithms for Minimization without Derivatives", 
#     Prentice-Hall.
#
      zmin = 1.e-6
      zmax = 20.00
      distmod = lambda z: distance_modulus(z) - dmod
#      print('zdistmodn: zmin, zmax, dmod', zmin, zmax, dmod)
      zz = brentq(distmod, zmin, zmax, xtol=1.e-12, rtol=1.e-11, maxiter=100, full_output=False, disp=True)
#      print("zdistmodn: zz", zz, dmod)
      return zz 
#
#----------------------------------------------------------------------
#
#    Find z from r
#
# Given a proper distance in Mpc, find the corresponding redshift
def z_from_r(r):
    zmin =  0.0
    zmax = 16.0
    zfromr = lambda z : dl(z) - r
    zz = brentq(zfromr, zmin, zmax, xtol=2e-12, rtol=8.881784197001252e-16, maxiter=100, full_output=False, disp=True)
    return zz
#
#
#----------------------------------------------------------------------
#
#    Find z from formation time
#
# Given a formation time, find the corresponding redshift
def z_from_tform(t_in_hubbles):
    zmin =  0.0
    zmax = 16.0
    zfromt = lambda z : tform(z) - time
    zz = brentq(zfromt, zmin, zmax, xtol=1.0e-6, rtol=1.0e-6, maxiter=100, full_output=False, disp=True)
    return zz
#
#----------------------------------------------------------------------
#
def zplus1ezi(z):
    zpezi = 1.0/(ez(z)*(1+z)) 
    return zpezi
#
#----------------------------------------------------------------------
#
def zvolume(solid_angle, zmin, zmax):
    
    if(zmax == 0 or zmin == zmax):
        volume = 0.0
        return volume
        
    dh = cee/h0

# lower limit of integration
    
    if(zmin <= 1.e-12):
        dvol0 = 0.0
    else :
        rm    = dm(zmin)
        if(omega_k == 0):
            dvol0 = (solid_angle/3.0) * rm**3
        else:
            rm_over_dh = rm/dh
            term1 = rm_over_dh*np.sqrt(1.0 + omega_k * rm_over_dh**2)
            term3 = rm_over_dh*np.sqrt(abs(omega_k))

            if(omega_k > 0):
                term2 = (1.0/np.sqrt(abs(omega_k))) * arcsinh(term3)
            else:
                term2 = (1.0/np.sqrt(abs(omega_k))) * np.arcsin(term3)
            
            dvol0 = (solid_angle*dh**3)/(2.00*omega_k) * (term1-term2)

# upper limit
    rm  = dm (zmax)

    if(omega_k ==  0.0):
        volume = (solid_angle/3.0) * rm**3 - dvol0
        return volume

    rm_over_dh = rm/dh
    term1 = rm_over_dh*np.sqrt(1.0 + omega_k * rm_over_dh**2)
    term3 = rm_over_dh*np.sqrt(abs(omega_k))

    if(omega_k > 0):
        term2 = (1.0/np.sqrt(abs(omega_k))) * arcsinh(term3)
    else:
        term2 = (1.0/np.sqrt(abs(omega_k))) * np.arcsin(term3)
        
    volume = ((solid_angle*dh**3)/(2.0*omega_k)) *(term1-term2)
    volume = volume - dvol0
    return volume
