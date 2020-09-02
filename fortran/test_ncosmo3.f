c
c     Example on how to use the code to calculate cosmological corrections.
c     The lines enclosed by  "c*****"  must be included in the 
c     main programme. The latter must minimally include values for 
c     H0, omega_m and omega_l. 
c
c     The code calculates several paramaters for different values of 
c     H0, Omega_m, Omega_nu, Omega_r, Omega_lambda, W, and W' and can
c     be compared to David Hogg's 1999 astro-ph/990511v4 paper
c     and Ned Wright's Cosmology calculator at
c     http://www.astro.ucla.edu/~wright/CosmoCalc.html
c     http://www.astro.ucla.edu/~wright/ACC.html
c     
c     Constants used in the routines adopt values published in 2010 in
c     published in 2010, Reviews of Modern Physics, 84, 1527
c     and are available at 
c     http://www.nist.gov/pml/div684/fcdc/codata.cfm
c     http://www.nist.gov/pml/div684/fcdc/upload/Wall-chart-3.pdf
c     Other physical/astronomical constants come from
c     Review of Particle Physics 2014 edition: 
c     Olive, K. A. et al. 2014 Chin. Phys. C, 38, 090001
c     and are available at
c     http://pdg.lbl.gov/2014/reviews/rpp2014-rev-astrophysical-constants.pdf
c     http://pdg.lbl.gov/2014/reviews/contents_sports.html
c
c     In addition to the cosmology routines, two functions to calculate 
c     the numerical integrations and root-finding are included. To the best of
c     my knowledge there are no associated royalties with the use and 
c     dissemination of these routines.
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-03-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision az, a, scale, tf, tl
      double precision z, result, solid_angle, h, t0, tnow, th, age
      double precision hh, omm, oll, orr, onn, oww, opp
      double precision th_in_gyr,
     *     omega_neutrino, rho_de
      double precision dc, dm, da, dl, tform, tlookback, zvolumen
      integer i
c
      dimension hh(6), omm(6), oll(6), orr(6),onn(6), oww(6), opp(6)
c*****
      double precision  h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
      common /cosmology/ h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
c*****
c
      tnow    = 2.72528d0       ! fiducial used by Wright 2006
      tnow    = 2.7255d0        !± 0.0006 K (Fixsen 2009, ApJ, 707, 916)
c
      t0      = tnow
c
c     Set a solid angle value
      solid_angle = 4.d0 * dacos(-1.0d0)
c
c     Values used by Hogg
c     Einstein-de Sitter
      hh(1)  = 100.0d0          ! Hubble Constant (km/s/Mpc)
c     hh(1)  =  69.6d0
      omm(1) =   1.0d0          ! Omega_matter
      oll(1) =   0.0d0          ! Omega_lambda
      orr(1) =   0.0d0          ! Omega_radiation
      onn(1) =   0.0d0          ! Omega_neutrino
      oww(1) =  -1.0d0          ! w  in the dark energy equation of state
      opp(1) =   0.0d0          ! w' in the dark energy equation of state
c
c     No dark matter
c
      hh(2)  = 100.0d0          
c     hh(2)  =  69.6d0           
      omm(2) =   0.05d0         ! Omega_matter                           
      oll(2) =   0.0d0          ! Omega_lambda                         
      orr(2) =   0.0d0          ! Omega_radiation                        
      onn(2) =   0.0d0          ! Omega_neutrino                         
      oww(2) =  -1.0d0          ! w  in the dark energy equation of state
      opp(2) =   0.0d0          ! w' in the dark energy equation of state
c
c     Omega_m = 0.2, Omega_l=0.8
      hh(3)  = 100.0d0          ! Value used by D. Hogg in his paper
      omm(3) =   0.2d0          ! idem
      oll(3) =   0.8d0          ! ibdem
      hh(3)  =  70.0d0          
      omm(3) =   0.3d0          ! Omega_matter                           
      oll(3) =   0.7d0          ! Omega_lambda                           
      orr(3) =   0.0d0          ! Omega_radiation                        
      onn(3) =   0.0d0          ! Omega_neutrino                         
      oww(3) =  -1.0d0          ! w  in the dark energy equation of state
      opp(3) =   0.0d0          ! w' in the dark energy equation of state
c
c     Ned Wright's calculator uses these values in the case of
c     massless neutrinos
c
      hh(4)  =  69.6d0
      omm(4) =   0.286d0
      oll(4) =   0.714d0
c      hh(4)  =  70.0d0
c      omm(4) =   0.3d0
c      oll(4) =   0.7d0
      h      =  hh(4)/100.d0
      orr(4) =   4.165d-05/(h*h)
      onn(4) =   0.0d0
      oww(4) =  -1.0d0
      opp(4) =   0.0d0
c     
c     or of Neutrinos with masses
c
      hh(5)  =  69.6d0
      omm(5) =   0.286d0
      oll(5) =   0.714d0
      h      =  hh(5)/100.d0
      orr(5) =   2.477d-05/(h*h)*(t0/tnow)**4 
      z      = 0.d0
      h0     = hh(5) 
      onn(5) =   omega_neutrino(z)
      oww(5) =  -1.0d0
      opp(5) =   0.d0
c
c     Planck 2015
c
      hh(6)  =  67.8d0
      omm(6) =   0.308d0
      oll(6) =   0.714d0
      h      =  hh(6)/100.d0
      orr(6) =  2.473d-05/(h*h)*(t0/tnow)**4 
      z      =  0.d0
      h0     =  hh(6) 
      onn(6) =  omega_neutrino(z)
      oww(6) =  -1.006d0
      opp(6) =   0.d0
c
      z = 3.0d0
      do i = 1, 6
         h0       = hh(i)
         omega_m  = omm(i)
         omega_l  = oll(i)
         omega_r  = orr(i)
         omega_nu = onn(i)
         w        = oww(i)
         wprime   = opp(i)
         th       = 1.d0/h0
         print 110, h0, omega_m, omega_l, omega_r, omega_nu, omega_k,
     *        w, wprime
 110     format('H0 ',f6.2,' M, L, R, Nu, K, W Wprime:', 5(1x, f8.5),
     *        2(1x,f4.1))
         print 120
 120     format('parameter  redshift                     value')
c     Calculate the age since formation at redshift z
         tf       = tform(z)
c     Calculate the look-back time from the present to z
         tl       = tlookback(z)
c     Calculate the present age of the Universe
         age      = tf+tl
         print *,'Age ', z, age,  th_in_gyr(age/th), ' Th, Gyr'
         print *,'Tf  ', z, tf, th_in_gyr(tf/th), ' Th, Gyr'
         print *,'TL  ', z, tl, th_in_gyr(tl/th), ' Th, Gyr'
c     Calculate the co-moving radial distance
         result   = dc(z)
         print *,'DC  ', z, result,' Mpc'
c     Calculate the co-moving volume between zmin (=0 in this case) and z
         result = zvolumen(solid_angle,0.0d0,z)
         print *,'Vol ', z, result/1.d09,' Gpc**3'
c     Calculate the angular size distance
         result   = da(z)
         print *,'DA  ', z, result,' Mpc'
c     Calculate the projected scale on the sky
         scale    = result*1000.d0/(3600.d0*180.d0/dacos(-1.0d0))
         print *,'scale', z, scale ,' kpc/arc sec'
c     Calculate the luminosity distance
         result   = dl(z)
         print *,'DL  ', z, result,' Mpc'
c     Calculate the proper motion distance
         result   = dm(z)
         print *,'DM ', z, result,' Mpc'
      end do
      stop
      end
