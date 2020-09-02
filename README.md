# ncosmo
Series of Fortran and Python routines that allow calculating cosmological
corrections. The first lines of code were contributed by undergrad David Reiss
c. 1992 at UC Santa Cruz, on a project supervised by Prof. David C. Koo.
Major modifications were made following the publication in 1999 by 
David Hogg of "Distance Measures in Cosmology":
https://arxiv.org/abs/astro-ph/9905116). 

Additional changes follow Ned Wright's Cosmology calculator described in
Wright, E. L. 2006, PASP, 118, 1711

The python code was tested using the anaconda astroconda environment.

The FORTRAN code uses two routines that were not developed by me:

havie_integration.f, by G. Wallick, 1970, Comm. ACM,13, p. 622 (Algorithm 400)

zeroin.f - which is available from netlib, and is based on R. Brent's algol 60,
in "Algorithms for Minimization without Derivatives", 1973, Prentice-Hall.
This version used here was modified for use under Starlink.

To the best of my knowledge there are no associated royalties with the use and 
dissemination of these routines.

2020-09-01


