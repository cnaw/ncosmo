      subroutine read_constants
      implicit none
      double precision h
      double precision h0, m_nu, mpc_in_m, tnow, w, year,
     &     omega_b, omega_cm, omega_l, omega_m, omega_nu, omega_r,
     &     omega_k, wprime
      common /cosmology/ h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
      common /constants/ year, tnow
      open(1,file='constants.dat')
      read(1,*)
      read(1,*)
      read(1,*) omega_b
      read(1,*) omega_cm
      read(1,*) omega_m
      read(1,*) omega_l
      read(1,*) m_nu
      read(1,*) h0
      read(1,*) w
      read(1,*) Tnow
      read(1,*) omega_r
      read(1,*) omega_nu
      read(1,*) year
      read(1,*) mpc_in_m
      close(1)
      h = h0/100.d0
      omega_r = omega_r /(h*h)
      return
      end
