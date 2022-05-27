      real bad               ! bad/missing data flag in ascii files
      real sbad              ! bad/missing data flag in sweep files
      parameter (bad=9.9E9, sbad=-32768.)
      integer maxflds        ! maximum no. of fields in sweep file
      integer maxrays        ! maximum no. of rays in sweep file

c
c.... CLZ (1/9/09): increase maxflds to 16 to handle xpol on 6/19/02 IHOP
c
c.... CLZ (11/2/10): increase maxflds to 18 to handle noxp in VORTEX2
c
c.... CLZ (11/4/10): increase maxrays to 2300 to handle dow6 in VORTEX2
c
c
c.... CLZ (1/18/10): increase maxflds to 20 to handle noxp in VORTEX2
c
c
c.... CLZ (6/15/17): increase maxflds to 30 to handle DOWs in PECAN
c
c.... CLZ (11/5/19): increase maxflds to 36 to handle new Pyart edit fields in TORUS and VORTEX-SE
c
c.... DS (2/13/20): increase maxflds to 50 to account for potential increase in fields from CIMMS-RadarQC
c

      parameter (maxflds = 50, maxrays = 4000)

      real pi                ! pi
      parameter (pi=3.14159265)
      real dtor              ! conversion from degrees to radians
      real rtod              ! conversion from radians to degrees
      parameter (dtor=0.0174533, rtod=57.29578)
      integer pd(12)         ! cumulative previous days (for julian day computation)
      data pd /0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/
      real rearth            ! actual radius of earth (km)
      real eer               ! effective earth radius (km)
      parameter (rearth=6367.0)
      parameter (eer=rearth*4.0/3.0)
c      parameter (eer=rearth)
      real omega             ! rotation rate of earth (rad/s)
      parameter (omega=0.0000729)