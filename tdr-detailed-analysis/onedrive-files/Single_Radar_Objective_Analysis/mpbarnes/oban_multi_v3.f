      program oban_multi_v3

c############################################################################
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                     PROGRAM OBAN                     ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
c
c############################################################################
c
c     PURPOSE:
c
c     This program interpolates data from radar coordinates to an
c     unstaggered Cartesian coordinate system.  The input files should
c     be in the "sweep file" format (related to dorade format) that
c     is used by the solo radar data editing package.  The output can
c     be standard Cartesian or PPI (sweep-by-sweep).
c
c     inputs:
c     (1) parameter file
c     (2) radar data (dorade sweep files)
c
c     outputs:
c     (1) oban.output
c     (2) objectively analyzed data in text format
c     (3) "                          " vis5d format
c     (4) "                          " netcdf format
c
c############################################################################
c
c     Author of Original 1-pass Oban code: David Dowell
c
c     December 2001: Original 1-pass Oban code creation date
c
c     August 2005: Last DD modifications
c
c     Jan-April 2007: Modified for multipass Barnes by Penn State (see below)
c
c     May 2007: Multipass Barnes converted to PowerPC Mac G5 by Conrad ZIegler
c
c     Jan 2009: Multipass Barnes converted to Intel Mac Pro by Conrad Ziegler
c
c     April 2016: Multi-analysis vesion for Intel Mac developed by Conrad Ziegler
c
c     September 2018: Conrad Ziegler adds airborne radar analysis option
c
c############################################################################
c

c
c  Run this job in background with the following command line syntax:
c
c        x.oban_multi_v3 oban.input >&! oban.output &
c
c   where oban.input is the input namelist and oban.output is the printer output.
c
c

c      use ifport
      implicit none


      integer system
      integer ierrno
      integer(4) iflag
      integer(4) errnum

      include 'dow.inc'

      integer iradfilt
      integer isigfilt
      integer ifilter
      integer nuvfpass
      integer nwfpass
      real dbzmiss
      real dbzmiss0
      real amsg
      integer ivsmooth
      integer ismooth
      integer nuvsmooth
      integer nwsmooth
      
      integer nanal               ! nanal = number of analyses to perform
      integer noban               ! analysis step counter

      integer nx, ny, nz          ! no. of grid points in x, y, and z directions
      real xmin, ymin, zmin       ! coordinates of lower southwest corner
                                  !   of grid, relative to origin (km)
      real dx, dy, dz             ! grid spacing in x, y, and z directions (km)
      real rlat, rlon             ! radar latitude and longitude (deg)
      real ralt                   ! radar altitude (km MSL)
      real glat, glon             ! latitude and longitude of grid origin (deg)
      real galt                   ! altitude of grid origin (km MSL)
      integer map_proj            ! map projection (for relating lat, lon to x, y):
                                  !   0 = flat earth
                                  !   1 = oblique azimuthal
                                  !   2 = Lambert conformal
      integer nfld                ! number of data fields to be gridded
      integer ngrdtot             ! total number of gridded fields including data and factors
      integer(kind=2) cyr,cmo,cda ! central date
      integer(kind=2) chr,cmn,cse ! central time
      real ut, vt                 ! translation velocity (m/s)
      integer yrcor               ! correction to year
      integer mocor               ! correction to month
      integer dacor               ! correction to day
      integer minswpdc            ! begin sweep for dacor
      integer maxswpdc            ! end sweep for dacor
      integer umass_flag          ! flag for UMass radar data (1=yes, 0=no)
      integer az_corr_flag        ! method of additional azimuth-angle correction
                                  !   0 = none
                                  !   1 = average current and previous angle
      real elcor                  ! elevation angle correction (deg)
      real azcor                  ! azimuth-angle offset (deg)
      integer maxswp              ! dimension of n-vector sweep arrays.  Maximum number of sweeps.
      integer nswp                ! number of sweep files
      integer output_format       ! format for output:
                                  !   1=Cartesian
                                  !   2=PPI (sweep-by-sweep)
                                  !   3=DART (also sweep-by-sweep)
      integer iradtype            ! radar type iradtype = 0 (ground-based), = 1 (airborne)
      integer method              ! interpolation method:
                                  !   1=Cressman
                                  !   2=Barnes
                                  !   3=linear least squares
      real hsp, vsp               ! horizontal and vertical smoothing parameters:
                                  !   if method=1, then hsp and vsp are the Cressman radii of influence (km)
                                  !   if method=2, then hsp and vsp are the Barnes smoothing parameters (km*km)
                                  !   if method=3, then hsp and vsp define the influence region dimensions (km)
      integer i, j, k, n,mlv      ! loop variables
      integer ii,jj               ! more loop variables
      real minrange               ! minimum-range threshold (data closer to radar
                                  !   are discarded)

      real elmax                  ! CLZ input maximum elevation angle
      real rmax                   ! CLZ input maximum slant range
      real rhmax                  ! CLZ maximum allowable horizontal range to oban a datum.
      integer maxrh               ! CLZ input switch to control range thresholding

      integer mincount            ! threshold for minimum number of data gates required
                                  !   to produce an objectively-analyzed observation
      real minsum                 ! threshold for minimum sum of weights required to produce
                                  !   an objectively-analyzed observation
      integer extrapolate_flag    ! should extrapolation be allowed?
                                  !   1=yes (standard objective analysis)
                                  !   0=no (interpolation only)
      character(len=150) infile   ! input parameter file name
      character(len=150) outfile  ! output test file name
      character(len=150) v5dfile  ! output VIS5D file
      character(len=150) ncfile   ! output netcdf file

c
c.... date/time control parameters
c

c      read(lunr,*) yybeg, mmbeg, ddbeg, hrbeg, mnbeg, scbeg
c      read(lunr,*) yyend, mmend, ddend, hrend, mnend, scend

      integer*2 yerswp
      integer*2 monswp
      integer*2 dayswp
      integer*2 hrswp
      integer*2 minswp
      integer*2 secswp
      integer*2 yybeg
      integer*2 mmbeg
      integer*2 ddbeg
      integer*2 hrbeg
      integer*2 mnbeg
      integer*2 scbeg
      integer*2 yyend
      integer*2 mmend
      integer*2 ddend
      integer*2 hrend
      integer*2 mnend
      integer*2 scend

      real timeswp_fps
      real timeswp_beg
      real timeswp_end

      integer ncd
      integer nprefix
      integer ncm
      integer nc
      integer ncf
      integer iswp
      integer iread
      integer :: firstchar

      character*400 cmd           ! command line
      character*200 cmda          ! half command line "A"
      character*200 cmdb           ! half command line "B"
      character*200 swpfile       ! path and file of sweep files
      character*100 swp_pref      ! path of subdirectory containing sweeps comprising volume
      character*20 swp_mask       ! sweep file mask to identify candidate input sweep

      character*400 line          ! input record line from swp_files.txt

      character(len=200), allocatable :: sfname(:) ! sweep file names


      real, allocatable :: f(:,:,:,:)              ! data fields
      real, allocatable :: fld1d(:)
      real, allocatable :: fld2d(:,:)
      real, allocatable :: f_mpass(:,:,:,:,:)      ! multi-pass data fields (CLZ 2/26/08)
      character(len=8), allocatable :: fname(:)    ! field name array
      character(len=8) velfld                      ! name of velocity field to be analyzed
      character(len=8) reffld                      ! name of reflectivity field to be analyzed
      real, allocatable :: fsan(:)                 ! field sanity check value
      real, allocatable :: signf(:)                ! sign scale for fields
      real, allocatable :: sf(:)                   ! scaling factors for fields

      real, allocatable :: az(:,:,:)               ! interpolated azimuth angle (deg)
      real, allocatable :: el(:,:,:)               ! interpolated elevation angle (deg)
      real, allocatable :: time(:,:,:)             ! time, relative to central time (sec)

      real, allocatable :: xrd(:,:,:)             ! interpolated airborne radar x-coordinate (km)
      real, allocatable :: yrd(:,:,:)             ! interpolated airborne radar y-coordinate (km)
      real, allocatable :: zrd(:,:,:)             ! interpolated airborne radar z-coordinate (km)

      real, allocatable :: height(:,:,:)           ! interpolated height of obs (km)
      integer, allocatable :: count(:,:,:,:)       ! # gates used in interpolation
      integer, allocatable :: days(:)              ! Gregorian day (since beginning of base year) corresponding to sweep
      integer, allocatable :: secs(:)              ! seconds (total for current day) corresponding sweep
      real hsp0, vsp0                              ! hsp and vsp for multipass
      real gamma                                   ! gamma parameter
      integer pass,npass                           ! number of passes 

      integer indx, iargc, lunr

      integer istorm  ! istorm: = 0 (constant ut = ut0, vt = vt0); = 1 (time-varying ut, vt)
      real ut0, vt0                 ! translation velocity (m/s) at relative time 0
      real ut1, vt1                 ! translation velocity (m/s) at relative time 1
      real smfpt0, smfpt1   ! time (fp hrs) of initial and final (ut,vt)

      integer ikappa  ! ikappa: = 0 (constant hsp, vsp); = 1 (time-varying hsp, vsp)
      real hsp00, vsp00             ! hsp and vsp for multipass at relative time 0
      real hsp01, vsp01             ! hsp and vsp for multipass at relative time 1
      real spfpt0, spfpt1 ! time (fp hrs) of initial and final (hsp,vsp)

      real delfph
      real wgt0
      real wgt1
      real fhour
      real fmin
      real fsec
      real fphour
      real earthr
      real erad
      real olat
      real olon
      real oalt
      real zz
      real xpnt
      real ypnt
      real plat
      real plon
      real palt

      integer nwhalf
      real diff_thres
      real rhmax_near
      real hmax_near


      data earthr/6366.8056/
      parameter(amsg = 999.0)








c
c.... begin executable code
c









c############################################################################
c
c     Read the input parameters, and allocate space for data.
c
c############################################################################




      write(6,*)
      write(6,*) 'program oban_multi_v3'
      write(6,*)


c      write(6,*) 'Enter the name of the input parameter file:'
c
c      read(*,99) infile
c
c 99   format(A80)

c
c.... CLZ (4/4/16):  Specify the input file lunr and open file.  Lunr must
c       NOT be set to 10 or 11 which are used in output subroutines.  CAVEAT EMPTOR!
c


      indx = iargc( )
      
      IF ( indx .ge. 1 ) THEN
        call getarg ( 1, infile )
        write(0,'(a,a)') 'opening file: ',infile

c
c.... input unit number lunr = 9 is allowed
c

        lunr = 9

        open(unit = lunr,file = infile,status = 'unknown')
        IF ( indx .gt. 1 ) THEN
         write(0,*) 'Only reads first file! One argument, please!'
        ENDIF
      ELSE
        write(0,*) 'If nothing happens, remember to use ',
     :  ' x.oban < oban.in  OR  x.oban  namelist',
     :  ' on the command line because the namelist is not read by ',
     :  'default.'
      ENDIF
    
      GOTO 49
 39   CONTINUE
      write(0,*) 'where am i?'
      lunr = 9
      OPEN(unit=lunr,file='oban.in',status='old')

 49   CONTINUE

c
c.... now read in the namelist input data
c

c      open(unit=10,file=infile,status='old')

      write(6,*) 'about to read in 1st input namelist record'

c
c.... CLZ (9/10/18): new input parameter iradtype = 0 (ground-based), = 1 (airborne)
c

      read(lunr,*) iradtype

      if (iradtype .eq. 0) then
      write(6,*) 'iradtype=', iradtype,' (ground-based)'
      elseif (iradtype .eq. 1) then
      write(6,*) 'iradtype=', iradtype,' (airborne)'
      endif


      read(lunr,*) rlon, rlat, ralt  ! values not used in calculations if iradtype = 1
      write(6,*) 'rlon = ', rlon
      write(6,*) 'rlat = ', rlat
      write(6,*) 'ralt = ', ralt

      read(lunr,*) glon, glat, galt
      write(6,*) 'glon = ', glon
      write(6,*) 'glat = ', glat
      write(6,*) 'galt = ', galt

      read(lunr,*) map_proj
      write(6,*) 'map_proj = ', map_proj

      read(lunr,*) nx, ny, nz
      write(6,*) 'nx = ', nx
      write(6,*) 'ny = ', ny
      write(6,*) 'nz = ', nz

      read(lunr,*) dx, dy, dz
      write(6,*) 'dx = ', dx
      write(6,*) 'dy = ', dy
      write(6,*) 'dz = ', dz

      read(lunr,*) xmin, ymin, zmin
      write(6,*) 'xmin = ', xmin
      write(6,*) 'ymin = ', ymin
      write(6,*) 'zmin = ', zmin

c
c*****************************************************************************
c.... COMPUTE AND OUTPUT (LAT, LON) OF LOWER-LEFT AND UPPER RIGHT GRID CORNERS
c*****************************************************************************
c

c
c.... add new code to compute and print out (lat,lon) of origin and upper
c.... right corners of grid.  useful for plotting mobile mesonet data on grid.
c

c
c    CLZ addition:  compute lat/lon of grid points
c

      erad = earthr
      olat = glat
      olon = glon
      oalt = 0.0
      zz = 0.0

c
c.... check the input (lat,lon) against the (lat,lon) for input (x0,y0) = (0,0)
c          values, which should essentially verify the input grid-origin (lat,lon).
c

      xpnt = 0.0
      ypnt = 0.0

      write(6,*)
      write(6,*) 'calling grid_xy2ll for (lat,lon) at input (0,0)'
      write(6,*)'xpnt=',xpnt,' ypnt=',ypnt
      
      call grid_xy2ll(plat,plon,palt,xpnt,ypnt,zz,olat,olon,oalt,erad)

c
c..... CLZ (7/26/17): calculate lower left corner (lat,lon) coordinates of radar analysis grid
c

      xpnt = xmin
      ypnt = ymin

      write(6,*) 'calling grid_xy2ll for (lat,lon) at input (xmin,ymin)'
      write(6,*)'xpnt=',xpnt,' ypnt=',ypnt
      
      call grid_xy2ll(plat,plon,palt,xpnt,ypnt,zz,olat,olon,oalt,erad)

c
c..... CLZ (7/26/17): calculate upper right corner (lat,lon) coordinates of radar analysis grid
c

      xpnt = xmin + (nx - 1) * dx
      ypnt = ymin + (ny - 1) * dy

      write(6,*) 'calling grid_xy2ll for upper right corner (lat,lon)'
      write(6,*)'xpnt=',xpnt,' ypnt=',ypnt
      
      call grid_xy2ll(plat,plon,palt,xpnt,ypnt,zz,olat,olon,oalt,erad)

      write(6,*)

c
c.... done computing lower-left and upper-right grid corner (lat,lon)
c

c*****************************************************************************


c      integer istorm  ! istorm: = 0 (constant ut = ut0, vt = vt0); = 1 (time-varying ut, vt)
c      real ut0, vt0                 ! translation velocity (m/s) at relative time 0
c      real ut1, vt1                 ! translation velocity (m/s) at relative time 1
c      real smfpt0, smfpt1   ! time (fp hrs) of initial and final (ut,vt)

      read(lunr,*) istorm, ut0, vt0, smfpt0, ut1, vt1, smfpt1
c      read(lunr,*) ut, vt

      if (istorm .eq. 0) then
      write(6,*) 'constant (ut=ut0,vt=vt0)'

      elseif (istorm .eq. 1) then
      write(6,*) '(ut,vt) changing from (ut0,vt0) to (ut1,vt1)'
      write(6,*) 'ut0 = ', ut0
      write(6,*) 'vt0 = ', vt0
      write(6,*) 'ut1 = ', ut1
      write(6,*) 'vt1 = ', vt1
      write(6,*) 'smfpt0 = ', smfpt0
      write(6,*) 'smfpt1 = ', smfpt1

      endif

c
c.... CLZ (11/6/19): nfld is the number of observational data fields (e.g., velocity and dbz)
c

      read(lunr,*) nfld
      write(6,*) 'nfld = ', nfld

c
c.... CLZ (11/6/19): ngrdtot is the total number of gridded fields (i.e., data and 3 factors)
c

      ngrdtot = nfld + 3
      write(6,*) 'ngrdtot = ', ngrdtot

      allocate(fname(nfld+3))
      allocate(fsan(nfld))
      allocate(signf(nfld))

      read(lunr,*) velfld
      write(6,*)'velocity field name to be analyzed=',velfld
      read(lunr,*) reffld
      write(6,*)'reflectivity field name to be analyzed=',reffld

      do n=1, nfld
        read(lunr,*) fname(n), fsan(n), signf(n)
      enddo

      do n=1, nfld
        write(6,*) 'fname = ', fname(n)
        write(6,*) ' sanity check value fsan =',fsan(n)
        write(6,*) ' sign scale signf =',signf(n)
      enddo

c
c.... scale factor for fields sf
c

      allocate(sf(nfld+3))

      do n=1, nfld
      sf(n) = 1.0
        write(6,*) 'scale factor sf = ',sf(n),' for', fname(n)
      enddo

      sf(:) = 1.0


      read(lunr,*) method
      write(6,*) 'method = ', method

c      integer ikappa  ! ikappa: = 0 (constant hsp, vsp); = 1 (time-varying hsp, vsp)
c      real hsp00, vsp00             ! hsp and vsp for multipass at relative time 0
c      real hsp01, vsp01             ! hsp and vsp for multipass at relative time 1
c      real spfpt0, spfpt1 ! time (fp hrs) of initial and final (hsp,vsp)

      read(lunr,*) ikappa, hsp00, vsp00, spfpt0, hsp01, vsp01, spfpt1

c      read(lunr,*) hsp0
c      read(lunr,*) vsp0

      if (method.eq.1) then
        write(6,*) 'horizontal radius of influence = ', hsp00
        write(6,*) 'vertical radius of influence = ', vsp00

      elseif (method .eq. 2) then
        write(6,*) 'hsp00 = ', hsp00
        write(6,*) 'vsp00 = ', vsp00
        write(6,*) 'hsp01 = ', hsp01
        write(6,*) 'vsp01 = ', vsp01
        write(6,*) 'spfpt0 = ', spfpt0
        write(6,*) 'spfpt1 = ', spfpt1

      elseif (method.eq.3) then
        write(6,*) 'horizontal radius of influence = ', hsp00
        write(6,*) 'vertical radius of influence = ', vsp00
      endif

      read(lunr,*) extrapolate_flag
      write(6,*) 'extrapolate_flag = ', extrapolate_flag

      read(lunr,*) mincount, minsum
      write(6,*) 'mincount = ', mincount
      write(6,*) 'minsum = ', minsum

      read(lunr,*) minrange
      write(6,*) 'minrange = ', minrange

c
c.... CLZ (8/28/18): parameter maxrh controls oban to max range rmax on all sweeps
c

      read(lunr,*) elmax, rmax, maxrh
      write(6,*) 'elmax=',elmax
      write(6,*) 'rmax=',rmax
      write(6,*) 'maxrh=',maxrh

      if (maxrh .eq. 0) then
      write(6,*) 'doing oban to max range rmax on all sweeps'
      rhmax = rmax
      write(6,*) ' rhmax=',rhmax
      
      elseif (maxrh .eq. 1) then
      write(6,*) 'doing oban only to horiz range rhmax on all sweeps'
      rhmax = cos(dtor * elmax) * rmax
      write(6,*) ' rhmax=',rhmax

      endif

      write(6,*)' rmax=',rmax,' rhmax=',rhmax
      

      read(lunr,*) yrcor
      read(lunr,*) mocor
      write(6,*) 'yrcor = ', yrcor
      write(6,*) 'mocor = ', mocor

      read(lunr,*) dacor
c      read(lunr,*) dacor, minswpdc, maxswpdc
      minswpdc = 0
      maxswpdc = 0
      write(6,*) 'dacor = ', dacor

      read(lunr,*) umass_flag
      write(6,*) 'umass_flag = ', umass_flag

      read(lunr,*) az_corr_flag
      write(6,*) 'az_corr_flag = ', az_corr_flag

      read(lunr,*) elcor
      write(6,*) 'elcor = ', elcor

      read(lunr,*) azcor
      write(6,*) 'azcor=',azcor


c***********************************************************************


c
c.... CLZ (12/11/19): new control parameters for  detecting/flagging bad velocity
c

c      parameter (nwhalf = 5)
c      parameter (diff_thres = 10.0)
c      parameter (rhmax_near = 15.0)
c      parameter (hmax_near = 6.0)

      read(lunr,*) nwhalf, diff_thres, rhmax_near, hmax_near

      write(6,*)'nwhalf=',nwhalf,' diff_thres=',diff_thres
      write(6,*)'rhmax_near=',rhmax_near,' hmax_near=',hmax_near

c
c.... CLZ (9/8/17): new control parameters for cosmetic single-radar hole-filling
c

      read(lunr,*) iradfilt, isigfilt, ifilter, nuvfpass, nwfpass,
     * dbzmiss, dbzmiss0

      write(6,*)'iradfilt=',iradfilt,' ifilter=',ifilter
      write(6,*)'isigfilt=',isigfilt

      write(6,*)'nuvfpass=',nuvfpass
      write(6,*)'nwfpass=',nwfpass
      write(6,*)'dbzmiss=',dbzmiss
      write(6,*)'dbzmiss0=',dbzmiss0

c
c.... CLZ (9/8/17): new control parameters for cosmetic single-radar filtering
c

      read(lunr,*) ivsmooth, ismooth, nuvsmooth, nwsmooth

      write(6,*)'ivsmooth=',ivsmooth,' ismooth=',ismooth
      write(6,*)'nuvsmooth=',nuvsmooth
      write(6,*)'nwsmooth=',nwsmooth

c***********************************************************************


      read(lunr,*) npass
      write(6,*) 'npass=',npass

      read(lunr,*) gamma
      write(6,*) 'gamma=',gamma

      read(lunr,*) output_format
      write(6,*) 'output_format = ', output_format

c
c.... CLZ (8/21/18): introduce new control parameters for inputting sweep files
c


      read(lunr,*) maxswp
      write(6,*) 'maxswp=',maxswp

      read(lunr,*) swp_pref
      write(6,*) 'input prefix=',swp_pref

      read(lunr,*) swp_mask
      write(6,*) 'swp_mask=',swp_mask


c
c.... CLZ (4-1-16): new control parameter nanal = number of objective analyses to perform
c

      read(lunr,*) nanal
      write(6,*) 'nanal = ', nanal

      write(6,*)


c      open(unit=11, file='oban.output', status='unknown')


      write(6,*) 'nx= ',nx     
      write(6,*) 'ny= ',ny   
      write(6,*) 'nz= ',nz   
      write(6,*) 'nfld= ',nfld
      write(6,*) 'gamma= ',gamma
      write(6,*) 'hsp00= ',hsp00
      write(6,*) 'vsp00= ',vsp00
      write(6,*) 'npass= ',npass
      write(6,*)
      write(6,*)
      write(6,*)


c************************************************

c
c.... CLZ (8/21/18): create file list of sweep files based on DPJ's swps2cor.f90
c

c      character*400 cmd           ! command line
c      character*200 cmda          ! half command line "A"
c      character*200 cmdb           ! half command line "B"
c      character*200 swpfile       ! path and file of sweep files
c      character*100 swp_pref      ! path of subdirectory containing sweeps comprising volume
c      character*20 swp_mask       ! sweep file mask to identify candidate input sweep
c      character*400 line          ! input record line from swp_files.txt
c      character(len=200), allocatable :: sfname(:) ! sweep file names


  ! Read in the directory for the sweep files and the file name mask
  
c  Read (10,'(a)') Dir
      ncd = firstchar(swp_pref,' ') - 1

c
c.... length of prefix (including last forward slash) is ncd + 1
c

      nprefix = ncd + 1
      write(6,*) 'swp directory prefix length=',nprefix

c  ncd = FirstChar(Dir,' ') - 1

c  Read (10,'(a)') Mask
      ncm = firstchar(swp_mask,' ') - 1
c  ncm = FirstChar(Mask,' ') - 1

      write (6,*)
      write (6,*) 'Input DORADE sweep dir=', swp_pref(1:ncd)

c
c.... CLZ (8/27/19): Make the command to create a list of the sweep files *in the directory path prescribed on input*
c

      cmda = 'find ' // swp_pref(1:ncd) //
     * ' -maxdepth 1 -type f -name "swp.*"'

c      cmda = 'find ' // swp_pref(1:ncd) // ' -type f -name "swp.*"'
c      cmda = 'find ' // swp_pref(1:ncd) // ' -type f -name "swp.*_v1"'
      cmdb = ' > ' // swp_pref(1:ncd) // '/swp_files.txt'
      cmd = trim(cmda) // cmdb

c      cmda = 'ls ' // swp_pref(1:ncd) // '/' // swp_mask(1:ncm)
c      cmdb = ' > ' // swp_pref(1:ncd) // '/swp_files.txt'
c      cmd = trim(cmda) // cmdb

      nc = len_trim(cmd)
      write(6,*) 'cmd string length=',nc
      write(6,*) 'cmd=',cmd

c      nc = Len_Trim(Cmd)

c
c.... execute command to create list of sweep files in dir
c

c      use ifport
c      integer(4) iflag
c      integer(4) errnum

      iflag = system(cmd(1:nc))

c      call execute_command_line(cmd(1:nc))

      if (iflag .eq. -1) then
      errnum = ierrno()
      write(6,*) 'error=',errnum
      endif

      swpfile = swp_pref(1:ncd) // '/swp_files.txt'
c  File = Dir(1:ncd) // '/files.txt'

      ncf = firstchar(swpfile,' ') - 1
c  ncf = FirstChar(File,' ') - 1

      write (6,*) 'File of sweep files=', swpfile(1:ncf)

c************************************************




c
c.... CLZ (4-1-16): manage a sequence of nanal single-radar objective analyses
c


      do noban = 1, nanal

c
c############################################################################
c
c     Interpolate the sweep file data to a grid.
c
c############################################################################
c

      write(6,*)
      write(6,*) 'objective analysis noban=',noban

c
c.... read in control parameters and allocate arrays for this analysis time
c

      read(lunr,*) cyr, cmo, cda, chr, cmn, cse
      write(6,*) 'cyr = ', cyr
      write(6,*) 'cmo = ', cmo
      write(6,*) 'cda = ', cda
      write(6,*) 'chr = ', chr
      write(6,*) 'cmn = ', cmn
      write(6,*) 'cse = ', cse

      read(lunr,*) yybeg, mmbeg, ddbeg, hrbeg, mnbeg, scbeg
      write(6,*) 'begdate=',yybeg, mmbeg, ddbeg
      write(6,*) 'begtime=',hrbeg, mnbeg, scbeg
      
      timeswp_beg = real(hrbeg)*3600.0 + real(mnbeg)*60.0 + real(scbeg)
      write(6,*)'timeswp_beg=',timeswp_beg

      read(lunr,*) yyend, mmend, ddend, hrend, mnend, scend
      write(6,*) 'enddate=',yyend, mmend, ddend
      write(6,*) 'endtime=',hrend, mnend, scend
      
      timeswp_end = real(hrend)*3600.0 + real(mnend)*60.0 + real(scend)
      write(6,*)'timeswp_end=',timeswp_end


c      read(lunr,*) nswp
c      write(6,*) 'nswp = ', nswp

c
c.... allocate statements match the following deallocate statements at end of noban loop
c


c      deallocate(sfname)
c      deallocate(days)
c      deallocate(secs)
c      deallocate(fld2d)
c      deallocate(f)
c      deallocate(f_mpass)
c      deallocate(az)
c      deallocate(el)
c      deallocate(height)
c      deallocate(time)
c      deallocate(count)

c      allocate(sfname(nswp))

      write(6,*)' about to allocate sfname(maxswp)'
      write(6,*) 'maxswp=',maxswp
      allocate(sfname(maxswp))
      write(6,*)' allocated sfname(maxswp)'

c      do n=1, nswp
c        read(lunr,*) sfname(n)
c      enddo

      allocate(fld1d(nz))
      allocate(fld2d(nx,ny))
      allocate(f(nx,ny,nz,nfld+3))
      allocate(f_mpass(nx,ny,nz,nfld+3,npass))
      allocate(height(nx,ny,nz))

c.... ground-based radar
      allocate(az(nx,ny,nz))
      allocate(el(nx,ny,nz))
      allocate(time(nx,ny,nz))

c.... airborne radar
      allocate(xrd(nx,ny,nz))
      allocate(yrd(nx,ny,nz))
      allocate(zrd(nx,ny,nz))

      allocate(count(nx,ny,nz,nfld))



c************************************************

  !     Open file of sweep files

      open (20, file=swpfile(1:ncf), status='unknown', action='read')
c  Open (20, File=File(1:ncf), Status='Old', Action='Read')

c
c.... initialize sweep counter before scanning through file of sweep files
c

c      integer*2 yerswp
c      integer*2 monswp
c      integer*2 dayswp
c      integer*2 hrswp
c      integer*2 minswp
c      integer*2 secswp
c      integer*2 yybeg
c      integer*2 mmbeg
c      integer*2 ddbeg
c      integer*2 hrbeg
c      integer*2 mnbeg
c      integer*2 scbeg
c      integer*2 yyend
c      integer*2 mmend
c      integer*2 ddend
c      integer*2 hrend
c      integer*2 mnend
c      integer*2 scend

      iswp = 0
      iread = 0

      write(6,*)
      write(6,*) 'entering read(20) loop to scan swp_files.txt'

    1 read (20,'(a)',end=2) line
      nc = firstchar(line,' ') - 1

      iread = iread + 1
      
      read(line( (nprefix+10) : (nprefix+11) ),'(i2)') dayswp
      read(line( (nprefix+12) : (nprefix+13) ),'(i2)') hrswp
      read(line( (nprefix+14) : (nprefix+15) ),'(i2)') minswp
      read(line( (nprefix+16) : (nprefix+17) ),'(i2)') secswp
      
      timeswp_fps = real(hrswp)*3600.0 + real(minswp)*60.0 + real(secswp)

c
c.... check whether sweep is within date/time limits for current radar volume
c

      if (dayswp .ge. ddbeg .and. dayswp .le. ddend .and.
     * timeswp_fps .ge. timeswp_beg .and.
     * timeswp_fps .le. timeswp_end) then

c      if (dayswp .ge. ddbeg .and. hrswp .ge. hrbeg .and.
c     * minswp .ge. mnbeg .and. secswp .ge. scbeg .and.
c     * dayswp .le. ddend .and. hrswp .le. hrend .and.
c     * minswp .le. mnend .and. secswp .le. scend .and.
c     * hrswp .ge. 0 .and. minswp .ge. 0 .and. secswp .ge. 0)  then

      iswp = iswp + 1

      write(6,*)
      write(6,*) 'swp=',iswp,' file=',line(1:nc)
      write(6,*)'dy=',dayswp,' hr=',hrswp,' mn=',minswp,' sc=',secswp
      write(6,*)'timeswp_fps=',timeswp_fps

c
c....  scan file times from the file of files and populate sfname sweep array list
c

      sfname(iswp) = line(1:nc)
      write(6,*) 'sfname(iswp)=',sfname(iswp)

c
c.... endif (dayswp .ge. ddbeg .and. dayswp .le. ddend .and.
c     * timeswp_fph .ge. timeswp_beg .and.
c     * timeswp_fph .le. timeswp_end) then
c

      endif

      goto 1
  
    2 close (20)

      write(6,*)
      write(6,*) 'reached end of file of ', swpfile(1:ncf)

      nswp = iswp
      if ( (output_format.eq.2) .or. (output_format.eq.3) ) nz = nswp
      write(6,*) 'nswp = ', nswp

c************************************************



c************************************************

      allocate(days(nswp))
      allocate(secs(nswp))

c************************************************


      if ( (output_format.eq.2) .or. (output_format.eq.3) ) nz = nswp

      read(lunr,*) outfile
      write(6,*) 'outfile = ', outfile

      read(lunr,*) v5dfile
      write(6,*) 'v5dfile = ', v5dfile

      read(lunr,*) ncfile
      write(6,*) 'ncfile = ', ncfile


c
c.... CLZ (7/24/17): compute the analysis time in units of floating-point hours.
c

      fhour = chr
      fmin = cmn
      fsec = cse
      
      if (fhour .ge. 10.0 .and. fhour .le. 23.0) then
      fphour = fhour + (fmin/60.0) + (fsec/3600.0)
      write(6,*) 'analysis time (floating-point hours) =',fphour
      elseif (fhour .ge. 0.0 .and. fhour .lt. 10.0) then
      fphour = 24.0 + fhour + (fmin/60.0) + (fsec/3600.0)
      write(6,*) 'analysis time (floating-point hours) =',fphour
      endif


c
c.... CLZ (7/24/17): prescribe the (possibly time-dependent) storm motion
c

c      integer istorm  ! istorm: = 0 (constant ut = ut0, vt = vt0); = 1 (time-varying ut, vt)
c      real ut0, vt0                 ! translation velocity (m/s) at relative time 0
c      real ut1, vt1                 ! translation velocity (m/s) at relative time 1
c      real smfpt0, smfpt1   ! time (fp hrs) of initial and final (ut,vt)

      if (istorm .eq. 0 ) then
      ut = ut0
      vt = vt0
      
      elseif (istorm .eq. 1 ) then

      if (fphour .le. smfpt0) then
      ut = ut0
      vt = vt0
      elseif (fphour .ge. smfpt1) then
      ut = ut1
      vt = vt1
      elseif (fphour .gt. smfpt0 .and. fphour .lt. smfpt1) then
      delfph = smfpt1 - smfpt0
      wgt0 = 1.0 - ((fphour - smfpt0)/delfph)
      wgt1 = 1.0 - wgt0
      ut = (wgt0 * ut0) + (wgt1 * ut1)
      vt = (wgt0 * vt0) + (wgt1 * vt1)
      endif

      write(6,*) 'ut/vt = ', ut,vt,' at fphour =',fphour

      endif

c
c.... CLZ (7/24/17): prescribe the (possibly time-dependent) smoothing parameter Kappa
c


c      integer ikappa  ! ikappa: = 0 (constant hsp, vsp); = 1 (time-varying hsp, vsp)
c      real hsp00, vsp00             ! hsp and vsp for multipass at relative time 0
c      real hsp01, vsp01             ! hsp and vsp for multipass at relative time 1
c      real spfpt0, spfpt1 ! time (fp hrs) of initial and final (hsp,vsp)

      if (ikappa .eq. 0 ) then
      hsp0 = hsp00
      vsp0 = vsp00
      
      elseif (ikappa .eq. 1 ) then

      if (fphour .le. spfpt0) then
      hsp0 = hsp00
      vsp0 = vsp00
      elseif (fphour .ge. spfpt1) then
      hsp0 = hsp01
      vsp0 = vsp01
      elseif (fphour .gt. spfpt0 .and. fphour .lt. spfpt1) then
      delfph = spfpt1 - spfpt0
      wgt0 = 1.0 - ((fphour - spfpt0)/delfph)
      wgt1 = 1.0 - wgt0
      hsp0 = (wgt0 * hsp00) + (wgt1 * hsp01)
      vsp0 = (wgt0 * vsp00) + (wgt1 * vsp01)
      endif

      write(6,*) 'hsp0/vsp0 = ', hsp0,vsp0,' at fphour =',fphour

      endif



c
c.... CLZ (5/18/07):  Preferred method is multipass Barnes oban.
c      Multipass Barnes Originators: Mario Majcen, Paul Markowski, Yvette Richardson (PSU)
c      Creation Date:  January-April 2007.
c


c
c.... CLZ (9/6/2018):  Insert mincount into argument list of call to multi_oban
c


c
c.... CLZ (11/23/10): Added rhmax to argument list of call to multi_oban
c


      if (output_format .eq. 1) then

c
c.... CLZ (12/11/19): added bad data detection/rejection to new airborne radar option
c

      call multi_oban(method, gamma, npass, hsp0, vsp0, mincount,
     $   f, f_mpass, nfld, fname, velfld, reffld,
     $     fsan, signf, nswp, maxswp, sfname,
     $        iradtype, extrapolate_flag, minsum, minrange, rhmax,
     $          yrcor, mocor, dacor, minswpdc, maxswpdc,
     $            azcor, elcor, az_corr_flag, umass_flag,
     $              az, el, time, xrd, yrd, zrd, count,
     $                map_proj, glat, glon, galt, rlat, rlon, ralt,
     $                  nx, ny, nz, dx, dy, dz, xmin, ymin, zmin,
     $                    cyr, cmo, cda, chr, cmn, cse, ut, vt,
     $                      nwhalf, diff_thres, rhmax_near, hmax_near)



      elseif ( (output_format.eq.2) .or. (output_format.eq.3) ) then
        call ppi_oban(method, hsp,
     $                f, nfld, fname, sfname,
     $                minsum, minrange,
     $                yrcor, mocor, dacor,
     $                azcor, elcor, az_corr_flag, umass_flag,
     $                az, el, height, time, count,
     $                map_proj, glat, glon, galt, rlat, rlon, ralt,
     $                nx, ny, nswp, dx, dy, xmin, ymin,
     $                cyr, cmo, cda, chr, cmn, cse, ut, vt,
     $                days, secs)



      elseif (output_format .eq. 0) then

c
c.... CLZ (2/7/11): new option to output gate data for scatterplotting
c

        call output_gates(method, gamma, npass, hsp0, vsp0,
     $                   f, f_mpass, nfld, fname, nswp, sfname,
     $                   extrapolate_flag,minsum,minrange,rhmax,
     $                   yrcor, mocor, dacor,
     $                   azcor, elcor, az_corr_flag, umass_flag,
     $                   az, el, time, count,
     $                   map_proj, glat, glon, galt, rlat, rlon, ralt,
     $                   nx, ny, nz, dx, dy, dz, xmin, ymin, zmin,
     $                   cyr, cmo, cda, chr, cmn, cse, ut, vt)




      else
        write(6,*) 'unknown method:  method = ', method,
     $             ', output_format = ', output_format
        stop
      endif


      if (output_format .gt. 0) then



c      write(6,*) 'aggregating f_mpass into f from ',npass,'-pass oban'


      do pass = 1, npass

      do n=1,nfld
c      write(6,*) 'printing f_mpass(ii,jj,1,n,pass) at pass= ',pass
c      write(6,*) 'field nfld (1=VE, 2=DZ) = ',n
      do jj=ny,1,-1
c      write(6,9877) jj,(f_mpass(ii,jj,1,n,pass),ii=35,39)
 9877 format(i3,5e14.6)
      enddo
      enddo

      do n=1, nfld
        do k=1, nz
          do j=1, ny
            do i=1, nx

c
c.... printout suggests missing values = 9.9e9
c

            if(abs(f_mpass(i,j,k,n,pass)) .lt. 900.0) then
                f(i,j,k,n) = f_mpass(i,j,k,n,pass)
            endif

            enddo
          enddo
        enddo
      enddo


      do n=1,nfld
c      write(6,*) 'printing f(ii,jj,1,n) at pass= ',pass
c      write(6,*) 'field nfld (1=VE, 2=DZ) = ',n
      do jj=ny,1,-1
c      write(6,9877) jj,(f(ii,jj,1,n),ii=35,39)
      enddo
      enddo


      enddo



c
c############################################################################
c
c     Eliminate gridded values that were determined by only a few
c     raw observations.
c
c############################################################################
c


c      write(6,*) 'threshold for count = ', mincount
      do n=1, nfld
        do k=1, nz
          do j=1, ny
            do i=1, nx
              if (count(i,j,k,n).lt.mincount) then
                f(i,j,k,n) = bad
              endif
            enddo
          enddo
        enddo
      enddo



      if (iradtype .eq. 0) then

c
c.... ground-based
c

      f(:,:,:,nfld+1) = az(:,:,:)
      fname(nfld+1) = 'AZ'
      f(:,:,:,nfld+2) = el(:,:,:)
      fname(nfld+2) = 'EL'
      f(:,:,:,nfld+3) = time(:,:,:)
      fname(nfld+3) = 'TIME'

      elseif (iradtype .eq. 1) then

c
c.... airborne
c

      f(:,:,:,nfld+1) = xrd(:,:,:)
      fname(nfld+1) = 'XR'
      f(:,:,:,nfld+2) = yrd(:,:,:)
      fname(nfld+2) = 'YR'
      f(:,:,:,nfld+3) = zrd(:,:,:)
      fname(nfld+3) = 'ZR'

      endif




c
c***************************************************************************
c***************************************************************************
c

c
c.... CLZ (9/8/17): perform single-radar hole-fill between closely spaced rays
c

      if(iradfilt .eq. 1 .and. ifilter .le. 3) then

c      write(6,*) ' about to call subroutine radfill3d'
      write(6,*) ' about to call subroutine radfill1d'
c      write(6,*) ' about to call subroutine radfill'

c
c.... now call radfill1d
c

c      call radfill3d(f, fld2d, nfld, nx, ny, nz,
c     * dx, dy, dz, xmin, ymin, zmin,
c     * bad, ifilter)

      call radfill1d(f, fld1d, ngrdtot, nx, ny, nz,
     * dx, dy, dz, xmin, ymin, zmin,
     * bad, ifilter)

c      call radfill1d(f, fld1d, nfld, nx, ny, nz,
c     * dx, dy, dz, xmin, ymin, zmin,
c     * bad, ifilter)

c      call radfill(f, fld2d, nfld, nx, ny, nz,
c     * dx, dy, dz, xmin, ymin, zmin,
c     * bad, ifilter)

      write(6,*)'returned to oban_multi from subroutine radfill1d'
c      write(6,*)'returned to oban_multi from subroutine radfill'

c
c.... endif(iradfilt .eq. 1 .and. ifilter .le. 3) then
c

      endif


c
c.... done with optional single-radar hole-filling
c


c
c***************************************************************************
c***************************************************************************
c




c
c############################################################################
c
c     Output the results.
c
c############################################################################
c



      write(6,*) 'creating output files...'


      if ( (outfile.ne.'none') .and. (outfile.ne.'NONE') ) then

        write(6,*) 'text file...'
        if (output_format.eq.1) then
          call sdout(glon, glat, galt, cyr, cmo, cda, chr, cmn, cse,
     $               nx, ny, nz, dx, dy, dz, xmin, ymin, zmin, ut, vt,
     $               nfld+3, fname, f, count, az, el, time, outfile)

        else if (output_format.eq.2) then
          call vdrasout(glon, glat, galt,
     $                  nx, ny, nswp, dx, dy, xmin, ymin, zmin,
     $                  f(1,1,1,2), f(1,1,1,1),
     $                  count, az, el, height, outfile)

        else if (output_format.eq.3) then
          call DARTout(glon, glat, galt, rlon, rlat, ralt, map_proj,
     $                 nx, ny, nswp, dx, dy, xmin, ymin,
     $                 f, nfld,
     $                 az, el, height, outfile, days, secs)
        endif

      endif

c
c.... vis5d output
c

      if ( (v5dfile.ne.'none') .and. (v5dfile.ne.'NONE') ) then

        if (nz.gt.100) then
          write(6,*) 'not creating vis5d file because nz = ', nz
        else
          write(6,*) 'creating vis5d file...'
          call writev5d(v5dfile, nx, ny, nz, dx, dy, dz,
     $                  xmin, ymin, zmin, nfld+3, fname, f, 
     $                  glat, glon, cyr, cmo, cda, chr, cmn, cse)

        endif

      endif



c
c.... netcdf output
c

      if ( (ncfile.ne.'none') .and. (ncfile.ne.'NONE') ) then

        write(6,*) 'creating netcdf file...'
        call writenetcdf(ncfile, nfld+3, fname, f, sf,
     *       glat, glon, galt, rlat, rlon, ralt,
     *       nx, ny, nz, dx, dy, dz,
     *       xmin, ymin, zmin, ut, vt,
     *       int(cyr), int(cmo), int(cda),
     *       int(chr), int(cmn), int(cse))

      endif



c
c....endif (output_format .gt. 0) then
c

      endif

c
c... deallocate selected arrays for this time level
c

      deallocate(sfname)
      deallocate(days)
      deallocate(secs)
      deallocate(fld2d)
      deallocate(f)
      deallocate(f_mpass)
      deallocate(height)
      deallocate(az)
      deallocate(el)
      deallocate(time)
      deallocate(xrd)
      deallocate(yrd)
      deallocate(zrd)
      deallocate(count)

c
c.... enddo noban = 1, nanal
c

      enddo

c
c.... done looping through all analysis times
c


      close(lunr)
      deallocate(sf)


c
c############################################################################
c
c     Clean up.
c
c############################################################################
c

      if (output_format .gt. 0) then

      close(11)

      deallocate(fname)

c      deallocate(f)
c      deallocate(az)
c      deallocate(el)
c      deallocate(height)
c      deallocate(time)
c      deallocate(count)
c      deallocate(sfname)
c      deallocate(days)
c      deallocate(secs)

      endif





      write(6,*) 'finished oban_multi_v3'

      stop
      end




      Integer Function FirstChar(String,Char)

c
c      Find first occurence of character in character string.               
c      Ex.     FirstChar('001234','2') will equal 4
c      Ex.     FirstChar('0012','2') will equal 4   
c      Ex.     FirstChar('223456','2') will equal 1  
c      Ex.     FirstChar('123','4') will equal 0,
c        i.e., zero will be returned if character is not present in string
c        (although one might have logically guessed the string_length + 1
c        would be returned).
c

      Implicit none                                                             
      Character String*(*),Char*1
      Integer*4 LenString
  
      LenString = Len(String)   
      FirstChar = 1

      Do While (FirstChar .le. LenString)
         If (String(FirstChar:FirstChar) .eq. Char) Return
         FirstChar = FirstChar +1
      End Do
  
      FirstChar = 0    ! character not found!!
  
      Return                                                                    
      End Function FirstChar ! Integer*4 Function FirstChar ends







      subroutine grid_xy2ll(plat, plon, palt, x, y, z, olat, olon,
     * oalt, R_earth)

c
c.... calculate (plat,plon) of a point at (x,y) relative to origin (olat,olon).
c.... all dimensions in km.
c

c    /* transform to earth coordinates and then to lat/lon/alt */
c
c    /* These calculations are from the book
c     * "Aerospace Coordinate Systems and Transformations"
c     * by G. Minkler/J. Minkler
c     * these are the ECEF/ENU point transformations
c     */


      double precision delta_o, lambda_o
      double precision delta_p, lambda_p
      double precision R, rr
      double precision sinLambda, cosLambda, sinDelta, cosDelta
      double precision xe, ye, ze, h
      double precision rad_cnvt

      data rad_cnvt / .017453292 /
      data deg_cnvt / 57.29577951 /



c
c.... begin executable code
c




      h = R_earth + oalt
      delta_o = rad_cnvt * ( olat )	
      lambda_o = rad_cnvt * ( olon )
      sinLambda = sin( lambda_o )
      cosLambda = cos( lambda_o )
      sinDelta = sin( delta_o )
      cosDelta = cos( delta_o )

c      write(6,*) 
c      write(6,*) 'cosdelta =',cosdelta
c      write(6,*) 'sinlambda =',sinlambda
c      write(6,*) 'coslambda =',coslambda
c      write(6,*) 'x =',x
c      write(6,*) 'sindelta =',sindelta
c      write(6,*) 'y =',y
c      write(6,*) 'z =',z
c      write(6,*) 

c    
c	/* transform to earth coordinates */
c

      xe = h * sinDelta + cosDelta * y + sinDelta * z

      ye = -h * cosDelta * sinLambda   -cosLambda * x
     * + sinLambda * sinDelta * y -sinLambda * cosDelta * z

      ze = h * cosDelta * cosLambda   -sinLambda * x
     * - cosLambda * sinDelta * y + cosLambda * cosDelta * z

	  lambda_p = datan2( -ye, ze )
	  delta_p = datan2( xe, sqrt( ye * ye + ze * ze ))

	  plat = deg_cnvt * ( delta_p )
	  plon = deg_cnvt * ( lambda_p )
	  palt = dsqrt( xe * xe + ye * ye + ze * ze ) - R_earth

c
c.... write out the returned (lat,lon) coordinate for input point (x,y)
c

      write(6,*) 'in xy2ll',' plat=',plat,' plon=',plon,' palt=',palt,
     * ' x=',x,' y=',y,' z=',z


c
c.... done with grid_xy2ll
c

      return
      end