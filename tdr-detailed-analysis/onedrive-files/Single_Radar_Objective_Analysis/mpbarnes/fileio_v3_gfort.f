c############################################################################
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######               SUBROUTINE WRITEDDTEXT                 ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
c
c############################################################################
c
c     PURPOSE:
c
c     This subroutine writes out an ascii text file containing the
c     synthesized dual-Doppler wind field.
c
c############################################################################


      subroutine writeddtext(synthfile, nx, ny, nz, dx, dy, dz,
     $                       xmin, ymin, zmin, rho,
     $                       u, v, w, rf)

      implicit none

      integer nx, ny, nz          ! no. of grid points in x, y, and z directions
      real xmin, ymin, zmin       ! coordinates of lower southwest corner
                                  !   of grid, relative to origin (km)
      real dx, dy, dz             ! grid spacing in x, y, and z directions (km)
      character(len=80) synthfile ! output text synthesis file
      integer i, j, k             ! grid indices
      real rf(nx,ny,nz)           ! reflectivity (dBZ)
      real u(nx,ny,nz)            ! u component of velocity (m/s)
      real v(nx,ny,nz)            ! v component of velocity (m/s)
      real w(nx,ny,nz)            ! w component of velocity (m/s)
      real rho                    ! density at lowest grid level (kg m**-3)


      open(unit=11, file=synthfile, status='unknown')
      write(11,*) nx, ny, nz
      write(11,*) dx*1000.0, dy*1000.0, dz*1000.0
      write(11,*) xmin*1000.0, ymin*1000.0, zmin*1000.0
      write(11,*) rho
 88   format(E11.4)
c 88   format(F8.2)
      do k=1, nz
        do j=1, ny
          write(11,88) (u(i,j,k), i=1, nx)
        enddo
      enddo
      do k=1, nz
        do j=1, ny
          write(11,88) (v(i,j,k), i=1, nx)
        enddo
      enddo
      do k=1, nz
        do j=1, ny
          write(11,88) (w(i,j,k), i=1, nx)
        enddo
      enddo
      do k=1, nz
        do j=1, ny
          write(11,88) (rf(i,j,k), i=1, nx)
        enddo
      enddo
      close(11)

      return
      end





c############################################################################
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                 SUBROUTINE WRITEV5D                  ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
c
c############################################################################
c
c     PURPOSE:
c
c     This subroutine writes out a VIS5D format file containing the
c     synthesized dual-Doppler fields.
c
c############################################################################

      subroutine writev5d(v5dfile, nx, ny, nz, dx, dy, dz,
     $                    xmin, ymin, zmin, nfld, fname, f,
     $                    glat, glon, cyr, cmo, cda, chr, cmn, cse)

      implicit none

      include 'v5df.h'
      include 'dow.inc'

c---- input parameters

      character(len=150) v5dfile   ! output file
      integer nx, ny, nz           ! no. of grid points in x, y, and z directions
      real dx, dy, dz              ! grid spacing in x, y, and z directions (km)
      real xmin, ymin, zmin        ! coordinates of lower southwest corner
                                   !   of grid, relative to origin (km)
      integer nfld                 ! number of fields
      character(len=8) fname(nfld) ! field names
      real f(nx,ny,nz,nfld)        ! data

      real glat, glon              ! grid origin

      integer(kind=2) cyr,cmo,cda ! central date
      integer(kind=2) chr,cmn,cse ! central time

c---- local variables

      integer i, j, k, s           ! grid indices

c---- Vis5d variables

      integer n
      integer it, iv
      real(kind=4) G(ny, nx, nz)

      integer nr, nc, nl
      integer numtimes
      integer numvars
      character(len=10) varname(MAXVARS)
      integer dates(MAXTIMES)
      integer times(MAXTIMES)
      real northlat
      real latinc
      real westlon
      real loninc
      real bottomhgt
      real hgtinc

      integer ndays(12)

c---- Vis5d variables initialized to missing values

      data nr,nc,nl / IMISSING, IMISSING, IMISSING /
      data numtimes,numvars / IMISSING, IMISSING /
      data (varname(i),i=1,MAXVARS) / MAXVARS*"          " /
      data (dates(i),i=1,MAXTIMES) / MAXTIMES*IMISSING /
      data (times(i),i=1,MAXTIMES) / MAXTIMES*IMISSING /
      data northlat, latinc / MISSING, MISSING /
      data westlon, loninc / MISSING, MISSING /
      data bottomhgt, hgtinc / MISSING, MISSING /
      data ndays /31,29,31,30,31,30,31,31,30,31,30,31/






c
c.... begin executable code
c





      write(6,*)'now inside writev5d'



c     initialize variables

      nr = ny
      nc = nx
      nl = nz
      numtimes = 1
      numvars = nfld

      write(6,*)'numvars (nfld)=',numvars

      do s=1, nfld
        varname(s) = fname(s)
      enddo
      dates(1) = 0
      times(1) = 0
c      northlat = ymin + (ny-1)*dy
      northlat = glat + (ny-1)*dy*(1./111.) + ymin*(1./111.) 
c      latinc = dy
      latinc = dy*(1./111.)
c      westlon = -xmin
      westlon = -glon-xmin*(180./3.14159)/(6371.*cos(glat*3.14159/180.))
c      loninc = dx
      loninc = dx*(180./3.14159)/(6371.*cos(glat*3.14159/180.))
c      do while ( (northlat.gt.90.0) .or. ((northlat-(ny-1)*latinc).lt.-90.0) .or.
c     $           (westlon.lt.-90.0) .or. ((westlon+(nx-1)*loninc).gt.90.0) )
c        northlat = 0.1*northlat
c        latinc = 0.1*latinc
c        westlon = 0.1*westlon
c        loninc = 0.1*loninc
c      enddo
      bottomhgt = zmin
      hgtinc = dz

! Compute Julian day for vis5d.

      dates(1) = 2004000
      if (cmo.gt.1) then
         do i = 1, cmo-1
            dates(1) = dates(1) + ndays(i)
         enddo
      endif
      dates(1) = dates(1) + cda
      times(1) = chr*10000 + cmn*100 + cse



c
c.... create the v5d file.
c

      n = v5dcreatesimple(v5dfile, numtimes, numvars, nr, nc, nl,
     *                 varname, times, dates,
     *                 northlat, latinc, westlon, loninc,
     *                 bottomhgt, hgtinc )



      write(6,*)'just after v5dcreatesimple'



      if (n .eq. 0) then
        write(*,*) '*** Error creating v5d file.'
        stop
      endif





      do it=1,numtimes
         do iv=1,numvars

           do k=1, nl
             do j=1, nr
               do i=1, nc

                 G(nr-j+1,i,k) = f(i,j,k,iv)

                 if (G(nr-j+1,i,k).eq.bad) then
                   G(nr-j+1,i,k) = 9.9E30
                 endif

               enddo
             enddo
           enddo



      write(6,*)'just after G array write-out with iv=',iv



c          write the 3-D grid to the v5d file

           n = v5dwrite( it, iv, G )


      write(6,*)'just after call to v5dwrite with returned n=',n

           if (n .eq. 0) then
             write(*,*) '*** Error writing to v5d file.'
             stop
           endif


c
c.... enddo iv=1,numvars
c

         enddo

c
c.... enddo it=1,numtimes
c

      enddo








c     close the v5d file and exit

      n = v5dclose()

      if (n .eq. 0) then
        write(*,*) '*** Error closing v5d file.'
        stop
      endif



      return
      end





c############################################################################
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                SUBROUTINE WRITENETCDF                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
c
c############################################################################
c
c     PURPOSE:
c
c     This subroutine writes out a netcdf file containing the gridded fields.
c
c############################################################################

c
c.... DS - 2/28/20 - Added a few lines of code in this subroutine to
c       make the ncgen input file name more unique, and then to remove
c       that file after the ncgen command completes. Previously, if running
c       several analyses at once in the same directory, it was possible for
c       the analyses to share whichever ncgen.input file was newest, thus
c       producing duplicate netcdf output files

      subroutine writenetcdf(ncfile, nfld, fname, f, sf,
     $                       glat, glon, galt, rlat, rlon, ralt,
     $                       nx, ny, nz, dx, dy, dz, xmin,
     $                       ymin, zmin, ut, vt, yr, mo, da, hr, mn, se)

c
c.... CLZ (5/17/07):  comment out "use iflport" -- not available
c

c      use iflport

      implicit none

      include 'dow.inc'

c---- input parameters

      character(len=150) ncfile    ! netcdf file name
      integer nfld                 ! number of fields
      character(len=8) fname(nfld) ! field names
      integer nx, ny, nz           ! no. of grid points in x, y, and z directions
      real f(nx,ny,nz,nfld)        ! data
      real sf(nfld)                ! scaling factor
      real glat, glon              ! latitude and longitude of grid origin (deg)
      real galt                    ! altitude of grid origin (km MSL)
      real rlat, rlon              ! latitude and longitude of radar (deg)
      real ralt                    ! altitude of radar (km MSL)
      real dx, dy, dz              ! grid spacing in x, y, and z directions (km)
      real xmin, ymin, zmin        ! coordinates of lower southwest corner
                                   !   of grid, relative to origin (km)
      real ut, vt                  ! storm translation velocity (m/s)
      integer yr, mo, da           ! year, month, and day
      integer hr, mn, se           ! hour, minute, and second


c---- local variables

      integer i, j, k, n, s
      integer ls                   ! string length
      character(len=320) command, command2
      character(len=170) ncgenf
      integer status









c     change bad data flag

      do n=1, nfld
        do k=1, nz
          do j=1, ny
            do i=1, nx
              if (f(i,j,k,n).eq.bad) f(i,j,k,n)=sbad
            enddo
          enddo
        enddo
      enddo

c     open ascii file that will be converted to netcdf

      ls = index(ncfile, ' ') - 1
      ncgenf = ncfile(1:ls-3) // '.ncgen'
      open(unit=11, file=ncgenf, status='unknown')

      
      write(11,'(3A)') 'netcdf ', ncfile(1:ls), ' {'

c     "dimensions" section

      write(11,'(A)') 'dimensions:'
      write(11,'(T9,A,I12,A)') 'fields = ', nfld, ' ;'
      write(11,'(T9,A)') 'long_string = 80 ;'
      write(11,'(T9,A)') 'short_string = 8 ;'
      write(11,'(T9,A)') 'date_string = 10 ;'
      write(11,'(T9,A)') 'time = UNLIMITED ; // (1 currently)'
      write(11,'(T9,A,I12,A)') 'x = ', nx, ' ;'
      write(11,'(T9,A,I12,A)') 'y = ', ny, ' ;'
      write(11,'(T9,A,I12,A)') 'z = ', nz, ' ;'
      write(11,'(T9,A,I12,A)') 'el = ', nz, ' ;'

c     "variables" section

      write(11,'(A)') 'variables:'
      write(11,'(T9,A)') 'char start_date(date_string) ;'
      write(11,'(T9,A)') 'char end_date(date_string) ;'
      write(11,'(T9,A)') 'char start_time(short_string) ;'
      write(11,'(T9,A)') 'char end_time(short_string) ;'
      write(11,'(T9,A)') 'int bad_data_flag ;'
      write(11,'(T9,A)') 'float grid_latitude ;'
      write(11,'(T17,A)') 'grid_latitude:value="Grid origin latitude" ;'
      write(11,'(T17,A)') 'grid_latitude:units = "deg" ;'
      write(11,'(T9,A)') 'float grid_longitude ;'
      write(11,'(T17,A)') 'grid_longitude:value="Gridorigin longitude";'
      write(11,'(T17,A)') 'grid_longitude:units = "deg" ;'
      write(11,'(T9,A)') 'float grid_altitude ;'
      write(11,'(T17,A)') 'grid_altitude:value="Altitudeof gridorigin";'
      write(11,'(T17,A)') 'grid_altitude:units = "km MSL" ;'
      write(11,'(T9,A)') 'float radar_latitude ;'
      write(11,'(T17,A)') 'radar_latitude:value = "Radar latitude" ;'
      write(11,'(T17,A)') 'radar_latitude:units = "deg" ;'
      write(11,'(T9,A)') 'float radar_longitude ;'
      write(11,'(T17,A)') 'radar_longitude:value = "Radar longitude" ;'
      write(11,'(T17,A)') 'radar_longitude:units = "deg" ;'
      write(11,'(T9,A)') 'float radar_altitude ;'
      write(11,'(T17,A)') 'radar_altitude:value = "Altitude of radar" ;'
      write(11,'(T17,A)') 'radar_altitude:units = "km MSL" ;'
      write(11,'(T9,A)') 'float x_min ;'
      write(11,'(T17,A)') 'x_min:value = "The minimum x grid value" ;'
      write(11,'(T17,A)') 'x_min:units = "km" ;'
      write(11,'(T9,A)') 'float y_min ;'
      write(11,'(T17,A)') 'y_min:value = "The minimum y grid value" ;'
      write(11,'(T17,A)') 'y_min:units = "km" ;'
      write(11,'(T9,A)') 'float z_min ;'
      write(11,'(T17,A)') 'z_min:value = "The minimum z grid value" ;'
      write(11,'(T17,A)') 'z_min:units = "km" ;'
      write(11,'(T9,A)') 'float x_max ;'
      write(11,'(T17,A)') 'x_max:value = "The maximum x grid value" ;'
      write(11,'(T17,A)') 'x_max:units = "km" ;'
      write(11,'(T9,A)') 'float y_max ;'
      write(11,'(T17,A)') 'y_max:value = "The maximum y grid value" ;'
      write(11,'(T17,A)') 'y_max:units = "km" ;'
      write(11,'(T9,A)') 'float z_max ;'
      write(11,'(T17,A)') 'z_max:value = "The maximum z grid value" ;'
      write(11,'(T17,A)') 'z_max:units = "km" ;'
      write(11,'(T9,A)') 'float x_spacing ;'
      write(11,'(T17,A)') 'x_spacing:units = "km" ;'
      write(11,'(T9,A)') 'float y_spacing ;'
      write(11,'(T17,A)') 'y_spacing:units = "km" ;'
      write(11,'(T9,A)') 'float z_spacing ;'
      write(11,'(T17,A)') 'z_spacing:units = "km" ;'

      write(11,'(T9,A)') 'float x(x) ;'
      write(11,'(T9,A)') 'float y(y) ;'
      write(11,'(T9,A)') 'float z(z) ;'
      write(11,'(T9,A)') 'float el(z) ;'
      write(11,'(T9,A)') 'float u_translation ;'
      write(11,'(T17,A)') 'u_translation:value = "storm motion u component" ;'
      write(11,'(T17,A)') 'u_translation:units = "m/s" ;'
      write(11,'(T9,A)') 'float v_translation ;'
      write(11,'(T17,A)') 'v_translation:value = "storm motion v component" ;'
      write(11,'(T17,A)') 'v_translation:units = "m/s" ;'
      write(11,'(T9,A)') 'char field_names(fields, short_string) ;'

      do n=1, nfld
        ls = index(fname(n), ' ') - 1
        write(11,'(T9,A,A,A)') 'float ', fname(n)(1:ls),
     $                         '(time, z, y, x) ;'
        write(11,'(T17,A,A,F12.9,A)') fname(n)(1:ls), ':scale_factor = '
     $                           , sf(n), ' ;'
        write(11,'(T17,A,A)') fname(n)(1:ls), ':add_offset = 0.0 ;'
        write(11,'(T17,A,A,F13.5,A,F13.5,A)') fname(n)(1:ls),
     $                       ':missing_value = ', sbad, ', ', sbad, ' ;'
      enddo

c     "data" section

      write(11,'(A)') 'data:'
      write(11,*)
      write(11,'(A,I2.2,A,I2.2,A,I4.4,A)') ' start_date = "', 
     $                                     mo, '/', da, '/', yr, '" ;'
      write(11,*)
      write(11,'(A,I2.2,A,I2.2,A,I4.4,A)') ' end_date = "', 
     $                                     mo, '/', da, '/', yr, '" ;'
      write(11,*)
      write(11,'(A,I2.2,A,I2.2,A,I2.2,A)') ' start_time = "', 
     $                                     hr, ':', mn, ':', se, '" ;'
      write(11,*)
      write(11,'(A,I2.2,A,I2.2,A,I2.2,A)') ' end_time = "', 
     $                                     hr, ':', mn, ':', se, '" ;'
      write(11,*)
      write(11,'(A,I7,A)') ' bad_data_flag = ', int(sbad), ' ;'
      write(11,*)
      write(11,'(A,F12.7,A)') ' grid_latitude = ', glat, ' ;'
      write(11,*)
      write(11,'(A,F12.7,A)') ' grid_longitude = ', glon, ' ;'
      write(11,*)
      write(11,'(A,F12.7,A)') ' grid_altitude = ', galt, ' ;'
      write(11,*)
      write(11,'(A,F12.7,A)') ' radar_latitude = ', rlat, ' ;'
      write(11,*)
      write(11,'(A,F12.7,A)') ' radar_longitude = ', rlon, ' ;'
      write(11,*)
      write(11,'(A,F12.7,A)') ' radar_altitude = ', ralt, ' ;'
      write(11,*)
      write(11,'(A,F12.7,A)') ' x_min = ', xmin, ' ;'
      write(11,*)
      write(11,'(A,F12.7,A)') ' y_min = ', ymin, ' ;'
      write(11,*)
      write(11,'(A,F12.7,A)') ' z_min = ', zmin, ' ;'
      write(11,*)
      write(11,'(A,F12.7,A)') ' x_max = ', xmin+(nx-1)*dx, ' ;'
      write(11,*)
      write(11,'(A,F12.7,A)') ' y_max = ', ymin+(ny-1)*dy, ' ;'
      write(11,*)
      write(11,'(A,F12.7,A)') ' z_max = ', zmin+(nz-1)*dz, ' ;'
      write(11,*)
      write(11,'(A,F12.7,A)') ' x_spacing = ', dx, ' ;'
      write(11,*)
      write(11,'(A,F12.7,A)') ' y_spacing = ', dy, ' ;'
      write(11,*)
      write(11,'(A,F12.7,A)') ' z_spacing = ', dz, ' ;'
      write(11,*)
      write(11,'(A,F12.7,A)') ' u_translation = ', ut, ' ;'
      write(11,*)
      write(11,'(A,F12.7,A)') ' v_translation = ', vt, ' ;'
      write(11,*)

      write(11,'(A)', advance='no') ' x = '
      do i=1, nx
        write(11,'(F15.7)', advance='no') xmin+(i-1)*dx
        if (i.eq.nx) then
          write(11,'(A)') ' ;'
        else
          write(11,'(A)', advance='no') ', '
          if (mod(i,5).eq.4) write(11,*)
        endif
      enddo
      write(11,*)

      write(11,'(A)', advance='no') ' y = '
      do j=1, ny
        write(11,'(F15.7)', advance='no') ymin+(j-1)*dy
        if (j.eq.ny) then
          write(11,'(A)') ' ;'
        else
          write(11,'(A)', advance='no') ', '
          if (mod(j,5).eq.4) write(11,*)
        endif
      enddo
      write(11,*)

      write(11,'(A)', advance='no') ' z = '
      do k=1, nz
        write(11,'(F15.7)', advance='no') zmin+(k-1)*dz
        if (k.eq.nz) then
          write(11,'(A)') ' ;'
        else
          write(11,'(A)', advance='no') ', '
          if (mod(k,5).eq.4) write(11,*)
        endif
      enddo
      write(11,*)

      write(11,'(A)', advance='no') ' el = '
      do k=1, nz
        write(11,'(F15.7)', advance='no') zmin+(k-1)*dz
        if (k.eq.nz) then
          write(11,'(A)') ' ;'
        else
          write(11,'(A)', advance='no') ', '
          if (mod(k,5).eq.4) write(11,*)
        endif
      enddo
      write(11,*)

      write(11,'(A)') ' field_names = '
      do n=1, nfld
        write(11,'(A,A8,A)', advance='no') '  "', fname(n), '"'
        if (n.eq.nfld) then
          write(11,'(A)') ' ;'
        else
          write(11,'(A)') ','
        endif
      enddo

c     Write 3D fields.

      do n=1, nfld
        write(11,*)
        ls = index(fname(n), ' ') - 1
        write(11,'(1X,A,A,A)') ' ', fname(n)(1:ls), ' ='
        s = 0
        do k=1, nz
          do j=1, ny
            do i=1, nx
              write(11,'(F15.7)', advance='no') f(i,j,k,n)
              s = s + 1
              if (s.eq.(nx*ny*nz)) then
                write(11,'(A)') ' ;'
              else
                write(11,'(A)', advance='no') ', '
                if (mod(s,5).eq.4) write(11,*)
              endif
            enddo
          enddo
        enddo
      enddo

c     Write final character and then close file.

      write(11,'(A)') '}'
      close(11)

c     Convert ascii file to netcdf.

      ls = index(ncfile, ' ') - 1
      write(command,*) 'ncgen -o ', ncfile(1:ls), ' ', ncgenf
      write(*,*) 'command = ', command
        
      write(command2,*) 'rm ', ncgenf
      

c
c.... CLZ (5/17/07):  call subroutine system and comment out status test
c
      call system(command)
      
      write(*,*) 'command2 = ', command2
      call system(command2)

c      status = system(command)
c      if (status.ne.0) then
c        write(*,*) 'unable to execute: ', command
c        stop
c      endif

      return
      end










c############################################################################
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                SUBROUTINE SDIN_NETCDF                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
c
c############################################################################
c
c     PURPOSE:
c
c     This subroutine reads gridded single-Doppler fields from the netcdf file.
c
c############################################################################

      subroutine sdin_netcdf(ncfile, nx, ny, nz, nfld, fname, f,
     $                       glon, glat, galt,
     $                       cyr, cmo, cda, chr, cmn, cse,
     $                       dx, dy, dz,
     $                       xmin, ymin, zmin, ut, vt)

      implicit none

      include 'dow.inc'
      include 'netcdf.inc'

c---- input parameters

      character(len=150) ncfile    ! netcdf file name
      integer nx, ny, nz           ! no. of grid points in x, y, and z directions
      integer nfld                 ! number of data fields needed
      character(len=8) fname(nfld) ! field names

c---- returned variables

      real f(nx,ny,nz,nfld)        ! data fields
      real glat, glon              ! latitude and longitude of grid origin (deg)
      real galt                    ! altitude of grid origin (km MSL)
      integer cyr,cmo,cda          ! central date
      integer chr,cmn,cse          ! central time
      real dx, dy, dz              ! grid spacing in x, y, and z directions (km)
      real xmin, ymin, zmin        ! coordinates of lower southwest corner
                                   !   of grid, relative to origin (km)
      real ut, vt                  ! translation velocity (m/s)

c---- local variables

      integer ncid, status, id
      character(len=10) date
      character(len=8) time
      integer i, j, k, n
      real z(nz)


      real dsbad
      parameter (dsbad = 10.0)











c     Open netcdf file.

      status = NF_OPEN(ncfile, NF_NOWRITE, ncid)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem opening file: ', NF_STRERROR(status), ncid
        stop
      endif

c     Read scalar variables.

      status = NF_INQ_VARID(ncid, 'grid_latitude', id)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining id for glat: ',NF_STRERROR(status)
        stop
      endif
      status = NF_GET_VAR_REAL(ncid, id, glat)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining glat: ', NF_STRERROR(status)
        stop
      endif

      status = NF_INQ_VARID(ncid, 'grid_longitude', id)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining id for glon: ',NF_STRERROR(status)
        stop
      endif
      status = NF_GET_VAR_REAL(ncid, id, glon)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining glon: ', NF_STRERROR(status)
        stop
      endif

      status = NF_INQ_VARID(ncid, 'grid_altitude', id)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining id for galt: ',NF_STRERROR(status)
        stop
      endif
      status = NF_GET_VAR_REAL(ncid, id, galt)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining galt: ', NF_STRERROR(status)
        stop
      endif

      status = NF_INQ_VARID(ncid, 'start_date', id)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining id for start_date: ', 
     $             NF_STRERROR(status)
        stop
      endif
      status = NF_GET_VAR_TEXT(ncid, id, date)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining start_date: ', NF_STRERROR(status)
        stop
      endif
      read(date(7:10),*) cyr
      read(date(1:2),*) cmo
      read(date(4:5),*) cda

      status = NF_INQ_VARID(ncid, 'start_time', id)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining id for start_time: ', 
     $             NF_STRERROR(status)
        stop
      endif
      status = NF_GET_VAR_TEXT(ncid, id, time)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining start_time: ', NF_STRERROR(status)
        stop
      endif
      read(time(1:2),*) chr
      read(time(4:5),*) cmn
      read(time(7:8),*) cse

      status = NF_INQ_VARID(ncid, 'x_spacing', id)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining id for dx: ', NF_STRERROR(status)
        stop
      endif
      status = NF_GET_VAR_REAL(ncid, id, dx)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining dx: ', NF_STRERROR(status)
        stop
      endif

      status = NF_INQ_VARID(ncid, 'y_spacing', id)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining id for dy: ', NF_STRERROR(status)
        stop
      endif
      status = NF_GET_VAR_REAL(ncid, id, dy)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining dy: ', NF_STRERROR(status)
        stop
      endif

      status = NF_INQ_VARID(ncid, 'z_spacing', id)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining id for dz: ', NF_STRERROR(status)
        stop
      endif
      status = NF_GET_VAR_REAL(ncid, id, dz)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining dz: ', NF_STRERROR(status)
        stop
      endif

      status = NF_INQ_VARID(ncid, 'z_spacing', id)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining id for dz: ', NF_STRERROR(status)
        stop
      endif
      status = NF_GET_VAR_REAL(ncid, id, dz)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining dz: ', NF_STRERROR(status)
        stop
      endif

      status = NF_INQ_VARID(ncid, 'x_min', id)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining id for xmin: ',NF_STRERROR(status)
        stop
      endif
      status = NF_GET_VAR_REAL(ncid, id, xmin)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining xmin: ', NF_STRERROR(status)
        stop
      endif

      status = NF_INQ_VARID(ncid, 'y_min', id)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining id for ymin: ',NF_STRERROR(status)
        stop
      endif
      status = NF_GET_VAR_REAL(ncid, id, ymin)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining ymin: ', NF_STRERROR(status)
        stop
      endif

      status = NF_INQ_VARID(ncid, 'z', id)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining id for z: ', NF_STRERROR(status)
        stop
      endif
      status = NF_GET_VAR_REAL(ncid, id, z)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining zmin: ', NF_STRERROR(status)
        stop
      endif
      zmin = z(1)

      status = NF_INQ_VARID(ncid, 'u_translation', id)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining id for ut: ', NF_STRERROR(status)
        stop
      endif
      status = NF_GET_VAR_REAL(ncid, id, ut)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining ut: ', NF_STRERROR(status)
        stop
      endif

      status = NF_INQ_VARID(ncid, 'v_translation', id)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining id for vt: ', NF_STRERROR(status)
        stop
      endif
      status = NF_GET_VAR_REAL(ncid, id, vt)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining vt: ', NF_STRERROR(status)
        stop
      endif

c     Read 3D arrays.

      do n=1, nfld
        status = NF_INQ_VARID(ncid, fname(n), id)
        if (status .ne. NF_NOERR) then
          write(*,*) 'Problem obtaining id for ', fname(n),
     $               ': ', NF_STRERROR(status)
          stop
        endif
        status = NF_GET_VAR_REAL(ncid, id, f(1,1,1,n))
        if (status .ne. NF_NOERR) then
          write(*,*) 'Problem obtaining ', fname(n),
     $               ': ', NF_STRERROR(status)
          stop
        endif

        do k=1, nz
          do j=1, ny
            do i=1, nx

c           if (f(i,j,k,n) .lt. (sbad + dsbad) .and.
c     *       f(i,j,k,n) .gt. (sbad - dsbad) ) f(i,j,k,n)=bad

              if (f(i,j,k,n).eq.sbad) f(i,j,k,n)=bad

            enddo
          enddo
        enddo

      enddo

c     Close  netcdf file

      status=NF_CLOSE(ncid)
      if (status.ne.NF_NOERR) then
        write(*,*) 'Error closing file: ', NF_STRERROR(status)
      endif


      return
      end









c############################################################################
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######               SUBROUTINE GET_DIMS_NETCDF             ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
c
c############################################################################
c
c     PURPOSE:
c
c     This subroutine reads gridded fields from the netcdf file.
c
c############################################################################

      subroutine get_dims_netcdf(ncfile, nx, ny, nz)

      implicit none

      include 'netcdf.inc'

c---- input parameters

      character(len=150) ncfile    ! netcdf file name

c---- returned variables

      integer nx, ny, nz           ! no. of grid points in x, y, and z directions

c---- local variables

      integer ncid, status, id

c     Open netcdf file.

      status = NF_OPEN(ncfile, NF_NOWRITE, ncid)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem opening file: ', NF_STRERROR(status), ncid
        stop
      endif

c     Read dimensions.

      status = NF_INQ_DIMID(ncid, 'x', id)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining id for nx: ', NF_STRERROR(status)
        stop
      endif
      status = NF_INQ_DIMLEN(ncid, id, nx)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining nx: ', NF_STRERROR(status)
        stop
      endif

      status = NF_INQ_DIMID(ncid, 'y', id)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining id for ny: ', NF_STRERROR(status)
        stop
      endif
      status = NF_INQ_DIMLEN(ncid, id, ny)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining ny: ', NF_STRERROR(status)
        stop
      endif

      status = NF_INQ_DIMID(ncid, 'z', id)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining id for nz: ', NF_STRERROR(status)
        stop
      endif
      status = NF_INQ_DIMLEN(ncid, id, nz)
      if (status .ne. NF_NOERR) then
        write(*,*) 'Problem obtaining nz: ', NF_STRERROR(status)
        stop
      endif

c     Close  netcdf file

      status=NF_CLOSE(ncid)
      if (status.ne.NF_NOERR) then
        write(*,*) 'Error closing file: ', NF_STRERROR(status)
      endif


      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                   SUBROUTINE SDOUT                   ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine outputs a Cartesian single-Doppler file in
c     ascii format.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  14 December 2001
c
c     Modifications:
c
c
c############################################################################

      subroutine sdout(glon, glat, galt, cyr, cmo, cda, chr, cmn, cse,
     $                 nx, ny, nz, dx, dy, dz, xmin, ymin, zmin, ut, vt,
     $                 nfld, fname, f, count, az, el, time, outfile)

      implicit none

      integer nx, ny, nz            ! no. of grid points in x, y, and z directions
      real xmin, ymin, zmin         ! coordinates of lower southwest corner
                                    !   of grid, relative to origin (km)
      real dx, dy, dz               ! grid spacing in x, y, and z directions (km)
      real glat, glon               ! latitude and longitude of grid origin (deg)
      real galt                     ! altitude of grid origin (km MSL)
      real ut, vt                   ! translation velocity (m/s)
      integer nfld                  ! number of data fields to be gridded
      integer(kind=2) cyr,cmo,cda   ! central date
      integer(kind=2) chr,cmn,cse   ! central time
      character(len=8) fname(nfld)  ! field names
      character(len=8) fnameo
      real f(nx,ny,nz,nfld)         ! data fields
      real az(nx,ny,nz)             ! interpolated azimuth angle (deg)
      real el(nx,ny,nz)             ! interpolated elevation angle (deg)
      real time(nx,ny,nz)           ! time, relative to central time (sec)
      integer count(nx,ny,nz,nfld)  ! no. gates used in interpolation
      integer i, j, k, n            ! loop variables
      character(len=150) outfile     ! output file name


c
c.... begin executable code
c



      open(unit=10, file=outfile, status='unknown')

      write(10,20) glon, glat, galt
 20   format(F9.4, 2x, F9.4, 2x, F6.3)
      write(10,30) cyr, cmo, cda, chr, cmn, cse
 30   format(i4, 2x, i2, 2x, i2, 2x, i2, 2x, i2, 2x, i2)
      write(10,40) nx, ny, nz
 40   format(i4, 2x, i4, 2x, i4)
      write(10,50) dx, dy, dz
 50   format(F7.3, 2x, F7.3, 2x, F7.3)
      write(10,60) xmin, ymin, zmin
 60   format(F9.3, 2x, F9.3, 2x, F9.3)
      write(10,60) ut, vt

c
c.... CLZ (10/7/19):  NOTE that single-radar ascii format is changed here from Dowell version
c

      write(10,*) nfld

c      do n=1, nfld

        write(10,99) (fname(n),n=1,nfld)
 99     format(5(A8))

        do k=1, nz
          do j=1, ny
            do i=1, nx

            write(10,*) (f(i,j,k,n),n=1,nfld)
  70  format(5(2x,f15.6))
c              write(10,*) f(i,j,k,n), '  ', count(i,j,k,n)
c              write(10,70) f(i,j,k,n), count(i,j,k,n)

            enddo
          enddo
        enddo

c      enddo


      close(10)

      return
      end

c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                   SUBROUTINE VDRASOUT                ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine outputs a PPI single-Doppler file in
c     the ascii format that can be read by VDRAS.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  6 May 2002
c
c     Modifications:
c
c
c############################################################################

      subroutine vdrasout(glon, glat, galt,
     $              nx, ny, nswp, dx, dy, xmin, ymin, zmin,
     $              dbz, vr, count, az, el, height, outfile)

      implicit none

      include 'dow.inc'

      integer nx, ny            ! no. of grid points in x and y directions
      integer nswp              ! number of radar sweeps
      real xmin, ymin, zmin     ! coordinates of lower southwest corner
                                !   of grid, relative to origin (km)
      real dx, dy               ! grid spacing in x and y directions (km)
      real glat, glon           ! latitude and longitude of grid origin (deg)
      real galt                 ! altitude of grid origin (km MSL)
      real dbz(nx,ny,nswp)      ! interpolated reflectivity (dBZ)
      real vr(nx,ny,nswp)       ! interpolated radial velocity (m/s)
      real az(nx,ny,nswp)       ! interpolated azimuth angle (deg)
      real el(nx,ny,nswp)       ! interpolated elevation angle (deg)
      real height(nx,ny,nswp)   ! interpolated height of obs (km)
      integer count(nx,ny,nswp) ! no. gates used in interpolation
      real q(nx,ny,nswp)        ! output data values
      integer i, j, k           ! loop variables
      character(len=80) outfile ! output file name
      real outbad               ! bad data flag for output file
      parameter(outbad=-1000.0)


      open(unit=10, file=outfile, status='unknown')

      write(10,30) glon, glat, galt, xmin, ymin
c      write(10,30) glon, glat, 0.0, xmin, ymin
      write(10,31) zmin, dx, dy, 0.0, nx
      write(10,32) ny, nswp, outbad

 30   format(f9.3, 2x, f9.3, 2x, f9.3, 2x, f9.3, 2x, f9.3)
 31   format(f9.3, 2x, f9.3, 2x, f9.3, 2x, f9.3, 2x, i9)
 32   format(i9, 2x, i9, 2x, e11.4)

 70   format(e11.4, 2x, e11.4, 2x, e11.4, 2x, e11.4, 2x, e11.4)

      do k=1, nswp
        do j=1, ny
          do i=1, nx
            if (dbz(i,j,k).eq.bad) then
              q(i,j,k)=outbad
            else
              q(i,j,k)=dbz(i,j,k)
            endif
          enddo
        enddo
      enddo
      write(10,70) (((q(i,j,k),i=1,nx),j=1,ny),k=1,nswp)

      do k=1, nswp
        do j=1, ny
          do i=1, nx
            if (vr(i,j,k).eq.bad) then
              q(i,j,k)=outbad
            else
              q(i,j,k)=vr(i,j,k)
            endif
          enddo
        enddo
      enddo
      write(10,70) (((q(i,j,k),i=1,nx),j=1,ny),k=1,nswp)

      do k=1, nswp
        do j=1, ny
          do i=1, nx
            q(i,j,k)=count(i,j,k)
          enddo
        enddo
      enddo
      write(10,70) (((q(i,j,k),i=1,nx),j=1,ny),k=1,nswp)

      do k=1, nswp
        do j=1, ny
          do i=1, nx
            if (az(i,j,k).eq.bad) then
              q(i,j,k)=outbad
            else
              q(i,j,k)=az(i,j,k)
            endif
          enddo
        enddo
      enddo
      write(10,70) (((q(i,j,k),i=1,nx),j=1,ny),k=1,nswp)

      do k=1, nswp
        do j=1, ny
          do i=1, nx
            if (el(i,j,k).eq.bad) then
              q(i,j,k)=outbad
            else
              q(i,j,k)=el(i,j,k)
            endif
          enddo
        enddo
      enddo
      write(10,70) (((q(i,j,k),i=1,nx),j=1,ny),k=1,nswp)

      do k=1, nswp
        do j=1, ny
          do i=1, nx
            if (height(i,j,k).eq.bad) then
              q(i,j,k)=outbad
            else
              q(i,j,k)=height(i,j,k)
            endif
          enddo
        enddo
      enddo
      write(10,70) (((q(i,j,k),i=1,nx),j=1,ny),k=1,nswp)

      close(10)

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                   SUBROUTINE DARTOUT                 ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine outputs a PPI single-Doppler file in
c     the ascii format that can be read by DART.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  23 March 2005
c
c     Modifications:
c
c
c############################################################################

      subroutine DARTout(glond, glatd, galt, rlond, rlatd, ralt,
     $                   map_proj, nx, ny, nswp, dx, dy, xmin, ymin,
     $                   f, nfld, az, el, height, outfile, days, secs)

      implicit none

      include 'dow.inc'

c---- Passed variables

      real glatd, glond         ! latitude and longitude of grid origin (deg)
      real galt                 ! altitude of grid origin (km MSL)
      real rlatd, rlond         ! latitude and longitude of radar (deg)
      real ralt                 ! altitude of radar (km MSL)
      integer map_proj          ! map projection (for relating lat, lon to x, y):
                                !   0 = flat earth
                                !   1 = oblique azimuthal
                                !   2 = Lambert conformal
      integer nx, ny            ! no. of grid points in x and y directions
      integer nswp              ! number of radar sweeps
      real dx, dy               ! grid spacing in x and y directions (km)
      real xmin, ymin           ! coordinates of southwest corner
                                !   of grid, relative to origin (km)
      real f(nx,ny,nswp,2)      ! interpolated reflectivity (dBZ) and radial velocity (m/s)
      integer nfld              ! number of valid fields in f
      real az(nx,ny,nswp)       ! interpolated azimuth angle (deg)
      real el(nx,ny,nswp)       ! interpolated elevation angle (deg)
      real height(nx,ny,nswp)   ! interpolated height of obs (km)
      character(len=80) outfile ! output file name
      integer secs(nswp), days(nswp)

c---- Local variables

      integer i, j, k, n                               ! loop variables
      integer fi; parameter(fi=10)                     ! file unit number
      integer dbz_index, vr_index
      integer num_copies; parameter(num_copies=1)
      integer num_qc; parameter(num_qc=0)
      integer num_obs
      integer max_num_obs
      character(len=129) copy_meta_data(num_copies)
      character(len=129) qc_meta_data(num_qc)
      integer first_time
      integer last_time
      real(kind=8) values(num_copies)
      real(kind=8) qc(num_qc)
      integer obs_kind                                 ! DART observation type:
                                                       !   100 = Doppler velocity
                                                       !   101 = radar reflectivity
      real glat, glon                                  ! grid lat and lon (rad)
      real(kind=8) rlat, rlon, rheight                 ! radar lat, lon (rad) and height (m)
      real lat, lon                                    ! lat and lon (rad)
      real(kind=8) olat, olon, oheight                 ! observation lat, lon (rad) and height (m)
      real x, y                                        ! observation coordinates (km)
      real(kind=8) error_variance


      glat = dtor*glatd
      glon = dtor*glond
      rlat = dtor*rlatd
      rlon = dtor*rlond
      rheight = 1000.0*ralt

      if (nfld.le.0) then
        write(*,*) 'error:  nfld = ', nfld
        stop
      else if (nfld.eq.1) then
        dbz_index = 1
        vr_index = 0
      else
        dbz_index = 2
        vr_index = 1
      endif

c     Count number of valid observations.

      num_obs = 0

      do k=1, nswp
        do j=1, ny
          do i=1, nx
            if ( (f(i,j,k,dbz_index).ne.bad) .and. (az(i,j,k).ne.bad)
     $          .and. (el(i,j,k).ne.bad) .and. (height(i,j,k).ne.bad) 
     $         ) then
              num_obs = num_obs + 1
            endif
            if ( (vr_index.ne.0) .and.
     $           (f(i,j,k,vr_index).ne.bad) .and. (az(i,j,k).ne.bad)
     $          .and. (el(i,j,k).ne.bad) .and. (height(i,j,k).ne.bad) 
     $           ) then
              num_obs = num_obs + 1
            endif
          enddo
        enddo
      enddo

      max_num_obs = num_obs

c     Write header information.

      open(unit=fi, file=outfile, status='unknown')

      write(fi,*) ' num_copies: ', num_copies, ' num_qc: ', num_qc
      write(fi,*) ' num_obs: ', num_obs,    ' max_num_obs: ',max_num_obs

      copy_meta_data(1) = 'observations'
      do n=1, num_copies
        write(fi, '(a129)') copy_meta_data(n)
      enddo

      do n=1, num_qc
        write(fi, '(a129)') qc_meta_data(n)
      enddo

      first_time = 1
      last_time = max_num_obs
      write(fi,*) ' first: ', first_time, ' last: ', last_time

c     Write observations.

      num_obs = 0

      do k=1, nswp
        do j=1, ny
          do i=1, nx
            x = xmin + dx*(i-1.0)
            y = ymin + dy*(j-1.0)
            call xy_to_ll(lat, lon, map_proj, x, y, glat, glon)
            olat = lat
            olon = lon
            oheight = 1000.0 * (galt + height(i,j,k))
            if ( (f(i,j,k,dbz_index).ne.bad) .and. (az(i,j,k).ne.bad)
     $           .and. (el(i,j,k).ne.bad) .and. (height(i,j,k).ne.bad) 
     $           ) then
              num_obs = num_obs + 1
              values(1) = f(i,j,k,dbz_index)
              obs_kind = 101
              error_variance = 2.0*2.0
              call write_DART_ob(fi, num_obs, max_num_obs, values, 
     $                           num_copies, qc, num_qc, olat, olon, 
     $                           oheight, az(i,j,k), el(i,j,k), rlat, 
     $                           rlon, rheight, obs_kind, secs(k), 
     $                           days(k), error_variance)
            endif
            if ( (vr_index.ne.0) .and.
     $           (f(i,j,k,vr_index).ne.bad) .and. (az(i,j,k).ne.bad)
     $           .and. (el(i,j,k).ne.bad) .and. (height(i,j,k).ne.bad) 
     $            ) then
              num_obs = num_obs + 1
              values(1) = f(i,j,k,vr_index)
              obs_kind = 100
              error_variance = 2.0*2.0
              call write_DART_ob(fi, num_obs, max_num_obs, values, 
     $                           num_copies, qc, num_qc, olat, olon, 
     $                           oheight, az(i,j,k), el(i,j,k), rlat, 
     $                           rlon, rheight, obs_kind, secs(k), 
     $                           days(k), error_variance)
            endif
          enddo
        enddo
      enddo

      close(fi)

      return
      end

c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                SUBROUTINE WRITE_DART_OB              ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine writes a single observation in DART format to the
c     output file
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  23 March 2005
c
c     Modifications:
c
c
c############################################################################

      subroutine write_DART_ob(fi, o, max_num_obs, values, num_copies, 
     $                         qc, num_qc, olat, olon, oheight, az, el, 
     $                         rlat, rlon, rheight, obs_kind, secs, 
     $                         days, error_variance)

      implicit none

      include 'dow.inc'

c---- Passed variables

      integer fi                            ! file unit number
      integer o                             ! current observation number
      integer max_num_obs
      integer num_copies
      real(kind=8) values(num_copies)
      integer num_qc
      real(kind=8) qc(num_qc)
      real(kind=8) olat, olon               ! observation lat and lon (rad)
      real(kind=8) oheight                  ! observation height (m)
      real az                               ! azimuth angle (deg)
      real el                               ! elevation angle (deg)
      real(kind=8) rlat, rlon               ! radar lat and lon (rad)
      real(kind=8) rheight                  ! radar height (m)
      integer obs_kind                      ! DART observation type:
                                            !   100 = Doppler velocity
                                            !   101 = radar reflectivity
      integer secs, days
      real(kind=8) error_variance

c---- Local variables

      integer n
      integer cov_group; parameter(cov_group=-1)
      integer prev_time, next_time
      integer which_vert; parameter(which_vert=3)
      real(kind=8) dir(3)


      write(fi,*) 'OBS ', o
      do n=1, num_copies
        write(fi,*) values(n)
      enddo
      do n=1, num_qc
        write(fi,*) qc(n)
      enddo

      if (o.eq.1) then
        prev_time = -1
      else
        prev_time = o-1
      endif
      if (o.eq.max_num_obs) then
        next_time = -1
      else
        next_time = o+1
      endif
      write(fi,*) prev_time, next_time, cov_group

      write(fi,11)
 11   format('obdef')

c     write_location

      write(fi, '(''loc3d'')' ) 
      write(fi,*) olon, olat, oheight, which_vert

c     write_kind

      write(fi, '(''kind'')' )
      write(fi,*) obs_kind

c     write_time

      write(fi,*) secs, days

c     error_variance

      write(fi,*) error_variance

c     write_platform

      dir(1) = sin(dtor*az)*cos(dtor*el)
      dir(2) = cos(dtor*az)*cos(dtor*el)
      dir(3) = sin(dtor*el)
      write(fi, '(''platform'')' )
      write(fi, '(''loc3d'')' ) 
      write(fi,*) rlon, rlat, rheight, which_vert
      write(fi, '(''dir3d'')' )
      write(fi,*) dir(1), dir(2), dir(3)

      return
      end










c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                   SUBROUTINE NCASCII                 ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine outputs an ascii file containing a particular
c     field in a raw sweep.  The output file is in a format that
c     can be converted to netcdf with the "ncgen" command.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  16 January 2003
c
c     Modifications:
c
c
c############################################################################

      subroutine ncascii(outfile, fieldnum, timestamp,
     $               vold,radd,celv,cfac,parm,swib,ryib,asib,rdat,
     $               rlon, rlat, ralt, firstray, lastray)

      implicit none

      include 'dow.inc'
      include 'structures.inc'

      character(len=150) outfile    ! output netcdf ascii file name
      integer fieldnum             ! type of output field:
                                   !    1 = radial velocity (m/s)
                                   !    2 = reflectivity (dBZ)
      character(len=10) timestamp  ! time stamp (seconds since 1/1/70)
      type(vold_info)                     :: vold
      type(radd_info)                     :: radd
      type(celv_info)                     :: celv
      type(cfac_info)                     :: cfac
      type(parm_info), dimension(maxflds) :: parm
      type(swib_info)                     :: swib
      type(ryib_info), dimension(maxrays) :: ryib
      type(asib_info), dimension(maxrays) :: asib
      type(rdat_info), dimension(maxrays) :: rdat
      real rlat, rlon          ! radar latitude and longitude (deg)
      real ralt                ! radar altitude (km MSL)
      integer r                ! ray number
      integer g                ! gate number
      real elavg               ! average elevation angle
      real bw                  ! beam width
      integer gl               ! gate length
      integer firstray         ! number of first ray in output ascii file
      integer lastray          ! number of last ray in output ascii file
      integer nr               ! number of rays
      integer dir              ! preferred order of data:
                               !   1=forward,  -1=backward 
      real el(maxrays)         ! elevation angles
      real az(maxrays)         ! azimuth angles
      real data(maxrays, 3000) ! radar data
      real MISSING             ! missing data flag in netcdf file
      parameter(MISSING=-99900.)
      real RFOLD               ! range folding flag in netcdf file
      parameter(RFOLD=-99901.)



      real dsbad
      parameter (dsbad = 10.0)













c     do a little processing of the data

      if (lastray.gt.swib%num_rays) then
        lastray = swib%num_rays
      endif
      nr = 1 + iabs(lastray-firstray)
      dir = isign(1, lastray-firstray)

      elavg = 0.0
      do r=1, nr

        el(r) = ryib(firstray + dir*(r-1))%elevation
        if (el(r) .gt. 180.0) then
          el(r) = el(r) - 360.0
        endif
        elavg = elavg + el(r)

        az(r) = ryib(firstray + dir*(r-1))%azimuth
        if (az(r) .gt. 360.0) then
          az(r) = az(r) - 360.0
        endif
        if (az(r) .lt. 0.0) then
          az(r) = az(r) + 360.0
        endif


        do g=1, celv%total_gates

          data(r,g) = rdat(firstray + dir*(r-1))%data(g)

c          if (data(r,g) .lt. (sbad + dsbad) .and.
c     *       data(r,g) .gt. (sbad - dsbad) ) then

          if (data(r,g) .eq. sbad) then

            data(r,g) = MISSING
          endif

        enddo


      enddo



      elavg = elavg / nr

c     create the ascii file

      open(unit=11, file=outfile, status='unknown')

c     header

      write(11,100) outfile(1:15)
 100  format('netcdf ', A15, ' {')
      write(11,101)
 101  format('dimensions:')
      write(11,102) nr
 102  format(T9, 'Azimuth =', I4, ' ;')
      write(11,103) celv%total_gates
 103  format(T9, 'Gate =', I4, ' ;')
      write(11,104)
 104  format('variables:')
      write(11,105)
 105  format(T9, 'float Azimuth(Azimuth) ;')
      write(11,106)
 106  format(T17, 'Azimuth:Units = "Degrees" ;')
      write(11,107)
 107  format(T9, 'float BeamWidth(Azimuth) ;')
      write(11,108)
 108  format(T17, 'BeamWidth:Units = "Degrees" ;')
      write(11,109)
 109  format(T9, 'float GateWidth(Azimuth) ;')
      write(11,110)
 110  format(T17, 'GateWidth:Units = "Meters" ;')
      if (fieldnum.eq.1) then
        write(11,111)
        write(11,112)
 111    format(T9, 'float Velocity(Azimuth, Gate) ;')
 112    format(T17, 'Velocity:Units = "MetersPerSecond" ;')
      else if (fieldnum.eq.2) then
        write(11,113)
        write(11,114)
 113    format(T9, 'float Reflectivity(Azimuth, Gate) ;')
 114    format(T17, 'Reflectivity:Units = "dBZ" ;')
      else
        write(*,*) 'invalid field number'
        write(*,*) 'fieldnum = ', fieldnum
        stop
      endif
      write(11,*)

      write(11,200)
 200  format('// global attributes:')
      if (fieldnum.eq.1) then
        write(11,201)
 201    format(T17, ':TypeName = "Velocity" ;')
      else if (fieldnum.eq.2) then
        write(11,202)
 202    format(T17, ':TypeName = "Reflectivity" ;')
      endif
      write(11,203)
 203  format(T17, ':DataType = "RadialSet" ;')
      write(11,204) rlat
 204  format(T17, ':Latitude = ', F7.4, ' ;')
      write(11,205) rlon
 205  format(T17, ':Longitude = ', F8.4, ' ;')
      write(11,206) ralt*1000.0
 206  format(T17, ':Height =', F5.0, ' ;')
      write(11,207) timestamp
 207  format(T17, ':Time = ', A10, ' ;')
      write(11,208)
 208  format(T17, ':FractionalTime = 0.0 ;')
      write(11,209)
 209  format(T17, ':attributes = " Nyquist_Vel radarName vcp" ;')
      write(11,210)
 210  format(T17, ':Nyquist_Vel-unit = "MetersPerSecond" ;')
      write(11,211) radd%unambig_vel
 211  format(T17, ':Nyquist_Vel-value = "', F5.2, '" ;')
      write(11,212)
 212  format(T17, ':radarName-unit = "dimensionless" ;')
      write(11,213)
 213  format(T17, ':radarName-value = "DOW" ;')
      write(11,214)
 214  format(T17, ':vcp-unit = "dimensionless" ;')
      write(11,215)
 215  format(T17, ':vcp-value = "11" ;')

c      elavg = 0.8

      write(11,216) elavg
 216  format(T17, ':Elevation = ', F4.1, ' ;')
      write(11,217)
 217  format(T17, ':ElevationUnits = "Degrees" ;')
      write(11,218) celv%gate_spacing(1)   ! range to near side of first gate
 218  format(T17, ':RangeToFirstGate = ', F6.1 , ' ;')
      write(11,219)
 219  format(T17, ':RangeToFirstGateUnits = "Meters" ;')
      write(11,220) MISSING
 220  format(T17, ':MissingData = ', F7.0, ' ;')
      write(11,221) RFOLD
 221  format(T17, ':RangeFolded = ', F7.0, ' ;')
      write(11,222)
 222  format('data:')
      write(11,*)

c     azimuth angles

      write(11,223) az(1), az(2), az(3), az(4), az(5)
 223  format(' Azimuth = ', F8.4, ', ', F8.4, ', ', F8.4, ', ',
     $       F8.4, ', ', F8.4, ',')

      do r=6, 5 * int((nr-1)/5) - 4, 5
        write(11,224) az(r), az(r+1), az(r+2), az(r+3), az(r+4)
      enddo
 224  format('    ', F8.4, ', ', F8.4, ', ', F8.4, ', ',
     $       F8.4, ', ', F8.4, ',')

      r = 5 + 5 * int((nr-1)/5) - 4
      if ( (nr-r) .eq. 0 ) then
        write(11,225) az(r)
      else if ( (nr-r) .eq. 1 ) then
        write(11,226) az(r), az(r+1)
      else if ( (nr-r) .eq. 2 ) then
        write(11,227) az(r), az(r+1), az(r+2)
      else if ( (nr-r) .eq. 3 ) then
        write(11,228) az(r), az(r+1), az(r+2), az(r+3)
      else
        write(11,229) az(r), az(r+1), az(r+2), az(r+3), az(r+4)
      endif
 225  format('    ', F8.4, ' ;')
 226  format('    ', F8.4, ', ', F8.4, ' ;')
 227  format('    ', F8.4, ', ', F8.4, ', ', F8.4, ' ;')
 228  format('    ', F8.4, ', ', F8.4, ', ', F8.4, ', ',
     $       F8.4, ' ;')
 229  format('    ', F8.4, ', ', F8.4, ', ', F8.4, ', ',
     $       F8.4, ', ', F8.4, ' ;')
      write(11,*)

c     beamwidths

      bw = radd%horiz_beam_width

      write(11,233) bw, bw, bw, bw, bw
 233  format(' BeamWidth = ', F8.6, ', ', F8.6, ', ', F8.6, ', ',
     $       F8.6, ', ', F8.6, ',')

      do r=6, 5 * int((nr-1)/5) - 4, 5
        write(11,234) bw, bw, bw, bw, bw
      enddo
 234  format('    ', F8.6, ', ', F8.6, ', ', F8.6, ', ',
     $       F8.6, ', ', F8.6, ',')

      r = 5 + 5 * int((nr-1)/5) - 4
      if ( (nr-r) .eq. 0 ) then
        write(11,235) bw
      else if ( (nr-r) .eq. 1 ) then
        write(11,236) bw, bw
      else if ( (nr-r) .eq. 2 ) then
        write(11,237) bw, bw, bw
      else if ( (nr-r) .eq. 3 ) then
        write(11,238) bw, bw, bw, bw
      else
        write(11,239) bw, bw, bw, bw, bw
      endif
 235  format('    ', F8.6, ' ;')
 236  format('    ', F8.6, ', ', F8.6, ' ;')
 237  format('    ', F8.6, ', ', F8.6, ', ', F8.6, ' ;')
 238  format('    ', F8.6, ', ', F8.6, ', ', F8.6, ', ',
     $       F8.6, ' ;')
 239  format('    ', F8.6, ', ', F8.6, ', ', F8.6, ', ',
     $       F8.6, ', ', F8.6, ' ;')
      write(11,*)

c     gate length

      gl = nint(celv%gate_spacing(2)-celv%gate_spacing(1))

      write(11,243) gl, gl, gl, gl, gl
 243  format(' GateWidth = ', I4, ', ', I4, ', ', I4, ', ',
     $       I4, ', ', I4, ',')

      do r=6, 5 * int((nr-1)/5) - 4, 5
        write(11,244) gl, gl, gl, gl, gl
      enddo
 244  format('    ', I4, ', ', I4, ', ', I4, ', ',
     $       I4, ', ', I4, ',')

      r = 5 + 5 * int((nr-1)/5) - 4
      if ( (nr-r) .eq. 0 ) then
        write(11,245) gl
      else if ( (nr-r) .eq. 1 ) then
        write(11,246) gl, gl
      else if ( (nr-r) .eq. 2 ) then
        write(11,247) gl, gl, gl
      else if ( (nr-r) .eq. 3 ) then
        write(11,248) gl, gl, gl, gl
      else
        write(11,249) gl, gl, gl, gl, gl
      endif
 245  format('    ', I4, ' ;')
 246  format('    ', I4, ', ', I4, ' ;')
 247  format('    ', I4, ', ', I4, ', ', I4, ' ;')
 248  format('    ', I4, ', ', I4, ', ', I4, ', ',
     $       I4, ' ;')
 249  format('    ', I4, ', ', I4, ', ', I4, ', ',
     $       I4, ', ', I4, ' ;')
      write(11,*)

c     finally, the data

      if (fieldnum.eq.1) then
        write(11,300)
      else if (fieldnum.eq.2) then
        write(11,301)
      endif
 300  format(' Velocity =')
 301  format(' Reflectivity =')

      do r=1, nr

        do g=1, 5 * int((celv%total_gates-1)/5) - 4, 5
          write(11,310) data(r,g), data(r,g+1),
     $                  data(r,g+2), data(r,g+3),
     $                  data(r,g+4)
        enddo
 310    format('  ', F8.1, ', ', F8.1, ', ', F8.1, ', ',
     $         F8.1, ', ', F8.1, ',')

        g = 5 + 5 * int((celv%total_gates-1)/5) - 4
        if (r.eq.nr) then
          if ( (celv%total_gates-g) .eq. 0 ) then
            write(11,330) data(r,g)
          else if ( (celv%total_gates-g) .eq. 1 ) then
            write(11,331) data(r,g), data(r,g+1)
          else if ( (celv%total_gates-g) .eq. 2 ) then
            write(11,332) data(r,g), data(r,g+1),
     $                    data(r,g+2)
          else if ( (celv%total_gates-g) .eq. 3 ) then
            write(11,333) data(r,g), data(r,g+1),
     $                    data(r,g+2), data(r,g+3)
          else
            write(11,334) data(r,g), data(r,g+1),
     $                    data(r,g+2), data(r,g+3),
     $                    data(r,g+4)
          endif
        else
          if ( (celv%total_gates-g) .eq. 0 ) then
            write(11,320) data(r,g)
          else if ( (celv%total_gates-g) .eq. 1 ) then
            write(11,321) data(r,g), data(r,g+1)
          else if ( (celv%total_gates-g) .eq. 2 ) then
            write(11,322) data(r,g), data(r,g+1),
     $                    data(r,g+2)
          else if ( (celv%total_gates-g) .eq. 3 ) then
            write(11,323) data(r,g), data(r,g+1),
     $                    data(r,g+2), data(r,g+3)
          else
            write(11,324) data(r,g), data(r,g+1),
     $                    data(r,g+2), data(r,g+3),
     $                    data(r,g+4)
          endif
        endif
 320    format('  ', F8.1, ',')
 321    format('  ', F8.1, ', ', F8.1, ',')
 322    format('  ', F8.1, ', ', F8.1, ', ', F8.1, ',')
 323    format('  ', F8.1, ', ', F8.1, ', ', F8.1, ', ',
     $         F8.1, ',')
 324    format('  ', F8.1, ', ', F8.1, ', ', F8.1, ', ',
     $         F8.1, ', ', F8.1, ',')
 330    format('  ', F8.1, ' ;')
 331    format('  ', F8.1, ', ', F8.1, ' ;')
 332    format('  ', F8.1, ', ', F8.1, ', ', F8.1, ' ;')
 333    format('  ', F8.1, ', ', F8.1, ', ', F8.1, ', ',
     $         F8.1, ' ;')
 334    format('  ', F8.1, ', ', F8.1, ', ', F8.1, ', ',
     $         F8.1, ', ', F8.1, ' ;')

      enddo

      write(11,400)
 400  format('}')

      close(11)

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                   SUBROUTINE SDIN1                   ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine reads in a Cartesian single-Doppler volume
c     (ascii format).
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  14 December 2001
c
c     Modifications:
c
c
c############################################################################

      subroutine sdin1(infile, glon, glat, galt,
     $                 cyr, cmo, cda, chr, cmn, cse,
     $                 nx, ny, nz, dx, dy, dz,
     $                 xmin, ymin, zmin, ut, vt,
     $                 nfld, fname, f)


      implicit none

      integer nx, ny, nz       ! no. of grid points in x, y, and z directions
      real xmin, ymin, zmin    ! coordinates of lower southwest corner
                               !   of grid, relative to origin (km)
      real dx, dy, dz          ! grid spacing in x, y, and z directions (km)
      real glat, glon          ! latitude and longitude of grid origin (deg)
      real galt                ! altitude of grid origin (km MSL)
      real ut, vt              ! translation velocity (m/s)
      integer nfldi            ! number of data fields in input file
      character(len=8) fnamei  ! field name in input file
      integer cyr,cmo,cda      ! central date
      integer chr,cmn,cse      ! central time
      integer nfld             ! number of data fields needed
      character(len=8) fname(nfld)  ! field names
      real f(nx,ny,nz,nfld)    ! data fields
      integer fi               ! index of data field
      integer i, j, k, n       ! loop variables
      character(len=150) infile ! input file name
      real dummy

      write(*,*) 'opening ', infile

      open(unit=12, file=infile, status='old')

      read(12,*) glon, glat, galt
      read(12,*) cyr, cmo, cda, chr, cmn, cse
      read(12,*) nx, ny, nz
      read(12,*) dx, dy, dz
      read(12,*) xmin, ymin, zmin
      read(12,*) ut, vt
      read(12,*) nfldi

      do n=1, nfldi

        read(12,*) fnamei
        fi = 0
        do i=1, nfld
          if (fnamei .eq. fname(i)) then
            fi=i
            write(*,*) 'reading field ', fname(i)
          endif
        enddo

        do k=1, nz
          do j=1, ny
            do i=1, nx
              if (fi.eq.0) then
                read(12,*) dummy
              else
                read(12,*) f(i,j,k,fi)
              endif
            enddo
          enddo
        enddo

      enddo

      close(12)

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                   SUBROUTINE SDIN2                   ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine reads in a Cartesian single-Doppler volume
c     (VDRAS ascii format).
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  21 February 2002
c
c     Modifications:
c
c
c############################################################################

      subroutine sdin2(infile, nx, ny, nz, dx, dy, dz,
     $                 xminr, yminr, zminr, vr, rf, az, el)


      implicit none

      include 'dow.inc'

      integer nx, ny, nz       ! no. of grid points in x, y, and z directions
      real xminr, yminr, zminr ! coordinates of lower southwest corner
                               !   of grid, relative to radar (km)
      real dx, dy, dz          ! grid spacing in x, y, and z directions (km)
      integer i, j, k          ! loop variables
      character(len=150) infile ! input file name

      real vr(nx,ny,nz)        ! radial velocity (m/s)
      real rf(nx,ny,nz)        ! reflectivity (dBZ)
      real az(nx,ny,nz)        ! azimuth angle (deg)
      real el(nx,ny,nz)        ! elevation angle (deg)
      real dummy               ! dummy value
      real spval               ! bad data flag
      real x, y, z             ! coords. of grid point with respect to radar
      real r                   ! range to grid point
      real comptoaz            ! function used by this subroutine

      write(*,*) 'opening ', infile

      open(unit=12, file=infile, status='old')

      read(12,*) dummy, dummy, dummy, xminr, yminr
      read(12,*) zminr, dx, dy, dz, nx
      read(12,*) ny, nz, spval

      read(12,99) (((rf(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      read(12,99) (((vr(i,j,k),i=1,nx),j=1,ny),k=1,nz)
 99   format(e11.4, 2x, e11.4, 2x, e11.4, 2x, e11.4, 2x, e11.4)

      close(12)

      do k=1, nz
        do j=1, ny
          do i=1, nx
            if (rf(i,j,k).eq.spval) then
              rf(i,j,k) = bad
            endif
            if (vr(i,j,k).eq.spval) then
              vr(i,j,k) = bad
            endif
          enddo
        enddo
      enddo

c     Compute azimuth and elevation angles.

      do k=1, nz
        z = zminr + (k-1.0)*dz
        do j=1, ny
          y = yminr + (j-1.0)*dy
          do i=1, nx
            x = xminr + (i-1.0)*dx
            r = sqrt(x*x + y*y + z*z)
            if (r.eq.0.0) then
              az(i,j,k) = bad
            else
              az(i,j,k) = comptoaz(x, y)
            endif

            el(i,j,k) = rtod*asin(z/r)

          enddo
        enddo
      enddo

      return
      end

