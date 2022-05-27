      subroutine radfill1d(f, fld1d, ngrdtot, nx, ny, nz,
     * dx, dy, dz, xmin, ymin, zmin,
     * bad, ifilter)

      real bad
      real envavgtop
      real f(nx,ny,nz,ngrdtot)
      real fld1d(nz)

c
c.... work arrays for 1-D scanning fill values
c

      real fld1d_up(nz)
      real fld1d_dn(nz)
      integer kup(nz)
      integer kdn(nz)

      integer idb, jdb, kdb
      common /debug/ idb, jdb, kdb 




c
c.... begin executable code
c




c
c.... CLZ (11/5/2019): maxhole is the maximum allowable levels to fill a data hole vertically
c

      maxhole = 10

      nxm1 = nx - 1
      nym1 = ny - 1
      nzm1 = nz - 1
      nxm2 = nx - 2
      nym2 = ny - 2


      rdx = 1.0 /(2000.0 * dx)
      rdy = 1.0 /(2000.0 * dy)
      rdz = 1.0 /(2000.0 * dz)
      r1dz = 1.0 /(1000.0 * dz)
      
      write(6,*)
      write(6,*) 'entered radfill1d'
      write(6,*)


      nxmid = real(nx) / 2.0
      write(6,*)
      write(6,*) 'nx=',nx,' nxmid=',nxmid
      write(6,*)



c
c.... CLZ (11/6/19): scan through the total of all single-radar data and factor fields
c

      do n = 1, ngrdtot

      write(6,*)'in radfill1d, field n=',n

c
c.... first set the grid row
c

      do j = 1, ny

c
c.... then set the column
c

      do i = 1, nx

c
c.... copy the current field column into a 1-D work array
c

      do k = 1, nz
      fld1d(k) = f(i,j,k,n)
      enddo

c
c.... set work arrays for the current column to missing value
c

      do k = 1, nz
      fld1d_up(k) = bad
      fld1d_dn(k) = bad
      kup(k) = 0
      kdn(k) = 0
      enddo

c
c.... scan upward, placing closest good value in work array at missing location
c

      fldup_save = bad

      do k = 2, nz

      if (fld1d(k) .eq. bad .and. fld1d(k-1) .ne. bad) then
      fldup_save = fld1d(k-1)
      kl_save = k-1
      endif

      if (fld1d(k) .eq. bad .and. fldup_save .ne. bad) then
      fld1d_up(k) = fldup_save
      kup(k) = kl_save
      endif

      enddo

c
c.... scan downward, placing closest good value in work array at missing location
c

      flddown_save = bad

      do k = nz - 1, 1, -1

      if (fld1d(k) .eq. bad .and. fld1d(k+1) .ne. bad) then
      flddown_save = fld1d(k+1)
      kh_save = k+1
      endif

      if (fld1d(k) .eq. bad .and. flddown_save .ne. bad) then
      fld1d_dn(k) = flddown_save
      kdn(k) = kh_save
      endif

      enddo

c
c.... fill holes in this column using non-missing values from upward and downward passes
c

      do k = 1, nz

      if (fld1d(k) .eq. bad) then

      if (fld1d_up(k) .ne. bad .and. fld1d_dn(k) .ne. bad) then
      rkup = real(kup(k))
      rkdn = real(kdn(k))
      delk = rkdn - rkup

c
c.... CLZ (11/5/2019): only take fill action if missing layer is less than maxhole deep
c

      if ((kdn(k) - kup(k)) .le. maxhole) then

      rk = real(k)
      wgtup = (rkdn - rk) / delk

      if (wgtup .gt. 1.0) then
      wgtup = 1.0
      elseif (wgtup .lt. 0.0) then
      wgtup = 0.0
      endif

      wgtdn = 1.0 - wgtup
      fld1d(k) = (wgtup * fld1d_up(k)) + (wgtdn * fld1d_dn(k))

c
c.... endif ((kdn(k) - kup(k)) .le. maxhole) then
c

      endif

c
c.... endif (fld1d_up(k) .ne. bad .and. fld1d_dn(k) .ne. bad) then
c

      endif

c
c.... endif (fld1d(k) .eq. bad) then
c

      endif

c
c.... enddo k = 1, nz
c

      enddo

c
c.... finally, copy 2-D work arrays back to 4-D field array f
c

      do k = 1, nz
      f(i,j,k,n) = fld1d(k)
      enddo

c
c.... enddo 1 = 1, nx
c

      enddo

c
c.... enddo j = 1, ny
c

      enddo

c
c.... enddo n = 1, nfld
c

      enddo

c
c.... done with radfill
c

      return
      end