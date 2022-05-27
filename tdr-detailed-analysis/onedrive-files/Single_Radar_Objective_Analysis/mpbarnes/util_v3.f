c############################################################################
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                 SUBROUTINE XDERIV                    ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
c
c############################################################################
c
c     PURPOSE:
c
c     This subroutine computes the derivative in the x direction of
c     a gridded quantity at each grid point.  At the sides, derivatives
c     are first order, uncentered.  Within the interior, derivatives
c     are second order, centered when possible.  If the value of
c     q at i-1 or i+1 is bad/missing, then a first order, uncentered
c     derivative is considered.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  17 December 2001
c
c############################################################################

      subroutine xderiv(q, dqdx, nx, imax, ny, jmax, nz, kmax, dx)

      implicit none

      include 'dow.inc'

      integer nx, ny, nz       ! grid dimensions
      integer imax, jmax, kmax ! no. of valid grid points in each direction 
      real dx                  ! grid spacing in x direction (km)
      real q(nx, ny, nz)       ! gridded variable
      real dqdx(nx, ny, nz)    ! x derivative of q
      integer i, j, k          ! grid indices
      real dx2                 ! 2.0*dx

      dx2 = 2.0*dx

      do k=1, kmax
        do j=1, jmax

          if ( (q(1,j,k).eq.bad) .or. (q(2,j,k).eq.bad) ) then
            dqdx(1,j,k) = bad
          else
            dqdx(1,j,k) = ( (q(2,j,k)) - (q(1,j,k)) ) / dx
          endif

          if ( (q(imax-1,j,k).eq.bad) .or. (q(imax,j,k).eq.bad) ) then
            dqdx(imax,j,k) = bad
          else
            dqdx(imax,j,k) = ( (q(imax,j,k)) - (q(imax-1,j,k)) ) / dx
          endif

          do i=2, imax-1
            if (q(i+1,j,k).eq.bad) then
              if ( (q(i-1,j,k).eq.bad) .or. (q(i,j,k).eq.bad) ) then
                dqdx(i,j,k) = bad
              else
                dqdx(i,j,k) = ( q(i,j,k) - q(i-1,j,k) ) / dx
              endif
            else if (q(i-1,j,k).eq.bad) then
              if ( (q(i,j,k).eq.bad) .or. (q(i+1,j,k).eq.bad) ) then
                dqdx(i,j,k) = bad
              else
                dqdx(i,j,k) = ( q(i+1,j,k) - q(i,j,k) ) / dx
              endif
            else
              dqdx(i,j,k) = ( q(i+1,j,k) - q(i-1,j,k) ) / dx2
            endif
          enddo

        enddo
      enddo

      return
      end


c############################################################################
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                 SUBROUTINE YDERIV                    ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
c
c############################################################################
c
c     PURPOSE:
c
c     This subroutine computes the derivative in the y direction of
c     a gridded quantity at each grid point.  At the sides, derivatives
c     are first order, uncentered.  Within the interior, derivatives
c     are second order, centered when possible.  If the value of
c     q at j-1 or j+1 is bad/missing, then a first order, uncentered
c     derivative is considered.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  17 December 2001
c
c############################################################################

      subroutine yderiv(q, dqdy, nx, imax, ny, jmax, nz, kmax, dy)

      implicit none

      include 'dow.inc'

      integer nx, ny, nz       ! grid dimensions
      integer imax, jmax, kmax ! no. of valid grid points in each direction 
      real dy                  ! grid spacing in y direction (km)
      real q(nx, ny, nz)       ! gridded variable
      real dqdy(nx, ny, nz)    ! y derivative of q
      integer i, j, k          ! grid indices
      real dy2                 ! 2.0*dy

      dy2 = 2.0*dy

      do k=1, kmax
        do i=1, imax

          if ( (q(i,1,k).eq.bad) .or. (q(i,2,k).eq.bad) ) then
            dqdy(i,1,k) = bad
          else
            dqdy(i,1,k) = ( (q(i,2,k)) - (q(i,1,k)) ) / dy
          endif

          if ( (q(i,jmax-1,k).eq.bad) .or. (q(i,jmax,k).eq.bad) ) then
            dqdy(i,jmax,k) = bad
          else
            dqdy(i,jmax,k) = ( (q(i,jmax,k)) - (q(i,jmax-1,k)) ) / dy
          endif

          do j=2, jmax-1
            if (q(i,j+1,k).eq.bad) then
              if ( (q(i,j-1,k).eq.bad) .or. (q(i,j,k).eq.bad) ) then
                dqdy(i,j,k) = bad
              else
                dqdy(i,j,k) = ( q(i,j,k) - q(i,j-1,k) ) / dy
              endif
            else if (q(i,j-1,k).eq.bad) then
              if ( (q(i,j,k).eq.bad) .or. (q(i,j+1,k).eq.bad) ) then
                dqdy(i,j,k) = bad
              else
                dqdy(i,j,k) = ( q(i,j+1,k) - q(i,j,k) ) / dy
              endif
            else
              dqdy(i,j,k) = ( q(i,j+1,k) - q(i,j-1,k) ) / dy2
            endif
          enddo

        enddo
      enddo

      return
      end


c############################################################################
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                 SUBROUTINE ZDERIV                    ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
c
c############################################################################
c
c     PURPOSE:
c
c     This subroutine computes the derivative in the z direction of
c     a gridded quantity at each grid point.  At the top/bottom, derivatives
c     are first order, uncentered.  Within the interior, derivatives
c     are second order, centered when possible.  If the value of
c     q at k-1 or k+1 is bad/missing, then a first order, uncentered
c     derivative is considered.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  17 December 2001
c
c############################################################################

      subroutine zderiv(q, dqdz, nx, imax, ny, jmax, nz, kmax, dz)

      implicit none

      include 'dow.inc'

      integer nx, ny, nz       ! grid dimensions
      integer imax, jmax, kmax ! no. of valid grid points in each direction 
      real dz                  ! grid spacing in z direction (km)
      real q(nx, ny, nz)       ! gridded variable
      real dqdz(nx, ny, nz)    ! z derivative of q
      integer i, j, k          ! grid indices
      real dz2                 ! 2.0*dz

      dz2 = 2.0*dz

      do j=1, jmax
        do i=1, imax

          if ( (q(i,j,1).eq.bad) .or. (q(i,j,2).eq.bad) ) then
            dqdz(i,j,1) = bad
          else
            dqdz(i,j,1) = ( (q(i,j,2)) - (q(i,j,1)) ) / dz
          endif

          if ( (q(i,j,kmax-1).eq.bad) .or. (q(i,j,kmax).eq.bad) ) then
            dqdz(i,j,kmax) = bad
          else
            dqdz(i,j,kmax) = ( (q(i,j,kmax)) - (q(i,j,kmax-1)) ) / dz
          endif

          do k=2, kmax-1
            if (q(i,j,k+1).eq.bad) then
              if ( (q(i,j,k-1).eq.bad) .or. (q(i,j,k).eq.bad) ) then
                dqdz(i,j,k) = bad
              else
                dqdz(i,j,k) = ( q(i,j,k) - q(i,j,k-1) ) / dz
              endif
            else if (q(i,j,k-1).eq.bad) then
              if ( (q(i,j,k).eq.bad) .or. (q(i,j,k+1).eq.bad) ) then
                dqdz(i,j,k) = bad
              else
                dqdz(i,j,k) = ( q(i,j,k+1) - q(i,j,k) ) / dz
              endif
            else
              dqdz(i,j,k) = ( q(i,j,k+1) - q(i,j,k-1) ) / dz2
            endif
          enddo

        enddo
      enddo

      return
      end


c############################################################################
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                 SUBROUTINE XDERIVC                   ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
c
c############################################################################
c
c     PURPOSE:
c
c     This subroutine computes the derivative in the x direction of
c     a gridded quantity at each grid point.  Only second order,
c     centered derivatives are computed.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  17 December 2001
c
c############################################################################

      subroutine xderivc(q, dqdx, nx, imax, ny, jmax, nz, kmax, dx)

      implicit none

      include 'dow.inc'

      integer nx, ny, nz       ! grid dimensions
      integer imax, jmax, kmax ! no. of valid grid points in each direction 
      real dx                  ! grid spacing in x direction (km)
      real q(nx, ny, nz)       ! gridded variable
      real dqdx(nx, ny, nz)    ! x derivative of q
      integer i, j, k          ! grid indices
      real dx2                 ! 2.0*dx

      dx2 = 2.0*dx

      do k=1, kmax
        do j=1, jmax

          dqdx(1,j,k) = bad
          dqdx(imax,j,k) = bad

          do i=2, imax-1
            if ( (q(i-1,j,k).eq.bad) .or. (q(i+1,j,k).eq.bad) ) then
              dqdx(i,j,k) = bad
            else
              dqdx(i,j,k) = ( q(i+1,j,k) - q(i-1,j,k) ) / dx2
            endif
          enddo

        enddo
      enddo

      return
      end


c############################################################################
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                 SUBROUTINE YDERIVC                   ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
c
c############################################################################
c
c     PURPOSE:
c
c     This subroutine computes the derivative in the y direction of
c     a gridded quantity at each grid point.  Only second order,
c     centered derivatives are computed.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  17 December 2001
c
c############################################################################

      subroutine yderivc(q, dqdy, nx, imax, ny, jmax, nz, kmax, dy)

      implicit none

      include 'dow.inc'

      integer nx, ny, nz       ! grid dimensions
      integer imax, jmax, kmax ! no. of valid grid points in each direction 
      real dy                  ! grid spacing in y direction (km)
      real q(nx, ny, nz)       ! gridded variable
      real dqdy(nx, ny, nz)    ! y derivative of q
      integer i, j, k          ! grid indices
      real dy2                 ! 2.0*dy

      dy2 = 2.0*dy

      do k=1, kmax
        do i=1, imax

          dqdy(i,1,k) = bad
          dqdy(i,jmax,k) = bad

          do j=2, jmax-1
            if ( (q(i,j-1,k).eq.bad) .or. (q(i,j+1,k).eq.bad) ) then
              dqdy(i,j,k) = bad
            else
              dqdy(i,j,k) = ( q(i,j+1,k) - q(i,j-1,k) ) / dy2
            endif
          enddo

        enddo
      enddo

      return
      end

c############################################################################
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                 SUBROUTINE ZDERIVC                   ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
c
c############################################################################
c
c     PURPOSE:
c
c     This subroutine computes the derivative in the z direction of
c     a gridded quantity at each grid point.  Only second order,
c     centered derivatives are computed.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  17 December 2001
c
c############################################################################

      subroutine zderivc(q, dqdz, nx, imax, ny, jmax, nz, kmax, dz)

      implicit none

      include 'dow.inc'

      integer nx, ny, nz       ! grid dimensions
      integer imax, jmax, kmax ! no. of valid grid points in each direction 
      real dz                  ! grid spacing in z direction (km)
      real q(nx, ny, nz)       ! gridded variable
      real dqdz(nx, ny, nz)    ! z derivative of q
      integer i, j, k          ! grid indices
      real dz2                 ! 2.0*dz

      dz2 = 2.0*dz

      do j=1, jmax
        do i=1, imax

          dqdz(i,j,1) = bad
          dqdz(i,j,kmax) = bad

          do k=2, kmax-1
            if ( (q(i,j,k-1).eq.bad) .or. (q(i,j,k+1).eq.bad) ) then
              dqdz(i,j,k) = bad
            else
              dqdz(i,j,k) = ( q(i,j,k+1) - q(i,j,k-1) ) / dz2
            endif
          enddo

        enddo
      enddo

      return
      end


c############################################################################
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                   SUBROUTINE HSMTH                   ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
c
c############################################################################
c
c     PURPOSE:
c
c     This program applies a 5-point horizontal smoother
c     at each horizontal level to the input 3D scalar field
c     
c
c############################################################################
c
c     Author: Steve Weygandt  Dec 1999
c
c
c############################################################################

      subroutine hsmth(fld,nx,ny,nz,npass)

      implicit none

      include 'dow.inc'

      integer numfld, numdvg
      integer i,j,k,np,npass
      integer nx,ny,nz
      real fld(nx,ny,nz)
      real temfld(nx,ny)
      real temdvg(nx,ny)
      real avgfld, avgdvg


      do np = 1,npass

        do k = 1,nz
          do i = 2,nx-1
            do j = 2,ny-1

              numfld = 0
              avgfld = 0.0

              if(fld(i,j,k).ne.bad) then
                numfld = numfld + 1
                avgfld = avgfld + fld(i,j,k)
              end if

              if(fld(i+1,j,k).ne.bad) then
                numfld = numfld + 1
                avgfld = avgfld + fld(i+1,j,k)
              end if

              if(fld(i-1,j,k).ne.bad) then
                numfld = numfld + 1
                avgfld = avgfld + fld(i-1,j,k)
              end if

              if(fld(i,j+1,k).ne.bad) then
                numfld = numfld + 1
                avgfld = avgfld + fld(i,j+1,k)
              end if

              if(fld(i,j-1,k).ne.bad) then
                numfld = numfld + 1
                avgfld = avgfld + fld(i,j-1,k)
              end if

              if(numfld.gt.0) then
                temfld(i,j) = avgfld/float(numfld)
              else
                temfld(i,j) = bad
              end if

            enddo
          enddo
        enddo

        do i = 2,nx-1
          do j = 2,ny-1
            fld(i,j,k) = temfld(i,j)
          enddo
        enddo

      enddo     ! do np = 1,npass

      return
      end



c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                   SUBROUTINE FILL2DD                 ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine fills data voids in a 2-D array.  Dirichlet
c     boundary conditions are used for the poisson equation.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  17 December 2001
c
c############################################################################

      subroutine fill2dd(p,missing,nx,ny,dx,dy,concr)

      implicit none

      include 'dow.inc'

      real alpha                      ! relaxation parameter
      parameter(alpha=1.35)
      integer nx, ny                  ! number of gridpoints in each direction
      real p(nx, ny)                  ! array to be filled
      integer missing(nx, ny)         ! 1 (0) if the data value is (not) missing
      real dx, dy                     ! grid spacing in each direction (m)
      real concr                      ! convergence criterion
      real pbar                       ! value obtained from Gauss-Seidel soln.
      integer i, j                    ! loop variables
      real cx, cy, c                  ! coefficients in residual equation
      real diff                       ! difference between old and new p
      real maxdiff                    ! maximum difference in the domain
      integer done                    ! 1 (0) if (not) finished with solution
      integer iter                    ! no. of iterations at particular level


      c = 0.5*dx*dx*dy*dy / (dx*dx + dy*dy)
      cx = c / (dx*dx)
      cy = c / (dy*dy)

c     Check boundary values.

      do i=1, nx
        if (p(i,1).eq.bad) then
          write(6,*) 'fill2dd:  bad value at ', i, 1
          stop
        endif
        if (p(i,ny).eq.bad) then
          write(6,*) 'fill2dd:  bad value at ', i, ny
          stop
        endif
      enddo
      do j=1, ny
        if (p(1,j).eq.bad) then
          write(6,*) 'fill2dd:  bad value at ', 1, j
          stop
        endif
        if (p(nx,j).eq.bad) then
          write(6,*) 'fill2dd:  bad value at ', nx, j
          stop
        endif
      enddo



c     Here is the main loop.

      done=0
      iter=0

      do while (done.eq.0)

        iter = iter + 1
        maxdiff = 0.0
        do i=2, nx-1
          do j=2, ny-1
            if (missing(i,j).eq.1) then
              pbar = cx * ( p(i+1,j)+p(i-1,j) )
     $             + cy * ( p(i,j+1)+p(i,j-1) )
              if (abs(pbar).lt. 1.0E-30) then
                pbar = 0.0         ! try to avoid "underflow" error messages
              endif
              diff = abs(pbar - p(i,j))
              if (diff.gt.maxdiff) maxdiff=diff
              p(i,j) = (1.0-alpha)*p(i,j) + alpha*pbar
            endif
          enddo
        enddo

        if (iter.eq.500) done=1
        if (maxdiff.le.concr) done=1

      enddo

      write(6,*) 'fill2dd: ', iter, ' iterations'
      write(6,*) 'concr=', concr

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                   SUBROUTINE FILLN                   ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine fills data voids in a 3-D array.  The filling
c     procedure is a least squares method in 2-D (horizontal).  Neumann
c     (zero gradient) boundary conditions are used for the poisson equation.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  17 December 2001
c
c############################################################################

      subroutine filln(pf,nx,ny,nz,dx,dy,concr)

      implicit none

      include 'dow.inc'

      real alpha                      ! relaxation parameter
      parameter(alpha=1.35)
      integer nx, ny, nz              ! number of gridpoints in each direction
      real pf(nx, ny, nz)             ! array to be filled
      real p(0:nx+1, 0:ny+1, nz)      ! temporary array
      integer missing(nx, ny)         ! 1 (0) if the data value is (not) missing
      real dx, dy                     ! grid spacing in each direction (m)
      real concr                      ! convergence criterion
      real pbar                       ! value obtained from Gauss-Seidel soln.
      integer i, j, k                 ! loop variables
      real cx, cy, c                  ! coefficients in residual equation
      real diff                       ! difference between old and new p
      real maxdiff                    ! maximum difference in the domain
      integer done                    ! 1 (0) if (not) finished with solution
      integer iter                    ! no. of iterations at particular level


      c = 0.5*dx*dx*dy*dy / (dx*dx + dy*dy)
      cx = c / (dx*dx)
      cy = c / (dy*dy)

      do k=1, nz

c       Initialize temporary array.

        do j=1, ny
          do i=1, nx
            if (pf(i,j,k).eq.bad) then
              p(i,j,k) = 0.0
              missing(i,j) = 1
            else
              p(i,j,k) = pf(i,j,k)
              missing(i,j) = 0
            endif
          enddo
        enddo

c       Here is the main loop.

        done=0
        iter=0

        do while (done.eq.0)

c         Boundary conditions.

          do i=1, nx
            p(i,0,k) = p(i,1,k)
            p(i,ny+1,k) = p(i,ny,k)
          enddo
          do j=1, ny
            p(0,j,k) = p(1,j,k)
            p(nx+1,j,k) = p(nx,j,k)
          enddo

c         Interior of array.

          iter = iter + 1
          maxdiff = 0.0
          do i=1, nx
            do j=1, ny
              if (missing(i,j).eq.1) then
                pbar = cx * ( p(i+1,j,k)+p(i-1,j,k) )
     $               + cy * ( p(i,j+1,k)+p(i,j-1,k) )
                if (abs(pbar).lt. 1.0E-30) then
                  pbar = 0.0         ! try to avoid "underflow" error messages
                endif
                diff = abs(pbar - p(i,j,k))
                if (diff.gt.maxdiff) maxdiff=diff
                p(i,j,k) = (1.0-alpha)*p(i,j,k) + alpha*pbar
              endif
            enddo
          enddo

          if (iter.eq.500) done=1
          if (maxdiff.le.concr) done=1

        enddo

        write(6,*) 'filln: ', iter, ' iterations, concr=', concr

        do j=1, ny
          do i=1, nx
            pf(i,j,k) = p(i,j,k)
          enddo
        enddo

      enddo

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                 SUBROUTINE POISSON2DN                ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine solves the 2-D Poisson equation
c
c                         d2p/dx2 + d2p/dy2 = f
c
c     with Neumann boundary conditions.  All values of f and the boundary
c     conditions must be valid.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  18 August 2003
c
c############################################################################

      subroutine poisson2dn(f, pfinal,
     $                      bcwest, bceast, bcsouth, bcnorth,
     $                      nx, ny, dx, dy, concr)

      implicit none

      include 'dow.inc'

      real alpha                      ! relaxation parameter
      parameter(alpha=1.35)
      integer nx, ny                  ! number of gridpoints in each direction
      real pfinal(nx, ny)             ! solution of poisson equation
      real p(0:nx+1, 0:ny+1)          ! temporary array
      real f(nx, ny)                  ! forcing function
      real bcwest(ny)                 ! west boundary condition for dp/dx
      real bceast(ny)                 ! east boundary condition for dp/dx
      real bcnorth(nx)                ! north boundary condition for dp/dy
      real bcsouth(nx)                ! south boundary condition for dp/dy
      real dx, dy                     ! grid spacing in each direction (m)
      real concr                      ! convergence criterion
      real pbar                       ! value obtained from Gauss-Seidel soln.
      integer i, j                    ! loop variables
      real cx, cy, c                  ! coefficients in residual equation
      real diff                       ! difference between old and new p
      real maxdiff                    ! maximum difference in the domain
      integer done                    ! 1 (0) if (not) finished with solution
      integer iter                    ! no. of iterations at particular level


      c = 0.5*dx*dx*dy*dy / (dx*dx + dy*dy)
      cx = c / (dx*dx)
      cy = c / (dy*dy)

c     Initialize temporary array.

      p(:,:) = 0.0

c     Here is the main loop.

      done=0
      iter=0

      do while (done.eq.0)

c       Neumann boundary conditions.

        do i=1, nx
          p(i,0) = p(i,2) - 2.0*dy*bcsouth(i)
          p(i,ny+1) = p(i,ny-1) + 2.0*dy*bcnorth(i)
        enddo
        do j=1, ny
          p(0,j) = p(2,j) - 2.0*dx*bcwest(j)
          p(nx+1,j) = p(nx-1,j) + 2.0*dx*bceast(j)
        enddo

c       Interior of array.

        iter = iter + 1
        maxdiff = 0.0
        do i=1, nx
          do j=1, ny
            pbar = cx * ( p(i+1,j)+p(i-1,j) )
     $           + cy * ( p(i,j+1)+p(i,j-1) )
     $           - c*f(i,j)
            if (abs(pbar).lt. 1.0E-30) then
              pbar = 0.0         ! try to avoid "underflow" error messages
            endif
            diff = abs(pbar - p(i,j))
            if (diff.gt.maxdiff) maxdiff=diff
            p(i,j) = (1.0-alpha)*p(i,j) + alpha*pbar
          enddo
        enddo

        if (iter.eq.5000) done=1
        if (maxdiff.le.concr) done=1

      enddo

      write(6,*) 'poisson2dn: iter=', iter, ', concr=', concr,
     $           ', maxdiff=', maxdiff

      do j=1, ny
        do i=1, nx
          pfinal(i,j) = p(i,j)
        enddo
      enddo

      return
      end


c############################################################################
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                 SUBROUTINE REMOVEHAVG                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
c
c############################################################################
c
c     PURPOSE:
c
c     This subroutine subtracts the horizontal average of q from q
c     at each level in a 3D array.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  18 August 2003
c
c############################################################################

      subroutine removehavg(q, nx, ny, nz)

      implicit none

      integer nx, ny, nz       ! grid dimensions
      real q(nx, ny, nz)       ! gridded variable
      integer i, j, k          ! grid indices
      real qsum                ! sum of values at current level
      real qavg                ! average of q at current level


      do k=1, nz
        qsum = 0.0
        do j=1, ny
          do i=1, nx
            qsum = qsum + q(i,j,k)
          enddo
        enddo
        qavg = qsum / (nx*ny)
        do j=1, ny
          do i=1, nx
            q(i,j,k) = q(i,j,k) - qavg
          enddo
        enddo
      enddo

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                  SUBROUTINE DEFORMATION              ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine computes the following:
c     1.  magnitude of the resultant horizontal deformation
c     2.  x component of horizontal deformation "vector"
c     3.  y component of horizontal deformation "vector" (/s)
c     4.  Dave Schultz' diagnostic quantity:  (def)**2 - (vorz)**2
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  24 November 2004
c
c############################################################################

      subroutine deformation(u, v, w, def, xdef, ydef, defvort,
     $                       nx, ny, nz, dx, dy)

      implicit none

      include 'dow.inc'

      integer nx, ny, nz       ! grid dimensions 
      real dx, dy              ! grid spacing in each direction (km)
      real dxm, dym            ! grid spacing (m)
      real u(nx,ny,nz)         ! x component of velocity (u) (m/s)
      real v(nx,ny,nz)         ! y component of velocity (v) (m/s)
      real w(nx,ny,nz)         ! z component of velocity (w) (m/s)
      real def(nx,ny,nz)       ! magnitude of the resultant horizontal deformation (/s)
      real xdef(nx,ny,nz)      ! x component of horizontal deformation "vector" (/s)
      real ydef(nx,ny,nz)      ! y component of horizontal deformation "vector" (/s)
      real defvort(nx,ny,nz)   ! Dave Schultz' diagnostic quantity:
                               !   (def)**2 - (vorz)**2  (/s**2)
      real d1, d2, vorz
      integer i, j, k          ! grid indices

      real, allocatable :: dvdx(:,:,:)  ! velocity derivatives
      real, allocatable :: dudy(:,:,:)  ! "                  "
      real, allocatable :: dudx(:,:,:)  ! "                  "
      real, allocatable :: dvdy(:,:,:)  ! "                  "


      dxm = dx * 1000.0
      dym = dy * 1000.0

      allocate(dvdx(nx,ny,nz))
      allocate(dudy(nx,ny,nz))
      allocate(dudx(nx,ny,nz))
      allocate(dvdy(nx,ny,nz))

      call xderivc(v, dvdx, nx, nx, ny, ny, nz, nz, dxm)
      call xderivc(u, dudx, nx, nx, ny, ny, nz, nz, dxm)
      call yderivc(u, dudy, nx, nx, ny, ny, nz, nz, dym)
      call yderivc(v, dvdy, nx, nx, ny, ny, nz, nz, dym)

      do k=1, nz
        do j=1, ny
          do i=1, nx
            if ( (dudx(i,j,k).eq.bad) .or. (dudy(i,j,k).eq.bad) .or.
     $           (dvdx(i,j,k).eq.bad) .or. (dvdy(i,j,k).eq.bad) ) then
              def(i,j,k) = bad
              xdef(i,j,k) = bad
              ydef(i,j,k) = bad
              defvort(i,j,k) = bad
            else
              d1 = dudx(i,j,k) - dvdy(i,j,k)
              d2 = dvdx(i,j,k) + dudy(i,j,k)
              vorz = dvdx(i,j,k) - dudy(i,j,k)
              def(i,j,k) = sqrt(d1*d1 + d2*d2)
              xdef(i,j,k) = sqrt(0.5*(def(i,j,k)**2 + d1*def(i,j,k)))
              ydef(i,j,k) = sqrt(0.5*(def(i,j,k)**2 - d1*def(i,j,k)))
              if (d2.lt.0.0) ydef(i,j,k)=-ydef(i,j,k)
              defvort(i,j,k) = def(i,j,k)**2 - vorz**2
            endif
          enddo
        enddo
      enddo

      deallocate(dvdx)
      deallocate(dudx)
      deallocate(dudy)
      deallocate(dvdy)

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                  SUBROUTINE VORTICITY                ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine computes the 3-D vorticity.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  17 December 2001
c
c############################################################################

      subroutine vorticity(u, v, w, vorx, vory, vorz,
     $                     nx, ny, nz, dx, dy, dz)

      implicit none

      include 'dow.inc'

      integer nx, ny, nz       ! grid dimensions 
      real dx, dy, dz          ! grid spacing in each direction (km)
      real dxm, dym, dzm       ! grid spacing (m)
      real u(nx,ny,nz)         ! x component of velocity (u) (m/s)
      real v(nx,ny,nz)         ! y component of velocity (v) (m/s)
      real w(nx,ny,nz)         ! z component of velocity (w) (m/s)
      real vorx(nx,ny,nz)      ! x component of vorticity (/s)
      real vory(nx,ny,nz)      ! y component of vorticity (/s)
      real vorz(nx,ny,nz)      ! z component of vorticity (/s)
      integer i, j, k          ! grid indices

      real, allocatable :: dvdx(:,:,:)  ! velocity derivatives
      real, allocatable :: dudy(:,:,:)  ! "                  "
      real, allocatable :: dwdy(:,:,:)  ! "                  "
      real, allocatable :: dvdz(:,:,:)  ! "                  "
      real, allocatable :: dudz(:,:,:)  ! "                  "
      real, allocatable :: dwdx(:,:,:)  ! "                  "


      dxm = dx * 1000.0
      dym = dy * 1000.0
      dzm = dz * 1000.0

      allocate(dvdx(nx,ny,nz))
      allocate(dudy(nx,ny,nz))
      allocate(dwdy(nx,ny,nz))
      allocate(dvdz(nx,ny,nz))
      allocate(dudz(nx,ny,nz))
      allocate(dwdx(nx,ny,nz))

      call xderivc(v, dvdx, nx, nx, ny, ny, nz, nz, dxm)
c      call xderiv(v, dvdx, nx, nx, ny, ny, nz, nz, dxm)
      call yderivc(u, dudy, nx, nx, ny, ny, nz, nz, dym)
c      call yderiv(u, dudy, nx, nx, ny, ny, nz, nz, dym)
      call xderivc(w, dwdx, nx, nx, ny, ny, nz, nz, dxm)
      call yderivc(w, dwdy, nx, nx, ny, ny, nz, nz, dym)
      call zderivc(v, dvdz, nx, nx, ny, ny, nz, nz, dzm)
      call zderivc(u, dudz, nx, nx, ny, ny, nz, nz, dzm)

c      write(6,*) 'IGNORING W SHEAR IN VORTICITY'
c      dwdx(:,:,:) = 0.0
c      dwdy(:,:,:) = 0.0

      do k=1, nz
        do j=1, ny
          do i=1, nx
            if ( (dwdy(i,j,k).eq.bad) .or.
     $           (dvdz(i,j,k).eq.bad) ) then
              vorx(i,j,k) = bad
            else
              vorx(i,j,k) = dwdy(i,j,k) - dvdz(i,j,k)
            endif
            if ( (dudz(i,j,k).eq.bad) .or.
     $           (dwdx(i,j,k).eq.bad) ) then
              vory(i,j,k) = bad
            else
              vory(i,j,k) = dudz(i,j,k) - dwdx(i,j,k)
            endif
            if ( (dvdx(i,j,k).eq.bad) .or.
     $           (dudy(i,j,k).eq.bad) ) then
              vorz(i,j,k) = bad
            else
              vorz(i,j,k) = dvdx(i,j,k) - dudy(i,j,k)
            endif
          enddo
        enddo
      enddo

      deallocate(dvdx)
      deallocate(dudy)
      deallocate(dwdy)
      deallocate(dvdz)
      deallocate(dudz)
      deallocate(dwdx)

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                   SUBROUTINE TILTING                 ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine computes tilting of horizontal vorticity into
c     the vertical.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  17 December 2001
c
c############################################################################

      subroutine tilting(u, v, w, tilt, nx, ny, nz, dx, dy, dz)

      implicit none

      include 'dow.inc'

      integer nx, ny, nz       ! grid dimensions 
      real dx, dy, dz          ! grid spacing in each direction (km)
      real dxm, dym, dzm       ! grid spacing (m)
      real u(nx,ny,nz)         ! u component of velocity (m/s)
      real v(nx,ny,nz)         ! v component of velocity (m/s)
      real w(nx,ny,nz)         ! w component of velocity (m/s)
      real tilt(nx,ny,nz)      ! tilting (/s**2)
      integer i, j, k          ! grid indices

      real, allocatable :: dudz(:,:,:)  ! velocity derivatives
      real, allocatable :: dwdy(:,:,:)  ! "                  "
      real, allocatable :: dvdz(:,:,:)  ! "                  "
      real, allocatable :: dwdx(:,:,:)  ! "                  "


      dxm = dx * 1000.0
      dym = dy * 1000.0
      dzm = dz * 1000.0

      allocate(dudz(nx,ny,nz))
      allocate(dwdy(nx,ny,nz))
      allocate(dvdz(nx,ny,nz))
      allocate(dwdx(nx,ny,nz))

      call zderivc(u, dudz, nx, nx, ny, ny, nz, nz, dzm)
      call yderivc(w, dwdy, nx, nx, ny, ny, nz, nz, dym)
      call zderivc(v, dvdz, nx, nx, ny, ny, nz, nz, dzm)
      call xderivc(w, dwdx, nx, nx, ny, ny, nz, nz, dxm)

      do k=1, nz
        do j=1, ny
          do i=1, nx
            if ( (dudz(i,j,k).eq.bad) .or.
     $           (dwdy(i,j,k).eq.bad) .or.
     $           (dvdz(i,j,k).eq.bad) .or.
     $           (dwdx(i,j,k).eq.bad) ) then
              tilt(i,j,k) = bad
            else
              tilt(i,j,k) = dudz(i,j,k)*dwdy(i,j,k)
     $                    - dvdz(i,j,k)*dwdx(i,j,k)
            endif
          enddo
        enddo
      enddo

      deallocate(dudz)
      deallocate(dwdy)
      deallocate(dvdz)
      deallocate(dwdx)

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                SUBROUTINE STRETCHING                 ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine computes stretching of vertical vorticity.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  17 December 2001
c
c############################################################################

      subroutine stretching(u, v, stretch, nx, ny, nz, dx, dy)

      implicit none

      include 'dow.inc'

      integer nx, ny, nz       ! grid dimensions 
      real dx, dy              ! grid spacing in each direction (km)
      real dxm, dym            ! grid spacing (m)
      real u(nx,ny,nz)         ! u component of velocity (m/s)
      real v(nx,ny,nz)         ! v component of velocity (m/s)
      real stretch(nx,ny,nz)   ! stretching of vertical vorticity (/s**2)
      integer i, j, k          ! grid indices

      real, allocatable :: dvdx(:,:,:)  ! velocity derivatives
      real, allocatable :: dudy(:,:,:)  ! "                  "
      real, allocatable :: dudx(:,:,:)  ! "                  "
      real, allocatable :: dvdy(:,:,:)  ! "                  "


      dxm = dx * 1000.0
      dym = dy * 1000.0

      allocate(dvdx(nx,ny,nz))
      allocate(dudy(nx,ny,nz))
      allocate(dudx(nx,ny,nz))
      allocate(dvdy(nx,ny,nz))

      call xderivc(u, dudx, nx, nx, ny, ny, nz, nz, dxm)
      call xderivc(v, dvdx, nx, nx, ny, ny, nz, nz, dxm)
      call yderivc(u, dudy, nx, nx, ny, ny, nz, nz, dym)
      call yderivc(v, dvdy, nx, nx, ny, ny, nz, nz, dym)

      do k=1, nz
        do j=1, ny
          do i=1, nx
            if ( (dudx(i,j,k).eq.bad) .or.
     $           (dvdx(i,j,k).eq.bad) .or.
     $           (dudy(i,j,k).eq.bad) .or.
     $           (dvdy(i,j,k).eq.bad) ) then
              stretch(i,j,k) = bad
            else
              stretch(i,j,k) = - (dvdx(i,j,k)-dudy(i,j,k))
     $                          *(dudx(i,j,k)+dvdy(i,j,k))
            endif
          enddo
        enddo
      enddo

      deallocate(dvdx)
      deallocate(dudy)
      deallocate(dudx)
      deallocate(dvdy)

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                  SUBROUTINE ADVECTION                ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine computes the 3-D advection of q.  A regular
c     (not staggered) grid is assumed.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  15 August 2003
c
c############################################################################

      subroutine advection(q, qadv, u, v, w, nx, ny, nz, dx, dy, dz)

      implicit none

      include 'dow.inc'

      integer nx, ny, nz       ! grid dimensions 
      real dx, dy, dz          ! grid spacing in each direction (m)
      real u(nx,ny,nz)         ! x component of velocity (u) (m/s)
      real v(nx,ny,nz)         ! y component of velocity (v) (m/s)
      real w(nx,ny,nz)         ! z component of velocity (w) (m/s)
      real q(nx,ny,nz)         ! advected quantity
      real qadv(nx,ny,nz)      ! advection of q
      real dqdx(nx,ny,nz)      ! spatial derivative of q
      real dqdy(nx,ny,nz)      ! "                     "
      real dqdz(nx,ny,nz)      ! "                     "
      integer i, j, k          ! grid indices


      call xderiv(q, dqdx, nx, nx, ny, ny, nz, nz, dx)
      call yderiv(q, dqdy, nx, nx, ny, ny, nz, nz, dy)
      call zderiv(q, dqdz, nx, nx, ny, ny, nz, nz, dz)

      do k=1, nz
        do j=1, ny
          do i=1, nx
            if ( (dqdx(i,j,k).eq.bad) .or.
     $           (dqdy(i,j,k).eq.bad) .or.
     $           (dqdz(i,j,k).eq.bad) .or.
     $           (u(i,j,k).eq.bad) .or.
     $           (v(i,j,k).eq.bad) .or.
     $           (w(i,j,k).eq.bad) ) then
              qadv(i,j,k) = bad
            else
              qadv(i,j,k) = -u(i,j,k)*dqdx(i,j,k)
     $                     - v(i,j,k)*dqdy(i,j,k)
     $                     - w(i,j,k)*dqdz(i,j,k)
            endif
          enddo
        enddo
      enddo

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                  SUBROUTINE MEANWIND                 ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine outputs the mean horizontal wind at each level.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  17 December 2001
c
c############################################################################

      subroutine meanwind(u, v, nx, ny, nz, dz, zmin)

      implicit none

      include 'dow.inc'

      integer nx, ny, nz       ! grid dimensions 
      real dz                  ! vertical grid spacing (km)
      real zmin                ! height of lowest grid level (km)
      real u(nx,ny,nz)         ! u component of velocity (m/s)
      real v(nx,ny,nz)         ! v component of velocity (m/s)
      integer i, j, k          ! grid indices
      integer i1, i2, j1, j2   ! horizontal limits
      real usum, vsum          ! sums of u and v at each level
      integer nsum             ! number of summed values


      i1 = 1
      i2 = nx
      j1 = 1
      j2 = ny

      write(20,*)
      write(20,*) 'MEAN WINDS AT EACH LEVEL'
      write(20,*)
      write(20,*) 'i1, i2 = ', i1, i2
      write(20,*) 'j1, j2 = ', j1, j2
      write(20,*)

      do k=1, nz

        nsum = 0
        usum = 0.0
        vsum = 0.0

        do j=j1, j2
          do i=i1, i2
            if ( (u(i,j,k).ne.bad) .and. (v(i,j,k).ne.bad) ) then
              usum = usum + u(i,j,k)
              vsum = vsum + v(i,j,k)
              nsum = nsum + 1
            endif
          enddo
        enddo

        if (nsum.eq.0) then
          write(20,98) zmin+(k-1.0)*dz
 98       format('z = ', F7.2, '   NO DATA')
        else
          write(20,99) zmin+(k-1.0)*dz, usum/nsum, vsum/nsum
 99       format('z = ', F7.2, '     U = ', F7.2, ', V = ', F7.2)
        endif

      enddo

      write(20,*)

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                   SUBROUTINE LL_TO_XY                ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine computes the projected (x, y) coordinates of the
c     point (lat2, lon2) relative to (lat1, lon1).  Various map projections
c     are possible.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  25 February 2005
c
c############################################################################

      subroutine ll_to_xy(x, y, map_proj, lat1, lon1, lat2, lon2)

      implicit none

      include 'dow.inc'

c---- Passed variables

      real lat1, lon1             ! coordinates of first point (rad)
      real lat2, lon2             ! coordinates of second point (rad)
      integer map_proj            ! map projection:
                                  !   0 = flat earth
                                  !   1 = oblique azimuthal
                                  !   2 = Lambert conformal

c---- Returned variables

      real x, y                   ! distance (km)


      if (map_proj.eq.0) then
        x = rearth * cos(0.5*(lat1+lat2)) * (lon2-lon1)
        y = rearth * (lat2-lat1)
      else
        write(6,*) 'map projection unavailable:  ', map_proj
        stop
      endif

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                   SUBROUTINE XY_TO_LL                ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine computes the projected (lat, lon) coordinates of the
c     point (x, y) relative to (lat0, lon0).  Various map projections
c     are possible.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  25 February 2005
c
c############################################################################

      subroutine xy_to_ll(lat, lon, map_proj, x, y, lat0, lon0)

      implicit none

      include 'dow.inc'

c---- Passed variables

      integer map_proj            ! map projection:
                                  !   0 = flat earth
                                  !   1 = oblique azimuthal
                                  !   2 = Lambert conformal
      real x, y                   ! distance (km)
      real lat0, lon0             ! coordinates (rad) of origin (where x=0, y=0)

c---- Returned variables

      real lat, lon               ! coordinates (rad) of point


      if (map_proj.eq.0) then
        lat = lat0 + y / rearth
        lon = lon0 + x / ( rearth * cos(0.5*(lat0+lat)) )
      else
        write(6,*) 'map projection unavailable:  ', map_proj
        stop
      endif

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                   SUBROUTINE DIST                    ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine computes the approximate distance (along a curved
c     path following the earth's surface) from (lon1, lat1, alt1)
c     to (lon2, lat2, alt2).
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  17 December 2001
c
c############################################################################

      subroutine dist(x, y, z, lon1, lat1, alt1, lon2, lat2, alt2)

      implicit none

      include 'dow.inc'

      real x, y, z                     ! distance (km)
      real lon1, lon2                  ! longitude (deg)
      real lat1, lat2                  ! latitude (deg)
      real alt1, alt2                  ! altitude (km MSL)

      x = 2.0*pi*rearth * cos(0.5*(lat1+lat2)*dtor)
     $                  * (lon2-lon1)/360.0
      y = 2.0*pi*rearth * (lat2-lat1)/360.0
      z = alt2 - alt1

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######               SUBROUTINE CORRECT_AZ_EL               ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine accomplishes three things:
c     1. adds a user specified correction to the elevation angles
c     2. adds a user specified correction to the azimuth angles
c     3. adjusts the beam locations to correspond to the middle of the
c        integration period rather than the end
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  5 March 2003
c
c############################################################################

      subroutine correct_az_el(ryib, num_rays, azcor, elcor,
     * az_corr_flag)

      implicit none

      include 'dow.inc'
      include 'structures.inc'

      type(ryib_info), dimension(maxrays) :: ryib    ! ray info. block
      integer num_rays      ! number of rays
      real azcor            ! user specified azimuth angle correction (deg)
      real elcor            ! user specified elevation angle correction (deg)
      integer r             ! ray number
      real azcompx          ! x component of azimuth direction
      real azcompy          ! y component of azimuth direction
      real comptoaz         ! function used by this subroutine
      real oldaz, oldel     ! previous values of az and el
      integer az_corr_flag  ! method of additional azimuth-angle correction
                            !   0 = none
                            !   1 = average current and previous angle


      if (num_rays.lt.3) then
        write(6,*) 'correct_az_el:  not enough rays'
        write(6,*) 'num_rays = ', num_rays
        stop
      endif

c     Add the user specified corrections.

      do r=1, num_rays
        ryib(r)%azimuth = ryib(r)%azimuth + azcor
        if (ryib(r)%azimuth .gt. 360.0) then
          ryib(r)%azimuth = ryib(r)%azimuth - 360.0
        endif
        if (ryib(r)%azimuth .lt. 0.0) then
          ryib(r)%azimuth = ryib(r)%azimuth + 360.0
        endif

c        if (r .eq. 10) then
c        write(6,*)
c        write(6,*) 'ryib(r)%elevation=',ryib(r)%elevation
c        endif

        ryib(r)%elevation = ryib(r)%elevation + elcor
        if (ryib(r)%elevation .gt. 180.0) then
          ryib(r)%elevation = ryib(r)%elevation - 360.0
        endif
      enddo

      if (az_corr_flag.eq.1) then

c       To determine the middle location of beam r, average
c       the end locations of beams r and r-1.

        do r=num_rays, 2, -1
          azcompx = 0.5*sin(dtor*ryib(r)%azimuth)
     $            + 0.5*sin(dtor*ryib(r-1)%azimuth)
          azcompy = 0.5*cos(dtor*ryib(r)%azimuth)
     $            + 0.5*cos(dtor*ryib(r-1)%azimuth)
          ryib(r)%azimuth = comptoaz(azcompx, azcompy)
          ryib(r)%elevation = 0.5*ryib(r)%elevation
     $                      + 0.5*ryib(r-1)%elevation
        enddo

c       We can't determine exactly the middle location of the first beam.
c       Try extrapolating from the locations of the second and
c       third beams.

        oldaz = ryib(1)%azimuth
        oldel = ryib(1)%elevation

        azcompx = 2.0*sin(dtor*ryib(2)%azimuth)
     $          - 1.0*sin(dtor*ryib(3)%azimuth)
        azcompy = 2.0*cos(dtor*ryib(2)%azimuth)
     $          - 1.0*cos(dtor*ryib(3)%azimuth)
        ryib(1)%azimuth = comptoaz(azcompx, azcompy)
        ryib(1)%elevation = 2.0*ryib(2)%elevation
     $                    - 1.0*ryib(3)%elevation

        if ( (abs(oldaz-ryib(1)%azimuth).gt.1.0) .or.
     $       (abs(oldel-ryib(1)%elevation).gt.1.0) ) then
          write(6,*) '*** warning from correct_az_el about ',
     $               'extrapolated beam location'
          write(6,*) 'old azimuth and elevation: ', oldaz, oldel
          write(6,*) 'new azimuth and elevation: ',
     $                ryib(1)%azimuth, ryib(1)%elevation
        endif

      endif        ! if (az_corr_flag.eq.1)

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                 REAL FUNCTION COMPTOAZ               ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This function computes the azimuth angle from the two components
c     of the beam direction.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  5 March 2003
c
c############################################################################

      real function comptoaz(azcompx, azcompy)

      implicit none

      include 'dow.inc'

      real azcompx        ! x component of azimuth direction
      real azcompy        ! y component of azimuth direction

      if ( (azcompx.eq.0.0) .and.
     $     (azcompy.eq.0.0) ) then
        comptoaz = bad
      else if (azcompy.eq.0.0) then
        if (azcompx.gt.0.0) then
          comptoaz = 90.0
        else
          comptoaz = 270.0
        endif
      else if ( (azcompx.ge.0.0) .and.
     $          (azcompy.gt.0.0) ) then
        comptoaz = atan(azcompx/azcompy)
     $           * rtod
      else if ( (azcompx.ge.0.0) .and.
     $          (azcompy.lt.0.0) ) then
        comptoaz = -atan(azcompx
     $                   /abs(azcompy))
     $              * rtod + 180.0
      else if ( (azcompx.lt.0.0) .and.
     $          (azcompy.lt.0.0) ) then
        comptoaz = atan(azcompx/azcompy)
     $             * rtod + 180.0
      else
        comptoaz = -atan(abs(azcompx)
     $                   /azcompy)
     $              * rtod + 360.0
      endif

      return
      end








c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######               SUBROUTINE CORRECT_DATE                ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     Correct the date fields in the input radar data.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  February 2005
c
c############################################################################

      subroutine correct_date(year, yrcor, month, mocor, day, dacor)

      implicit none

      integer(kind=2) year, month, day  ! date
      integer yrcor                     ! correction to year
      integer mocor                     ! correction to month
      integer dacor                     ! correction to day

      year = year + yrcor
      month = month + mocor
      day = day + dacor
      write(6,*) 'date:  ', year, month, day

      return
      end










c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######             SUBROUTINE CORRECT_UMASS_DATA            ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutines tries to correct some of the errors in UMass radar data.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  15 February 2005
c
c############################################################################

      subroutine correct_umass_data(num_rays, ryib, total_gates, rdat)

      implicit none

      include 'dow.inc'
      include 'structures.inc'

      integer num_rays
      integer total_gates
      type(ryib_info), dimension(maxrays) :: ryib
      type(rdat_info), dimension(maxrays) :: rdat
      integer r, g

c     Corrections:
c     1. Elevation angles are not known, and the values that are stored in
c        the radar data are meaningless.  For now, assume the elevation
c        angle should be 0.5 degrees.
c     2. Bad/missing data are represented by many different flags in UMass radar
c        data.  For now, assume any number <= -60.0 represents bad data.

      write(6,*) 'correcting UMass radar data...'

      do r=1, num_rays
        ryib(r)%elevation = 0.5
        do g=1, total_gates
          if (rdat(r)%data(g) .le. -60.00) then
            rdat(r)%data(g) = sbad
          endif
        enddo
      enddo

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######               SUBROUTINE GRID_COORDINATES            ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     Fill the array "x" with the Cartesian coordinates corresponding
c     to the input parameters.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  February 2005
c
c############################################################################

      subroutine grid_coordinates(x, nx, dx, xmin)

      implicit none

      integer nx           ! number of grid points
      real x(nx)           ! grid coordinates
      real dx              ! distance between grid points
      real xmin            ! x(1)
      integer i

      do i=1, nx
        x(i) = xmin+(i-1.0)*dx
      enddo

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######               SUBROUTINE CHECK_SWEEP_SIZE            ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     Determine whether array bounds "maxflds" and "maxrays" are large
c     enough to store the input sweep file.  If they are not, then halt.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  February 2005
c
c############################################################################

      subroutine check_sweep_size(num_param_desc, num_rays)

      implicit none

      include 'dow.inc'

      integer(kind=2) num_param_desc
      integer(kind=4) num_rays

      if (num_param_desc .gt. maxflds) then
        write(6,*) 'too many fields in sweep file'
        write(6,*) 'num_param_desc = ', num_param_desc
        write(6,*) 'maxflds = ', maxflds
        stop
      endif
      if (num_rays .gt. maxrays) then
        write(6,*) 'too many rays in sweep file'
        write(6,*) 'num_rays = ', num_rays
        write(6,*) 'maxrays = ', maxrays
        stop
      endif

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######              REAL FUNCTION RANGEKM_TO_GATE           ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     Return the range, in km, to the center of the current gate.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  February 2005
c
c############################################################################

      real function rangekm_to_gate(celv, g)

      implicit none

      include 'dow.inc'
      include 'structures.inc'

      type(celv_info) :: celv    ! information about ranges, in m, to the
                                 !   near edges of the gates
      integer g                  ! current gate number


      if ( (celv%total_gates.le.1) .or. (g.lt.1) .or.
     * (g.gt.celv%total_gates) ) then

        rangekm_to_gate = bad
      else if (g.eq.celv%total_gates) then
        rangekm_to_gate = ( celv%gate_spacing(g)
     $                     +0.5*( celv%gate_spacing(g)
     $                           -celv%gate_spacing(g-1) ) )
     $                  / 1000.0
      else
        rangekm_to_gate = 0.5 * ( celv%gate_spacing(g)
     $                           +celv%gate_spacing(g+1) )
     $                  / 1000.0
      endif

      return
      end






c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######               SUBROUTINE UPDATE_DIR_BINS             ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     Update the runnings sums of weights for each directional bin.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Creation Date:  February 2005
c
c############################################################################

      subroutine update_dir_bins(sum, wgt, x, y, z, xg, yg, zg)

      implicit none

c
c.... CLZ (9/15/08): test option to prevent horizontal but allow vertical extrapolation
c

c
c.... CLZ (3/9/10): set option to allow/prevent horizontal/vertical extrapolation
c


c      real x, y, z                    ! location of gate relative to grid origin (km)
c      real wgt                        ! datum weight
c      sumdir(:,:,:,:,:) --> sum(8 - or - 4) which is weights in directional bins

c      xg(:) are the x coordinates (km) of grid points
c      yg(:) are the y coordinates (km) of grid points
c      zg(:) are the z coordinates (km) of grid points


c
c.... sum(8) if not allowing extrapolation (extrapolate_flag = 0)
c
      real sum(8)

c
c.... sum(4) if allowing vertical but not horizontal extrapolation (extrapolate_flag = 1)
c
c      real sum(4)

      real wgt
      real x, y, z
      real xg, yg, zg
      integer ihorextraponly


c
c.... begin executable code
c



      ihorextraponly = 0

c
c.... ihorextraponly = 0 is default option
c


      if (ihorextraponly .eq. 0) then

c
c.... CLZ (3/9/10): must have 1st dimension of sumdir in calling routine set = 8
c

      if ( (x.le.xg) .and. (y.le.yg) .and. (z.le.zg) ) sum(1) = sum(1) + wgt
      if ( (x.le.xg) .and. (y.le.yg) .and. (z.ge.zg) ) sum(2) = sum(2) + wgt
      if ( (x.le.xg) .and. (y.ge.yg) .and. (z.le.zg) ) sum(3) = sum(3) + wgt
      if ( (x.le.xg) .and. (y.ge.yg) .and. (z.ge.zg) ) sum(4) = sum(4) + wgt
      if ( (x.ge.xg) .and. (y.le.yg) .and. (z.le.zg) ) sum(5) = sum(5) + wgt
      if ( (x.ge.xg) .and. (y.le.yg) .and. (z.ge.zg) ) sum(6) = sum(6) + wgt
      if ( (x.ge.xg) .and. (y.ge.yg) .and. (z.le.zg) ) sum(7) = sum(7) + wgt
      if ( (x.ge.xg) .and. (y.ge.yg) .and. (z.ge.zg) ) sum(8) = sum(8) + wgt

      endif


      if (ihorextraponly .eq. 1) then

c
c.... CLZ (3/9/10): must have 1st-dimension of sumdir in calling routine set = 4
c

c      if ( (x.le.xg) .and. (y.le.yg) .and. (z.le.zg) ) sum(1) = sum(1) + wgt
      if ( (x.le.xg) .and. (y.le.yg) ) sum(1) = sum(1) + wgt

c      if ( (x.le.xg) .and. (y.le.yg) .and. (z.ge.zg) ) sum(2) = sum(2) + wgt
c      if ( (x.le.xg) .and. (y.ge.yg) .and. (z.le.zg) ) sum(3) = sum(3) + wgt
      if ( (x.le.xg) .and. (y.ge.yg) ) sum(2) = sum(2) + wgt
c      if ( (x.le.xg) .and. (y.ge.yg) .and. (z.ge.zg) ) sum(4) = sum(4) + wgt
c      if ( (x.ge.xg) .and. (y.le.yg) .and. (z.le.zg) ) sum(5) = sum(5) + wgt
      if ( (x.ge.xg) .and. (y.le.yg) ) sum(3) = sum(3) + wgt
c      if ( (x.ge.xg) .and. (y.le.yg) .and. (z.ge.zg) ) sum(6) = sum(6) + wgt
c      if ( (x.ge.xg) .and. (y.ge.yg) .and. (z.le.zg) ) sum(7) = sum(7) + wgt
      if ( (x.ge.xg) .and. (y.ge.yg) ) sum(4) = sum(4) + wgt
c      if ( (x.ge.xg) .and. (y.ge.yg) .and. (z.ge.zg) ) sum(8) = sum(8) + wgt

      endif




c
c.... done with subroutine update_dir_bins
c

      return
      end









c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######              SUBROUTINE GET_BEAM_INFO                ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     For all beams in a set of sweep files, this subroutine catalogs
c     each beam's time offset, azimuth angle, and elevation angle.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Created:  August 2005
c
c     Modified:  
c
c############################################################################

      subroutine get_beam_info(nswp, sfname, fname,
     $                         yrcor, mocor, dacor,
     $                         azcor, elcor, az_corr_flag, umass_flag,
     $                         cyr, cmo, cda, chr, cmn, cse,
     $                         beam_info, num_beams)

      implicit none

      include 'dow.inc'
      include 'structures.inc'

c---- Passed variables

      integer nswp                    ! number of sweep files
      character(len=200) sfname(nswp) ! sweep file names
      character(len=8) fname          ! field name
      integer yrcor                   ! correction to year
      integer mocor                   ! correction to month
      integer dacor                   ! correction to day
      real elcor                      ! elevation angle correction (deg)
      real azcor                      ! azimuth-angle offset (deg)
      integer az_corr_flag            ! method of additional azimuth-angle correction
                                      !   0 = none
                                      !   1 = average current and previous angle
      integer umass_flag              ! apply UMass data corrections? (1=yes, 0=no)
      integer(kind=2) cyr,cmo,cda     ! central date
      integer(kind=2) chr,cmn,cse     ! central time

c---- Returned variables

      integer num_beams(nswp)         ! number of beams in each sweep
      real beam_info(maxrays,nswp,3)  ! beam information:
                                      !   (1) time offset (s) from central time
                                      !   (2) azimuth angle (deg)
                                      !   (3) elevation angle (deg)

c---- Local variables

      integer(kind=2) cms; parameter(cms=0.0)
      integer s                       ! sweep file number
      integer r                       ! ray number
      real ti                         ! time (s) relative to central time
      real timediff                   ! functions used by this subroutine

      type(vold_info)                     :: vold
      type(radd_info)                     :: radd
      type(celv_info)                     :: celv
      type(cfac_info)                     :: cfac
      type(parm_info), dimension(maxflds) :: parm
      type(swib_info)                     :: swib
      type(ryib_info), dimension(maxrays) :: ryib
      type(asib_info), dimension(maxrays) :: asib
      type(rdat_info), dimension(maxrays) :: rdat


      do s=1, nswp

      write(6,*) 'reading from ', sfname(s)
      call sweepread(sfname(s),vold,radd,celv,cfac,parm,swib,ryib,
     * asib,rdat,fname)

      call check_sweep_size(radd%num_param_desc, swib%num_rays)
      call correct_date(vold%year, yrcor, vold%mon, mocor, vold%day,
     * dacor)
      call correct_az_el(ryib, swib%num_rays, azcor, elcor,
     * az_corr_flag)

      if (umass_flag.eq.1) then
      call correct_umass_data(swib%num_rays, ryib, celv%total_gates,
     * rdat)
      endif

      num_beams(s) = swib%num_rays

      do r=1, swib%num_rays

      ti = timediff(vold%year, vold%mon, vold%day,
     * ryib(r)%hour, ryib(r)%min, ryib(r)%sec, ryib(r)%msec,
     * cyr, cmo, cda, chr, cmn, cse, cms)

            beam_info(r,s,1) = ti
            beam_info(r,s,2) = ryib(r)%azimuth
            beam_info(r,s,3) = ryib(r)%elevation

          enddo       ! r=1, swib%num_rays
      enddo       ! s=1, nswp

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######              SUBROUTINE CARTESIAN_OBAN               ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine reads in radar data from sweep files and
c     interpolates the desired fields to a Cartesian grid.  Multiple options
c     are available for the 3D weighting function used by the interpolation scheme:
c
c     1.  Cressman weighting function:
c
c         wgt = ( 1.0 - deltax**2/hsp**2 - deltay**2/hsp**2 - deltaz**2/vsp**2 )
c             / ( 1.0 + deltax**2/hsp**2 + deltay**2/hsp**2 + deltaz**2/vsp**2 )
c
c
c     2.  Barnes weighting function:
c
c         wgt = exp ( -deltax**2/hsp - deltay**2/hsp - deltaz**2/vsp )
c
c############################################################################
c
c     Author:  David Dowell
c
c     Created:  December 2001
c
c     Modified:  February 2005
c
c############################################################################

      subroutine cartesian_oban(method, hsp, vsp,
     $                          f, nfld, fname, nswp, sfname,
     $                          extrapolate_flag, minsum, minrange,
     $                          yrcor, mocor, dacor,
     $                          azcor, elcor, az_corr_flag, umass_flag,
     $                          az, el, time, count,
     $                          map_proj, glat, glon, galt, rlat, rlon, ralt,
     $                          nx, ny, nz, dx, dy, dz, xmin, ymin, zmin,
     $                          cyr, cmo, cda, chr, cmn, cse, ut, vt)

      implicit none

      include 'dow.inc'
      include 'structures.inc'

c---- Passed variables

      integer method       ! interpolation method:
                           !   1=Cressman
                           !   2=Barnes
      real hsp, vsp        ! horizontal and smoothing parameters:
                           !   method=1:  hsp and vsp are the Cressman radii of influence (km)
                           !   method=2:  hsp and vsp are the Barnes smoothing parameters (km*km)
      integer nfld                    ! number of data fields to be gridded
      character(len=8) fname(nfld)    ! field names
      integer nswp                    ! number of sweep files
      character(len=200) sfname(nswp) ! sweep file names
      integer extrapolate_flag        ! should extrapolation be allowed?
                                      !   1=yes (standard objective analysis)
                                      !   0=no (interpolation only)
      real minsum                     ! threshold for minimum sum of weights required to produce
                                      !   an objectively-analyzed observation
      real minrange                   ! minimum-range threshold (data closer to radar
                                      !   are discarded)
      integer yrcor                   ! correction to year
      integer mocor                   ! correction to month
      integer dacor                   ! correction to day
      real elcor                      ! elevation angle correction (deg)
      real azcor                      ! azimuth-angle offset (deg)
      integer az_corr_flag            ! method of additional azimuth-angle correction
                                      !   0 = none
                                      !   1 = average current and previous angle
      integer umass_flag              ! apply UMass data corrections? (1=yes, 0=no)
      integer map_proj                ! map projection (for relating lat, lon to x, y):
                                      !   0 = flat earth
                                      !   1 = oblique azimuthal
                                      !   2 = Lambert conformal
      real glat, glon                 ! latitude and longitude of grid origin (deg)
      real galt                       ! altitude of grid origin (km MSL)
      real rlat, rlon                 ! radar latitude and longitude (deg)
      real ralt                       ! radar altitude (km MSL)
      integer nx, ny, nz              ! no. of grid points in x, y, and z directions
      real dx, dy, dz                 ! grid spacing in x, y, and z directions (km)
      real xmin, ymin, zmin           ! coordinates of lower southwest corner
                                      !   of grid, relative to origin (km)
      integer(kind=2) cyr,cmo,cda     ! central date
      integer(kind=2) chr,cmn,cse     ! central time
      real ut, vt                     ! storm translation velocity (m/s)

c---- Returned variables

      real f(nx,ny,nz,nfld)           ! data fields
      real az(nx,ny,nz)               ! interpolated azimuth angle (deg)
      real el(nx,ny,nz)               ! interpolated elevation angle (deg)
      real time(nx,ny,nz)             ! time, relative to central time (sec)
      integer count(nx,ny,nz,nfld)    ! no. of gates used in interpolation

c---- Local variables

      real glatr, glonr               ! latitude and longitude of grid origin (rad)
      real rlatr, rlonr               ! radar latitude and longitude (rad)
      integer s                       ! sweep file number
      integer n                       ! field number
      integer r                       ! ray number
      integer g                       ! gate number
      integer d                       ! directional bin number
      real x, y, z                    ! location of gate relative to grid origin (km)
      real ti                         ! time (s) relative to central time
      real timediff, comptoaz         ! functions used by this subroutine
      real rangekm_to_gate            ! "                               "
      integer(kind=2) cms; parameter(cms=0.0)
      real xrad, yrad, zrad           ! coords. (km) of radar with respect to origin
      integer i, j, k                 ! grid indices
      integer imin, imax              ! minimimum and maximum i for search
      integer jmin, jmax              ! minimimum and maximum j for search
      integer kmin, kmax              ! minimimum and maximum k for search
      integer irad, jrad, krad        ! radius of search area
      integer ii, jj, kk              ! nearest grid point to gate location
      real dh2                        ! square of horizontal distance (km*km)
      real dv2                        ! square of vertical distance (km*km)
      real wgt                        ! Cressman weight
      real c1                         ! hsp*hsp
      real c2                         ! vsp*vsp
      real range                      ! distance from radar to gate (km)
      integer valid_sum               ! 1 (0) if the sum of weights is (not) sufficient
                                      !   to produce an objectively analyzed observations

      type(vold_info)                     :: vold
      type(radd_info)                     :: radd
      type(celv_info)                     :: celv
      type(cfac_info)                     :: cfac
      type(parm_info), dimension(maxflds) :: parm
      type(swib_info)                     :: swib
      type(ryib_info), dimension(maxrays) :: ryib
      type(asib_info), dimension(maxrays) :: asib
      type(rdat_info), dimension(maxrays) :: rdat

      real, allocatable :: sumf(:,:,:,:)     ! sum of weights for each data field
      real, allocatable :: sumdir(:,:,:,:,:) ! sum of weights in directional bins
      real, allocatable :: azcomp(:,:,:,:)   ! azimuth angle components
      real, allocatable :: xg(:)             ! x coordinates (km) of grid points
      real, allocatable :: yg(:)             ! y coordinates (km) of grid points
      real, allocatable :: zg(:)             ! z coordinates (km) of grid points

c############################################################################
c
c     Initialize variables.
c
c############################################################################

      allocate(sumf(nx,ny,nz,0:nfld))
      if (extrapolate_flag.eq.0) then
        allocate(sumdir(8,nx,ny,nz,nfld))
      endif
      allocate(azcomp(nx,ny,nz,2))
      allocate(xg(nx))
      allocate(yg(ny))
      allocate(zg(nz))

      c1 = hsp*hsp
      c2 = vsp*vsp
      f(:,:,:,:) = 0.0
      azcomp(:,:,:,:) = 0.0
      el(:,:,:) = 0.0
      time(:,:,:) = 0.0
      count(:,:,:,:) = 0
      sumf(:,:,:,:) = 0.0
      if (extrapolate_flag.eq.0) then
        sumdir(:,:,:,:,:) = 0.0
      endif

      call grid_coordinates(xg, nx, dx, xmin)
      call grid_coordinates(yg, ny, dy, ymin)
      call grid_coordinates(zg, nz, dz, zmin)

      if (method.eq.1) then
        irad = 1+nint(hsp/dx)
        jrad = 1+nint(hsp/dy)
        krad = 1+nint(vsp/dz)
      else if (method.eq.2) then
        irad = 1+nint(sqrt(7.0*hsp)/dx)
        jrad = 1+nint(sqrt(7.0*hsp)/dy)
        krad = 1+nint(sqrt(7.0*vsp)/dz)
      else
        write(6,*) 'unknown method:  ', method
        stop
      endif

      glatr = dtor*glat
      glonr = dtor*glon
      rlatr = dtor*rlat
      rlonr = dtor*rlon
      call ll_to_xy(xrad, yrad, map_proj, glatr, glonr, rlatr, rlonr)
      zrad = ralt - galt
      write(6,*) 'location of radar relative to grid origin:'
      write(6,*) xrad, yrad, zrad

c############################################################################
c
c     Read in sweep files and compute weighted sums.
c
c############################################################################

      do s=1, nswp

        do n=1, nfld

          write(6,*) 'reading from ', sfname(s)
          call sweepread(sfname(s),vold,radd,celv,cfac,parm,swib,
     $                   ryib,asib,rdat,fname(n))

          call check_sweep_size(radd%num_param_desc, swib%num_rays)
          call correct_date(vold%year, yrcor, vold%mon, mocor, vold%day, dacor)
          call correct_az_el(ryib, swib%num_rays, azcor, elcor, az_corr_flag)

          if (umass_flag.eq.1) then
            call correct_umass_data(swib%num_rays, ryib, celv%total_gates, rdat)
          endif

          do r=1, swib%num_rays

            ti = timediff(vold%year, vold%mon, vold%day,
     $                    ryib(r)%hour, ryib(r)%min, ryib(r)%sec, ryib(r)%msec,
     $                    cyr, cmo, cda, chr, cmn, cse, cms)

c            write(6,*) ryib(r)%azimuth, ryib(r)%elevation
c            write(6,*) vold%year, vold%mon, vold%day
c            write(6,*) ryib(r)%hour, ryib(r)%min, ryib(r)%sec, ryib(r)%msec
c            write(6,*)

c            write(6,*) 'gate_spacing 1, 2, 3:  ',
c     $                 celv%gate_spacing(1),
c     $                 celv%gate_spacing(2),
c     $                 celv%gate_spacing(3)
c            write(6,*) 'elevation angle: ', ryib(r)%elevation

            do g=1, celv%total_gates

              if (rdat(r)%data(g) .ne. sbad) then

                range = rangekm_to_gate(celv, g)
                call xyzloc(x, y, z, range, dtor*ryib(r)%azimuth, dtor*ryib(r)%elevation,
     $                      map_proj, glatr, glonr, galt, rlatr, rlonr, ralt,
     $                      ut, vt, ti)

                ii = 1 + nint((x-xmin)/dx)
                jj = 1 + nint((y-ymin)/dy)
                kk = 1 + nint((z-zmin)/dz)

                imin = max(ii-irad, 1)
                imax = min(ii+irad, nx)
                jmin = max(jj-jrad, 1)
                jmax = min(jj+jrad, ny)
                kmin = max(kk-krad, 1)
                kmax = min(kk+krad, nz)

                if (range.ge.minrange) then

                  do k=kmin, kmax
                    do j=jmin, jmax
                      do i=imin, imax

                        if (method.eq.1) then
                          dh2 = (x-xg(i))*(x-xg(i)) + (y-yg(j))*(y-yg(j))
                          dv2 = (z-zg(k))*(z-zg(k))
                          wgt = (1.0 - dh2/c1 - dv2/c2) / (1.0 + dh2/c1 + dv2/c2)
                        else if (method.eq.2) then
                          wgt = exp( -(x-xg(i))*(x-xg(i))/hsp
     $                               -(y-yg(j))*(y-yg(j))/hsp
     $                               -(z-zg(k))*(z-zg(k))/vsp )
                        endif

                        if (wgt.gt.0.0) then
                          f(i,j,k,n) = f(i,j,k,n)
     $                               + wgt*rdat(r)%data(g)
                          sumf(i,j,k,n) = sumf(i,j,k,n) + wgt
                          count(i,j,k,n) = count(i,j,k,n) + 1
                          if (extrapolate_flag.eq.0) then
                            call update_dir_bins(sumdir(1,i,j,k,n), wgt,
     $                                           x, y, z, xg(i), yg(j), zg(k))
                          endif
                          azcomp(i,j,k,1) = azcomp(i,j,k,1)
     $                          + wgt*sin(dtor*ryib(r)%azimuth)
                          azcomp(i,j,k,2) = azcomp(i,j,k,2)
     $                          + wgt*cos(dtor*ryib(r)%azimuth)
                          el(i,j,k) = el(i,j,k)
     $                          + wgt*ryib(r)%elevation
                          time(i,j,k) = time(i,j,k) + wgt*ti
                          sumf(i,j,k,0) = sumf(i,j,k,0) + wgt
                        endif   ! if (wgt.gt.0.0)

                      enddo
                    enddo
                  enddo

                endif       ! (range.ge.minrange)

              endif       ! (rdat(r)%data(g) .ne. sbad)

            enddo       ! g=1, celv%total_gates
          enddo       ! r=1, swib%num_rays
        enddo       ! n=1, nfld

      enddo       ! s=1, nswp

c############################################################################
c
c     Compute interpolated values from weighted sums.
c
c############################################################################

      do k=1, nz
        do j=1, ny
          do i=1, nx

            do n=1, nfld

              valid_sum = 1
              if (sumf(i,j,k,n).le.minsum) then
                valid_sum = 0
              endif
              if (extrapolate_flag.eq.0) then
                do d=1, 8
                  if (sumdir(d,i,j,k,n).le.(0.1*minsum) ) then
                    valid_sum = 0
                  endif
                enddo
              endif

              if (valid_sum.eq.1) then
                f(i,j,k,n) = f(i,j,k,n) / sumf(i,j,k,n)
              else 
                f(i,j,k,n) = bad
              endif

            enddo

            if (sumf(i,j,k,0) .gt. minsum) then
              azcomp(i,j,k,1) = azcomp(i,j,k,1) / sumf(i,j,k,0)
              azcomp(i,j,k,2) = azcomp(i,j,k,2) / sumf(i,j,k,0)
              az(i,j,k) = comptoaz(azcomp(i,j,k,1),azcomp(i,j,k,2))
              el(i,j,k) = el(i,j,k) / sumf(i,j,k,0)
              time(i,j,k) = time(i,j,k) / sumf(i,j,k,0)
            else 
              az(i,j,k) = bad
              el(i,j,k) = bad
              time(i,j,k) = bad
            endif

          enddo
        enddo
      enddo

      deallocate(sumf)
      if (extrapolate_flag.eq.0) then
        deallocate(sumdir)
      endif
      deallocate(azcomp)
      deallocate(xg)
      deallocate(yg)
      deallocate(zg)

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                 SUBROUTINE PPI_OBAN                  ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine reads in radar data from sweep files and
c     interpolates the desired fields to a "semi-Cartesian" grid;
c     i.e., the horizontal coordinates of the grid are Cartesian, and
c     the vertical coordinate is the scan number.  Height is one
c     of the computed fields.
c
c     Multiple options are available for the 2D weighting function
c     used by the interpolation scheme:
c
c     1.  Cressman weighting function:
c
c         wgt = ( 1.0 - deltax**2/hsp**2 - deltay**2/hsp**2 )
c             / ( 1.0 + deltax**2/hsp**2 + deltay**2/hsp**2 )
c
c     2.  Barnes weighting function:
c
c         wgt = exp ( -deltax**2/hsp - deltay**2/hsp )
c
c############################################################################
c
c     Author:  David Dowell
c
c     Created:  December 2001
c
c     Modified:  February 2005
c
c############################################################################

      subroutine ppi_oban(method, hsp,
     $                    f, nfld, fname, sfname,
     $                    minsum, minrange,
     $                    yrcor, mocor, dacor,
     $                    azcor, elcor, az_corr_flag, umass_flag,
     $                    az, el, height, time, count,
     $                    map_proj, glat, glon, galt, rlat, rlon, ralt,
     $                    nx, ny, nswp, dx, dy, xmin, ymin,
     $                    cyr, cmo, cda, chr, cmn, cse, ut, vt,
     $                    days, secs)

      implicit none

      include 'dow.inc'
      include 'structures.inc'

c---- Passed variables

      integer method       ! interpolation method:
                           !   1=Cressman
                           !   2=Barnes
      real hsp             ! horizontal smoothing parameter:
                           !   method=1:  hsp is the Cressman radii of influence (km)
                           !   method=2:  hsp is the Barnes smoothing parameters (km*km)
      integer nfld                    ! number of data fields to be gridded
      character(len=8) fname(nfld)    ! field names
      integer nswp                    ! number of sweep files
      character(len=200) sfname(nswp) ! sweep file names
      real minsum                     ! threshold for minimum sum of weights required to produce
                                      !   an objectively-analyzed observation
      real minrange                   ! minimum-range threshold (data closer to radar
                                      !   are discarded)
      integer yrcor                   ! correction to year
      integer mocor                   ! correction to month
      integer dacor                   ! correction to day
      real azcor                      ! azimuth-angle offset (deg)
      real elcor                      ! elevation angle correction (deg)
      integer az_corr_flag            ! method of additional azimuth-angle correction
                                      !   0 = none
                                      !   1 = average current and previous angle
      integer umass_flag              ! apply UMass data corrections? (1=yes, 0=no)
      integer map_proj                ! map projection:
                                      !   0 = flat earth
                                      !   1 = oblique azimuthal
                                      !   2 = Lambert conformal
      real glat, glon                 ! latitude and longitude of grid origin (deg)
      real galt                       ! altitude of grid origin (km MSL)
      real rlat, rlon                 ! radar latitude and longitude (deg)
      real ralt                       ! radar altitude (km MSL)
      integer nx, ny                  ! no. of grid points in x, y, and z directions
      real dx, dy                     ! grid spacing in x and y directions (km)
      real xmin, ymin                 ! coordinates of lower southwest corner
                                      !   of grid, relative to origin (km)
      integer(kind=2) cyr,cmo,cda     ! central date
      integer(kind=2) chr,cmn,cse     ! central time
      real ut, vt                     ! storm translation velocity (m/s)

c---- Returned variables

      real f(nx,ny,nswp,nfld)         ! data fields
      real az(nx,ny,nswp)             ! interpolated azimuth angle (deg)
      real el(nx,ny,nswp)             ! interpolated elevation angle (deg)
      real time(nx,ny,nswp)           ! time, relative to central time (sec)
      real height(nx,ny,nswp)         ! height of ob, relative to grid origin (km)
      integer count(nx,ny,nswp,nfld)  ! no. of gates used in interpolation
      integer days(nswp)              ! Gregorian day (since beginning of base year) of first beam in sweep
      integer secs(nswp)              ! seconds of first beam in sweep

c---- Local variables

      real glatr, glonr               ! latitude and longitude of grid origin (rad)
      real rlatr, rlonr               ! radar latitude and longitude (rad)
      integer s                       ! sweep file number
      integer n                       ! field number
      integer r                       ! ray number
      integer g                       ! gate number
      real x, y, z                    ! location of gate relative to grid origin (km)
      real ti                         ! time (s) relative to central time
      real timediff, comptoaz         ! functions used by this subroutine
      real rangekm_to_gate            ! "                               "
      integer(kind=2) cms; parameter(cms=0.0)
      real xrad, yrad, zrad           ! coords. of radar with respect to origin (km)
      integer i, j                    ! grid indices
      integer imin, imax              ! minimimum and maximum i for search
      integer jmin, jmax              ! minimimum and maximum j for search
      integer irad, jrad              ! radius of search area
      integer ii, jj                  ! nearest grid point to gate location
      real dh2                        ! square of horizontal distance (km*km)
      real wgt                        ! Cressman weight
      real c1                         ! hsp*hsp
      real range                      ! distance from radar to gate (km)
      integer num_valid_obs           ! total number of valid raw observations

      type(vold_info)                     :: vold
      type(radd_info)                     :: radd
      type(celv_info)                     :: celv
      type(cfac_info)                     :: cfac
      type(parm_info), dimension(maxflds) :: parm
      type(swib_info)                     :: swib
      type(ryib_info), dimension(maxrays) :: ryib
      type(asib_info), dimension(maxrays) :: asib
      type(rdat_info), dimension(maxrays) :: rdat

      real, allocatable :: sumf(:,:,:,:)   ! sum of weights for each data field
      real, allocatable :: azcomp(:,:,:,:) ! azimuth angle components
      real, allocatable :: xg(:)           ! x coordinates (km) of grid points
      real, allocatable :: yg(:)           ! y coordinates (km) of grid points

c############################################################################
c
c     Initialize variables.
c
c############################################################################

      allocate(sumf(nx,ny,nswp,0:nfld))
      allocate(azcomp(nx,ny,nswp,2))
      allocate(xg(nx))
      allocate(yg(ny))

      c1 = hsp*hsp
      f(:,:,:,:) = 0.0
      azcomp(:,:,:,:) = 0.0
      el(:,:,:) = 0.0
      height(:,:,:) = 0.0
      time(:,:,:) = 0.0
      count(:,:,:,:) = 0
      sumf(:,:,:,:) = 0.0
      num_valid_obs = 0

      call grid_coordinates(xg, nx, dx, xmin)
      call grid_coordinates(yg, ny, dy, ymin)

      if (method.eq.1) then
        irad = 1+nint(hsp/dx)
        jrad = 1+nint(hsp/dy)
      else if (method.eq.2) then
        irad = 1+nint(sqrt(7.0*hsp)/dx)
        jrad = 1+nint(sqrt(7.0*hsp)/dy)
      else
        write(6,*) 'unknown method:  ', method
        stop
      endif

      glatr = dtor*glat
      glonr = dtor*glon
      rlatr = dtor*rlat
      rlonr = dtor*rlon
      call ll_to_xy(xrad, yrad, map_proj, glatr, glonr, rlatr, rlonr)
      zrad = ralt - galt
      write(6,*) 'location of radar relative to grid origin:'
      write(6,*) xrad, yrad, zrad

c############################################################################
c
c     Read in sweep files and compute weighted sums.
c
c############################################################################

      do s=1, nswp

        do n=1, nfld

          write(6,*) 'reading from ', sfname(s)
          call sweepread(sfname(s),vold,radd,celv,cfac,parm,swib,
     $                   ryib,asib,rdat,fname(n))

          call check_sweep_size(radd%num_param_desc, swib%num_rays)
          call correct_date(vold%year, yrcor, vold%mon, mocor, vold%day, dacor)
          call correct_az_el(ryib, swib%num_rays, azcor, elcor, az_corr_flag)

          if (umass_flag.eq.1) then
            call correct_umass_data(swib%num_rays, ryib, celv%total_gates, rdat)
          endif

          if (n.eq.1) then
            call set_date_gregorian(days(s), secs(s),
     $                              int(vold%year), int(vold%mon), int(vold%day),
     $                              int(ryib(1)%hour), int(ryib(1)%min), int(ryib(1)%sec))
          endif

          do r=1, swib%num_rays

            ti = timediff(vold%year, vold%mon, vold%day,
     $                    ryib(r)%hour, ryib(r)%min, ryib(r)%sec, ryib(r)%msec,
     $                    cyr, cmo, cda, chr, cmn, cse, cms)

c            write(6,*) 'gate_spacing 1, 2, 3:  ',
c     $                 celv%gate_spacing(1),
c     $                 celv%gate_spacing(2),
c     $                 celv%gate_spacing(3)
c            write(6,*) 'elevation angle: ', ryib(r)%elevation

            do g=1, celv%total_gates

              if (rdat(r)%data(g) .ne. sbad) then
                num_valid_obs = num_valid_obs + 1
              endif

              range = rangekm_to_gate(celv, g)
              call xyzloc(x, y, z, range, dtor*ryib(r)%azimuth, dtor*ryib(r)%elevation,
     $                    map_proj, glatr, glonr, galt, rlatr, rlonr, ralt,
     $                    ut, vt, ti)

              ii = 1 + nint((x-xmin)/dx)
              jj = 1 + nint((y-ymin)/dy)

              imin = max(ii-irad, 1)
              imax = min(ii+irad, nx)
              jmin = max(jj-jrad, 1)
              jmax = min(jj+jrad, ny)

              if (range.ge.minrange) then

                do j=jmin, jmax
                  do i=imin, imax

                    if (method.eq.1) then
                      dh2 = (x-xg(i))*(x-xg(i)) + (y-yg(j))*(y-yg(j))
                      wgt = (1.0 - dh2/c1) / (1.0 + dh2/c1)
                    else if (method.eq.2) then
                      wgt = exp( -(x-xg(i))*(x-xg(i))/hsp - (y-yg(j))*(y-yg(j))/hsp )
                    endif

                    if (wgt.gt.0.0) then
                      if (rdat(r)%data(g) .ne. sbad) then
                        f(i,j,s,n) = f(i,j,s,n)
     $                             + wgt*rdat(r)%data(g)
                        sumf(i,j,s,n) = sumf(i,j,s,n) + wgt
                        count(i,j,s,n) = count(i,j,s,n) + 1
                      endif
                      azcomp(i,j,s,1) = azcomp(i,j,s,1)
     $                       + wgt*sin(dtor*ryib(r)%azimuth)
                      azcomp(i,j,s,2) = azcomp(i,j,s,2)
     $                       + wgt*cos(dtor*ryib(r)%azimuth)
                      el(i,j,s) = el(i,j,s) + wgt*ryib(r)%elevation
                      height(i,j,s) = height(i,j,s) + wgt*z
                      time(i,j,s) = time(i,j,s) + wgt*ti
                      sumf(i,j,s,0) = sumf(i,j,s,0) + wgt
                    endif    ! if (wgt.gt.0.0)

                  enddo
                enddo

              endif       ! (range.ge.minrange)

            enddo       ! g=1, celv%total_gates
          enddo       ! r=1, swib%num_rays
        enddo       ! n=1, nfld

      enddo       ! s=1, nswp


c############################################################################
c
c     Compute interpolated values from weighted sums.
c
c############################################################################

      do s=1, nswp
        do j=1, ny
          do i=1, nx

            do n=1, nfld
              if (sumf(i,j,s,n) .gt. minsum) then
                f(i,j,s,n) = f(i,j,s,n) / sumf(i,j,s,n)
              else 
                f(i,j,s,n) = bad
              endif
            enddo

            if (sumf(i,j,s,0) .gt. minsum) then
              azcomp(i,j,s,1) = azcomp(i,j,s,1) / sumf(i,j,s,0)
              azcomp(i,j,s,2) = azcomp(i,j,s,2) / sumf(i,j,s,0)
              az(i,j,s) = comptoaz(azcomp(i,j,s,1),azcomp(i,j,s,2))
              el(i,j,s) = el(i,j,s) / sumf(i,j,s,0)
              height(i,j,s) = height(i,j,s) / sumf(i,j,s,0)
              time(i,j,s) = time(i,j,s) / sumf(i,j,s,0)
            else
              az(i,j,s) = bad
              el(i,j,s) = bad
              height(i,j,s) = bad
              time(i,j,s) = bad
            endif

          enddo
        enddo
      enddo

      write(6,*)
      write(6,*) 'total number of valid raw observations processed = ', num_valid_obs
      write(6,*)

      deallocate(sumf)
      deallocate(azcomp)
      deallocate(xg)
      deallocate(yg)

      return
      end


c###########################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                 SUBROUTINE XYZLOC                    ######
c     ######                                                      ######
c     ##################################################################
c
c     PURPOSE:
c
c     This subroutine computes the projected (x, y, z) coordinates of a
c     radar observation, relative to the grid origin.  The observation
c     location is adjusted for storm motion.
c
c############################################################################
c
c     Author:  David Dowell
c
c     Created:  December 2001
c
c     Modified:  February 2005
c
c############################################################################

      subroutine xyzloc(x, y, z, r, az0, el0, 
     $                  map_proj, glat, glon, galt, rlat, rlon, ralt,
     $                  ut, vt, ti)

      implicit none

      include 'dow.inc'

c---- Passed variables

      real r                   ! slant-path range (km) to observation
      real az0                 ! azimuth angle (rad) at radar
      real el0                 ! elevation angle (rad) at radar
      integer map_proj         ! map projection:
                               !   0 = flat earth
                               !   1 = oblique azimuthal
                               !   2 = Lambert conformal
      real glat, glon          ! latitude and longitude of grid origin (rad)
      real galt                ! altitude of grid origin (km MSL)
      real rlat, rlon          ! radar latitude and longitude (rad)
      real ralt                ! radar altitude (km MSL)
      real ti                  ! time (s) relative to central time
      real ut, vt              ! storm translation velocity (m/s)

c---- Returned variables

      real x, y, z             ! location of gate relative to grid origin (km)

c---- Local variables

      real h                   ! height (km) of observation relative to radar
      real s                   ! great-circle distance (km) along spherical earth
      real d                   ! great-circle distance, expressed in radians
      real lat, lon            ! observation latitude and longitude (rad)
      real dlon                ! longitude difference (rad)


c     observation height and great-circle distance [Doviak and Zrnic 1993, p. 21]

      h = sqrt( r*r + eer*eer + 2.0*r*eer*sin(el0) ) - eer
      z = h + ralt - galt
      s = eer * asin( (r*cos(el0)) / (eer+h) )
      d = s / rearth

c     observation (lat, lon) [E. Williams, "Aviation Formulary V1.42"]

      lat = asin(sin(rlat)*cos(d) + cos(rlat)*sin(d)*cos(az0))
      dlon = atan2(sin(az0)*sin(d)*cos(rlat), cos(d)-sin(rlat)*sin(lat))
      lon = mod(rlon+dlon+pi, 2.0*pi) - pi

c     observation (x, y)

      call ll_to_xy(x, y, map_proj, glat, glon, lat, lon)

c     storm-motion adjustment

      x = x - ti*(ut/1000.0)
      y = y - ti*(vt/1000.0)

      return
      end


c###########################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                 SUBROUTINE XYZLOC_OLD                ######
c     ######                                                      ######
c     ##################################################################
c
c     PURPOSE:
c
c     This subroutine computes the (x, y, z) location of a gate,
c     relative to the grid origin.  The location is adjusted
c     for storm motion.
c
c     THIS SUBROUTINE IS NOW OBSOLETE.
c
c############################################################################

      subroutine xyzloc_old(x, y, z, xrad, yrad, zrad,
     $                      az, el, d, ut, vt, ti)

      implicit none

      include 'dow.inc'

      real xrad, yrad, zrad    ! coords. of radar with respect to origin (km)
      real ut, vt              ! storm translation velocity (m/s)
      real az                  ! azimuth angle (deg)
      real el                  ! elevation angle (deg)
      real x, y, z             ! location of gate relative to grid origin (km)
      real ti                  ! time (s) relative to central time
      real d                   ! distance from radar to gate (km)

      x = xrad + d*sin(az*dtor)*cos(el*dtor)
      y = yrad + d*cos(az*dtor)*cos(el*dtor)
      z = zrad + sqrt( eer*eer + d*d 
     $                 + 2.0*d*eer*sin(el*dtor) )
     $         - eer

      x = x - ti*(ut/1000.0)
      y = y - ti*(vt/1000.0)

      return
      end


c###########################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                REAL FUNCTION TLINT                   ######
c     ######                                                      ######
c     ##################################################################
c
c     PURPOSE:
c
c     This function returns the tri-linearly interpolated value of a scalar
c     quantity at the given (x,y,z) location.
c
c############################################################################

      real function tlint(v, x, y, z, nx, ny, nz, dx, dy, dz,
     $                    xlsw, ylsw, zlsw)

      implicit none

      include 'dow.inc'

      integer nx, ny, nz                 ! number of gridpoints in x, y, and
                                         !   z directions
      real v(nx,ny,nz)                   ! scalar array
      real x, y, z                       ! parcel location
      real dx, dy, dz                    ! grid spacings in x/y/z
      real xlsw, ylsw, zlsw              ! coordinates of lower southwest
                                         !   corner of grid
      integer i, j, k                    ! indices of lower southwest corner of
                                         !   gridbox in which parcel is located
      real xi, yj, zk                    ! coordinates corresponding to these
                                         !   indices
      real xb, yb, zb                    ! coordinates of parcel with respect
                                         !   to grid box corner
      real q1, q2, vb, vt                ! temporary quantities

      i = aint(1.0 + (x-xlsw)/dx )
      j = aint(1.0 + (y-ylsw)/dy )
      k = aint(1.0 + (z-zlsw)/dz )
      xi = xlsw + (i-1)*dx
      yj = ylsw + (j-1)*dy
      zk = zlsw + (k-1)*dz
      xb = x - xi
      yb = y - yj
      zb = z - zk

      if ((i.ge.1).and.(j.ge.1).and.(k.ge.1).and.
     $    (i.lt.nx).and.(j.lt.ny).and.(k.lt.nz)) then
        if ((v(i,j,k).eq.BAD) .or. (v(i+1,j,k).eq.BAD) .or.
     $      (v(i+1,j+1,k).eq.BAD) .or. (v(i,j+1,k).eq.BAD) .or.
     $      (v(i,j,k+1).eq.BAD) .or. (v(i+1,j,k+1).eq.BAD) .or.
     $      (v(i+1,j+1,k+1).eq.BAD) .or. (v(i,j+1,k+1).eq.BAD)) then
          tlint = BAD
        else
          q1 = ((dx-xb)*v(i,j,k)+xb*v(i+1,j,k)) / (dx*dy)
          q2 = ((dx-xb)*v(i,j+1,k)+xb*v(i+1,j+1,k)) / (dx*dy)
          vb = q1*(dy-yb) + q2*yb
          q1 = ((dx-xb)*v(i,j,k+1)+xb*v(i+1,j,k+1)) / (dx*dy)
          q2 = ((dx-xb)*v(i,j+1,k+1)+xb*v(i+1,j+1,k+1)) / (dx*dy)
          vt = q1*(dy-yb) + q2*yb
          tlint = (vb*(dz-zb)+vt*zb)/dz
        endif
      else
        tlint = BAD
      endif

      return
      end


c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######              REAL FUNCTION TIMEDIFF                  ######
c     ######                                                      ######
c     ##################################################################
c
c     PURPOSE:
c
c     This function computes the difference (in seconds) between
c     two times (time1 - time2).
c
c     Note:  the possiblity of a leap year is not considered.
c
c############################################################################

      real function timediff(y1, mo1, d1, h1, mi1, s1, ms1,
     $                       y2, mo2, d2, h2, mi2, s2, ms2)

      implicit none

      include 'dow.inc'

      integer(kind=2) y1, y2                ! year
      integer(kind=2) mo1, mo2              ! month
      integer(kind=2) d1, d2                ! day
      integer(kind=2) h1, h2                ! hour
      integer(kind=2) mi1, mi2              ! minute
      integer(kind=2) s1, s2                ! second
      integer(kind=2) ms1, ms2              ! milliseconds
      real(kind=8) deltat
      integer jd1, jd2                      ! julian days


      jd1 = d1 + pd(mo1)
      jd2 = d2 + pd(mo2)

      deltat = 31536000.0 * (y1-y2)
      deltat = deltat + 86400.0 * (jd1-jd2)
      deltat = deltat + 3600.0 * (h1-h2)
      deltat = deltat + 60.0 * (mi1-mi2)
      deltat = deltat + (s1-s2)
      deltat = deltat + 0.001*(ms1-ms2)

      timediff = deltat

      return
      end

  
c###########################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######           SUBROUTINE SET_DATE_GREGORIAN              ######
c     ######                                                      ######
c     ##################################################################
c
c     PURPOSE:
c
c     Computes time corresponding to date for gregorian calendar.
c
c############################################################################
c
c     Author:  Data Assimilation Research Testbed -- DART
C              Data Assimilation Initiative, University Corporation for Atmospheric Research
c
c     Created:  September 9, 2004
c
c     Modified:  May 31, 2005 (David Dowell)
c
c############################################################################

      subroutine set_date_gregorian(days, secs, year, month, day, hours, minutes, seconds)

      implicit none

!---- returned variables

      integer :: days                  ! Gregorian day since beginning of base year
      integer :: secs                  ! seconds since beginning of day

!---- input variables

      integer :: day, month, year
      integer :: seconds, minutes, hours

!---- local variables

      integer :: ndays, m, nleapyr
      integer :: base_year = 1601
      logical :: leap
      integer :: days_per_month(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)

      ! Need to check for bogus dates

      if ( seconds > 59 .or. seconds  < 0 .or.
     $     minutes > 59 .or. minutes  < 0 .or.
     $     hours   > 23 .or. hours    < 0 .or.
     $                       day      < 1 .or.
     $     month   > 12 .or. month    < 1 .or.
     $                       year     < base_year ) then
        write(6,*) 'set_date_gregorian:  s,m,h,d,mn,y',
     $             seconds,minutes,hours,day,month,year,
     $             ' not a valid date'
        stop
      endif

      if (month /= 2 .and. day > days_per_month(month)) then
         write(6,*) 'month (',month,') does not have ',day,' days'
         stop
      endif

      ! Is this a leap year? Gregorian calandar assigns each year evenly
      ! divisible by 4 that is not a century year unevenly divisible by 400
      ! as a leap-year. (i.e. 1700,1800,1900 are not leap-years, 2000 is)

      leap=(modulo(year,4).eq.0)
      if((modulo(year,100).eq.0).and.(modulo(year,400).ne.0))then
       leap=.false.
      endif

      ! Finish checking for day specification errors

      if (month == 2 .and. (day > 29 .or. ((.not. leap) .and. day > 28))) then
        write(6,*) 'month (',month,') does not have ',
     $             day,' days in a lon-leap year'
        stop
      endif

      ! Compute number of leap years fully past since base_year

      nleapyr = (year - base_year) / 4 - (year - base_year) / 100 + (year - base_year) / 400

      ! Count up days in this year

      ndays = 0
      do m=1,month-1
       ndays = ndays + days_per_month(m)
       if(leap .and. m == 2) ndays = ndays + 1
      enddo

      secs = seconds + 60*(minutes + 60 * hours)
      days = day - 1 +  ndays + 365*(year - base_year - nleapyr) + 366*(nleapyr)

      return
      end subroutine set_date_gregorian  



c###########################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######                 SUBROUTINE FALLSPEED                 ######
c     ######                                                      ######
c     ##################################################################
c
c     PURPOSE:
c
c     This subroutine computes an estimate of the precipitation fallspeed
c     from an empirial relationship involving reflectivity and density.
c
c############################################################################

      subroutine fallspeed(wt, rf, rho, nx, ny, nz)

      implicit none

      include 'dow.inc'

      integer nx, ny, nz        ! no. of grid points in x, y, and z directions
      real rf(nx,ny,nz)         ! reflectivity (dBZ)
      real rho(nz)              ! density (kg m**-3)
      real wt(nx,ny,nz)         ! precipitation terminal fall speed (m/s)
      integer i, j, k           ! grid indices

      do k=1, nz
        do j=1, ny
          do i=1, nx
            if (rf(i,j,k).eq.bad) then
              wt(i,j,k) = bad
            else
              wt(i,j,k) = -2.6 * (10.0**(0.1*rf(i,j,k)))**0.107
     $                  * (1.2/rho(k))**0.4
            endif
          enddo
        enddo
      enddo

      return
      end



























      subroutine multi_oban(method, gamma, npass, hsp0, vsp0, mincount,
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



c        call multi_oban(method, gamma, npass, hsp0, vsp0, mincount,
c     $      f, f_mpass, nfld, fname, velfld, reffld,
c     $           fsan, signf, nswp, maxswp, sfname,
c     $             iradtype, extrapolate_flag, minsum, minrange, rhmax,
c     $                   yrcor, mocor, dacor, minswpdc, maxswpdc,
c     $                   azcor, elcor, az_corr_flag, umass_flag,
c     $                   az, el, time, xrd, yrd, zrd, count,
c     $                   map_proj, glat, glon, galt, rlat, rlon, ralt,
c     $                   nx, ny, nz, dx, dy, dz, xmin, ymin, zmin,
c     $                   cyr, cmo, cda, chr, cmn, cse, ut, vt)





c
c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######              SUBROUTINE multi_oban                   ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine reads in radar data from sweep files and
c     interpolates the desired fields to a Cartesian grid.  Multiple options
c     are available for the 3D weighting function used by the interpolation scheme:
c
c     1.  Cressman weighting function:
c
c         wgt = ( 1.0 - deltax**2/hsp**2 - deltay**2/hsp**2 - deltaz**2/vsp**2 )
c             / ( 1.0 + deltax**2/hsp**2 + deltay**2/hsp**2 + deltaz**2/vsp**2 )
c
c
c     2.  Barnes weighting function:
c
c         wgt = exp ( -deltax**2/hsp - deltay**2/hsp - deltaz**2/vsp )
c
c############################################################################
c
c     Author of initial code:  Mario Majcen
c
c     Created:  January 2007
c
c     CLZ (9/22/09):  fixed time difference function call bug for npass > 1
c
c      CLZ (2/17/10):  Markowski's bug fix, maxgate = ngate(1) --> maxgate = ngate(s)
c
c     CLZ (12/18/13): restored option to control extrapolation in all 3 coordinate directions
c                     (call update_dir_bins, located in util.f).  Must read in/set
c                       extrapolate_flag = 0 (do not extrapolate) or = 1 (extrapolate).
c
c     CLZ (10/30/2018): Conrad Ziegler comp[letes addition of airborne radar analysis option.
c
c############################################################################
c


      implicit none

      include 'dow.inc'
      include 'structures.inc'

c---- Passed variables

      integer iradtype

      integer method       ! interpolation method:
                           !   1=Cressman
                           !   2=Barnes
      
      
      integer nfld                    ! number of data fields to be gridded

      character(len=8) fname(nfld)    ! field names
      character(len=8) velfld         ! name of velocity field to be analyzed
      character(len=8) reffld         ! name of reflectivity field to be analyzed

      real fsan(nfld)                 ! field sanity check value
      real signf(nfld)                ! field sign scale value

      integer nswp                    ! number of sweep files
      integer maxswp                  ! dimension of n-vector sweep arrays
      character(len=200) sfname(maxswp) ! sweep file names
c      character(len=200) sfname(nswp) ! sweep file names
      integer extrapolate_flag        ! should extrapolation be allowed?
                                      !   1=yes (standard objective analysis)
                                      !   0=no (interpolation only)
      integer mincount            ! threshold for minimum number of data gates required
                                  !   to produce an objectively-analyzed observation
      real minsum                     ! threshold for minimum sum of weights required to produce
                                      !   an objectively-analyzed observation
      real minrange                   ! minimum-range threshold (data closer to radar
                                      !   are discarded)
  
c
c.... CLZ (11/23/10): rhmax is optional maximum allowable horizontal range to oban a datum
c
c.... CLZ (5/11/18): do no threshold if maxrh = 0; only threshold if maxrh = 1
c

      real rhmax                ! rhmax (km) = cos(dtor * elmax) * rmax
      real hrange               ! horizontal range to gate (km)

      real hpbwlen
      real dbzclut
      real vclut
      real angclut
      real angclut0
      real dangdr
      real dtr
      parameter (dtr = 0.0174532925)

      integer yrcor                   ! correction to year
      integer mocor                   ! correction to month
      integer dacor                   ! correction to day
      integer minswpdc            ! begin sweep for dacor
      integer maxswpdc            ! end sweep for dacor
      real elcor                      ! elevation angle correction (deg)
      real azcor                      ! azimuth-angle offset (deg)
      integer az_corr_flag            ! method of additional azimuth-angle correction
                                      !   0 = none
                                      !   1 = average current and previous angle
      integer umass_flag              ! apply UMass data corrections? (1=yes, 0=no)
      integer map_proj                ! map projection (for relating lat, lon to x, y):
                                      !   0 = flat earth
                                      !   1 = oblique azimuthal
                                      !   2 = Lambert conformal
      real glat, glon                 ! latitude and longitude of grid origin (deg)
      real galt                       ! altitude of grid origin (km MSL)
      real rlat, rlon                 ! radar latitude and longitude (deg)
      real ralt                       ! radar altitude (km MSL)
      integer nx, ny, nz              ! no. of grid points in x, y, and z directions
      real dx, dy, dz                 ! grid spacing in x, y, and z directions (km)
      real xmin, ymin, zmin           ! coordinates of lower southwest corner
                                      !   of grid, relative to origin (km)
      integer(kind=2) cyr,cmo,cda     ! central date
      integer(kind=2) chr,cmn,cse     ! central time
      real ut, vt                     ! storm translation velocity (m/s)

c---- Returned variables

      real f(nx,ny,nz,nfld+3)             ! data fields
      real f_mpass(nx,ny,nz,nfld+3,npass) ! multi-pass data fields

      real az(nx,ny,nz)                   ! interpolated azimuth angle (deg)
      real el(nx,ny,nz)                   ! interpolated elevation angle (deg)
      real time(nx,ny,nz)                 ! time, relative to central time (sec)

      real xrd(nx,ny,nz)                 ! interpolated airborne radar x-coordinate (km)
      real yrd(nx,ny,nz)                 ! interpolated airborne radar y-coordinate (km)
      real zrd(nx,ny,nz)                 ! interpolated airborne radar z-coordinate (km)

      integer count(nx,ny,nz,nfld)        ! no. of gates used in interpolation


c---- Local variables

      real hsp, vsp   ! horizontal and smoothing parameters:
                      !   method=1:  hsp and vsp are the Cressman radii of influence (km)
                      !   method=2:  hsp and vsp are the Barnes smoothing parameters (km*km)  
      real c1                         ! hsp*hsp
      real c2                         ! vsp*vsp



      real rot_raw
      real tilt
      real roll
      real rot_rc
      real pitc
      real drft
      real head
      real trck
      real azim_north
      real elev_hor
      real elev_zen
      real tilt_tr

      real azimang
      real elevang


      real glatr, glonr               ! latitude and longitude of grid origin (rad)
      real rlatr, rlonr               ! radar latitude and longitude (rad)
      integer s                       ! sweep file number
      integer n                       ! field number
      integer r                       ! ray number
      integer g                       ! gate number
      integer gg                      ! median filter window element number
      integer d                       ! directional bin number
      real x, y, z                    ! location of gate relative to grid origin (km)
      real ti                         ! time (s) relative to central time
      real timediff, comptoaz         ! functions used by this subroutine
      real rangekm_to_gate            ! "                               "
      integer(kind=2) cms; parameter(cms=0.0)
      real xrad, yrad, zrad           ! coords. (km) of radar with respect to origin
      integer i, j, k                 ! grid indices
      integer imin, imax              ! minimimum and maximum i for search
      integer jmin, jmax              ! minimimum and maximum j for search
      integer kmin, kmax              ! minimimum and maximum k for search
      integer irad, jrad, krad        ! radius of search area
      integer ii, jj, kk              ! nearest grid point to gate location
      real dh2                        ! square of horizontal distance (km*km)
      real dv2                        ! square of vertical distance (km*km)
      real wgt                        ! Cressman weight

      real range                      ! distance from radar to gate (km)
      integer valid_sum               ! 1 (0) if the sum of weights is (not) sufficient
                                      !   to produce an objectively analyzed observations

      type(vold_info)                     :: vold
      type(radd_info)                     :: radd
      type(celv_info)                     :: celv
      type(cfac_info)                     :: cfac
      type(parm_info), dimension(maxflds) :: parm
      type(swib_info)                     :: swib
      type(ryib_info), dimension(maxrays) :: ryib
      type(asib_info), dimension(maxrays) :: asib
      type(rdat_info), dimension(maxrays) :: rdat

      real, allocatable :: sumf(:,:,:,:)     ! sum of weights for each data field
      real, allocatable :: sumdir(:,:,:,:,:) ! sum of weights in directional bins
      real, allocatable :: azcomp(:,:,:,:)   ! azimuth angle components
      real, allocatable :: xg(:)             ! x coordinates (km) of grid points
      real, allocatable :: yg(:)             ! y coordinates (km) of grid points
      real, allocatable :: zg(:)             ! z coordinates (km) of grid points

c-----added variables (Mario)

      integer pass, npass
      integer maxray, maxgate
      real hsp0, vsp0, gamma, goodf
      real, allocatable :: baddata(:,:,:)
      real, allocatable :: f_old(:,:,:,:)
      integer, allocatable :: f_flag(:,:,:,:)      
      real, allocatable :: f_b(:,:,:,:)     ! backinterpolated data 
      real, allocatable :: sumf_b(:,:,:,:)  ! backinterpolated sum of wgt
      real, allocatable :: count_b(:,:,:,:) ! backinterpolated count

      integer, allocatable :: nray(:)
      integer, allocatable :: ngate(:)

c
c.... CLZ (7/27/18): added others variables previously
c

      real dataval
      real dsbad
      integer kfill
      parameter (dsbad = 10.0)

c
c.... CLZ (11/13/19): added arrays and parameters for median filter
c

      real, allocatable :: rayvel(:)
      real, allocatable :: work(:)
      real, allocatable :: med(:)
      real, allocatable :: window(:)

c      integer, allocatable :: rayvel(:)
c      integer, allocatable :: work(:)
c      integer, allocatable :: window(:)

      integer nrayout
      integer rayout
      integer maxgnwhalf
      integer nw
      integer nwhalf
      integer iw
      integer midwnd
      real diff_thres
      real absdiff
      real rhmax_near
      real hmax_near

c
c.... CLZ (12/10/19): filter parameters
c

      parameter (rayout = 0)     ! no diagnostic printout of ray data editing
c      parameter (rayout = 1)





c
c.... begin executable code
c






      write(6,*)
      write(6,*) 'entered multi_oban'

      write(6,*)'velfld=',velfld,' fname(1)=',fname(1)
      write(6,*)'reffld=',reffld,' fname(2)=',fname(2)

c
c.... open test output file containing radar ray data that was filtered
c

      if (rayout .eq. 1) then
      open(unit = 50, file = 'output_gates.dat',form='formatted',
     * status = 'unknown')
      endif



c
c############################################################################
c
c     Initialize variables.
c
c############################################################################ 
c


      goodf = 0.0

      allocate (f_flag(nx,ny,nz,nfld))
      f_flag(:,:,:,:) = 0

      allocate (f_old(nx,ny,nz,nfld))
      allocate (nray(nswp)) 
      allocate (ngate(nswp))
 
      allocate(sumf(nx,ny,nz,0:nfld))


      if (extrapolate_flag.eq.0) then

        allocate(sumdir(8,nx,ny,nz,nfld))
c        allocate(sumdir(4,nx,ny,nz,nfld))

      endif



      
      allocate(azcomp(nx,ny,nz,2))
      allocate(xg(nx))
      allocate(yg(ny))
      allocate(zg(nz))
      
c      write(6,*) 'error check'



c
c.... CLZ (5/18/07):  David Dowell's oban code modified for multipass.
c      Programmed by Mario Majcen, Paul Markowski, and Yvette Richardson (PSU)
c      during January-April 2007.



c-----add loop through passes (Mario)
c.... MM (1/1/07):  add multi-pass do-loop
c









c
c.... print out some control parameters for median filter
c

      write(6,*) 'diff_thres=',diff_thres
      write(6,*) 'nwhalf=',nwhalf




c
c------------------------------------------------------------------
c------------------------------------------------------------------
c
c.... do pass
c


      do pass = 1, npass






      azcomp(:,:,:,:) = 0.0

      if (pass .eq. 1) then
       f_old(:,:,:,:) = 0.0
      endif

      el(:,:,:) = 0.0
      time(:,:,:) = 0.0

      xrd(:,:,:) = 0.0
      yrd(:,:,:) = 0.0
      zrd(:,:,:) = 0.0

      count(:,:,:,:) = 0
      sumf(:,:,:,:) = 0.0

      if (extrapolate_flag .eq. 0) then
       sumdir(:,:,:,:,:) = 0.0
      endif







c
c.... CLZ (9/26/08):  "Base" h- and z-smoothing parameter values are computed below for
c         the current pass number.
c


      hsp = hsp0 * (gamma**(pass-1))
      write (6,*) 'on pass = ',pass,' hsp = ',hsp

      vsp = vsp0 * (gamma**(pass-1))
      write (6,*) 'on pass = ',pass,' vsp = ',vsp





      c1 = hsp * hsp
      c2 = vsp * vsp

      if (method.eq.1) then
        irad = 1+nint(hsp/dx)
        jrad = 1+nint(hsp/dy)
        krad = 1+nint(vsp/dz)
      else if (method.eq.2) then
        irad = 1+nint(sqrt(7.0*hsp)/dx)
        jrad = 1+nint(sqrt(7.0*hsp)/dy)
        krad = 1+nint(sqrt(7.0*vsp)/dz)
      else
        write(6,*) 'unknown method:  ', method
        stop
      endif ! if (method.eq.1) 









c
c.... CLZ (9/26/08):  Original location of following 4 lines was just below "c2 = ...".
c

      f(:,:,:,:) = 0.0

      call grid_coordinates(xg, nx, dx, xmin)
      call grid_coordinates(yg, ny, dy, ymin)
      call grid_coordinates(zg, nz, dz, zmin)



      glatr = dtor * glat
      glonr = dtor * glon

      if (iradtype .eq. 0) then

c
c.... fixed ground-based radar
c

      rlatr = dtor * rlat
      rlonr = dtor * rlon
      call ll_to_xy(xrad, yrad, map_proj, glatr, glonr, rlatr, rlonr)
      zrad = ralt - galt
      write(6,*) 'location of ground radar relative to grid origin:'
      write(6,*) xrad, yrad, zrad

c
c.... endif (iradtype .eq. 0) then
c

      endif









c
c############################################################################
c
c     Read in sweep files and compute weighted sums.
c
c############################################################################
c








c
c.... CLZ (5/18/07):  First scan the sweep file to obtain control parameters
c


       if (pass.eq.1) then

       do s=1, nswp

        do n=1, nfld

           write(6,*) 'Now on pass = ', pass
           write(6,*) '1st scan reading from ', sfname(s)


           call sweepread(sfname(s),vold,radd,celv,cfac,parm,swib,
     $                    ryib,asib,rdat,fname(n))


      write(6,*) 'finished 1st call to sweepread in multioban'

         call check_sweep_size(radd%num_param_desc, swib%num_rays)


        if (s .ge. minswpdc .and. s .le. maxswpdc) then
        
         call correct_date(vold%year, yrcor, vold%mon, mocor, vold%day, dacor)
         write(6,*)'1st sweepread on swp=',s,'corrected day=',vold%day
         
        endif


         call correct_az_el(ryib, swib%num_rays, azcor, elcor, az_corr_flag)

           if (umass_flag.eq.1) then
             call correct_umass_data(swib%num_rays, ryib, celv%total_gates, rdat)
           endif

           ngate(s)=celv%total_gates      
           nray(s)=swib%num_rays

         enddo
       enddo


c
c           find max number of rays and max number of gates in a ray (Mario)
c

c
c.... CLZ (2/17/10):  P. Markowski reported bug, should have maxgate = ngate(s) in if-test
c

            maxray = nray(1)
            maxgate = ngate(1)


            do s = 1, nswp


              if (nray(s).gt.maxray) then
                 maxray=  nray(s)
              endif

c
c.... must update maxgate as well as maxray (previous code mistakenly duplicated)
c

c              if (nray(s).gt.maxray) then
c                 maxray=  nray(s)
c              endif

              if (ngate(s) .gt. maxgate) then
                 maxgate =  ngate(s)
              endif

c
c.... enddo s=1, nswp
c

            enddo

c
c.... compute and store median filter parameters
c

      if (nwhalf .ge. 1) then

        write(6,*)'in multi_oban, nwhalf=',nwhalf
        write(6,*)'nwhalf > 0, so median filter velocity'

        nw = (2 * nwhalf) + 1
        write(6,*)'in multi_oban, nw=',nw

        midwnd = nwhalf + 1
        write(6,*)'in multi_oban, midwnd=',midwnd

        maxgnwhalf = maxgate - nwhalf

        write(6,*)'in multi_oban, maxgate=',maxgate
        write(6,*)'in multi_oban, maxgnwhalf=',maxgnwhalf

      elseif (nwhalf .eq. 0) then
        write(6,*)'in multi_oban, nwhalf=',nwhalf
        write(6,*)'nwhalf = 0, so do not median filter velocity'
      endif

      if (rayout .eq. 1) then
        write(50,*) maxgate
      endif

      if (nwhalf .ge. 1) then
        allocate (rayvel(maxgate))
        allocate (work(maxgate))
        allocate (med(maxgate))
        allocate (window(nw))
      endif


        allocate(baddata(nswp,maxray,maxgate))
        baddata(:,:,:)= 0

c
c.... endif pass.eq.1
c

      endif

     
c
c.... CLZ (5/18/07):  Scan through sweep to read and interpolate data
c    

      do s = 1, nswp

        do n=1, nfld
           write(6,*) 'Now on pass = ', pass
           write(6,*) '2nd scan reading from ', sfname(s)


           call sweepread(sfname(s),vold,radd,celv,cfac,parm,swib,
     $                    ryib,asib,rdat,fname(n))

           call check_sweep_size(radd%num_param_desc, swib%num_rays)



        if (s .ge. minswpdc .and. s .le. maxswpdc) then
        
         call correct_date(vold%year, yrcor, vold%mon, mocor, vold%day, dacor)
         write(6,*)'2nd sweepread on swp=',s,'corrected day=',vold%day

        endif




c           call correct_date(vold%year, yrcor, vold%mon, mocor, vold%day, dacor)
           call correct_az_el(ryib, swib%num_rays, azcor, elcor, az_corr_flag)

           if (umass_flag.eq.1) then
             call correct_umass_data(swib%num_rays, ryib, celv%total_gates, rdat)
           endif

            
           do r=1, swib%num_rays

            ti = timediff(vold%year, vold%mon, vold%day,
     $                    ryib(r)%hour, ryib(r)%min, ryib(r)%sec, ryib(r)%msec,
     $                    cyr, cmo, cda, chr, cmn, cse, cms)




c************************************************

      if (iradtype .eq. 0) then

c                call xyzloc(x, y, z, range, dtor*ryib(r)%azimuth, dtor*ryib(r)%elevation,
c     $                      map_proj, glatr, glonr, galt, rlatr, rlonr, ralt,
c     $                      ut, vt, ti)

c
c.... CLZ (9/3/19): azimuth angle for ground-based radar relative to north
c

      azimang = ryib(r)%azimuth

c
c.... CLZ CLZ (9/3/19): elevation angle for ground-based radar from horizontal (+ or - 90 deg)
c

      elevang = ryib(r)%elevation

c
c      if (s .eq. 1 .and. r .eq. 30) then
c      write(6,*)'azimang=',azimang,' elevang=',elevang
c      endif

c
c.... endif (iradtype .eq. 0) then
c
      endif


      if (iradtype .eq. 1) then

c
c.... CLZ (9/20/18): extract ray "s" TA radar rlat (deg), rlon (deg) , and ralt (m MSL)
c

c
c.... P-3 radar latitude
c

      rlat = asib(r)%lat
c      rlat = radd%radar_lat

c
c.... P-3 radar longitude
c

      rlon = asib(r)%lon
c      rlon = radd%radar_lon

c
c.... P-3 radar altitude (MSL)
c

      ralt = asib(r)%alt_msl
c     Alt_AGL = asib(i)%alt_agl
c      ralt = radd%radar_alt

c
c.... moving airborne radar at (rlat, rlon, ralt) where ralt is altitude (km MSL)
c

      rlatr = dtor * rlat
      rlonr = dtor * rlon
      call ll_to_xy(xrad, yrad, map_proj, glatr, glonr, rlatr, rlonr)
      zrad = ralt - galt
c      write(6,*) 'ray=',r,'xr=',xrad,'yr=',yrad,'zr=',zrad

c
c.... print angle parameters for this ray
c


      rot_raw = asib(r)%rot_ang      ! Raw rotation angle [deg]
      tilt = asib(r)%tilt_ang      ! Tilt angle [deg], raw (fore/aft direction)
      roll = asib(r)%roll       ! Roll angle [deg]
      rot_rc = Amod((rot_raw + roll + 360.0), 360.0)      ! Correct the rotation angle for roll
      pitc = asib(r)%pitch       ! Pitch angle [deg] (uncorrected for radar mounting error)
      drft = asib(r)%drift       ! Drift angle [deg]
      head = asib(r)%head       ! Heading [deg]
      trck = Amod((head + drft + 360.0), 360.0)                   !  Track angle [deg]

c
c.... CLZ (10/30/18): call subroutine trans_coord to obtain meteorological azimuth and elevation
c

      call trans_coords(rot_rc, tilt, drft, pitc, trck, azim_north,
     * elev_hor, elev_zen, tilt_tr)

c
c.... DPJ (personal communication, 2018): Trans_coords transforms heading from trigonometric 
c       units (+/- 180 deg) to meteorological units (0-360 deg).
c

c      head = amod(head+360.0, 360.0)


c
c.... CLZ (10/30/18): azimuth angle relative to north
c

      azimang = azim_north
      
c      azimang = head + 90.0 - asib(r)%tilt_ang
c      if (azimang .ge. 360.0) then
c      azimang = azimang - 360.0
c      elseif (azimang .lt. 0.0) then
c      azimang = azimang + 360.0
c      endif


c
c.... CLZ (10/30/18): elevation angle from horizontal (+ or - 90 deg)
c

      elevang = elev_hor

c      elevang = 90.0 - asib(r)%rot_ang

c      write(6,*) 'ray=',r,'az=',azimang,'el=',elevang
c      write(6,*) 'ray=',r,'rotangle=',rotangle,'tilt=',tilt
c      write(6,*) 'ray=',r,'head=',head,'roll=',roll
c      write(6,*) 'ray=',r,'pitc=',pitc,'drft=',drft

c
c.... endif (iradtype .eq. 1) then
c

      endif

c************************************************

c
c.... CLZ (12/11/19): perform median filtering of velocity prior to objectively
c                       analyzing the current ray data
c

          if(nwhalf .ge. 1 .and. n .eq. 1) then
c          if(n .le. 2) then

            do g = 1, maxgate

c
c.... compute the slant range (km) to the gate
c

            range = rangekm_to_gate(celv, g)

c
c.... CLZ (11/23/10): now compute the horizontal range (km) from the slant range
c

            hrange = cos(dtor * elevang) * range

c
c.... CLZ CLZ (9/3/19): corrected error to compute azimang and elevang for ground-based radar via
c                       azimuth and elevation angle arguments in commented xyzloc call below
c

            call xyzloc(x, y, z, range, dtor*azimang, dtor*elevang,
     *        map_proj, glatr, glonr, galt, rlatr, rlonr, ralt,
     *           ut, vt, ti)

c
c.... CLZ(12/10/19): eliminate close-range, low-altitude data points
c

            if (hrange .le. rhmax_near .and.
     *           (z - zmin) .le. hmax_near) then
            rdat(r)%data(g) = sbad
            endif

            work(g) = rdat(r)%data(g)

            if (fname(n) .eq. velfld) then
            rayvel(g) = work(g)
            endif

            enddo       ! g=1, celv%total_gates

c
c.... median filter pass on interior
c

            do g = nwhalf + 1, maxgnwhalf - 1
c            do g = 3, maxg5

c
c.... Generate and apply the median filter window
c

            iw = 0
            
            do gg = g - nwhalf, g + nwhalf
            
            iw = iw + 1
            window(iw) = work(gg)

c
c.... enddo gg = g - nwhalf, g + nwhalf
c

            enddo

c
c.... bubble-sort window elements (including missing values) to find median
c

            call bubble(window, nw)

c
c.... CLZ(12/10/19): identify and save the median value
c

            med(g) = window(midwnd)

c
c.... CLZ(12/11/19):  apply the median value if datum is significantly different from median
c

            if (abs(window(midwnd) - work(g)) .ge. diff_thres) then

            work(g) = med(g)

c
c.... endif (abs(window(midwnd) - work(g)) .ge. diff_thres) then
c

            endif

c
c.... enddo g = nwhalf + 1, maxgnwhalf - 1
c

            enddo

c
c.... CLZ (12/11/19): simple sequential difference threshold correction at beginning points
c

            do g = nwhalf, 1, -1

            absdiff = abs(work(g) - work(g+1))
            if (absdiff .ge. diff_thres) then
            work(g) = work(g+1)
            endif

            enddo

c
c.... CLZ (12/11/19): simple sequential difference threshold correction at ending points
c

            do g = maxgnwhalf, maxgate

            absdiff = abs(work(g) - work(g-1))
            if (absdiff .ge. diff_thres) then
            work(g) = work(g-1)
            endif

            enddo

c
c.... CLZ (12/11/19): copy the filtered ray field work(g) to rdat(r)%data(g) for interpolation
c

            do g = 1, maxgate
            rdat(r)%data(g) = work(g)
            enddo

c
c.... output the ray velocity data if rayout = 1 and fname(n) = 'VT'
c

      if (rayout .eq. 1 .and. (fname(n) .eq. velfld) .and.
     * rot_rc .ge. 250.0 .and. rot_rc .le. 300.0) then

      write(50,*)
      write(50,*)
      write(50,*) fname(n)
      write(50,*)'s=',s,' r=',r,' rot_rc=',rot_rc
      write(50,*)
      
      do g = 1, maxgate
      write(50,*) g,' dat=',rayvel(g),' wrk=',work(g),' med=',med(g)
      enddo

c
c.... endif (rayout .eq. 1 .and. (fname(n) .eq. velfld) .and.
c    * rot_rc .ge. 250.0 .and. rot_rc .le. 300.0) then
c

      endif


c
c.... endif(nwhalf .ge. 1 .and. n .eq. 1) then
c

          endif


c
c.... now weight the data
c

            do g=1, celv%total_gates

c           if (rdat(r)%data(g) .lt. (sbad + dsbad) .and.
c     *       rdat(r)%data(g) .gt. (sbad - dsbad) ) then

              if (rdat(r)%data(g) .eq. sbad ) then
                 baddata(s,r,g) = 1
              endif


              if (rdat(r)%data(g) .ne. sbad) then

c      if (s .eq. 1 .and. r .eq. 30 .and. (n .eq. 1 .or. n .eq. 2)) then
c
c      if (n .eq. 1) then
c      write(6,*)'g=',g,' field1=',rdat(r)%data(g)
c      elseif (n .eq. 2) then
c      write(6,*)'g=',g,' field2=',rdat(r)%data(g)
c      endif
c
c      endif


c
c.... compute the slant range (km) to the gate
c

                range = rangekm_to_gate(celv, g)

c
c.... CLZ (11/23/10): now compute the horizontal range (km) from the slant range
c

                hrange = cos(dtor * elevang) * range
c                hrange = cos(dtor * ryib(r)%elevation) * range

c
c.... CLZ CLZ (9/3/19): corrected error to compute azimang and elevang for ground-based radar via
c                       azimuth and elevation angle arguments in commented xyzloc call below
c

                call xyzloc(x, y, z, range, dtor*azimang, dtor*elevang,
     $                      map_proj, glatr, glonr, galt, rlatr, rlonr, ralt,
     $                      ut, vt, ti)

c                call xyzloc(x, y, z, range, dtor*ryib(r)%azimuth, dtor*ryib(r)%elevation,
c     $                      map_proj, glatr, glonr, galt, rlatr, rlonr, ralt,
c     $                      ut, vt, ti)



                ii = 1 + nint((x-xmin)/dx)
                jj = 1 + nint((y-ymin)/dy)
                kk = 1 + nint((z-zmin)/dz)

c
c.... CLZ (9/26/08):  current gate is at grid level = kk at this point in the code
c


                imin = max(ii-irad, 1)
                imax = min(ii+irad, nx)
                jmin = max(jj-jrad, 1)
                jmax = min(jj+jrad, ny)
                kmin = max(kk-krad, 1)
                kmax = min(kk+krad, nz)

c
c.... CLZ (9/18/19): vclut is hypothesized magnitude of a clutter gate radial velocity
c

                if ( iradtype .eq. 1) then

                vclut = 2.5

c                dbzclut = 30.0
c                dangdr = 0.1  !degrees per kilometer
c                angclut0 = 2.0
c                angclut = angclut0 + (dangdr * range)

c
c.... CLZ (9/23/19): compute the threshold half-power beamwidth length (km)
c

                hpbwlen = 0.0
c				hpbwlen = range * sin(angclut * dtr)

c
c.... CLZ (9/24/19): no clutter check for ground-based radar
c

                elseif ( iradtype .eq. 0) then
                hpbwlen = 0.0

c
c.... endif ( iradtype .eq. 1) then
c

                endif

c
c.... CLZ (9/23/19): impose prerequisite checks on gate to objectively weigh:
c                      (1) gate in allowable range interval;
c                      (2) lower half of main lobe (half-power beamwidth) about ground;
c                      (3) gate value within sanity limits;
c                      (4) gate velocity outside ground target velocity interval.
c

                if (range .ge. minrange .and. hrange .le. rhmax .and.
     *              z .gt. hpbwlen .and.
     *               rdat(r)%data(g) .gt. -fsan(n) .and.
     *                rdat(r)%data(g) .lt. fsan(n)) then
c                if (range .ge. minrange .and. hrange .le. rhmax) then

c
c.... CLZ (12/9/19): if velocity field velfld, do not weight abs(V) < vclut
c

                  if ( iradtype .eq. 0 .or.
     *             (iradtype .eq. 1 .and. fname(n) .ne. velfld) .or.
     *              (iradtype .eq. 1 .and.
     *               fname(n) .eq. velfld .and.
     *               (rdat(r)%data(g) .lt. -vclut .or.
     *                rdat(r)%data(g) .gt. vclut) ) ) then


                  do k=kmin, kmax
                    do j=jmin, jmax
                      do i=imin, imax

                        if (method.eq.1) then
                          dh2 = (x-xg(i))*(x-xg(i)) + (y-yg(j))*(y-yg(j))
                          dv2 = (z-zg(k))*(z-zg(k))
                          wgt = (1.0 - dh2/c1 - dv2/c2) / (1.0 + dh2/c1 + dv2/c2)
                        else if (method.eq.2) then
                          wgt = exp( -(x-xg(i))*(x-xg(i))/hsp
     $                               -(y-yg(j))*(y-yg(j))/hsp
     $                               -(z-zg(k))*(z-zg(k))/vsp )

                        endif

c                        write (6,*) 'wgt',wgt


c*********************

                        if (wgt .gt. 0.0) then


                          if (pass .gt. 1) then

                          if (baddata(s,r,g) .eq. 0) then
                           if (f_flag(i,j,k,n) .eq. 0) then

                            f(i,j,k,n) = f(i,j,k,n)
     $                        + wgt*((signf(n)*rdat(r)%data(g))-f_b(s,r,g,n))
                            sumf(i,j,k,n) = sumf(i,j,k,n) + wgt
                          
                            count(i,j,k,n) = count(i,j,k,n) + 1

c                            if (sumf(i,j,k,n).lt.1.0E-20) then
c                              write(6,*) 'bad sumf somehow got in at',i,j,k,n
c                              write(6,*) 'sumf',sumf(i,j,k,n)
c                            endif
c                            write(6,*) 'wgt',wgt
c                            write(6,*) 'rdat(r)%data(g)',rdat(r)%data(g)
c                            write(6,*) 'f_b(s,r,g,n)',f_b(s,r,g,n)
c                            write(6,*) 'f(i,j,k,n)',f(i,j,k,n)

                          endif ! if (f_flag(i,j,k,n).eq.0)
                          endif ! baddata(s,r,g).eq.0)

                          endif !(pass.gt.1)




                          if (pass.eq.1) then

                           f(i,j,k,n) = f(i,j,k,n)
     $                               + wgt*(signf(n)*rdat(r)%data(g))

                          sumf(i,j,k,n) = sumf(i,j,k,n) + wgt

                          count(i,j,k,n) = count(i,j,k,n) + 1
                          endif ! (pass.eq.1)




                          if (extrapolate_flag.eq.0) then

                            call update_dir_bins(sumdir(1,i,j,k,n), wgt,
     $                                           x, y, z, xg(i), yg(j), zg(k))

                          endif
                          



                          if (iradtype .eq. 0) then

c........................ case of ground-based radar (compute weighted sum of az, el, and time)

                          azcomp(i,j,k,1) = azcomp(i,j,k,1)
     $                          + wgt*sin(dtor*ryib(r)%azimuth)
                          azcomp(i,j,k,2) = azcomp(i,j,k,2)
     $                          + wgt*cos(dtor*ryib(r)%azimuth)
                          el(i,j,k) = el(i,j,k)
     $                          + wgt*ryib(r)%elevation
                          time(i,j,k) = time(i,j,k) + wgt*ti
                          sumf(i,j,k,0) = sumf(i,j,k,0) + wgt

                          elseif (iradtype .eq. 1) then

c........................ case of airborne radar (compute weighted sum of xrd, yrd, and zrd)

                          xrd(i,j,k) = xrd(i,j,k) + (wgt * xrad)
                          yrd(i,j,k) = yrd(i,j,k) + (wgt * yrad)
                          zrd(i,j,k) = zrd(i,j,k) + (wgt * zrad)
                          sumf(i,j,k,0) = sumf(i,j,k,0) + wgt

c
c.... endif (iradtype .eq. 0) then
c

                          endif



                        endif   ! if (wgt.gt.0.0)

c*********************

c
c.... enddo i=imin, imax
c

                      enddo

c
c.... enddo j=jmin, jmax
c

                    enddo

c
c.... enddo do k=kmin, kmax
c

                  enddo

c
c.... endif ( iradtype .eq. 0 .or.
c      (iradtype .eq. 1 .and. fname(n) .ne. 'VT') .or. 
c       (iradtype .eq. 1 .and. fname(n) .eq. 'VT' .and.
c        (rdat(r)%data(g) .lt. -vclut .or.
c         rdat(r)%data(g) .gt. vclut) ) ) then
c

                     endif

c
c.... endif (range .ge. minrange .and. hrange .le. rhmax .and.
c     *              z .gt. hpbwlen .and.
c     *               rdat(r)%data(g) .gt. -fsan(n) .and.
c     *                rdat(r)%data(g) .lt. fsan(n)) then
c

                endif       ! (range.ge.minrange)

              endif       ! (rdat(r)%data(g) .ne. sbad)
                     
            enddo       ! g=1, celv%total_gates
             
          enddo       ! r=1, swib%num_rays
            
        enddo       ! n=1, nfld
        
      enddo       ! s=1, nswp



c      write(6,*) 'Finished summing wgts and wgts*data on pass = ', pass

c      write(6,*) 'error check'

c
c############################################################################
c
c     Compute interpolated values from weighted sums.
c
c############################################################################
c

c      print*,'i  got here'

c      write (6,*) 'minsum', minsum

      do k=1, nz
        do j=1, ny
          do i=1, nx


            do n=1, nfld

             valid_sum = 1



             if (pass .eq. 1) then
 
                if (sumf(i,j,k,n) .le. minsum .or.
     *           count(i,j,k,n) .lt. mincount) then
c                if (sumf(i,j,k,n).le.minsum) then

                  valid_sum = 0
                  f_flag(i,j,k,n) = 1
                endif
 
                if (extrapolate_flag .eq. 0) then

c                  do d=1, 8
                  do d=1,4
                    if (sumdir(d,i,j,k,n).le.(0.1*minsum) ) then
                      valid_sum = 0
                      f_flag(i,j,k,n) = 1
                    endif
                  enddo
                endif

             endif ! (pass.eq.1) then



              if (pass .gt. 1) then

                 if (f_flag(i,j,k,n) .eq. 1) then
                   valid_sum = 0
                 endif

                if (sumf(i,j,k,n) .le. minsum .or.
     *           count(i,j,k,n) .lt. mincount) then
c                 if (sumf(i,j,k,n).le.minsum) then

                  valid_sum = 0
                  f_flag(i,j,k,n) = 1
                 endif

                 if (extrapolate_flag .eq. 0) then

c                  do d=1, 8
                  do d=1,4
                    if (sumdir(d,i,j,k,n).le.(0.1*minsum) ) then
                      valid_sum = 0
                      f_flag(i,j,k,n) = 1
                    endif
                  enddo
                 endif

c
c.... endif (pass.gt.1)
c

              endif

             

             if (valid_sum .eq. 1) then

                goodf = goodf + 1.0

                f(i,j,k,n) = f(i,j,k,n) / sumf(i,j,k,n) + f_old(i,j,k,n)
c                write (6,*) 'f(',i,',',j,',',k,',',n,')=',f(i,j,k,n)

             else 

                f(i,j,k,n) = bad

             endif

            enddo ! enddo n=1, nfld




            if (sumf(i,j,k,0) .gt. minsum) then


             if (iradtype .eq. 0) then

c............ case of ground-based radar (compute weighted sum of az, el, and time)

              azcomp(i,j,k,1) = azcomp(i,j,k,1) / sumf(i,j,k,0)
              azcomp(i,j,k,2) = azcomp(i,j,k,2) / sumf(i,j,k,0)
              az(i,j,k) = comptoaz(azcomp(i,j,k,1),azcomp(i,j,k,2))
              el(i,j,k) = el(i,j,k) / sumf(i,j,k,0) 
              time(i,j,k) = time(i,j,k) / sumf(i,j,k,0)

             elseif (iradtype .eq. 1) then

c............ case of airborne radar (compute weighted sum of xrd, yrd, and zrd)

              xrd(i,j,k) = xrd(i,j,k) / sumf(i,j,k,0)
              yrd(i,j,k) = yrd(i,j,k) / sumf(i,j,k,0)
              zrd(i,j,k) = zrd(i,j,k) / sumf(i,j,k,0)

             endif

            else

             if (iradtype .eq. 0) then

c............ case of ground-based radar (compute weighted sum of az, el, and time)

              az(i,j,k) = bad
              el(i,j,k) = bad
              time(i,j,k) = bad

             elseif (iradtype .eq. 1) then

c............ case of airborne radar (compute weighted sum of xrd, yrd, and zrd)

              xrd(i,j,k) = bad
              yrd(i,j,k) = bad
              zrd(i,j,k) = bad

             endif

            endif

          enddo
        enddo
      enddo
      


c        write(6,*) 'goodf/2 =',goodf/2.










c
c.... CLZ (5/18/07):  Now perform interpolation from nth-pass gridded fields
c       back to the gate locations.
c

c----  back-interpolation (mario)

     
c
c############################################################################
c
c     Read in sweep files and compute weighted sums.
c
c############################################################################
c
           
c             add allocation for new variables (Mario)

            if (pass.eq.1) then
             allocate(f_b(nswp,maxray,maxgate,nfld))
             allocate(sumf_b(nswp,maxray,maxgate,nfld))
             allocate(count_b(nswp,maxray,maxgate,nfld))
            endif ! (pass.eq.1)          

            
            f_b(:,:,:,:) = 0.0
            sumf_b(:,:,:,:) = 0.0
            count_b(:,:,:,:) = 0.0

c
c.... CLZ (5/18/07):  Back-interpolate on nth pass to prepare for pass n+1
c      So.......perform the back-interpolation only if pass > 1
c      (i.e., only if npass > 1 and pass < npass)
c
c.... E.g., Back-interpolate if npass = 2 and pass = 1, but do NOT
c         back-interpolate if npass = 1 and pass = 1.
c
c

        if (pass .ne. npass) then
      
        do s=1, nswp
         do n=1, nfld
c           write(6,*) 'back-interpolation on pass = ', pass,
c     *      ' of ', npass,' total passes'
c           write(6,*) 'this is pass=', pass
c           write(6,*) 'reading from ', sfname(s)
           call sweepread(sfname(s),vold,radd,celv,cfac,parm,swib,
     $                    ryib,asib,rdat,fname(n))

           call check_sweep_size(radd%num_param_desc, swib%num_rays)





        if (s .ge. minswpdc .and. s .le. maxswpdc) then
        
         call correct_date(vold%year, yrcor, vold%mon, mocor, vold%day, dacor)
c         write(6,*)'3rd sweepread on swp=',s,'corrected day=',vold%day

        endif






c           call correct_date(vold%year, yrcor, vold%mon, mocor, vold%day, dacor)
           call correct_az_el(ryib, swib%num_rays, azcor, elcor, az_corr_flag)

           if (umass_flag.eq.1) then
             call correct_umass_data(swib%num_rays, ryib, celv%total_gates, rdat)
           endif
          





          do r=1, nray(s)

c
c.... CLZ (9/22/09):  Add call to timediff here!!!  (It was missing in original code.)
c

            ti = timediff(vold%year, vold%mon, vold%day,
     $                    ryib(r)%hour, ryib(r)%min, ryib(r)%sec, ryib(r)%msec,
     $                    cyr, cmo, cda, chr, cmn, cse, cms)


             do g=1, ngate(s)



c           if (rdat(r)%data(g) .gt. (sbad + dsbad) .or.
c     *       rdat(r)%data(g) .lt. (sbad - dsbad) ) then

               if (rdat(r)%data(g) .ne. sbad) then



c
c.... compute the slant range (km) to the gate
c


                 range = rangekm_to_gate(celv, g)

c
c.... CLZ (11/23/10): now compute the horizontal range (km) from the slant range
c

                hrange = cos(dtor * ryib(r)%elevation) * range


                 call xyzloc(x, y, z, range, dtor*ryib(r)%azimuth, dtor*ryib(r)%elevation,
     $                      map_proj, glatr, glonr, galt, rlatr, rlonr, ralt,ut, vt, ti)



                if (range .ge. minrange .and. hrange .le. rhmax .and.
     *               rdat(r)%data(g) .gt. -fsan(n) .and.
     *                rdat(r)%data(g) .lt. fsan(n)) then
c                if (range .ge. minrange .and. hrange .le. rhmax) then

c                if (range.ge.minrange) then



                ii = 1 + nint((x-xmin)/dx)
                jj = 1 + nint((y-ymin)/dy)
                kk = 1 + nint((z-zmin)/dz)


c
c.... CLZ (9/26/08):  current gate is at grid level = kk at this point in the code
c



                imin = max(ii-irad, 1)
                imax = min(ii+irad, nx)
                jmin = max(jj-jrad, 1)
                jmax = min(jj+jrad, ny)
                kmin = max(kk-krad, 1)
                kmax = min(kk+krad, nz)




                  
                
                   
                  do k=kmin, kmax
                   do j=jmin, jmax
                    do i=imin, imax
                       if (f_flag(i,j,k,n).eq.0) then
                        if (method.eq.1) then
                          dh2 = (x-xg(i))*(x-xg(i)) + (y-yg(j))*(y-yg(j))
                          dv2 = (z-zg(k))*(z-zg(k))
                          wgt = (1.0 - dh2/c1 - dv2/c2) / (1.0 + dh2/c1 + dv2/c2)
                        else if (method.eq.2) then
                          wgt = exp( -(x-xg(i))*(x-xg(i))/hsp
     $                               -(y-yg(j))*(y-yg(j))/hsp
     $                               -(z-zg(k))*(z-zg(k))/vsp )
                        endif
                        if (wgt.gt.0.0) then
                          
                            f_b(s,r,g,n) = f_b(s,r,g,n) + wgt*f(i,j,k,n)
                            sumf_b(s,r,g,n) = sumf_b(s,r,g,n) + wgt
                            count_b(s,r,g,n) = count_b(s,r,g,n) + 1
c                           write (6,*) 'sumf_b(s,r,g,n)',sumf_b(s,r,g,n)
c                           write (6,*) 'sumf_b(s,r,g,n)',sumf_b(s,r,g,n)                                           
                        endif   ! if (wgt.gt.0.0)
                       endif ! (f_flag(i,j,k,n).eq.0)
                      
                 
                 enddo        ! k=kmin, kmax
                 enddo        ! j=jmin, jmax
                 enddo        ! i=imin, imax
               endif       ! (range.ge.minrange)
               endif       ! (rdat(r)%data(g) .ne. sbad)
             enddo       ! g=1, ngate(s)
           enddo       ! r=1, nray(s)
         enddo       ! n=1, nfld
       enddo       ! s=1, nswp
c      write(6,*) 'error check2'      
      

c
c.... CLZ (5/17/07):  Compute interpolated fields for final pass
c

c
c############################################################################
c
c     Compute interpolated values from weighted sums.
c
c############################################################################
c


        do s=1, nswp
c         write (6,*) 'error check3'
         do n=1, nfld
          do r=1, nray(s)
             do g=1, ngate(s)

         


c              range = rangekm_to_gate(celv, g)
c              if (range.ge.minrange) then

                if (baddata(s,r,g).eq.0) then
                  f_b(s,r,g,n) = f_b(s,r,g,n) / sumf_b(s,r,g,n) 

c                 else 
c                   f_b(s,r,g,n) = bad
                endif

c              endif ! if (range.ge.minrange)

            enddo ! g=
          enddo   ! r=
        enddo     ! n=
      enddo       ! s=
      


      

c----- variables f passed to f_old (Mario)

      do i=1,nx
       do j=1,ny
        do k=1,nz
         do n=1,nfld
               
          f_old(i,j,k,n) = f(i,j,k,n)
         enddo
        enddo
       enddo
      enddo

c
c.... endif (pass.ne.npass)
c      
      
      endif

c
c.... Have now completed the back-interpolation loop
c





c     minsum normalization
c      minsum=minsum*gamma

      minsum = minsum * exp(-1./gamma)







c      write(6,*) 'copying f to f_mpass after completing pass = ', pass


      do n=1, nfld
        do k=1, nz
          do j=1, ny
            do i=1, nx
                f_mpass(i,j,k,n,pass) = f(i,j,k,n)
            enddo
          enddo
        enddo
      enddo









c
c..... enddo (pass=1,npass)
c

      enddo




c
c.... CLZ (12/18/13):  constant profile from sfc to first data level if k < kfill
c

      if (nz. gt. 3) then
      kfill = 3
      elseif (nz .le. 3) then
      kfill = nz
      endif

      if (extrapolate_flag .eq. 1) then

      do n = 1, nfld
      do k = 1, kfill
      do j = 1, ny
      do i = 1, nx

      if (f(i, j, k, n) .eq. 9.9000003E+09 .and.
     * f(i, j, kfill + 1, n) .ne. 9.9000003E+09) then
      f(i, j, k, n) = f(i, j, kfill + 1, n)
      endif

      enddo
      enddo
      enddo
      enddo

      endif




c
c------------------------------------------------------------------
c------------------------------------------------------------------
c




c
c.... close the edited ray output file if rayout = 1
c

      if (rayout .eq. 1) then
      rewind(50)
      close(50)
      endif





c
c.... deallocate
c

      deallocate(f_old)
      deallocate(f_b)
      deallocate(sumf_b)
      deallocate(count_b)
      deallocate(nray)
      deallocate(ngate)
      
      deallocate(sumf)
      if (extrapolate_flag.eq.0) then
        deallocate(sumdir)
      endif
      deallocate(azcomp)
      deallocate(xg)
      deallocate(yg)
      deallocate(zg)
      deallocate(f_flag)
      deallocate(baddata)

      write(6,*)
      write(6,*) 'Have completed multi_oban, now returning to driver'
      write(6,*)

      return
      end





      subroutine bubble(window, nw)

c      call bubble(window, nw)


c
c.... CLZ (11/13/19): bubble sort values in a 1-D data window from largest to smallest
c

      real window(nw)
      real savewj
      integer nw
      integer nwm1
      integer ibub

      nwm1 = nw - 1

c
c.... seed the do while loop by setting ibub = 1
c

      ibub = 1

c
c.... resort values inside window from largest to smallest
c

      do while (ibub .eq. 1)

c
c.... will remain inside do while if & only if one or more pairs are resorted in j loop
c

      ibub = 0

c
c.... perform a pass scanning through window, resorting elements one at a time as needed
c

      do j = 1, nwm1

      if (window(j+1) .gt. window(j)) then
      savewj = window(j)
      window(j) = window(j+1)
      window(j+1) = savewj
      ibub = 1
      endif

c
c.... enddo j = 1, nwm1
c

      enddo

c
c.... enddo while (ibub .eq. 1)
c

      enddo

c
c.... done with bubble
c

      return
      end






      subroutine trans_coords(rotation, tilt_raw, drft, ptch, trck,
     * azimuth, elev_hor, elev_zen, tilt_tr)
  
c  Routine to translate the aircraft relative tail-radar measured angles
c  to track relative angles

c  Input parameters:
c       Rotation:  Rotation angle of the antenna from zenith
c                  (used to be called "azimuth") corrected for roll
c       Tilt_Raw:  Fore/Aft angle measured relative to a normal plane
c                  to the airframe
c       Drft:      Drift angle
c       Ptch:      Pitch angle
c       Trck:      Track angle

c  Output Parameters:
c       Azimuth:   Compass direction of the beam measured from North
c       Elev_Zen:  Elevation angle of the beam measured from zenith
c       Tilt_TR:   Track relative tilt angle (Fore/Aft)

      real dtr
      parameter (dtr = 0.0174532925)

c  Convert angles to radians

      rot = rotation * dtr

      If (tilt_raw .gt. 100.0) then
      tilt_raw = tilt_raw - 180.0     ! for RVP-8 tilt
      endif

      tlt = tilt_raw * dtr
      dri = drft * dtr
      pit = ptch * dtr
      trk = trck * dtr

c  Direction cosine for x (distance normal to the track)

      x = (cos(rot) * sin(dri) * sin(pit) * cos(tlt)) +
     * (cos(dri)*sin(rot)*cos(tlt)) -
     * (sin(dri)*cos(pit)*sin(tlt))

c  Direction cosine for y (distance along track due to fore/aft pointing)

      y = (-cos(rot)*cos(dri)*sin(pit)*cos(tlt)) +
     * (sin(dri)*sin(rot)*cos(tlt)) +
     * (cos(dri)*cos(pit)*sin(tlt))

c  Direction cosine for z (height)

      z = (cos(pit)*cos(rot)*cos(tlt)) + (sin(pit)*sin(tlt))

c  Track relative tilt angle (fore/aft looking angle)

      tilt_tr = asin(y)/dtr

c  Azimuth (x,y) angle (compass direction) relative to the track

      azm_tr = atan2(x,y)/dtr

c  Azimuth angle relative to North

      azimuth = amod((azm_tr + trck), 360.0)

c  Elevation angle from the horizontal (+ or - 90 deg)

      elev_hor = asin(z)/dtr

c  Elevation angle from zenith (0 to 180 deg)

      elev_zen = 90.0 - elev_hor

c
c.... done with trans_coords
c

      return
      end







      subroutine output_gates(method, gamma,npass, hsp0, vsp0,
     $                      f, f_mpass, nfld, fname, nswp, sfname,
     $                      extrapolate_flag,minsum,minrange,rhmax,
     $                      yrcor, mocor, dacor,
     $                      azcor, elcor, az_corr_flag, umass_flag,
     $                      az, el, time, count,
     $                      map_proj, glat, glon, galt, rlat, rlon, ralt,
     $                      nx, ny, nz, dx, dy, dz, xmin, ymin, zmin,
     $                      cyr, cmo, cda, chr, cmn, cse, ut, vt)







c
c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######              SUBROUTINE output_gates                 ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     This subroutine reads in radar data from sweep files and
c     interpolates the desired fields to a Cartesian grid.  Multiple options
c     are available for the 3D weighting function used by the interpolation scheme:
c
c     1.  Cressman weighting function:
c
c         wgt = ( 1.0 - deltax**2/hsp**2 - deltay**2/hsp**2 - deltaz**2/vsp**2 )
c             / ( 1.0 + deltax**2/hsp**2 + deltay**2/hsp**2 + deltaz**2/vsp**2 )
c
c
c     2.  Barnes weighting function:
c
c         wgt = exp ( -deltax**2/hsp - deltay**2/hsp - deltaz**2/vsp )
c
c############################################################################
c
c     Author:  Conrad Ziegler
c
c     Created:  February 2011 (derived from multi_oban.f)
c
c     CLZ (9/22/09):  fixed time difference function call bug for npass > 1
c
c      CLZ (2/17/10):  Markowski's bug fix, maxgate = ngate(1) --> maxgate = ngate(s)
c
c     CLZ (3/9/10): restored option for no extrapolation in all 3 coordinate directions
c                     (call update_dir_bins, located in util.f)
c                     read in/set extrapolate_flag = 0
c
c############################################################################
c


      implicit none

      include 'dow.inc'
      include 'structures.inc'

c---- Passed variables

      integer method       ! interpolation method:
                           !   1=Cressman
                           !   2=Barnes
      
      
      integer nfld                    ! number of data fields to be gridded
      character(len=8) fname(nfld)    ! field names
      integer nswp                    ! number of sweep files
      character(len=200) sfname(nswp) ! sweep file names
      integer extrapolate_flag        ! should extrapolation be allowed?
                                      !   1=yes (standard objective analysis)
                                      !   0=no (interpolation only)
      real minsum                     ! threshold for minimum sum of weights required to produce
                                      !   an objectively-analyzed observation
      real minrange                   ! minimum-range threshold (data closer to radar
                                      !   are discarded)
  
c
c.... CLZ (11/23/10): rhmax is NEW maximum allowable horizontal range to oban a datum
c
  
      real rhmax                ! rhmax (km) = cos(dtor * elmax) * rmax
      real hrange               ! horizontal range to gate (km)

      integer yrcor                   ! correction to year
      integer mocor                   ! correction to month
      integer dacor                   ! correction to day
      real elcor                      ! elevation angle correction (deg)
      real azcor                      ! azimuth-angle offset (deg)
      integer az_corr_flag            ! method of additional azimuth-angle correction
                                      !   0 = none
                                      !   1 = average current and previous angle
      integer umass_flag              ! apply UMass data corrections? (1=yes, 0=no)
      integer map_proj                ! map projection (for relating lat, lon to x, y):
                                      !   0 = flat earth
                                      !   1 = oblique azimuthal
                                      !   2 = Lambert conformal
      real glat, glon                 ! latitude and longitude of grid origin (deg)
      real galt                       ! altitude of grid origin (km MSL)
      real rlat, rlon                 ! radar latitude and longitude (deg)
      real ralt                       ! radar altitude (km MSL)
      integer nx, ny, nz              ! no. of grid points in x, y, and z directions
      real dx, dy, dz                 ! grid spacing in x, y, and z directions (km)
      real xmin, ymin, zmin           ! coordinates of lower southwest corner
                                      !   of grid, relative to origin (km)
      integer(kind=2) cyr,cmo,cda     ! central date
      integer(kind=2) chr,cmn,cse     ! central time
      real ut, vt                     ! storm translation velocity (m/s)

c---- Returned variables

      real f(nx,ny,nz,nfld+3)             ! data fields
      real f_mpass(nx,ny,nz,nfld+3,npass) ! multi-pass data fields
      real az(nx,ny,nz)                   ! interpolated azimuth angle (deg)
      real el(nx,ny,nz)                   ! interpolated elevation angle (deg)
      real time(nx,ny,nz)                 ! time, relative to central time (sec)
      integer count(nx,ny,nz,nfld)        ! no. of gates used in interpolation


c---- Local variables
      real hsp, vsp   ! horizontal and smoothing parameters:
                      !   method=1:  hsp and vsp are the Cressman radii of influence (km)
                      !   method=2:  hsp and vsp are the Barnes smoothing parameters (km*km)  
      real c1                         ! hsp*hsp
      real c2                         ! vsp*vsp




      real glatr, glonr               ! latitude and longitude of grid origin (rad)
      real rlatr, rlonr               ! radar latitude and longitude (rad)
      integer s                       ! sweep file number
      integer n                       ! field number
      integer r                       ! ray number
      integer g                       ! gate number
      integer d                       ! directional bin number
      real x, y, z                    ! location of gate relative to grid origin (km)
      real ti                         ! time (s) relative to central time
      real timediff, comptoaz         ! functions used by this subroutine
      real rangekm_to_gate            ! "                               "
      integer(kind=2) cms; parameter(cms=0.0)
      real xrad, yrad, zrad           ! coords. (km) of radar with respect to origin
      integer i, j, k                 ! grid indices
      integer imin, imax              ! minimimum and maximum i for search
      integer jmin, jmax              ! minimimum and maximum j for search
      integer kmin, kmax              ! minimimum and maximum k for search
      integer irad, jrad, krad        ! radius of search area
      integer ii, jj, kk              ! nearest grid point to gate location
      real dh2                        ! square of horizontal distance (km*km)
      real dv2                        ! square of vertical distance (km*km)
      real wgt                        ! Cressman weight

      real range                      ! distance from radar to gate (km)
      integer valid_sum               ! 1 (0) if the sum of weights is (not) sufficient
                                      !   to produce an objectively analyzed observations

      type(vold_info)                     :: vold
      type(radd_info)                     :: radd
      type(celv_info)                     :: celv
      type(cfac_info)                     :: cfac
      type(parm_info), dimension(maxflds) :: parm
      type(swib_info)                     :: swib
      type(ryib_info), dimension(maxrays) :: ryib
      type(asib_info), dimension(maxrays) :: asib
      type(rdat_info), dimension(maxrays) :: rdat

      real, allocatable :: sumf(:,:,:,:)     ! sum of weights for each data field
      real, allocatable :: sumdir(:,:,:,:,:) ! sum of weights in directional bins
      real, allocatable :: azcomp(:,:,:,:)   ! azimuth angle components
      real, allocatable :: xg(:)             ! x coordinates (km) of grid points
      real, allocatable :: yg(:)             ! y coordinates (km) of grid points
      real, allocatable :: zg(:)             ! z coordinates (km) of grid points

c-----added variables (Mario)

      integer pass, npass
      integer maxray, maxgate
      real hsp0, vsp0, gamma, goodf
      real, allocatable :: baddata(:,:,:)
      real, allocatable :: f_old(:,:,:,:)
      integer, allocatable :: f_flag(:,:,:,:)      
      real, allocatable :: f_b(:,:,:,:)     ! backinterpolated data 
      real, allocatable :: sumf_b(:,:,:,:)  ! backinterpolated sum of wgt
      real, allocatable :: count_b(:,:,:,:) ! backinterpolated count
      integer, allocatable :: nray(:)
      integer, allocatable :: ngate(:)



c
c.... add variables (CLZ)
c

      real, allocatable :: gate_data(:,:)
      real dsbad
      parameter (dsbad = 10.0)
      integer deci
      integer maxgate_out
      integer numgate
      integer iprint
      integer actual_gates
      real ndavg
      real sumgd
      integer gg













c
c.... begin executable code
c
















c
c############################################################################
c
c     Initialize variables.
c
c############################################################################ 
c



      goodf = 0.0

      allocate (f_flag(nx,ny,nz,nfld))
      f_flag(:,:,:,:) = 0

      allocate (f_old(nx,ny,nz,nfld))
      allocate (nray(nswp)) 
      allocate (ngate(nswp))
 
      allocate(sumf(nx,ny,nz,0:nfld))


      if (extrapolate_flag.eq.0) then

        allocate(sumdir(8,nx,ny,nz,nfld))
c        allocate(sumdir(4,nx,ny,nz,nfld))

      endif



      
      allocate(azcomp(nx,ny,nz,2))
      allocate(xg(nx))
      allocate(yg(ny))
      allocate(zg(nz))
      
c      write(6,*) 'error check'



c
c.... CLZ (5/18/07):  David Dowell's oban code modified for multipass.
c      Programmed by Mario Majcen, Paul Markowski, and Yvette Richardson (PSU)
c      during January-April 2007.



c-----add loop through passes (Mario)
c.... MM (1/1/07):  add multi-pass do-loop
c













c
c------------------------------------------------------------------
c------------------------------------------------------------------
c



      azcomp(:,:,:,:) = 0.0



      pass = 1

c      if (pass .eq. 1) then
       f_old(:,:,:,:) = 0.0
c      endif

      el(:,:,:) = 0.0
      time(:,:,:) = 0.0
      count(:,:,:,:) = 0
      sumf(:,:,:,:) = 0.0

      if (extrapolate_flag .eq. 0) then
       sumdir(:,:,:,:,:) = 0.0
      endif







c
c.... CLZ (9/26/08):  "Base" h- and z-smoothing parameter values are computed below for
c         the current pass number.
c


      hsp = hsp0 * (gamma**(pass-1))
      write (6,*) 'on pass = ',pass,' hsp = ',hsp

      vsp = vsp0 * (gamma**(pass-1))
      write (6,*) 'on pass = ',pass,' vsp = ',vsp





      c1 = hsp * hsp
      c2 = vsp * vsp

      if (method.eq.1) then
        irad = 1+nint(hsp/dx)
        jrad = 1+nint(hsp/dy)
        krad = 1+nint(vsp/dz)
      else if (method.eq.2) then
        irad = 1+nint(sqrt(7.0*hsp)/dx)
        jrad = 1+nint(sqrt(7.0*hsp)/dy)
        krad = 1+nint(sqrt(7.0*vsp)/dz)
      else
        write(6,*) 'unknown method:  ', method
        stop
      endif ! if (method.eq.1) 









c
c.... CLZ (9/26/08):  Original location of following 4 lines was just below "c2 = ...".
c

      f(:,:,:,:) = 0.0

      call grid_coordinates(xg, nx, dx, xmin)
      call grid_coordinates(yg, ny, dy, ymin)
      call grid_coordinates(zg, nz, dz, zmin)






      glatr = dtor*glat
      glonr = dtor*glon
      rlatr = dtor*rlat
      rlonr = dtor*rlon
      call ll_to_xy(xrad, yrad, map_proj, glatr, glonr, rlatr, rlonr)
      zrad = ralt - galt
      write(6,*) 'location of radar relative to grid origin:'
      write(6,*) xrad, yrad, zrad









c
c############################################################################
c
c     Read in sweep files and compute weighted sums.
c
c############################################################################
c


      open(unit = 52, file = 'output_gates.dat',form='formatted',
     * status = 'unknown')





c
c.... CLZ (5/18/07):  First scan the sweep file to obtain control parameters
c


       if (pass.eq.1) then

       do s=1, nswp

        do n=1, nfld
           write(6,*) 'Now on pass = ', pass
           write(6,*) '1st scan reading from ', sfname(s)
           call sweepread(sfname(s),vold,radd,celv,cfac,parm,swib,
     $                    ryib,asib,rdat,fname(n))

      write(6,*) 'finished 1st call to sweepread in multioban'

           call check_sweep_size(radd%num_param_desc, swib%num_rays)
           call correct_date(vold%year, yrcor, vold%mon, mocor, vold%day, dacor)
           call correct_az_el(ryib, swib%num_rays, azcor, elcor, az_corr_flag)

           if (umass_flag.eq.1) then
             call correct_umass_data(swib%num_rays, ryib, celv%total_gates, rdat)
           endif
           ngate(s)=celv%total_gates      
           nray(s)=swib%num_rays
         enddo
       enddo


c
c           find max number of rays and max number of gates in a ray (Mario)
c

c
c.... CLZ (2/17/10):  P. Markowski reported bug, should have maxgate = ngate(s) in if-test
c

            maxray = nray(1)
            maxgate = ngate(1)

            do s = 1, nswp


              if (nray(s).gt.maxray) then
                 maxray=  nray(s)
              endif

c
c.... must update maxgate as well as maxray (previous code mistakenly duplicated)
c

c              if (nray(s).gt.maxray) then
c                 maxray=  nray(s)
c              endif

              if (ngate(s) .gt. maxgate) then
                 maxgate =  ngate(s)
              endif

c
c.... enddo s=1, nswp
c

            enddo

c
c.... for gate print-out, decimate by ndec
c

      deci = 20

      write(6,*)
      write(6,*)'deci=',deci
      
      maxgate_out = maxray * ((maxgate/2*deci) + 10)
      write(6,*)
      write(6,*)'maxray=',maxray,'maxgate=',maxgate,'maxgate_out=',
     * maxgate_out



        allocate (gate_data(maxgate_out,nfld))

        allocate(baddata(nswp,maxray,maxgate))
        baddata(:,:,:)= 0

c
c.... endif pass.eq.1
c

      endif










     
c
c.... CLZ (2/7/11):  Scan through sweep to read and output selected gate data
c    





      do s = 1, nswp



        do n = 1, nfld

c
c.... zero the number of gate for this field
c

         numgate = 0


           write(6,*) 'Now on pass = ', pass
           write(6,*) '2nd scan reading from ', sfname(s)
           call sweepread(sfname(s),vold,radd,celv,cfac,parm,swib,
     $                    ryib,asib,rdat,fname(n))

           call check_sweep_size(radd%num_param_desc, swib%num_rays)
           call correct_date(vold%year, yrcor, vold%mon, mocor, vold%day, dacor)
           call correct_az_el(ryib, swib%num_rays, azcor, elcor, az_corr_flag)

           if (umass_flag.eq.1) then
             call correct_umass_data(swib%num_rays, ryib, celv%total_gates, rdat)
           endif
           
           
            
           do r = 1, swib%num_rays

            ti = timediff(vold%year, vold%mon, vold%day,
     $                    ryib(r)%hour, ryib(r)%min, ryib(r)%sec, ryib(r)%msec,
     $                    cyr, cmo, cda, chr, cmn, cse, cms)



            do g = deci + 1, celv%total_gates - deci - 1, 2 * deci + 1

c            do g = 1, celv%total_gates, 10



c           if (rdat(r)%data(g) .lt. (sbad + dsbad) .and.
c     *       rdat(r)%data(g) .gt. (sbad - dsbad) ) then

              if (rdat(r)%data(g) .eq. sbad ) then
                 baddata(s,r,g) = 1
              endif



c           if (rdat(r)%data(g) .gt. (sbad + dsbad) .or.
c     *       rdat(r)%data(g) .lt. (sbad - dsbad) ) then

c              if (rdat(r)%data(g) .ne. sbad) then


               numgate = numgate + 1


c
c.... compute the slant range (km) to the gate
c

                range = rangekm_to_gate(celv, g)

c
c.... CLZ (11/23/10): now compute the horizontal range (km) from the slant range
c

                hrange = cos(dtor * ryib(r)%elevation) * range

                call xyzloc(x, y, z, range, dtor*ryib(r)%azimuth, dtor*ryib(r)%elevation,
     $                      map_proj, glatr, glonr, galt, rlatr, rlonr, ralt,
     $                      ut, vt, ti)


                ii = 1 + nint((x-xmin)/dx)
                jj = 1 + nint((y-ymin)/dy)
                kk = 1 + nint((z-zmin)/dz)

c
c.... CLZ (9/26/08):  current gate is at grid level = kk at this point in the code
c

                imin = max(ii-irad, 1)
                imax = min(ii+irad, nx)
                jmin = max(jj-jrad, 1)
                jmax = min(jj+jrad, ny)
                kmin = max(kk-krad, 1)
                kmax = min(kk+krad, nz)

               if (range .ge. minrange .and. hrange .le. rhmax) then

               ndavg = 0
               sumgd = 0.0

              do gg = g - deci, g + deci
               
               if (rdat(r)%data(gg) .ne. sbad ) then
                 sumgd = sumgd + rdat(r)%data(gg)
                 ndavg = ndavg + 1.0
               endif

              enddo

               gate_data(numgate,n) = sumgd / ndavg
c               gate_data(numgate,n) = rdat(r)%data(g)


c               write(52,*) rdat(r)%data(g)

               endif       ! (range.ge.minrange)



c              endif       ! (rdat(r)%data(g) .ne. sbad)




            enddo       ! g=1, celv%total_gates
             
          enddo       ! r=1, swib%num_rays
            
        enddo       ! n=1, nfld
        
      enddo       ! s=1, nswp











c
c------------------------------------------------------------------
c------------------------------------------------------------------
c



      write(6,*)
      write(6,*)'number of gates averaged=',numgate

      write(52,*) (fname(n),n = 1, nfld)

      actual_gates = 0


      do j = 1, numgate

      iprint = 1
      do n = 1, nfld
      if (gate_data(j,n) .eq. sbad .or. gate_data(j,n)
     * .eq. 0.0) iprint = 0
      enddo
      
      if (iprint .eq. 1 .and. gate_data(j,2) .gt. 15.0) then
      write(52,*) (gate_data(j,n), n = 1, nfld)
      actual_gates = actual_gates + 1
      endif

      enddo




      close(52)


      write(6,*) 'Finished outputting',actual_gates,' gate values'







c
c.... deallocate
c



c      deallocate(f_old)
c      deallocate(f_b)
c      deallocate(sumf_b)
c      deallocate(count_b)
c      deallocate(nray)
c      deallocate(ngate)
      
c      deallocate(sumf)
c      if (extrapolate_flag.eq.0) then
c        deallocate(sumdir)
c      endif
c      deallocate(azcomp)
c      deallocate(xg)
c      deallocate(yg)
c      deallocate(zg)
c      deallocate(f_flag)
c      deallocate(baddata)




      write(6,*)
      write(6,*) 'Done with output_gates.f, now returning to oban.f'
      write(6,*)

      return
      end



