! Implementing the interpolation in fortran due to speed issues.
! The numpy implementation was only 3 times faster than calling
! slalib directly.
!
! I will try to keep the ordering memory efficient too, this time,
! with x: (incomp,nsamp),  y: (oncomp,nsamp)
!   xbox: (incomp,2),      n: (incomp)
!  ygrid: (oncomp,ngrid), dy: (incomp,oncomp,ngrid)
!
! When collapsing dimensions, use C ordering, as it has no
! impact on performance in this case.

function ipol(x, xbox, n, ygrid, dygrid) result(y)
  implicit none
  real*8       :: x(:,:), xbox(:,:), ygrid(:,:), dygrid(:,:,:)
  real*8       :: y(size(ygrid,1),size(x,2))
  real*8       :: x0(size(x,1)), idx(size(x,1))
  real*8       :: xrel(size(x,1))
  integer*4    :: n(:), xind(size(x,1))
  integer*4    :: incomp, oncomp, nsamp, ngrid, ic, oc, is, ig
  integer*4    :: steps(size(x,1))
  incomp = size(x,1); oncomp = size(ygrid,1)
  nsamp  = size(x,2); ngrid  = size(ygrid,2)

  ! First build the nD to 1D translation
  steps(incomp) = 1
  do ic = incomp-1, 1, -1
     steps(ic) = steps(ic+1)*n(ic+1)
  end do
  x0 = xbox(:,1); idx = (n-1)/(xbox(:,2)-xbox(:,1))
  ! Then do the actual lookup
  do is = 1, nsamp
     xrel = (x(:,is)-x0)*idx
     xind = floor(xrel+0.5)
     xrel = xrel - xind
     ig   = sum(xind*steps)+1
     do oc = 1, oncomp
        y(oc,is) = ygrid(oc,ig) + sum(dygrid(:,oc,ig)*xrel)
     end do
  end do
end function
