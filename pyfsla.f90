! These functions expect their inputs as right handed north polar coordinates!
function aomulti(mjd, icoord, ao, am) result(ocoord)
  implicit none
  real*8    :: ao(:), am(:), mjd(:), icoord(:,:), ocoord(size(icoord,1),size(icoord,2))
  real*8    :: ra, dec, pi
  integer*4 :: i, n
  n  = size(mjd)
  pi = 4*atan(1d0)
  do i = 1, n
     call sla_aoppat(mjd(i), ao)
     call sla_oapqk("A", -icoord(1,i), icoord(2,i), ao, ra, dec)
     call sla_ampqk(ra, dec, am, ocoord(1,i), ocoord(2,i))
     ocoord(2,i) = pi/2-ocoord(2,i)
  end do
end function

function equ2gal(icoord) result(ocoord)
  implicit none
  real*8    :: icoord(:,:), ocoord(size(icoord,1),size(icoord,2)), pi
  integer*4 :: i
  pi = 4*atan(1d0)
  do i = 1, size(icoord,2)
     call sla_eqgal(icoord(1,i), pi/2-icoord(2,i), ocoord(1,i), ocoord(2,i))
     ocoord(2,i) = pi/2-ocoord(2,i)
  end do
end function
