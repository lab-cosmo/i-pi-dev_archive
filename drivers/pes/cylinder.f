      subroutine apply_constrain(en, force, nat, norder, at)
!!! Must pass coordinates, check
      implicit none
      real*8 :: force(3, nat), en
      integer :: nat
      real*8 :: norder(nat), at(3,nat) ! must transform O to 1 and H to 0
      character(*), parameter :: func = 'apply_constrain'
      en=0.d0
      call cylinder(norder, nat, en, force)
      end subroutine

      subroutine cylinder(norder, n_atoms, pot, forces)
      ! Cylinder along x direction that mimics a 6,6 carbon nanotube
      ! Analytical form: V(r)=a2*r^2+a4*r^4+a6*r^6+a8*r^8
      ! and parameters following Dellago and Naor, CPC 169 (2005) 36--39
      implicit none
      real*8, intent(in) :: norder, n_atoms
      real*8 :: a2, a4, a6, a8
      real*8  :: pot, forces(3, n_atoms), r
      integer i_atom, i_dim
      a2=-0.000101790427
      a4=0.0001362104651
      a6=8.1919580588d-06
      a8=3.188645093e-06
      pot=0.d0
      forces=0.d0
      do i_atom=1, n_atoms
         if (norder(i_atom).eq.1) then
            call cylindrical_distance(coords(:,i_atom), r)
            pot=pot+a2*(r**2)+a4*(r**4)+a6*(r**6)+a8*(r**8)
            if (r.gt.0.d0) then
               do i_dim=2, 3
                  forces(i_dim,i_atom)=-(2*a2*r+4*a4*(r**3)+6*a6*(r**5)+8*a8*(r**7))*coords(i_dim, i_atom)/r
               enddo
            endif
         endif
      enddo
      end subroutine

      subroutine cylindrical_distance(q, r)
      implicit none
      real*8 :: q(3), r
         r=dsqrt(q(2)**2 + q(3)**2)
      end subroutine 
      END MODULE

