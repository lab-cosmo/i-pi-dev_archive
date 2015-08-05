      subroutine apply_constrain(en, force, nat, lab, coords)
      implicit none
      real*8, intent(inout) :: force(nat, 3), en
      integer, intent(in) :: nat
      character*4, intent(in) :: lab(nat)
      real*8, intent(in) :: coords(nat,3) 
      integer :: norder(nat), i
      character(*), parameter :: func = 'apply_constrain'
      do i=1, nat
         if (lab(i).eq.'O') then
            norder(i)=1
         else
            norder(i)=0
         endif
      enddo
 
      en=0.d0
      call cylinder(norder, nat, en, force, coords)
      end subroutine

      subroutine cylinder(norder, n_atoms, pot, forces, coords)
      ! Cylinder along x direction that mimics a 6,6 carbon nanotube
      ! Analytical form: V(r)=a2*r^2+a4*r^4+a6*r^6+a8*r^8
      ! and parameters following Dellago and Naor, CPC 169 (2005) 36--39
      implicit none
      integer :: n_atoms
      integer, intent(in) :: norder(n_atoms)
      real*8 :: a2, a4, a6, a8
      real*8  :: pot, forces(n_atoms, 3), coords(n_atoms, 3), r
      integer i_atom, i_dim
      a2=-0.000101790427
      a4=0.0001362104651
      a6=8.1919580588d-06
      a8=3.188645093e-06
      pot=0.d0
      forces=0.d0
      do i_atom=1, n_atoms
         if (norder(i_atom).eq.1) then
            call cylindrical_distance(coords(i_atom,:), r)
            pot=pot+a2*(r**2)+a4*(r**4)+a6*(r**6)+a8*(r**8)
            if (r.gt.0.d0) then
               do i_dim=2, 3
                  forces(i_atom,i_dim)=-(2*a2*r+4*a4*(r**3)+6*a6*(r**5)+8*a8*(r**7))*coords(i_atom, i_dim)/r
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

