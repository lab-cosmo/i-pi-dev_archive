      module LJ
         use sys_vars
      implicit none

      private
      public :: get_all, nearest_neighbours
      public :: get_all2, nearest_neighbours2

!      double precision, parameter :: sigma = 0.8d0
!      double precision, parameter :: rc = 2.5*sigma
!      double precision, parameter :: rn = 300*rc
!      double precision, parameter :: eps = 1.0d0
      double precision, parameter :: sigma = 6.43452d0
      double precision, parameter :: rc = 46.651d0  !2.5*sigma
      double precision, parameter :: rn = 1.2*rc
      double precision, parameter :: eps = 0.0003793865d0

       double precision, parameter :: correction =  
     c      4*eps*((sigma/rc)**12 - (sigma/rc)**6)

      contains

         subroutine LJ_functions(r, pot, force)
            double precision, intent(in) :: r
            double precision, intent(out) :: pot
            double precision, intent(out) :: force

            double precision sigma_by_r6

            if (r > rc) then
               pot = 0.0d0
               force = 0.0d0
               return
            end if
      
            sigma_by_r6 = (sigma/r)**6

            pot = 4*eps*(sigma_by_r6*(sigma_by_r6 - 1)) - correction

            force = 4*eps*(6/r*sigma_by_r6*(2*sigma_by_r6 - 1))
      
         end subroutine

         subroutine separation(atoms, i, j, cell, rij, r)
            type(Atom), dimension(:), intent(in) :: atoms
            integer, intent(in) :: i
            integer, intent(in) :: j
            type(Cell_vec), intent(in) :: cell
            double precision, dimension(3), intent(out) :: rij
            double precision, intent(out) :: r

            integer k
            double precision, dimension(3) :: s

            s = matmul(cell%ih, atoms(i)%pos - atoms(j)%pos)
            do k = 1, 3
               s(k) = s(k) - dnint(s(k))
            end do
            rij = matmul(cell%h, s)
            r = sqrt(dot_product(rij, rij))

         end subroutine

         subroutine LJ_fij(atoms, i, j, cell, fij, rij, pot)
            type(Atom), dimension(:), intent(in) :: atoms
            integer, intent(in) :: i
            integer, intent(in) :: j
            type(Cell_vec), intent(in) :: cell
            double precision, dimension(3), intent(out) :: fij
            double precision, dimension(3), intent(out) :: rij
            double precision, intent(out) :: pot

            double precision f_tot, r
            integer k

            call separation(atoms, i, j, cell, rij, r)
            call LJ_functions(r, pot, f_tot)

            do k = 1, 3
               fij(k) = f_tot*rij(k)/r
            end do

         end subroutine

         subroutine nearest_neighbours(atoms, cell, n_list)
            type(Atom), dimension(:), intent(in) :: atoms
            type(Cell_vec), intent(in) :: cell
            logical, dimension(size(atoms),size(atoms)), intent(out)
     1:: n_list

            double precision, dimension(3) :: rij
            double precision r
            integer i, j

            n_list = .false.
            do i = 1, size(atoms)-1
               do j = i+1, size(atoms)
                  call separation(atoms, i, j, cell, rij, r)
                  if (r < rn) then
                     n_list(i,j) = .true.
                     n_list(j,i) = .true.
                  end if
               end do
            end do

         end subroutine 

         subroutine get_all(atoms, cell, n_list, pot, f, vir)
            type(Atom), dimension(:), intent(in) :: atoms
            type(Cell_vec), intent(in) :: cell
            logical, dimension(size(atoms),size(atoms)), intent(in) ::
     c         n_list
            double precision, intent(out) :: pot
            double precision, dimension(3,size(atoms)), 
     c          intent(out) :: f
            double precision, dimension(3,3), intent(out) :: vir
            
            double precision volume
            integer i, j, k, l
            double precision, dimension(3) :: fij, rij
            double precision pot_ij

            vir = 0.0d0
            pot = 0.0d0
            f = 0.0d0
            volume = cell%h(1,1)*cell%h(2,2)*cell%h(3,3)

            do i = 1, size(atoms)-1
               do j = i+1, size(atoms)
                  !if (n_list(i,j)) then
                     call LJ_fij(atoms, i, j, cell, fij, rij, pot_ij)

                     f(:,i) = f(:,i) + fij
                     f(:,j) = f(:,j) - fij
                     pot = pot + pot_ij
                  
                     do k = 1, 3
                        do l = k, 3
                           vir(k,l) = vir(k,l) + fij(k)*rij(l)
                        end do
                     end do
                  !end if

               end do
            end do
            !vir = vir/volume

         end subroutine

         subroutine nearest_neighbours2(atoms, cell, n_list, index_list)
            type(Atom), dimension(:), intent(in) :: atoms
            type(Cell_vec), intent(in) :: cell
            integer, dimension(size(atoms)-1), intent(out)
     c                  :: n_list
            integer, dimension(size(atoms)*(size(atoms)-1)/2),
     c                     intent(out) :: index_list

            double precision, dimension(3) :: rij
            double precision r
            integer i, j, counter

            n_list = 0
            index_list = 0
            counter = 1

            do i = 1, size(atoms)-1
               do j = i+1, size(atoms)
                  call separation(atoms, i, j, cell, rij, r)
                  if (r < rn) then
                     n_list(i) = n_list(i) + 1
                     index_list(counter) = j
                     counter = counter + 1
                  end if
               end do
            end do

         end subroutine

         subroutine get_all2(atoms, cell, n_list, index_list, 
     1pot, f, vir)
            type(Atom), dimension(:), intent(in) :: atoms
            type(Cell_vec), intent(in) :: cell
            integer, dimension(size(atoms)-1), intent(in) :: n_list
            integer, dimension(size(atoms)*(size(atoms)-1)/2), 
     1intent(in) :: index_list
            double precision, intent(out) :: pot
            double precision, dimension(3,size(atoms)), 
     1intent(out) :: f
            double precision, dimension(3,3), intent(out) :: vir
            
            double precision volume
            integer i, j, k, l
            integer start, finish
            double precision, dimension(3) :: fij, rij
            double precision pot_ij

            vir = 0.0d0
            pot = 0.0d0
            f = 0.0d0
            volume = cell%h(1,1)*cell%h(2,2)*cell%h(3,3)

            start = 1

            do i = 1, size(atoms)-1
               finish = start + n_list(i) - 1
               do j = start, finish !if n_list(i) == 0, what happens?
                  call LJ_fij(atoms, i, index_list(j), cell, fij, 
     1rij, pot_ij)

                  f(:,i) = f(:,i) + fij
                  f(:,index_list(j)) = f(:,index_list(j)) - fij
                  pot = pot + pot_ij
               
                  do k = 1, 3
                     do l = k, 3
                        vir(k,l) = vir(k,l) + fij(k)*rij(l)
                     end do
                  end do

               end do
               start = finish + 1
            end do

         end subroutine

      end module
