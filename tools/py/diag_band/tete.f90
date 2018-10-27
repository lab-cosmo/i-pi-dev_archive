module math 
IMPLICIT NONE
contains
!----------------------------------------------------------------
!----------------------------------------------------------------

subroutine diag2(aa,eigv,n)
!FUNCTION diag2(a,eigv,n,verbose) RESULT(eigv)

	REAL(KIND=8),DIMENSION(n,n)  :: aa,aux_a 
	INTEGER ( kind =4 ) :: n,i
	INTEGER,DIMENSION(n) :: jndx
	REAL(KIND=8), DIMENSION(n) :: DP
	REAL(KIND=8),DIMENSION(n)     ::  eigv,aux_eigv

	CALL TRED2(aa,n  ,n,EIGV,DP)
	CALL TQLI (EIGV,DP,n,n,aa)
 
	aux_eigv = eigv
	aux_a = aa
	CALL SORTING (n,aux_eigv,JNDX)

	DO i=1,n
	  eigv(i) = aux_eigv(jndx(i))
	  aa(1:n,i) = aux_a(1:n,jndx(i))
	ENDDO
	!write(6,*) 'EIGV',(EIGV(I),I=1,n)


!end function
end subroutine
!----------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine diagband(A,lambda,N,KD)

      !JOBZ ='V' 
      !SET LOWER FORM
      !        3.1 3.2 3.3 *   * 
      ! AB =   2.1 2.2 2.3 2.4 *   
      !        1.1 1.2 1.3 1.3 1.4

      INTEGER      :: INFO                     !INFO =0 means no error 
      CHARACTER(1) :: JOBZ                     !'N' compute evalues only, 'V' compute evectors.
      CHARACTER(1) :: UPLO                     !'U' Upper triangle, 'L' lower triangle
      INTEGER      :: N                        !Length of main diagonal 
      INTEGER      :: KD                       !THe number of super-diagonals (sub-diagonals)
      INTEGER      :: I,J
      REAL(KIND=8), DIMENSION(KD+1,N)  :: A      !Banded matrix
      REAL(KIND=8), DIMENSION(N)       :: lambda !if INFO =0, W contains the eigenvalues in ascending order   
      REAL(KIND=8), DIMENSION(1,3*N-2) :: WORK
      REAL(KIND=8), DIMENSION(2,2)     :: Z      !if INFO = 0, Z contains the orthonormal eigenvectors of the matrix A,


      UPLO ='L'
      JOBZ ='N'
      !DO I=1,KD+1
      !write(6,*) (A(I,J), J=1,N)
      !ENDDO
      CALL DSBEV(JOBZ, UPLO, N, KD, A, KD+1, lambda, Z, N, WORK, INFO)
      
      !WRITE (6 ,99999) (lambda(J),J=1,N)
!99999 FORMAT (3X,(8F8.4))

      END
!------------------------------------------------------------------
!------------------------------------------------------------------
end module math
