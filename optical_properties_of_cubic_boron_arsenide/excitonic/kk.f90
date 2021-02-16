PROGRAM MAIN
  IMPLICIT NONE

  INTEGER, PARAMETER :: nomega = 2601
  REAL(16), PARAMETER :: pi = 3.14159265359
  REAL(16) :: omega(nomega), epsilon2(nomega), epsilon1(nomega)
  REAL(16) :: omega_p


  INTEGER :: i,j,k
  REAL(16) :: x,y,z, sum
  REAL(16) :: dummy
  OPEN(UNIT=1,FILE='eps2.dat')
  !READ(1,*)
  !READ(1,*)
  DO i = 1, nomega
     READ(1,*) omega(i), epsilon2(i)
  END DO
  CLOSE(1)

  DO i = 1, nomega

     sum = 0.0

     DO j = 1, nomega-1
        omega_p = 0.5*(omega(j)+omega(j+1))
        !sum = sum + 0.5*( epsilon2(j) + epsilon2(j+1))*(omega(j+1)-omega(j))/(omega_p - omega(i))
        sum = sum + omega_p*((epsilon2(j)+epsilon2(j+1))*0.5)/(omega_p**2 - omega(i)**2)*(omega(j+1)-omega(j))
     END DO
     !epsilon1(i) = 1.0 + sum/pi
     epsilon1(i) = 1.0 + 2*sum/pi
     WRITE(*,"(2e20.10)") omega(i), epsilon1(i)
  END DO

  


END PROGRAM MAIN

