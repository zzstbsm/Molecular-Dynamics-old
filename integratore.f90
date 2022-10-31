SUBROUTINE RungeKutta(func,x,t,step)

	USE parametri

	
	IMPLICIT NONE

	DOUBLE PRECISION, INTENT(INOUT)	:: x(1:6,1:nAtoms)
	DOUBLE PRECISION, INTENT(IN)	:: t, step

	DOUBLE PRECISION	:: c(0:3)=(/ 0., 0.5, 0.5, 1. /) ! Nodi
	DOUBLE PRECISION	:: b(0:3)=(/ 1./6, 1./3, 1./3, 1./6 /) ! Pesi
	DOUBLE PRECISION	:: k(0:4,1:6,1:nAtoms) ! k(0)=0 e serve per ottimizzare l'algoritmo

	DOUBLE PRECISION	:: x_t(1:6,1:nAtoms)
	DOUBLE PRECISION	:: tt

	INTEGER				:: j

	EXTERNAL			:: func

	! Valuta i coefficienti
	k = 0
	DO j=1,4
		x_t = x + k(j-1,:,:)*step*c(j-1)
		tt = t + step*c(j-1)
		CALL func(k(j,:,:),x_t(:,:),tt)
		
	END DO

	DO j=1,4
		x = x + k(j,:,:)*b(j-1)*step
	END DO
	
END SUBROUTINE RungeKutta

SUBROUTINE Midpoint(func,x,t,step)

	USE parametri

	
	IMPLICIT NONE

	DOUBLE PRECISION, INTENT(INOUT)	:: x(1:6,1:nAtoms)
	DOUBLE PRECISION, INTENT(IN)	:: t, step

	DOUBLE PRECISION	:: c(0:1)=(/ 0., 0.5 /) ! Nodi
	DOUBLE PRECISION	:: b(0:1)=(/ 0., 1. /) ! Pesi
	DOUBLE PRECISION	:: k(0:2,1:6,1:nAtoms) ! k(0)=0 e serve per ottimizzare l'algoritmo

	DOUBLE PRECISION	:: x_t(1:6,1:nAtoms)
	DOUBLE PRECISION	:: tt

	INTEGER				:: j

	EXTERNAL			:: func

	! Valuta i coefficienti
	k = 0
	DO j=1,2
		x_t = x + k(j-1,:,:)*step*c(j-1)
		tt = t + step*c(j-1)
		CALL func(k(j,:,:),x_t(:,:),tt)
		
	END DO

	DO j=1,2
		x = x + k(j,:,:)*b(j-1)*step
	END DO
	
END SUBROUTINE Midpoint

SUBROUTINE Eulero(func,x,t,step)
	

	USE parametri

	
	IMPLICIT NONE

	DOUBLE PRECISION, INTENT(INOUT)	:: x(1:6,1:nAtoms)
	DOUBLE PRECISION, INTENT(IN)	:: t, step

	DOUBLE PRECISION	:: k(1:6,1:nAtoms)


	! Valuta i coefficienti
	k = 0
	CALL func(k(:,:),x(:,:),t)

	x = x + k(:,:)*step
		
END SUBROUTINE Eulero

SUBROUTINE Verlet(func,x,t,step)
	
	USE parametri
	
	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(INOUT)	:: x(1:6,1:nAtoms)
	DOUBLE PRECISION, INTENT(IN)	:: t, step
	
	DOUBLE PRECISION	:: k1(1:6,1:nAtoms),k2(1:6,1:nAtoms)
	
	! Valuta i coefficienti
	k1 = 0
	k2 = 0
	CALL func(k1(:,:),x(:,:),t)
	
	
	! k1(4:6,:) è l'accelerazione
	x(1:3,:) = x(1:3,:) + x(4:6,:)*step + (k1(4:6,:)*step**2)/2
	CALL func(k2(:,:),x(:,:),t+step)
	x(4:6,:) = x(4:6,:) + (k1(4:6,:)+k2(4:6,:))*step/2
	
	
END SUBROUTINE Verlet
