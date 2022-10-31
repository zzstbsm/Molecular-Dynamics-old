MODULE physics
	
	USE parametri
	
	DOUBLE PRECISION	:: energy
	DOUBLE PRECISION	:: kineticEnergy
	DOUBLE PRECISION	:: potentialEnergy
	DOUBLE PRECISION	:: pressure = -1
	DOUBLE PRECISION	:: temperature
	DOUBLE PRECISION	:: constPot = -(((sigma/cutOff)**12)-((sigma/cutOff)**6))
	DOUBLE PRECISION	:: momentum
	DOUBLE PRECISION	:: momentumVec(1:3)
	
	! callingCalcRefresh serve per capire se è stata chiamata calcRefresh
	INTEGER, PRIVATE	:: callingCalcRefresh = 0
	! computingForce serve per aggiornare pressure in modo idoneo
	INTEGER, PRIVATE	:: computingForce = 0
	! Viene aggiornato durante il calcolo delle forze, la pressione finale pressure viene calcolata da qui
	DOUBLE PRECISION, PRIVATE	:: tempPressure = 0
	
	DOUBLE PRECISION, PRIVATE	:: force(1:3,1:nAtoms)
	
	LOGICAL	:: callingMontecarlo = .FALSE.
	
CONTAINS

	SUBROUTINE calcDistance(distance,direction,position1,position2)

		IMPLICIT NONE

		DOUBLE PRECISION, INTENT(IN)	:: position1(1:3)
		DOUBLE PRECISION, INTENT(IN)	:: position2(1:3)

		DOUBLE PRECISION, INTENT(OUT)	:: distance
		! Direction non è normalizzato
		DOUBLE PRECISION, INTENT(OUT)	:: direction(1:3)
		
		DOUBLE PRECISION	:: halfBox
		
		! n assume valori che indicano se il secondo atomo va spostato e dove (condizioni periodiche
		! a destra o a sinistra
		INTEGER		:: n(1:3)
		DOUBLE PRECISION	:: dist(1:3)

		! Contatori
		INTEGER		:: i

		halfBox = boxLength/2

		! Se la differenza delle coordinate tra le due particelle su una direzione è maggiore halfBox, sposta 
		! la particella 2 a seconda della posizione della particella 1 rispetto alla particella (lo spostamento
		! avviene impostando un valore idoneo sulla variabile n)

		! Verifica che il secondo atomo sia nel cubo
		DO i=1,3
			
			dist(i) = ABS(position1(i)-position2(i))
			
			! La distanza è minore di halfBox
			IF (dist(i)<halfBox) THEN
				n(i) = 0
			! La distanza è maggiore di halfBox, sposta gli atomi
			ELSE
				dist(i) = boxLength-dist(i)
				! La particella 2 va spostata a sinistra dato che si trova a destra
				IF (position2(i)>position1(i)) THEN
					n(i) = -1
				! La particella 2 va spostata a destra dato che si trova a sinistra
				ELSE IF (position2(i)<position1(i)) THEN
					n(i) = +1
				END IF
			END IF

		END DO

		! Calcola la distanza
		direction = (position2-position1+n*boxLength)
		distance = SQRT((dist(1))**2 + (dist(2))**2 + (dist(3))**2)

	END SUBROUTINE calcDistance

	SUBROUTINE calcForceTwoParticles(force,position1,position2)

		IMPLICIT NONE

		DOUBLE PRECISION, INTENT(IN)	:: position1(1:3)
		DOUBLE PRECISION, INTENT(IN)	:: position2(1:3)

		DOUBLE PRECISION, INTENT(OUT)	:: force(1:3)

		! distance è la distanza, direction è il versore della forza
		DOUBLE PRECISION	:: distance, dist(1:3)
		DOUBLE PRECISION	:: direction(1:3), temp
		
		LOGICAL :: flagComputeDistance
		INTEGER	:: i
		
		
		flagComputeDistance = .TRUE.
		DO i=1,3
			dist(i) = ABS(position1(i)-position2(i))
			IF (dist(i) > cutOff .AND. dist(i) < (boxLength-cutOff)) THEN
				flagComputeDistance = .FALSE.
			END IF
		END DO
		
		IF (flagComputeDistance .eqv. .TRUE.) THEN
			CALL calcDistance(distance,direction(:),position1(:),position2(:))
		ELSE
			force = 0
			RETURN
		END IF

		IF (distance>cutOff) THEN
			force = 0
		ELSE
			temp = sigma/distance
			force = 4 * &
				( ((temp)**12 * 12) &
				- ((temp)**6 * 6) ) &
				*direction/(distance**2)
			
			! Questa parte aggiorna l'energia potenziale
			potentialEnergy = potentialEnergy + 4 * eps * (temp**12 - temp**6 + constPot)
			
			! Questa parte serve per il calcolo della pressione
			tempPressure = tempPressure + DOT_PRODUCT(direction,force)
			
		END IF
		
	END SUBROUTINE calcForceTwoParticles
	
!~ 	DOUBLE PRECISION FUNCTION calcPotEnergy(atoms,i)
		
!~ 		USE parametri
!~ 		USE list_mod
		
!~ 		IMPLICIT NONE
		
!~ 		DOUBLE PRECISION, INTENT(IN)	:: atoms(1:6,1:nAtoms)
!~ 		INTEGER, INTENT(IN)				:: i
		
!~ 		DOUBLE PRECISION				:: tempPot,distance,tempInv,dist(1:3)
!~ 		INTEGER							:: j,k,ix,iy,iz
!~ 		INTEGER							:: indexArray(1:nAtoms)
!~ 		LOGICAL							:: inCutOff,near
		
!~ 		tempPot = 0
		
!~ 		! Niente celle
!~ 		IF (m==1) THEN
			
!~ 			DO j=1,nAtoms
!~ 				IF (i==j) THEN
!~ 					CYCLE
!~ 				END IF
				
!~ 				inCutOff = .TRUE.
!~ 				DO k=1,3
!~ 					dist(k) = ABS(atoms(k,i) - atoms(k,j))
!~ 					IF (dist(k) > cutOff) THEN
!~ 						IF (dist(k)<(boxLength-cutOff)) THEN
!~ 							inCutOff = .FALSE.
!~ 						ELSE
!~ 							dist(k) = boxLength - dist(k)
!~ 						END IF
!~ 					END IF
!~ 				END DO
				
!~ 				! Fuori dal cut off
!~ 				IF (.NOT. inCutOff) THEN
!~ 					CYCLE
				
!~ 				! Potrebbe essere entro il cut off
!~ 				ELSE
!~ 					distance = SQRT(dist(1)*dist(1)+dist(2)*dist(2)+dist(3)*dist(3))
!~ 				END IF
				
!~ 				IF (distance<cutOff) THEN
!~ 					tempInv = 1/distance
!~ 					tempPot = tempPot + 4 *(tempInv**12 - tempInv**6 + constPot)
!~ 				END IF
				
!~ 			END DO
			
!~ 		! Con le celle
!~ 		ELSE
		
!~ 			! Identifica dove si trova l'atomo
!~ 			ix = atomsInCellCoordinates(1,i)
!~ 			iy = atomsInCellCoordinates(2,i)
!~ 			iz = atomsInCellCoordinates(3,i)
			
!~ 			DO j=1,nAtoms
!~ 				IF (i==j) THEN
!~ 					CYCLE
!~ 				END IF
!~ 				near = .FALSE.
!~ 				IF ((ix == atomsInCellCoordinates(1,j)) .OR. &
!~ 					(ix == MOD(atomsInCellCoordinates(1,j) + 1,m)).OR. &
!~ 					(ix == MOD(atomsInCellCoordinates(1,j) - 1,m))) THEN
					
!~ 					IF ((iy == atomsInCellCoordinates(2,j)) .OR. &
!~ 						(iy == MOD(atomsInCellCoordinates(2,j) + 1,m)).OR. &
!~ 						(iy == MOD(atomsInCellCoordinates(2,j) - 1,m))) THEN
						
!~ 						IF ((iz == atomsInCellCoordinates(3,j)) .OR. &
!~ 							(iz == MOD(atomsInCellCoordinates(3,j) + 1,m)).OR. &
!~ 							(iz == MOD(atomsInCellCoordinates(3,j) - 1,m))) THEN
!~ 							near = .TRUE.
!~ 						END IF
!~ 					END IF
!~ 				END IF
				
!~ 				! In una cella lontana
!~ 				IF (.NOT. near) THEN
!~ 					CYCLE
				
!~ 				! In una cella vicina
!~ 				ELSE
!~ 					DO k=1,3
!~ 						dist(k) = ABS(atoms(k,i) - atoms(k,j))
!~ 						IF (dist(k)>boxLength/2) THEN
!~ 							dist(k) = boxLength - dist(k)
!~ 						END IF
!~ 					END DO
!~ 					distance = SQRT(dist(1)**2+dist(2)**2+dist(3)**2)
					
!~ 					IF (distance<cutOff) THEN
!~ 						tempInv = 1/distance
!~ 						tempPot = tempPot + 4 *(tempInv**12 - tempInv**6 + constPot)
!~ 					END IF
!~ 				END IF
				
!~ 			END DO
			
!~ 		END IF
		
!~ 		calcPotEnergy = tempPot
		
!~ 	END FUNCTION calcPotEnergy
	
!~ 	SUBROUTINE calcForce_BK(x)
		
!~ 		USE list_mod
		
!~ 		IMPLICIT NONE
		
!~ 		DOUBLE PRECISION, INTENT(IN)	:: x(1:6,1:nAtoms)
!~ 		DOUBLE PRECISION	:: tempVec(1:3)

!~ 		INTEGER	:: i, j, ix, iy, iz
!~ 		LOGICAL	:: near

!~ 		! Calcola la forza: decide quale forze vanno valutate e valuta quelle
!~ 		force = 0
!~ 		potentialEnergy = 0
		
!~ 		IF (m == 1) THEN
!~ 			DO i=1,nAtoms,+1
!~ 				DO j=i+1,nAtoms,+1
					
!~ 					CALL calcForceTwoParticles(tempVec(:),x(1:3,i),x(1:3,j))
!~ 					force(:,i) = force(:,i) - tempVec
!~ 					force(:,j) = force(:,j) + tempVec
					
!~ 				END DO
!~ 			END DO
		
!~ 		ELSE
			
!~ 			! Calcola la forza: con la struttura a celle
!~ 			DO i=1,nAtoms,+1
		
!~ 				! Identifica dove si trova l'atomo
!~ 				ix = atomsInCellCoordinates(1,i)
!~ 				iy = atomsInCellCoordinates(2,i)
!~ 				iz = atomsInCellCoordinates(3,i)
				
!~ 				DO j=i+1,nAtoms,+1
					
					
!~ 					near = .FALSE.
!~ 					IF ((ix == atomsInCellCoordinates(1,j)) .OR. &
!~ 						(ix == MOD(atomsInCellCoordinates(1,j) + 1,m)).OR. &
!~ 						(ix == MOD(atomsInCellCoordinates(1,j) - 1,m))) THEN
						
!~ 						IF ((iy == atomsInCellCoordinates(2,j)) .OR. &
!~ 							(iy == MOD(atomsInCellCoordinates(2,j) + 1,m)).OR. &
!~ 							(iy == MOD(atomsInCellCoordinates(2,j) - 1,m))) THEN
							
!~ 							IF ((iz == atomsInCellCoordinates(3,j)) .OR. &
!~ 								(iz == MOD(atomsInCellCoordinates(3,j) + 1,m)).OR. &
!~ 								(iz == MOD(atomsInCellCoordinates(3,j) - 1,m))) THEN
!~ 								near = .TRUE.
!~ 							END IF
!~ 						END IF
!~ 					END IF
					
!~ 					! In una cella lontana
!~ 					IF (.NOT. near) THEN
!~ 						CYCLE
					
!~ 					! In una cella vicina
!~ 					ELSE
!~ 						CALL calcForceTwoParticles(tempVec(:),x(1:3,i),x(1:3,j))
!~ 						force(:,i) = force(:,i) - tempVec
!~ 						force(:,j) = force(:,j) + tempVec
!~ 					END IF
					
!~ 				END DO
			
!~ 			END DO
			
!~ 		END IF
		
!~ 	END SUBROUTINE calcForce_BK
	
	SUBROUTINE calcForce(x)
		
		USE list_mod
		
		IMPLICIT NONE
		
		DOUBLE PRECISION, INTENT(IN)	:: x(1:6,1:nAtoms)
		DOUBLE PRECISION	:: tempVec(1:3)

		INTEGER	:: i, j, tx, ty, tz, atomI, atomJ
		! Array che contiene l'indici degli atomi che interagiscono
		INTEGER	:: interactAtomsIndex(1:nAtoms)
		! Intero che indica quanto e' pieno interactAtomsIndex
		INTEGER	:: interactAtomsIndexFillingThisCell
		INTEGER	:: interactAtomsIndexFillingAllCells

		
		! Questo metodo è più veloce dell'altro con BK
		
		! Calcola la forza: decide quale forze vanno valutate e valuta quelle
		force = 0
		potentialEnergy = 0
		
		IF (m == 1) THEN
			DO i=1,nAtoms,+1
				DO j=i+1,nAtoms,+1
					
					CALL calcForceTwoParticles(tempVec(:),x(1:3,i),x(1:3,j))
					force(:,i) = force(:,i) - tempVec
					force(:,j) = force(:,j) + tempVec
					
				END DO
			END DO
		
		ELSE
			
			! Calcola la forza: con la struttura a celle
			DO tx=1,m
				DO ty=1,m
					DO tz=1,m
						
						interactAtomsIndex = 0
						interactAtomsIndexFillingThisCell = 0
						interactAtomsIndexFillingAllCells = 0
						
						! Colleziona atomi nella stessa cella (Blocco 0)
						CALL listCountAndCollect(interactAtomsIndexFillingThisCell,interactAtomsIndex,tx,ty,tz)
						
						! Se in questa cella non ci sono atomi, skippa
						IF (interactAtomsIndexFillingThisCell == 0) THEN
							CYCLE
						END IF
						
						! Colleziona gli atomi di tutte le celle
						interactAtomsIndexFillingAllCells = interactAtomsIndexFillingThisCell
						CALL fillInteractingArray(interactAtomsIndexFillingAllCells, interactAtomsIndex(:),tx,ty,tz)
						
						! Fai interagire tra di loro gli atomi nella stessa cella
						DO i=1,interactAtomsIndexFillingThisCell,+1
							DO j=i+1,interactAtomsIndexFillingThisCell,+1
							
								atomI = interactAtomsIndex(i)
								atomJ = interactAtomsIndex(j)
														
								CALL calcForceTwoParticles(tempVec(:), x(1:3,atomI), x(1:3,atomJ))
								force(:,atomI) = force(:,atomI) - tempVec
								force(:,atomJ) = force(:,atomJ) + tempVec
								
							END DO
						END DO
						
						! Fai interagire tra di loro gli atomi della cella (tx,ty,tz) con gli atomi delle celle adiacenti
						DO i=1,interactAtomsIndexFillingThisCell,+1
							DO j=interactAtomsIndexFillingThisCell+1,interactAtomsIndexFillingAllCells,+1
								
								atomI = interactAtomsIndex(i)
								atomJ = interactAtomsIndex(j)
								
								CALL calcForceTwoParticles(tempVec(:), x(1:3,atomI), x(1:3,atomJ))
								force(:,atomI) = force(:,atomI) - tempVec
								force(:,atomJ) = force(:,atomJ) + tempVec
								
							END DO
						END DO
						
						
					END DO
				END DO
			END DO
			
		END IF
		
	END SUBROUTINE calcForce
	
	SUBROUTINE calcKineticEnergy(atoms)
		
		IMPLICIT NONE
		
		DOUBLE PRECISION, INTENT(IN)	:: atoms(1:6,1:nAtoms)
		
		INTEGER	:: i,j
		
		kineticEnergy = 0
		DO i=1,nAtoms
			DO j = 4,6
				kineticEnergy = kineticEnergy + atoms(j,i)**2
			END DO
		END DO
		kineticEnergy = kineticEnergy*mass/2
		
	END SUBROUTINE
	
	SUBROUTINE calcTemperature(atoms)
		
		IMPLICIT NONE
		
		DOUBLE PRECISION, INTENT(IN)	:: atoms(1:6,1:nAtoms)
		
		IF (callingCalcRefresh==0) THEN
			CALL calcKineticEnergy(atoms(:,:))
		END IF
		
		temperature	= 2*kineticEnergy/(3*kb*nAtoms)
		
	END SUBROUTINE calcTemperature
	
	SUBROUTINE calcMomentum(atoms)
		
		IMPLICIT NONE
		
		DOUBLE PRECISION, INTENT(IN)	:: atoms(1:6,1:nAtoms)
		
		INTEGER	:: i
		
		momentumVec(1:3) = 0
		DO i=1,nAtoms
			momentumVec(1:3) = momentumVec(1:3) + atoms(4:6,i)
		END DO
		momentumVec(:) = momentumVec(:)/nAtoms
		momentum = SQRT(momentumVec(1)**2 + momentumVec(2)**2 + momentumVec(3)**2)
		
	END SUBROUTINE calcMomentum
	
	SUBROUTINE calcPressure(atoms)
		
		IMPLICIT NONE
		
		DOUBLE PRECISION, INTENT(IN)	:: atoms(1:6,1:nAtoms)
		
		IF (callingCalcRefresh == 0 .AND. (.NOT. callingMontecarlo)) THEN
			CALL calcTemperature(atoms(:,:))
		END IF
		
		IF (computingForce == 0) THEN
			CALL calcForce(atoms(:,:))
		END IF
		
		pressure = (2*kineticEnergy+tempPressure)/(3*(boxLength**3))
		tempPressure = 0
		
	END SUBROUTINE calcPressure
	
	SUBROUTINE adjustTemp(atoms,actualTemp,targetTemp)
		
		IMPLICIT NONE
		
		DOUBLE PRECISION, INTENT(INOUT)	:: atoms(1:6,1:nAtoms)
		DOUBLE PRECISION, INTENT(IN)		:: actualTemp,targetTemp
		
		INTEGER			:: i
		
		
		atoms(4:6,:) = atoms(4:6,:) * SQRT(targetTemp/actualTemp)
		PRINT*,actualTemp,"=>",targetTemp
		
		DO i=1,nAtoms
			atoms(4:6,i) = atoms(4:6,i) - momentumVec(1:3)
		END DO
		!CALL calcTemperature(atoms)
		!PRINT*,"NewTemp:",temperature
		
	END SUBROUTINE adjustTemp
	
	SUBROUTINE calcRefresh(atoms)
		
		USE parametri
		
		IMPLICIT NONE
		
		DOUBLE PRECISION, INTENT(IN)	:: atoms(1:6,1:nAtoms)
		
		callingCalcRefresh = 1
		
		CALL calcKineticEnergy(atoms(:,:))
		CALL calcTemperature(atoms(:,:))
		CALL calcPressure(atoms(:,:))
		CALL calcMomentum(atoms(:,:))
		
		energy = kineticEnergy + potentialEnergy
		
		callingCalcRefresh = 0
		
	END SUBROUTINE calcRefresh
	
	SUBROUTINE refreshQuantities()
		
		IMPLICIT NONE
		
		energy = 0
		kineticEnergy = 0
		potentialEnergy = 0
		temperature = 0
		pressure = -1
		tempPressure = 0
		momentum = 0
		
	END SUBROUTINE refreshQuantities
	
	SUBROUTINE dynamic(xPrime,x,t)
		
		IMPLICIT NONE

		DOUBLE PRECISION, INTENT(OUT)	:: xPrime(1:6,1:nAtoms)
		DOUBLE PRECISION, INTENT(IN)	:: x(1:6,1:nAtoms)
		DOUBlE PRECISION, INTENT(IN)	:: t
		
		INTEGER	:: i
		
		computingForce = 1
		CALL refreshQuantities()
		CALL calcForce(x(:,:))
		CALL calcRefresh(x(:,:))
		computingForce = 0
		
		DO i=1,3
			xPrime(i,:) = x(i+3,:)
			xPrime(i+3,:) = force(i,:)
		END DO
		

	END SUBROUTINE dynamic

END MODULE physics
