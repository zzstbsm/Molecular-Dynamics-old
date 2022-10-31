SUBROUTINE initStructure(atoms, nCells)
	
	USE parametri
	USE physics
	!USE list_mod

	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(INOUT)	:: atoms(1:6,1:nAtoms)
	INTEGER, INTENT(IN)	:: nCells
	
	! For velocity
	DOUBLE PRECISION	:: v, phi, theta
	DOUBLE PRECISION	:: totVel(1:3)
	
	INTEGER	:: seed

	INTEGER	:: i
	
	!CALL simpleCubic(atoms(:,:))
	CALL FCC(atoms(:,:))
	
	! Genera la parte delle velocità
	DO i=1,nAtoms
		
		CALL RANDOM_NUMBER(v)
		CALL RANDOM_NUMBER(phi)
		CALL RANDOM_NUMBER(theta)
		
		phi = phi*2*PI
		theta = theta*PI
		
		atoms(4,i) = v*SIN(theta)*COS(phi)
		atoms(5,i) = v*SIN(theta)*SIN(phi)
		atoms(6,i) = v*COS(theta)

	END DO
	
	IF (nAtoms == 2) THEN
		atoms(:,1) = (/2., 0., 0., -1., 0., 0./)
		atoms(:,2) = (/8. , 0., 0., 1., 0., 0./)
		m = 1
		flagOneCell = .TRUE.
	END IF
	
	
	IF (flagOneCell .eqv. .FALSE.) THEN
		! Inserisci nella struttura a celle
		CALL setCells(atoms(:,:),nCells)
	END IF
	
	! Ci si porta nel sistema del centro di massa
	totVel = 0
	DO i=1,nAtoms
		
		totVel = totVel + atoms(4:6,i)
		
	END DO
	totVel = totVel/nAtoms
	DO i=1,nAtoms
		atoms(4:6,i) = atoms(4:6,i) - totVel
	END DO
	
	! Si riscala la temperatura per ottenere quella richiesta
	CALL calcTemperature(atoms(:,:))
	atoms(4:6,:) = atoms(4:6,:) * SQRT(requiredTemp/temperature)
	
	PRINT*, "Inizializzato"

END SUBROUTINE initStructure

SUBROUTINE periodicConditions(positions)

	USE parametri

	IMPLICIT NONE

	DOUBLE PRECISION, INTENT(INOUT)	:: positions(1:3,1:nAtoms)

	INTEGER				:: i, j, n

	DO i=1,nAtoms
		DO j=1,3
			n = FLOOR(positions(j,i)/boxLength)
			positions(j,i) = positions(j,i) - n*boxLength
		END DO
	END DO
	
	IF (m .ne. 1 ) THEN
		CALL orderAtomsInCell(positions(:,:))
	END IF

END SUBROUTINE periodicConditions

SUBROUTINE simpleCubic(atoms)
	
	USE parametri
	
	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(INOUT)	:: atoms(1:6,1:nAtoms)
	INTEGER	:: i, ix, iy, iz
	
	INTEGER	:: edge
	
	edge = CEILING(nAtoms**(1.0/3))
	
	DO i=1,nAtoms
		
		ix = i
		IF (ix > edge) THEN
			iy = FLOOR(REAL(ix)/edge)
			ix = ix - edge*FLOOR(REAL(ix)/edge)
			IF (iy > edge) THEN					
				iz = FLOOR(REAL(iy)/edge)
				iy = iy - edge*FLOOR(REAL(iy)/edge)
			ELSE
				iz = 0
			END IF
		ELSE
			iy = 0
			iz = 0
		END IF
		
		atoms(1,i) = ix*(boxLength*0.9)/edge
		atoms(2,i) = iy*(boxLength*0.9)/edge
		atoms(3,i) = iz*(boxLength*0.9)/edge
	
	END DO
	
	! Imposta il passo del reticolo
	r0 = (boxLength*0.9)/edge
	
END SUBROUTINE simpleCubic

SUBROUTINE FCC(atoms)
	
	USE parametri
	
	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(INOUT)	:: atoms(1:6,1:nAtoms)
	INTEGER	:: i, ix, iy, iz, insertedAtoms, check
	
	INTEGER	:: edge
	
	edge = CEILING(nAtoms**(1./3.))*2
	
	i = 0
	insertedAtoms = 0
	
	DO WHILE (insertedAtoms < nAtoms)
		
		i = i+1
		
		ix = i
		IF (ix > edge) THEN
			iy = FLOOR(REAL(ix)/edge)
			ix = ix - edge*FLOOR(REAL(ix)/edge)
			IF (iy > edge) THEN					
				iz = FLOOR(REAL(iy)/edge)
				iy = iy - edge*FLOOR(REAL(iy)/edge)
			ELSE
				iz = 0
			END IF
		ELSE
			iy = 0
			iz = 0
		END IF
		
		check = ix + iy + iz
		IF (MOD(check,2) == 0) THEN
			insertedAtoms = insertedAtoms + 1
			atoms(1,insertedAtoms) = ix*(boxLength*0.9)/edge
			atoms(2,insertedAtoms) = iy*(boxLength*0.9)/edge
			atoms(3,insertedAtoms) = iz*(boxLength*0.9)/edge
		END IF
		
	END DO

	! Imposta il passo del reticolo
	r0 = (boxLength*0.9)/edge*SQRT(2.)
	
	
END SUBROUTINE FCC

SUBROUTINE orderAtomsInCell(atoms)

	USE parametri
	USE list_mod

	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(IN)	:: atoms(1:3,1:nAtoms)
	
	TYPE(listAtomsIndexType), POINTER	:: previous
	TYPE(listAtomsIndexType), POINTER	:: actual
	TYPE(listAtomsIndexType), POINTER	:: next

	DOUBLE PRECISION	:: normalizedPositions(1:3,1:nAtoms)

	INTEGER				:: newX,newY,newZ,oldX,oldY,oldZ
	INTEGER				:: atomIndex
	
	! Normalizza le posizioni su m, in questo modo normalizedPosition ha già l'indice incorporato
	normalizedPositions = atoms(1:3,:)/boxLength
	normalizedPositions = normalizedPositions*m

	DO atomIndex = 1,nAtoms
	
		newX = CEILING(normalizedPositions(1,atomIndex))
		newY = CEILING(normalizedPositions(2,atomIndex))
		newZ = CEILING(normalizedPositions(3,atomIndex))
		
		IF (newX==0) THEN
			newX = 1
		END IF
		IF (newY==0) THEN
			newY = 1
		END IF
		IF (newZ==0) THEN
			newZ = 1
		END IF
		
		oldX = atomsInCellCoordinates(1,atomIndex)
		oldY = atomsInCellCoordinates(2,atomIndex)
		oldZ = atomsInCellCoordinates(3,atomIndex)
		
		! Controlla se l'atomo è ancora nella celletta
		! Se è ancora nella celletta, passa al prossimo atomo
		IF ((newX == oldX) .AND. (newY == oldY) .AND. (newZ == oldZ)) THEN
			CYCLE
		END IF
		
		! Se l'atomo ha cambiato celletta, rimuovi l'atomo e reinseriscilo
		CALL listRmElement(atomIndex)
		CALL listAddHead(atomIndex,newX,newY,newZ)
		
	END DO

END SUBROUTINE orderAtomsInCell

SUBROUTINE setCells(atoms,nCells)
	
	USE parametri
	USE list_mod
	
	IMPLICIT NONE
	
	INTEGER, INTENT(IN)	:: nCells
	DOUBLE PRECISION, INTENT(INOUT)	:: atoms(1:6,1:nAtoms)
	
	DOUBLE PRECISION	:: normalizedPositions(1:3,1:nAtoms)
	INTEGER	:: i, ix, iy, iz
	
	m = nCells
	IF (nCells == 1) THEN
		flagOneCell = .TRUE.
	ELSE
		flagOneCell = .FALSE.
		CALL initList()
		
		! Normalizza le posizioni su m, in questo modo normalizedPosition ha già l'indice incorporato
		normalizedPositions = atoms(1:3,:)/boxLength
		normalizedPositions = normalizedPositions*m
		
		DO i = 1,nAtoms
		
			! Inizializza listAtomsIndex
			listAtomsIndex(i)%atomIndex = i
			
			ix = CEILING(normalizedPositions(1,i))
			iy = CEILING(normalizedPositions(2,i))
			iz = CEILING(normalizedPositions(3,i))
			
			IF (ix==0) THEN
				ix = 1
			END IF
			IF (iy==0) THEN
				iy = 1
			END IF
			IF (iz==0) THEN
				iz = 1
			END IF
			
			! Inserisci di testa
			CALL listAddHead(i,ix,iy,iz)
			
		END DO
	
	END IF

	
END SUBROUTINE setCells
