PROGRAM DinamicaMolecolare
	
	USE parametri
	USE list_mod
	USE physics
	
	IMPLICIT NONE
	
	DOUBLE PRECISION	:: atoms(1:6,1:nAtoms)
	DOUBLE PRECISION	:: temp, minTemp, maxTemp, deltaTemp,startTemp

	DOUBLE PRECISION	:: tempVel(1:3,1:nAtoms)

	DOUBLE PRECISION	:: t_fin, step
	DOUBLE PRECISION	:: x0,xf,deltaX,xNow
	
	LOGICAL				:: justOne, ifVerlet, ifRun = .TRUE.,ifFirst, exists, tempLogical
	
	INTEGER	:: i, nCells, saveEach, mcSteps
	
	DOUBLE PRECISION	:: requiredError
	
	DOUBLE PRECISION	:: jumpChain,en,pre,tem,pot,wow
	
	CHARACTER(LEN = 100):: fileName
	CHARACTER(LEN = 100):: nAtomsStr
	CHARACTER(LEN = 100):: boxLengthStr
	CHARACTER(LEN = 100):: requiredTempStr
	
	! Impostazione parametri
	justOne = .FALSE.
	ifFirst = .TRUE.
	
	step = 1e-2
	
	! Per Van Deer Waals e Monte Carlo
	x0 = .9 ! Che con 1000 atomi corrisponde a boxLenght = 100
	xf = -0.2
	deltaX = -0.2
	xNow = x0
	
	boxLength = EXP((LOG(DBLE(nAtoms)) + x0*LOG(10.))/3.)
!~ 	boxLength = 100.
	
	! Per latro che non ricordo
	t_fin = 1000.
	saveEach = 100
	
	! Per la fusione
	minTemp = .01
	maxTemp = 5.
	deltaTemp = .01
	startTemp = .01
	
	CALL initStructure(atoms(:,:), nCells)
	
	WRITE( nAtomsStr, '(i10)' )  nAtoms
	WRITE( boxLengthStr, '(f10.2)' )  boxLength
	WRITE( requiredTempStr, '(f10.2)' )  requiredTemp
	
	DO WHILE (ifRun)
		
		IF (justOne .eqv. .TRUE.) THEN
			ifRun = .FALSE.
		END IF
		
		! 1 = Energy conservation: Verlet VS Runge Kutta
		! 3 = Draw Wan Deer Waals curve at fixed temperature
		! 4 = Fondi il cristallo
		! 5 = Monte Carlo
		! 6 = Salva tutti tutti gli step di integrazione per l'analisi dell'autocorrelazione
		! 7 = Salva le velocità ogni saveEach steps
		! 8 = Fusione con Monte Carlo
		
		nCells = MIN(FLOOR(boxLength/2),maxCells)
		IF (nCells < 4) THEN
			nCells = 1
		END IF
		PRINT*,"Setting cells: ",nCells
		CALL setCells(atoms(:,:),nCells)
		
		IF (typeRun == 1) THEN
		
			ifVerlet = .TRUE.
			
			step = 1e-2
			CALL energyConservation(atoms(:,:),ifVerlet,t_fin,step,saveEach)
			
			step = 5e-3
			CALL energyConservation(atoms(:,:),ifVerlet,t_fin,step,saveEach)
			
			step = 1e-3
			CALL energyConservation(atoms(:,:),ifVerlet,t_fin,step,saveEach)
			
			!step = 5e-4
			!CALL energyConservation(atoms(:,:),ifVerlet,t_fin,step,saveEach)
			
			ifRun = .FALSE.
			
		ELSE IF (typeRun == 3) THEN
			
			requiredError = 1e-2
			
			IF (xNow >= xf) THEN
			
				fileName = "dati/gasReale/" // TRIM(nAtomsStr) // "-" // TRIM(requiredTempStr) // ".csv"
				CALL onePointMeasure(atoms(:,:),step,fileName,requiredError)
				CALL rescaleBox(atoms(:,:),xNow,deltaX)
				
			ELSE
				ifRun = .FALSE.
			END IF
		
		ELSE IF (typeRun == 4) THEN
			fileName = "dati/fusione/changeVolume/" // TRIM(nAtomsStr) // ".csv"
			!fileName = "dati/fusione/small/" // TRIM(nAtomsStr) // ".csv"
			!fileName = "dati/fusione/regular/" // TRIM(nAtomsStr) // ".csv"
			!fileName = "dati/fusione/big/" // TRIM(nAtomsStr) // ".csv"
			!fileName = "dati/fusione/bigger/" // TRIM(nAtomsStr) // ".csv"
			
			
			PRINT*, "Fusione del cristallo"
			requiredError = 0.005
			ifMFP = .TRUE.
			
			! Imposta il volume per avere pressione nulla
			IF (ifFirst) THEN
				CALL adjustTemp(atoms(:,:),requiredTemp,startTemp)
				requiredTemp = startTemp
				
				tempVel = atoms(4:6,:)
				atoms(4:6,:) = 0
				CALL adjustVolumeZeroPressure(atoms(:,:))
				atoms(4:6,:) = tempVel
				
				ifFirst = .FALSE.
				xNow = LOG10(boxLength**3/nAtoms)
				x0 = xNow
				PRINT*,"Volume trovato:", xNow
				
!~ 				WRITE( boxLengthStr, '(f10.2)' )  boxLength
!~ 				fileName = "dati/fusione/meh/" // TRIM(boxLengthStr) // "-" // TRIM(nAtomsStr) // ".csv"
				
			END IF
			
			IF (requiredTemp > 1.5) THEN
				deltaTemp = .2
			ELSE IF (requiredTemp > .7) THEN
				deltaTemp = .05
			ELSE
				deltaTemp = .01
			END IF
			
			IF (requiredTemp <= maxTemp .AND. requiredTemp >= minTemp) THEN
				
				IF (requiredTemp > 3.) THEN
					step = 5e-3
				ELSE
					step = 1e-2
				END IF
				
				IF (.NOT. ifFirst) THEN
					requiredTemp = requiredTemp + deltaTemp
					CALL adjustVolumeZeroPressure(atoms(:,:))
					!CALL initStructure(atoms(:,:), nCells)
				END IF
				
				crystal(:,:) = atoms(1:3,:)
				CALL onePointMeasure(atoms(:,:),step,fileName,requiredError)
				
			ELSE
				ifRun = .FALSE.
			END IF
		
		ELSE IF (typeRun == 5) THEN
			
			IF (xNow >= xf) THEN
				fileName = "dati/Montecarlo/" // TRIM(nAtomsStr) // "-" // TRIM(requiredTempStr) // ".csv"
				PRINT*, "Metodo montecarlo"
				PRINT*, "Numero atomi:", nAtoms
				
				IF (requiredTemp > 1.) THEN
					mcSteps = 500
				ELSE IF (requiredTemp > .7) THEN
					mcSteps = 1000
				ELSE
					mcSteps = 1500
				END IF
				
				IF (xNow > 1.) THEN
					jumpChain = r0
				ELSE
					jumpChain = r0
				END IF
				
				CALL initStructure(atoms(:,:), nCells)
				CALL montecarlo(atoms(:,:),jumpChain,mcSteps,fileName)
				CALL rescaleBox(atoms(:,:),xNow,deltaX)
				
			ELSE
				ifRun = .FALSE.
			END IF
			
		ELSE IF (typeRun == 6) THEN
		
			IF (xNow >= xf) THEN
				
				CALL saveAllSteps(atoms(:,:),step)
				CALL rescaleBox(atoms(:,:),xNow,deltaX)
				
			ELSE
				ifRun = .FALSE.
			END IF
		
		ELSE IF (typeRun == 7) THEN
			
			fileName = "dati/velocita/" // TRIM(nAtomsStr) // "-" // TRIM(requiredTempStr) // ".dat"
			
			CALL saveVelocities(atoms(:,:),step,saveEach,fileName)
			ifRun = .FALSE.
			
		ELSE IF (typeRun == 8) THEN
		
			!fileName = "dati/fusione montecarlo/smaller/" // TRIM(nAtomsStr) // ".csv"
			!fileName = "dati/fusione montecarlo/small/" // TRIM(nAtomsStr) // ".csv"
			!fileName = "dati/fusione montecarlo/regular/" // TRIM(nAtomsStr) // ".csv"
			fileName = "dati/fusione montecarlo/big/" // TRIM(nAtomsStr) // ".csv"
			!fileName = "dati/fusione montecarlo/bigger/" // TRIM(nAtomsStr) // ".csv"
			
			PRINT*, "Fusione del cristallo con Monte Carlo"
			
			! Imposta il volume per avere pressione nulla
			IF (ifFirst) THEN
				
				ifEnergy = .TRUE.
				CALL adjustTemp(atoms(:,:),requiredTemp,startTemp)
				requiredTemp = startTemp
				CALL adjustVolumeZeroPressure(atoms(:,:))
				ifFirst = .FALSE.
				xNow = LOG10(boxLength**3/nAtoms)
				x0 = xNow
				PRINT*,"Volume trovato:", xNow
				PRINT*,"Atomi:", nAtoms
				
			END IF
			
			PRINT*,requiredTemp
			IF (requiredTemp > 1.) THEN
				deltaTemp = -.2
				mcSteps = 500
				jumpChain = r0
			ELSE IF (requiredTemp > .7) THEN
				deltaTemp = -.05
				mcSteps = 1000
				jumpChain = r0/10.
			ELSE
				deltaTemp = -.01
				mcSteps = 1500
				jumpChain = r0/20.
			END IF
			
			IF (requiredTemp <= maxTemp .AND. requiredTemp >= minTemp) THEN
				
				CALL montecarlo(atoms(:,:),jumpChain,mcSteps,fileName)
				
				CALL adjustTemp(atoms(:,:),requiredTemp,requiredTemp + deltaTemp)
				requiredTemp = requiredTemp + deltaTemp
				CALL calcKineticEnergy(atoms)
				
				!CALL initStructure(atoms(:,:), nCells)
			ELSE
				ifRun = .FALSE.
			END IF
			
		END IF
		
	END DO

END PROGRAM DinamicaMolecolare

SUBROUTINE adjustVolumeZeroPressure(atoms)
	
	USE parametri
	USE physics
	
	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(INOUT)	:: atoms(1:6,1:nAtoms)
	
	DOUBLE PRECISION	:: deltaX,xNow,temp
	DOUBLE PRECISION	:: step=1e-3,nSteps=400
	DOUBLE PRECISION	:: press,en,tmp
	
	INTEGER				:: nCells,pressSign
	INTEGER,PARAMETER	:: mode=2
	LOGICAL				:: ifFlip
	
	ifFlip = .FALSE.
	
	deltaX = -0.001
	xNow = LOG10(boxLength**3/nAtoms)
	
	! Imposta le condizioni iniziali
	
	CALL calcPressure(atoms(:,:))
	
	IF (mode==1) THEN
		! Azzera la pressione
		DO WHILE (pressure .ne. 0)
			
			CALL rescaleBox(atoms(:,:),xNow,deltaX)		
			CALL calcPressure(atoms(:,:))
			
		END DO
	
	ELSE IF (mode==2) THEN
		PRINT*,"Mode 2"
		press=pressure
		! Ferma quando la pressione cambia segno, ovvero il prodotto è negativo
		DO WHILE (press*pressure>0)
			press = pressure
			
			CALL rescaleBox(atoms(:,:),xNow,deltaX)
			CALL calcPressure(atoms(:,:))
			
			PRINT*,pressure, xNow,boxLength,r0
			
			IF ((ifFlip .eqv. .FALSE.) .AND. pressure<0. .AND. pressure>press) THEN
				deltaX = -deltaX
				ifFlip = .TRUE.
			ELSE IF (xNow<0.5) THEN
				deltaX = +0.001
			END IF
			
			! Se la pressione è molto vicina allo zero senza , senza cambiare segno, si ha un gas
			IF (pressure < 1e-2) THEN
				EXIT
			END IF
			
		END DO
		
		! Il volume con pressione minima era il precedente, reimposta
		IF (ifFlip .eqv. .FALSE.) THEN
			deltaX = -deltaX
			CALL rescaleBox(atoms(:,:),xNow,deltaX)
		END IF
	END IF
	! Mode 3 is to use x0
	
	CALL calcRefresh(atoms(:,:))
	
	nCells = MIN(FLOOR(boxLength/2),maxCells)
	IF (nCells < 4) THEN
		nCells = 1
	END IF
	PRINT*,"Setting cells: ",nCells
	CALL setCells(atoms(:,:),nCells)
	
	PRINT*,pressure, xNow,boxLength,r0
	
END SUBROUTINE adjustVolumeZeroPressure

SUBROUTINE rescaleBox(atoms,xNow,deltaX)
	
	USE parametri
	
	IMPLICIT NONE
	
	DOUBLE PRECISION,INTENT(INOUT)	:: atoms(1:6,1:nAtoms)
	DOUBLE PRECISION,INTENT(INOUT)	:: xNow
	DOUBLE PRECISION,INTENT(IN)		:: deltaX
	DOUBLE PRECISION				:: temp
	
	temp = EXP(LOG(10.)*deltaX/3.)
	boxLength = boxLength*temp
	xNow = xNow + deltaX
	
	atoms(1:3,:) = atoms(1:3,:)*temp
	r0 = r0*temp
	
END SUBROUTINE rescaleBox
