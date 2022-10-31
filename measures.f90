SUBROUTINE onePointMeasure(atoms,step,fileName,requiredError)
	
	USE parametri
	USE physics
	
	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(INOUT)	:: atoms(1:6,1:nAtoms)
	CHARACTER(LEN=100), INTENT(IN)	:: fileName
	
	DOUBLE PRECISION,INTENT(IN)	:: step
	
	DOUBLE PRECISION, INTENT(IN)	:: requiredError
	
	INTEGER :: kShow=0, kSave=0, showEach=1000, saveEach, nSteps=0
	
	LOGICAL	:: ifSaveRun, firstRun, outBound
	DOUBLE PRECISION	:: errBound, errActual
	
	DOUBLE PRECISION	:: tc, collisionDistance
	
	INTEGER				:: nMeasures, newnMeasures, MIndex, i
	
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE	:: arrayEnergy
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE	:: arrayTemperature
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE	:: arrayPressure
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE	:: temp
	DOUBLE PRECISION							:: pot
	
	DOUBLE PRECISION	:: temp1, temp2, error, minTemp, maxTemp, meanTemp
	
	DOUBLE PRECISION	:: stepHere
	
	DOUBLE PRECISION, EXTERNAL	:: calcMean
	
	! Impostazione parametri
	
	time=0.
	stepHere = step
	
	errBound = 5*requiredError
	
	ifSaveRun = .FALSE.
	firstRun = .TRUE.
	MIndex = 0
	
	! Stima tempo tra due collisioni (la velocità media dell'ordine della radice della temperatura)
	tc = r0/SQRT(2*requiredTemp)
	collisionDistance = (boxLength**2)/(nAtoms**(2./3)*4*PI)
	tc = collisionDistance/SQRT(2*requiredTemp)
	
	IF (tc < 2.) THEN
		tc = 2.
	END IF
	
	! Si effettua una misura ogni tc/step steps
	saveEach = FLOOR(tc/stepHere)
	
	!nMeasures = MAX(2*FLOOR(5*tc),minMeasures)
	nMeasures = minMeasures
	ALLOCATE(arrayTemperature(1:nMeasures))
	ALLOCATE(arrayPressure(1:nMeasures))
	ALLOCATE(arrayEnergy(1:nMeasures))
	
	! Simulazione del sistema
	PRINT *, "Numero particelle:", nAtoms
	PRINT *, "Lato cubo:", boxLength
	PRINT *, "Temperatura preimpostata:", requiredTemp
	PRINT *, "Passo del reticolo:", r0
	PRINT *, "Stima velocità media:", SQRT(2*requiredTemp)
	PRINT *, "Libero cammino medio:", collisionDistance
	PRINT *, "Tempo tra due collisioni:", tc
	PRINT *, "Numero di misure minimo:", nMeasures
	
	
	PRINT*, "Indice ","Tempo ","Energia totale ","Temperatura ","Pressione"
	
	DO WHILE (ifSaveRun .eqv. .FALSE.)
		
		MIndex = MIndex + 1
		CALL evolve(atoms(:,:),stepHere,saveEach,arrayEnergy(MIndex),pot,arrayTemperature(MIndex),arrayPressure(MIndex))
		
		nSteps = nSteps + saveEach
		
		meanTemp = arrayTemperature(MIndex)
		errActual = ABS(requiredTemp-meanTemp)/requiredTemp
		
		IF (errActual>errBound) THEN
			outBound = .TRUE.
		ELSE
			outBound = .FALSE.
		END IF
		
		! Rescale temperature after the first rescaling of the volume
		IF (firstRun .eqv. .TRUE.) THEN
		
			firstRun = .FALSE.
			!PRINT*, nSteps,time,energy,potentialEnergy,temperature,pressure,"Scarta prima misura"
			PRINT*, nSteps,time,arrayEnergy(MIndex),pot,arrayTemperature(MIndex),arrayPressure(MIndex),"Scarta prima misura"
			CALL adjustTemp(atoms(:,:),arrayTemperature(MIndex),requiredTemp)
			MIndex = Mindex - 1				
			crystal(:,:) = atoms(1:3,:)
		
		ELSE IF (outBound .eqv. .TRUE.) THEN
			
			! Scarta solo se non continua a fluttuare troppo, ma le fluttuazioni devono rimanere dentro il 20%
			IF (errActual>.20) THEN
				!PRINT*, nSteps,time,energy,potentialEnergy,temperature,pressure,"Scarta misura fuori dalla zona di accettazione"
				PRINT*, nSteps,time,arrayEnergy(MIndex),pot,arrayTemperature(MIndex),arrayPressure(MIndex),&
							"Scarta misura fuori dalla zona di accettazione"
				CALL adjustTemp(atoms(:,:),arrayTemperature(MIndex),requiredTemp)
								
				crystal(:,:) = atoms(1:3,:)
				
				MIndex = 0
			END IF
		ELSE
			! Stampa su schermo i dati fisici
			PRINT*, nSteps,time,arrayEnergy(MIndex),pot,arrayTemperature(MIndex),arrayPressure(MIndex),"Salvo: misura numero", MIndex
			
		END IF
		
		IF (MIndex == nMeasures) THEN
			
			temp1 = calcMean(arrayTemperature(1:nMeasures/2),nMeasures/2)
			temp2 = calcMean(arrayTemperature(1+nMeasures/2:nMeasures),nMeasures/2)
			
			error = ABS((0.5*(temp1+temp2)-requiredTemp)/requiredTemp)
			
			! La media delle due misure sono in accordo entro l'errore richiesto
			IF (error < requiredError) THEN
				! Se le due misure sono vicine alla temperatura richiesta entro l'errore richiesto, salva
				error = ABS((0.5*(temp1+temp2)-requiredTemp)/requiredTemp)
				IF (error < requiredError) THEN
					ifSaveRun = .TRUE.
				END IF
			END IF
			
			IF (ifSaveRun) THEN
				
				PRINT*, "Salvo il run con errore ", error
				
				CALL saveRun(atoms(:,:),arrayEnergy(1:nMeasures),arrayPressure(1:nMeasures),arrayTemperature(1:nMeasures),fileName,nMeasures)
				
				DEALLOCATE(arrayTemperature)
				DEALLOCATE(arrayPressure)
				DEALLOCATE(arrayEnergy)
				
			ELSE
				
				! Aumenta il numero di misure e butta via la prima metà
				IF (nMeasures < maxMeasures) THEN
					newnMeasures = nMeasures + 2
				ELSE
					newnMeasures = nMeasures
				END IF
				
				ALLOCATE(temp(1:nMeasures/2))
				
				! Rialloca la temperatura
				temp(:) = arrayTemperature(1+nMeasures/2:nMeasures)
				DEALLOCATE(arrayTemperature)
				ALLOCATE(arrayTemperature(1:newnMeasures))
				arrayTemperature(1:nMeasures/2) = temp(:)
				
				! Rialloca la pressione
				temp(:) = arrayPressure(1+nMeasures/2:nMeasures)
				DEALLOCATE(arrayPressure)
				ALLOCATE(arrayPressure(1:newnMeasures))
				arrayPressure(1:nMeasures/2) = temp(:)
				
				! Rialloca l'energia
				temp(:) = arrayEnergy(1+nMeasures/2:nMeasures)
				DEALLOCATE(arrayEnergy)
				ALLOCATE(arrayEnergy(1:newnMeasures))
				arrayEnergy(1:nMeasures/2) = temp(:)
				
				! Sistema le cose
				MIndex = nMeasures/2
				nMeasures = newnMeasures
				DEALLOCATE(temp)
				
				! temp2 è considerata la temperatura del sistema
				! Riscala le energie cinetiche se temp2 è fuori dalla soglia di accettazione
				minTemp = requiredTemp - requiredError*requiredTemp
				maxTemp = requiredTemp + requiredError*requiredTemp
				IF (temp2 < minTemp .OR. temp2 > maxTemp) THEN
					PRINT*, "Non all'equilibrio entro l'errore, riscalamento temperatura"
					CALL adjustTemp(atoms(:,:),temp2,requiredTemp)
				ELSE
					PRINT*, "Si lascia evolvere senza riscalare"
				END IF
				PRINT*, "Temp1", temp1
				PRINT*, "Temp2", temp2
				PRINT*, "Escursione permessa:"
				PRINT*, minTemp, " - ", maxTemp
				
			END IF
			
		END IF
			
	END DO
	
END SUBROUTINE onePointMeasure

SUBROUTINE evolve(atoms,step,nSteps,measureEnergy,measurePotEnergy,measureTemperature,measurePressure)
	
	USE parametri
	USE physics
	
	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(INOUT)	:: atoms(1:6,1:nAtoms)
	DOUBLE PRECISION, INTENT(IN)	:: step
	INTEGER, INTENT(IN)				:: nSteps
	DOUBLE PRECISION, INTENT(OUT)	:: measureEnergy,measurePotEnergy,measureTemperature,measurePressure
	
	DOUBLE PRECISION	:: arrayEnergy(1:nSteps)
	DOUBLE PRECISION	:: arrayPotEnergy(1:nSteps)
	DOUBLE PRECISION	:: arrayTemperature(1:nSteps)
	DOUBLE PRECISION	:: arrayPressure(1:nSteps)
	
	INTEGER	:: stepNow
	
	DOUBLE PRECISION, EXTERNAL	:: calcMean
	
	stepNow = 0
	
	DO WHILE (stepNow<nSteps)
	
		CALL Verlet(dynamic,atoms(:,:),time,step)
		stepNow = stepNow + 1
		time=time+step
		
		arrayEnergy(stepNow) = energy
		arrayPotEnergy(stepNow) = potentialEnergy
		arrayTemperature(stepNow) = temperature
		arrayPressure(stepNow) = pressure
		
		CALL periodicConditions(atoms(1:3,:))
		
		!PRINT*, stepNow,time,energy,temperature,pressure,"Evolve"
		
	END DO
	
	measureEnergy = calcMean(arrayEnergy(:),nSteps)
	measurePotEnergy = calcMean(arrayPotEnergy(:),nSteps)	
	measureTemperature = calcMean(arrayTemperature(:),nSteps)
	measurePressure = calcMean(arrayPressure(:),nSteps)	
	
END SUBROUTINE

DOUBLE PRECISION FUNCTION meanFreePath(positionFin,positionIn)
	
	USE parametri
	USE physics
	
	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(IN)	:: positionFin(1:3,1:nAtoms), positionIn(1:3,1:nAtoms)
	DOUBLE PRECISION	:: temp
	INTEGER	:: i
	
	DOUBLE PRECISION	:: distance, direction(1:3)
	
	
	temp = 0
	distance = 0
	direction = 0
	DO i=1,nAtoms
		CALL calcDistance(distance,direction(1:3),positionFin(1:3,i),positionIn(1:3,i))
		temp = temp + distance**2
		
	END DO
	
	temp = temp/nAtoms
	
	meanFreePath = SQRT(temp)
			
END FUNCTION meanFreePath

SUBROUTINE radialDistribution(atoms)
	
	USE parametri
	USE physics
	
	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(IN)	:: atoms(1:6,1:natoms)
	DOUBLE PRECISION				:: distances(1:nAtoms*(nAtoms-1))
	
	DOUBLE PRECISION	:: distance, direction(1:3)
	
	CHARACTER(LEN=100)			:: nAtomsStr, radial,requiredTempStr
	
	INTEGER	:: i,j,k
	
	WRITE( nAtomsStr, '(i10)' )  nAtoms
	WRITE( requiredTempStr, '(f10.2)' )  requiredTemp
	
	radial = "dati/Radial distribution/chVolume/" // TRIM(nAtomsStr) // "-" // TRIM(requiredTempStr) // ".dat"
	
	! Salva tutte le distanze tra gli atomi
	k = 0
	DO i=1,nAtoms
		
		DO j=1,nAtoms
			
			IF (i==j) THEN
				CYCLE
			END IF
			k = k+1
			
			CALL calcDistance(distance,direction(:),atoms(1:3,i),atoms(1:3,j))
			
			distances(k) = distance
			
		END DO
		
	END DO
	
	OPEN(3,file=radial,form="unformatted")
		
		WRITE(3) distances
		
	CLOSE(3)
	
END SUBROUTINE radialDistribution

SUBROUTINE saveRun(atoms,arrayEnergy,arrayPressure,arrayTemperature,fileName,nMeasures)
	
	USE parametri
	USE physics
	
	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(IN)	:: atoms(1:6,1:nAtoms)
	DOUBLE PRECISION, INTENT(INOUT)	:: arrayEnergy(1:nMeasures),arrayPressure(1:nMeasures),arrayTemperature(1:nMeasures)
	CHARACTER(LEN=100), INTENT(IN)	:: fileName
	INTEGER, INTENT(IN)				:: nMeasures
	
	DOUBLE PRECISION	:: toSaveEnergy, toSaveTemperature, toSavePressure
	DOUBLE PRECISION	:: toSaveEnergy_std, toSaveTemperature_std, toSavePressure_std
	
	LOGICAL	:: exists
	
	DOUBLE PRECISION, EXTERNAL	:: calcMean
	DOUBLE PRECISION, EXTERNAL	:: calcStd
	DOUBLE PRECISION, EXTERNAL	:: meanFreePath
	
	! Compute average
	toSaveEnergy = calcMean(arrayEnergy(1:nMeasures),nMeasures)
	toSavePressure = calcMean(arrayPressure(1:nMeasures),nMeasures)
	toSaveTemperature = calcMean(arrayTemperature(1:nMeasures),nMeasures)
	
	! Compute std
	toSaveEnergy_std = calcStd(toSaveEnergy,arrayEnergy(1:nMeasures),nMeasures)
	toSavePressure_std = calcStd(toSavePressure,arrayPressure(1:nMeasures),nMeasures)
	toSaveTemperature_std = calcStd(toSaveTemperature,arrayTemperature(1:nMeasures),nMeasures)
	
	IF (ifMFP) THEN
		INQUIRE(file=fileName,exist=exists)
		IF (exists) THEN
			OPEN(2, file=fileName, status="old", position="append", action="write")
		ELSE
			OPEN(2, file=fileName, status="new", action="write")
			WRITE(2,*) "time,boxLength,volume,&
					energy,energy_std,&
					temperature,temperature_std,&
					pressure,pressure_std,&
					n_measures,mean_free_path"
		END IF
		
		! Write in file
		WRITE(2,*) time, ",", boxLength, ",", boxLength**3, ",", &
					toSaveEnergy, ",",toSaveEnergy_std, ",", &
					toSaveTemperature, ",", toSaveTemperature_std, ",", &
					toSavePressure, ",", toSavePressure_std, ",", nMeasures, ",", &
					meanFreePath(atoms(1:3,:),crystal(:,:))
		
		CALL radialDistribution(atoms(:,:))
		
	ELSE
		INQUIRE(file=fileName,exist=exists)
		IF (exists) THEN
			OPEN(2, file=fileName, status="old", position="append", action="write")
		ELSE
			OPEN(2, file=fileName, status="new", action="write")
			WRITE(2,*) "time,boxLength,volume,&
					energy,energy_std,&
					temperature,temperature_std,&
					pressure,pressure_std,&
					n_measures"
		END IF
		! Write in file
		WRITE(2,*) time, ",", boxLength, ",", boxLength**3, ",", &
					toSaveEnergy, ",",toSaveEnergy_std, ",", &
					toSaveTemperature, ",", toSaveTemperature_std, ",", &
					toSavePressure, ",", toSavePressure_std, ",", nMeasures
	END IF
	CLOSE(2)
	
	PRINT*, "energy, energy_std"
	PRINT*, toSaveEnergy, ",",toSaveEnergy_std
	
	PRINT*, "temperature, temperature_std"
	PRINT*, toSaveTemperature, ",", toSaveTemperature_std
	
	PRINT*, "pressure, pressure_std"
	PRINT*, toSavePressure, ",", toSavePressure_std

	
END SUBROUTINE saveRun

SUBROUTINE saveBoltzmann(fileName,atoms,preamble)
	
	USE parametri
	
	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(IN)	:: atoms(1:6,1:nAtoms)
	LOGICAL, INTENT(IN)	:: preamble
	
	INTEGER				:: i
	
	INTEGER, PARAMETER	:: printArraySize = 50
	DOUBLE PRECISION	:: velMax = 10
	INTEGER	:: boltzmann(1:printArraySize)
	INTEGER				:: vel
	
	CHARACTER(LEN = 500), INTENT(IN)	:: fileName
	
	! Il preambolo va da 0 a 50 e va normalizzato su 10 su python
	IF (preamble) THEN
		OPEN(1,file=fileName)
		
		WRITE(1,"(A,i2,A2)",ADVANCE="no") '"', 0, '",'
		DO i=1,printArraySize-1
			WRITE(1,"(A,i2,A3,i2,A2)",ADVANCE="no") '"', i, '","', i,'",'
		END DO
		WRITE(1,*) '"',printArraySize,'"'
		
		CLOSE(1)
	ELSE
		OPEN(1,file=fileName, status="old", position="append", action="write")
		
		boltzmann = 0
		! Riempie l'array che conta le velocità
		DO i=1,nAtoms
			
			! Calcola in quale colonna dell'istogramma deve essere posto
			vel = CEILING(SQRT(atoms(4,i)**2+atoms(5,i)**2+atoms(6,i)**2)*printArraySize/velMax)
			
			! Inserisci colonna dell'istogramma
			IF ((vel <= printArraySize) .AND. (vel>=0)) THEN
				IF (vel == 0) THEN
					boltzmann(1) = boltzmann(1) + 1
				ELSE
					boltzmann(vel) = boltzmann(vel) + 1
				END IF				
			ELSE
				boltzmann(printArraySize) = boltzmann(printArraySize) + 1
			END IF
		
		END DO
		
		! Stampa nel file
		WRITE(1,"(i3,A)",ADVANCE="no") boltzmann(1),","
		DO i=1,printArraySize-1
			WRITE(1,"(i3,A,i3,A)",ADVANCE="no") boltzmann(i),",",boltzmann(i+1),","
		END DO
		WRITE(1,"(i3)") boltzmann(printArraySize)
		
		CLOSE(1)
	END IF
		
	
END SUBROUTINE saveBoltzmann

SUBROUTINE singleMeasure(atoms)
	
	! Do just the first measure
	
	USE parametri
	USE physics
	
	IMPLICIT NONE
	
	DOUBLE PRECISION	:: atoms(1:6,1:nAtoms)
	
	DOUBLE PRECISION	:: step = 1e-3
	
	CHARACTER(LEN = 100):: fileName
	CHARACTER(LEN = 100):: nAtomsStr
	CHARACTER(LEN = 100):: boxLengthStr
	CHARACTER(LEN = 100):: requiredTempStr
	LOGICAL	:: exists
	
	
	WRITE( nAtomsStr, '(i10)' )  nAtoms
	WRITE( boxLengthStr, '(f10.2)' )  boxLength
	WRITE( requiredTempStr, '(f10.2)' )  requiredTemp
	
	!CALL Verlet(dynamic,atoms(:,:),time,step)
	CALL calcRefresh(atoms(:,:))
	
	fileName = "dati/single_measure/" // TRIM(nAtomsStr) // "-" // TRIM(requiredTempStr) // ".csv"
	
	INQUIRE(file=fileName,exist=exists)
	IF (exists) THEN
		OPEN(2, file=fileName, status="old", position="append", action="write")
	ELSE
		OPEN(2, file=fileName, status="new", action="write")
		WRITE(2,*) "time,boxLength,volume,&
				energy,energy_std,&
				temperature,temperature_std,&
				pressure,pressure_std,&
				n_measures"
	END IF
	
	! Write in file
	WRITE(2,*) time, ",", boxLength, ",", boxLength**3, ",", &
				energy, ",",0, ",", &
				temperature, ",", 0, ",", &
				pressure, ",", 0, ",", 1	
	CLOSE(2)
	
	PRINT*, "Lato ","Tempo ","Energia totale ","Energia potenziale ","Temperatura ","Pressione"
	PRINT*, boxLength,time,energy,potentialEnergy,temperature,pressure
	
	
END SUBROUTINE singleMeasure

SUBROUTINE montecarlo(atoms,jumpInput,nCycle,fileName)
	
	USE parametri
	USE physics
	
	IMPLICIT NONE
	
	DOUBLE PRECISION,INTENT(IN)		:: jumpInput
	INTEGER, INTENT(IN)				:: nCycle
	CHARACTER(LEN = 100),INTENT(IN)	:: fileName
	
	INTEGER				:: accepted, acceptedNow
	
	DOUBLE PRECISION,INTENT(INOUT)	:: atoms(1:6,1:nAtoms)
	
	DOUBLE PRECISION				:: oldPosition(1:3)
	DOUBLE PRECISION				:: jump
	
	DOUBLE PRECISION	:: beta,probability,prob,savePressure,displacement(1:3)
	DOUBLE PRECISION	:: saveEnergy,saveEnergyStd,savePressureStd
	DOUBLE PRECISION	:: Ui_new,Ui_old,oldPressure
	DOUBLE PRECISION	:: savePressure1,savePressure2,error
	
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE	:: pressureArray, chunkPress
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE	:: energyArray, chunkEn
	INTEGER										:: pressureArraySize, pressureArrayIndex
	
	LOGICAL	:: exists,measureFlag,full
	
	INTEGER	:: i,j,show, k
	
	INTEGER	:: showEach = 10,n
	
	DOUBLE PRECISION, EXTERNAL	:: calcMean
	DOUBLE PRECISION, EXTERNAL	:: calcStd
	
	beta=1/(kb*requiredTemp)
	accepted = 0
	acceptedNow = 0
	jump = jumpInput
	
	CALL calcRefresh(atoms(:,:))
	
	pressureArraySize = nCycle
	ALLOCATE(pressureArray(1:pressureArraySize))
	ALLOCATE(energyArray(1:pressureArraySize))
	ALLOCATE(chunkPress(1:nAtoms))
	ALLOCATE(chunkEn(1:nAtoms))
	
	pressureArray = 0
	
	IF (ifEnergy) THEN
		PRINT*, "Tentativo ", "Lato ", "Energia ","Temperatura ","Pressione ","Accetati ","Frazione accettata"
	ELSE
		PRINT*, "Tentativo ", "Lato ", "Temperatura ","Pressione ","Accetati ","Frazione accettata"
	END IF
	
	n = 0
	show = 0
	pressureArrayIndex = 0
	
	callingMontecarlo = .TRUE.
	measureFlag=.FALSE.
	full=.FALSE.
	
	DO WHILE (.NOT. measureFlag)
		
		
		n = n+1
		show = show + 1
		pressureArrayIndex = MOD(pressureArrayIndex,nCycle)
		pressureArrayIndex = pressureArrayIndex + 1
		
		DO i=1,nAtoms
			oldPosition = atoms(1:3,i)
			
			! Calcola l'energia potenziale dell'i-esimo atomo
			Ui_old = potentialEnergy
			oldPressure = pressure
			
			! Sposta l'atomo
			CALL RANDOM_NUMBER(displacement)
			atoms(1:3,i) = atoms(1:3,i) + (displacement-.5)*jump
			
			! Condizioni al contorno periodiche
			CALL periodicConditions(atoms(1:3,:))
			
			! Calcola la nuova energia potenziale
			CALL calcPressure(atoms(:,:))
			Ui_new = potentialEnergy
			
			probability = EXP(-beta*(Ui_new-Ui_old))
			
			CALL RANDOM_NUMBER(prob)
			
			! Accetta il trial
			IF (prob < probability) THEN
				acceptedNow = acceptedNow + 1
				accepted = accepted + 1
			ELSE
				atoms(1:3,i) = oldPosition(:)
				CALL periodicConditions(atoms(1:3,:))
				potentialEnergy = Ui_old
				pressure = oldPressure
			END IF
			
			energy = kineticEnergy + potentialEnergy
			
			chunkPress(i) = pressure
			chunkEn(i) = energy
			
		END DO
		
		pressureArray(pressureArrayIndex) = calcMean(chunkPress(:),nAtoms)
		IF (ifEnergy) THEN
			energyArray(pressureArrayIndex) = calcMean(chunkEn(:),nAtoms)
		END IF
		
		IF (show == showEach) THEN
			
			show = 0
			
			IF (ifEnergy) THEN
				PRINT*, n,boxLength,energyArray(pressureArrayIndex),requiredTemp,pressureArray(pressureArrayIndex),&
							accepted,(accepted/DBLE(n*nAtoms)),"Misura:",pressureArrayIndex
			ELSE
				PRINT*, n,boxLength,requiredTemp,pressureArray(pressureArrayIndex),&
							accepted,(accepted/DBLE(n*nAtoms)),"Misura:",pressureArrayIndex
			END IF
			
			IF (MOD(n,nCycle) == 0 .AND. (acceptedNow/DBLE(showEach*nAtoms)) < 0.1) THEN
				PRINT*,"Dimezzamento dello step: ",jump,"=>",jump/2.
				jump = jump/2.
			END IF
			
			acceptedNow = 0
			
		END IF
		
		! Si possono cominciare a fare le medie delle misure quando l'array delle pressioni si riempie
		IF (n == nCycle) THEN
			full = .TRUE.
		END IF
		
		! Fai le medie
		IF ((full) .AND. MOD(pressureArrayIndex,50) == 0) THEN
			! Prendi le nCycle/2 misure più vecchie e mediale
			savePressure1 = 0
			savePressure2 = 0
			DO i=0,nCycle/2-1
				savePressure1 = savePressure1 + pressureArray(MOD(i+pressureArrayIndex,nCycle))
				savePressure2 = savePressure2 + pressureArray(MOD(i+nCycle/2+pressureArrayIndex,nCycle))
			END DO
			savePressure1 = savePressure1/(nCycle/2)
			savePressure2 = savePressure2/(nCycle/2)
			
			! Verifica che siano uguali entro l'1%
			savePressure = (savePressure1+savePressure2)/2
			error = ABS((savePressure1-savePressure2)/MIN(ABS(savePressure1),ABS(savePressure2)))
			PRINT*,"Misura:", savePressure1,savePressure2,savePressure,error
			
			! Entro l'errore, salva la misura
			IF (error<.01) THEN
				measureFlag = .TRUE.
			END IF
			
		END IF
	END DO
	
	IF (ifEnergy) THEN
		saveEnergy = calcMean(energyArray,nCycle)
		saveEnergyStd = calcStd(saveEnergy,energyArray,nCycle)
	END IF
	savePressureStd = calcStd(savePressure,pressureArray,nCycle)
	
	DEALLOCATE(pressureArray)
	DEALLOCATE(chunkPress)
	DEALLOCATE(energyArray)
	DEALLOCATE(chunkEn)
	
	PRINT*, "Salva "
	PRINT*, n,boxLength,saveEnergy,requiredTemp,savePressure
	
	IF (ifEnergy) THEN
		! Salva i dati	
		INQUIRE(file=fileName,exist=exists)
		IF (exists) THEN
			OPEN(2, file=fileName, status="old", position="append", action="write")
		ELSE
			OPEN(2, file=fileName, status="new", action="write")
			WRITE(2,*) "n,boxLength,energy,energy_std,temperature,pressure,pressure_std,accepted"
		END IF
		
		WRITE(2,*) n*nAtoms, ",", boxLength, ",", saveEnergy, ",", saveEnergyStd, ",", requiredTemp, ",", &
				savePressure, ",", savePressureStd, ",", accepted
		
		CLOSE(2)
	ELSE
		! Salva i dati	
		INQUIRE(file=fileName,exist=exists)
		IF (exists) THEN
			OPEN(2, file=fileName, status="old", position="append", action="write")
		ELSE
			OPEN(2, file=fileName, status="new", action="write")
			WRITE(2,*) "n,boxLength,temperature,pressure,pressure_std,accepted"
		END IF
		
		WRITE(2,*) n*nAtoms, ",", boxLength, ",", requiredTemp, ",", &
				savePressure, ",", savePressureStd, ",", accepted
		
		CLOSE(2)
	END IF
	
END SUBROUTINE montecarlo
SUBROUTINE energyConservation(atoms,ifVerlet,t_fin,step, saveEach)
	
	! Per poter studiare la conservazione dell'energia confrontando verlet e runge kutta
	
	USE parametri
	USE physics
	
	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(INOUT)	:: atoms(1:6,1:nAtoms)
	LOGICAL, INTENT(IN)				:: ifVerlet ! if true, use verlet, if false use RK instead
	DOUBLE PRECISION, INTENT(IN)	:: t_fin, step
	INTEGER, INTENT(IN)				:: saveEach
	
	INTEGER :: nSteps=0
	
	DOUBLE PRECISION	:: t_in = 0, tc
	DOUBLE PRECISION	:: collisionDistance
	
	INTEGER				:: showEach=1000, kShow=0, kSave=0
	
	CHARACTER(LEN = 1000)	:: fileName
	CHARACTER(LEN = 100)	:: nAtomsStr
	CHARACTER(LEN = 100)	:: boxLengthStr
	CHARACTER(LEN = 100)	:: requiredTempStr
	
	! Impostazione parametri
	! Stima tempo tra due collisioni (la velocità media dell'ordine della radice della temperatura)
	collisionDistance = (boxLength**2)/(nAtoms**(2./3)*4*PI)
	tc = collisionDistance/SQRT(2*requiredTemp)
	
	! Simulazione del sistema
	time = t_in
	
	IF (ifVerlet) THEN
		PRINT *, "Using Verlet"
	ELSE
		PRINT *, "Using RK"
	END IF
	PRINT *, "Numero particelle:", nAtoms
	PRINT *, "Lato cubo:", boxLength
	PRINT *, "Temperatura preimpostata:", requiredTemp
	PRINT *, "Passo del reticolo:", r0
	PRINT *, "Stima velocità media:", SQRT(2*requiredTemp)
	PRINT *, "Tempo tra due collisioni:", tc
	PRINT *, "Tempo finale:", t_fin
		
	PRINT*, "Indice ","Tempo ","Energia totale ","Energia potenziale ","Temperatura ","Pressione"
	
	WRITE( nAtomsStr, '(i10)' )  nAtoms
	WRITE( boxLengthStr, '(f10.2)' )  boxLength
	WRITE( requiredTempStr, '(f10.2)' )  requiredTemp
	
			
	IF (ifVerlet) THEN
		IF (step == 1e-2) THEN
			fileName = "dati/energia/Verlet/step_1e-2/" &
				// TRIM(nAtomsStr) // "-" // TRIM(requiredTempStr) // "-" // TRIM(boxLengthStr) // ".csv"
		ELSE IF (step == 5e-3) THEN
			fileName = "dati/energia/Verlet/step_5e-3/" &
				// TRIM(nAtomsStr) // "-" // TRIM(requiredTempStr) // "-" // TRIM(boxLengthStr) // ".csv"
		ELSE IF (step == 1e-3) THEN
			fileName = "dati/energia/Verlet/step_1e-3/" &
				// TRIM(nAtomsStr) // "-" // TRIM(requiredTempStr) // "-" // TRIM(boxLengthStr) // ".csv"
		ELSE IF (step == 5e-4) THEN
			fileName = "dati/energia/Verlet/step_5e-4/" &
				// TRIM(nAtomsStr) // "-" // TRIM(requiredTempStr) // "-" // TRIM(boxLengthStr) // ".csv"
		END IF
	ELSE
		IF (step == 1e-2) THEN
			fileName = "dati/energia/RK/step_1e-2/" &
				// TRIM(nAtomsStr) // "-" // TRIM(requiredTempStr) // "-" // TRIM(boxLengthStr) // ".csv"
		ELSE IF (step == 5e-3) THEN
			fileName = "dati/energia/RK/step_5e-3/" &
				// TRIM(nAtomsStr) // "-" // TRIM(requiredTempStr) // "-" // TRIM(boxLengthStr) // ".csv"
		ELSE IF (step == 1e-3) THEN
			fileName = "dati/energia/RK/step_1e-3/" &
				// TRIM(nAtomsStr) // "-" // TRIM(requiredTempStr) // "-" // TRIM(boxLengthStr) // ".csv"
		ELSE IF (step == 5e-4) THEN
			fileName = "dati/energia/RK/step_5e-4/" &
				// TRIM(nAtomsStr) // "-" // TRIM(requiredTempStr) // "-" // TRIM(boxLengthStr) // ".csv"
		END IF
	END IF

	OPEN(2, file=fileName)
	WRITE(2,*) "index,time,boxLength,volume,energy,temperature,pressure"
	
	DO WHILE (time<t_fin)
		
		IF (ifVerlet) THEN
			CALL Verlet(dynamic,atoms(:,:),time,step)
		ELSE
			CALL RungeKutta(dynamic,atoms(:,:),time,step)
		END IF
		
		time=time+step
		nSteps = nSteps + 1
		kShow = kShow+1
		kSave = kSave+1
		
		CALL periodicConditions(atoms(1:3,:))
		
		IF (kSave == saveEach) THEN
			WRITE(2,*) nSteps, ",", time, ",", boxLength, ",", boxLength**3, ",", energy,",", temperature,",",pressure
			kSave = 0
		END IF
		
		IF (kShow == showEach) THEN
			kShow = 0
			
			! Stampa su schermo i dati fisici
			PRINT*, nSteps,time,energy,potentialEnergy,temperature,pressure
		END IF
		
	END DO
	
	CLOSE(2)
	
END SUBROUTINE energyConservation

SUBROUTINE saveVelocities(atoms,step,saveEach,fileName)
	
	! Per poter studiare la conservazione dell'energia confrontando verlet e runge kutta
	
	USE parametri
	USE physics
	
	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(INOUT)	:: atoms(1:6,1:nAtoms), step
	INTEGER, INTENT(IN)				:: saveEach
	CHARACTER(LEN = 100),INTENT(IN)	:: fileName
	
	INTEGER :: nSteps=0
	
	DOUBLE PRECISION	:: velocities(1:nAtoms)
	
	INTEGER				:: showEach=1000, kShow=0, kSave=0
	
	
	! Simulazione del sistema
	
	PRINT *, "Numero particelle:", nAtoms
	PRINT *, "Lato cubo:", boxLength
	PRINT *, "Temperatura preimpostata:", requiredTemp
		
	PRINT*, "Indice ","Tempo ","Energia totale ","Energia potenziale ","Temperatura ","Pressione"
	
	OPEN(2, file=fileName,form="unformatted")
	
	velocities(:) = SQRT(atoms(4,:)**2+atoms(5,:)**2+atoms(6,:)**2)
	!PRINT*,velocities
	WRITE(2) velocities
	
	DO WHILE (nSteps<10000)
		
		CALL Verlet(dynamic,atoms(:,:),time,step)
		
		nSteps = nSteps + 1
		kShow = kShow+1
		kSave = kSave+1
		time = time+step
		
		CALL periodicConditions(atoms(1:3,:))
		
		IF (kSave == saveEach) THEN
			velocities(:) = SQRT(atoms(4,:)**2+atoms(5,:)**2+atoms(6,:)**2)
			!PRINT*,velocities
			WRITE(2) velocities
			kSave = 0
		END IF
		
		IF (kShow == showEach) THEN
			kShow = 0
			
			! Stampa su schermo i dati fisici
			PRINT*, nSteps,time,energy,potentialEnergy,temperature,pressure
		END IF
		
	END DO
	
	CLOSE(2)
	
END SUBROUTINE saveVelocities

SUBROUTINE saveAllSteps(atoms,step)
	
	! Per poter studiare l'autocorrelazione
	
	USE parametri
	USE physics
	
	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(INOUT)	:: atoms(1:6,1:nAtoms)
	DOUBLE PRECISION, INTENT(IN)	:: step
	
	INTEGER :: nSteps=0
	
	DOUBLE PRECISION	:: t_in = 0, t_fin = 100d0
	
	DOUBLE PRECISION	:: tc, collisionDistance
	
	INTEGER				:: showEach=1000, kShow=0, i,j
	
	CHARACTER(LEN = 500):: fileName, velocitiesFilename
	CHARACTER(LEN = 100):: nAtomsStr
	CHARACTER(LEN = 100):: boxLengthStr
	CHARACTER(LEN = 100):: requiredTempStr
	
	! Impostazione parametri
	! Stima tempo tra due collisioni (la velocità media dell'ordine della radice della temperatura)
	collisionDistance = (boxLength**2)/(nAtoms**(2./3)*4*PI)
	tc = collisionDistance/SQRT(2*requiredTemp)
	
	! Simulazione del sistema
	time = t_in
	
	PRINT *, "Salva tutto per autocorrelazione"
	PRINT *, "Numero particelle:", nAtoms
	PRINT *, "Lato cubo:", boxLength
	PRINT *, "Temperatura preimpostata:", requiredTemp
	PRINT *, "Passo del reticolo:", r0
	PRINT *, "Stima velocità media:", SQRT(2*requiredTemp)
	PRINT *, "Tempo tra due collisioni:", tc
	PRINT *, "Tempo finale:", t_fin
		
	PRINT *, "Indice ","Tempo ","Energia totale ","Energia potenziale ","Temperatura ","Pressione ","Momento"
	
	WRITE( nAtomsStr, '(i10)' )  nAtoms
	WRITE( boxLengthStr, '(f10.2)' )  boxLength
	WRITE( requiredTempStr, '(f10.2)' )  requiredTemp
	
	fileName = "dati/autocorrelazione/" // TRIM(nAtomsStr) // "-" // TRIM(requiredTempStr) // "-" // TRIM(boxLengthStr) // ".csv"
	!velocitiesFileName = "dati/autocorrelazione/" // "velocities_" &
	!	// TRIM(nAtomsStr) // "-" // TRIM(requiredTempStr) // "-" // TRIM(boxLengthStr) // ".csv"
	
	!CALL saveBoltzmann(velocitiesFileName,atoms(:,:),.TRUE.)
	OPEN(2, file=fileName)
	WRITE(2,*) "index,time,boxLength,volume,energy,temperature,pressure,momentum"
	
	DO WHILE (time<t_fin)
		
		CALL Verlet(dynamic,atoms(:,:),time,step)
		
		time=time+step
		nSteps = nSteps + 1
		kShow = kShow+1
		
		CALL periodicConditions(atoms(1:3,:))
		
		WRITE(2,*) nSteps, ",", time, ",", boxLength, ",", boxLength**3, ",", energy,",", temperature,",",pressure,",",momentum
		!CALL saveBoltzmann(velocitiesFileName,atoms(:,:),.FALSE.)
		
		IF (kShow == showEach) THEN
			kShow = 0
			
			! Stampa su schermo i dati fisici
			PRINT*, nSteps,time,energy,potentialEnergy,temperature,pressure, momentum
		END IF
		
	END DO
	
	CLOSE(2)
	
END SUBROUTINE saveAllSteps
