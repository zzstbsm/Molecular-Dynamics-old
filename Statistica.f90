DOUBLE PRECISION FUNCTION calcMean(inArray, nElements) RESULT(mean)
	
	IMPLICIT NONE
	
	INTEGER, INTENT(IN)				:: nElements
	DOUBLE PRECISION, INTENT(IN)	:: inArray(1:nElements)
	
	DOUBLE PRECISION	:: temp1
	INTEGER				:: i
	
	temp1 = 0
	
	DO i=1,nElements
		
		temp1 = temp1 + inArray(i)
		
	END DO
	
	temp1 = temp1/nElements
	mean = temp1
	
END FUNCTION calcMean

DOUBLE PRECISION FUNCTION calcStd(mean, inArray, nElements) RESULT(std)
	
	IMPLICIT NONE
	
	INTEGER, INTENT(IN)				:: nElements
	DOUBLE PRECISION, INTENT(IN)	:: inArray(1:nElements)
	DOUBLE PRECISION, INTENT(IN)	:: mean
	
	DOUBLE PRECISION	:: temp1
	INTEGER				:: i
	
	temp1 = 0
	
	DO i=1,nElements
		
		temp1 = temp1 + (inArray(i))**2
		
	END DO
	
	temp1 = temp1/nElements
	temp1 = temp1-mean**2
	temp1 = temp1/nElements
	std = SQRT(temp1)
	
END FUNCTION calcStd

SUBROUTINE calcAutocorrelation(tc, atoms, step)
	
	USE parametri
	USE physics
	
	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(OUT)	:: tc
	DOUBLE PRECISION, INTENT(INOUT)	:: atoms(1:6,1:nAtoms)
	DOUBLE PRECISION, INTENT(IN)	:: step
	
	DOUBLE PRECISION, PARAMETER		:: t_start = 100., considered_t = 50., t_end = 200.
	DOUBLE PRECISION				:: t
	
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE	:: measureList
	INTEGER										:: nMeasures, nSteps
	
	t = 0.
	
	nMeasures = FLOOR(t_end/step)
	nSteps = 0
	
	ALLOCATE(measureList(1:nMeasures))
	
	PRINT *, "Raccolta di dati per il calcolo dell'autocorrelazione"
	
	DO WHILE(nSteps<=nMeasures)
		
		CALL Verlet(dynamic,atoms(:,:),time,step)
		
		time=time+step
		nSteps = nSteps + 1
		
		CALL periodicConditions(atoms(1:3,:))
		
		measureList(nSteps) = temperature
		
	END DO
	
	!autocorrelation(simulazione["temperature"][0:nMeasures],nMeasures,maxIndex)
	
	DEALLOCATE(measureList)
	
END SUBROUTINE calcAutocorrelation

SUBROUTINE autocorrelation(output,measurelist,nMeasures,cutOff)
	
	IMPLICIT NONE
	
	INTEGER, INTENT(IN)				:: nMeasures
	INTEGER, INTENT(IN)				:: cutOff
	DOUBLE PRECISION, INTENT(IN)	:: measurelist(1:nMeasures)
	
	DOUBLE PRECISION, INTENT(OUT)	:: output(1:cutOff)
	
	INTEGER	:: i, j, kShow=0, showEach=1000
	DOUBLE PRECISION	:: mean
	
	DOUBLE PRECISION	:: autocorr(1:cutOff)
	
	mean = 0d0
	DO i = 1,nMeasures
		mean = mean + measurelist(i)/nMeasures
	END DO
	
	autocorr = 0
	DO i=1,cutOff
		DO j = 1,nMeasures-i
			autocorr(i) = autocorr(i) + (measurelist(j)-mean)*(measurelist(j+i)-mean)
		END DO
		autocorr(i) = autocorr(i)/(nMeasures-i)
		kShow = kShow+1
		IF (kShow == showEach) THEN
			WRITE(*,"(A12,i8,A,i8)") "Running => ",i,"/",cutOff
			kShow = 0
		END IF
	END DO
	
	output = autocorr
	
END SUBROUTINE autocorrelation

SUBROUTINE meanInBin(output,measurelist,nMeasures,binSize)
	
	IMPLICIT NONE
	
	INTEGER, INTENT(IN)				:: nMeasures
	INTEGER, INTENT(IN)				:: binSize
	DOUBLE PRECISION, INTENT(IN)	:: measurelist(1:nMeasures)
	
	DOUBLE PRECISION, INTENT(OUT)	:: output(1:nMeasures-binSize)
	
	INTEGER	:: i
	DOUBLE PRECISION	:: tempOut(1:nMeasures-binSize)
	
	tempOut = 0
	DO i=0,binSize-1
			tempOut(1) = tempOut(1) + measureList(binSize-i)
	END DO
	! Togli il primo elemento dal bin precedente e aggiungi l'ultimo elemento del bin attuale
	DO i=2,nMeasures-binSize
		tempOut(i) = tempOut(i-1) - measureList(i-1) + measureList(binSize+i)
	END DO
	tempOut = tempOut/binSize
	
	output = tempOut
	
END SUBROUTINE meanInBin
