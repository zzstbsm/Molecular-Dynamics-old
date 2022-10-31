MODULE parametri

	INTEGER, PARAMETER	:: nAtoms = 260
	DOUBLE PRECISION	:: boxLength ! Lato del cubo
	DOUBLE PRECISION	:: r0 ! Passo del reticolo
	DOUBLE PRECISION	:: time = 0
	DOUBLE PRECISION, PARAMETER	:: cutOff = 2
	INTEGER			:: m ! numero di celle
	LOGICAL	:: flagOneCell
	
	LOGICAL	:: ifMFP = .FALSE. ! Se bisogna salvare il cammino libero medio
	LOGICAL	:: ifEnergy = .FALSE. ! Se bisogna salvare l'energia in Monte Carlo
	DOUBLE PRECISION	:: crystal(1:3,1:nAtoms)
	
	INTEGER, PARAMETER	:: minMeasures = 60
	INTEGER, PARAMETER	:: maxMeasures = 60
	
	DOUBLE PRECISION, PARAMETER	:: eps = 1.
	DOUBLE PRECISION, PARAMETER	:: sigma = 1.
	DOUBLE PRECISION, PARAMETER	:: mass = 1.
	
	DOUBLE PRECISION, PARAMETER	:: kb = 1
	
	INTEGER, PARAMETER	:: PI=3.141592654
	
	DOUBLE PRECISION	:: requiredTemp = 2.
	INTEGER, PARAMETER	:: maxCells = 14
	INTEGER, PARAMETER	:: typeRun = 4
	! 1 = Energy conservation: Verlet VS Runge Kutta
	! 3 = Draw Wan Deer Waals curve at fixed temperature
	! 4 = Fondi il cristallo
	! 5 = Monte Carlo
	! 6 = Salva tutti tutti gli step di integrazione per l'analisi dell'autocorrelazione
	! 7 = Salva le velocità ogni saveEach steps
	! 8 = Fusione con Monte Carlo
	
END MODULE parametri
