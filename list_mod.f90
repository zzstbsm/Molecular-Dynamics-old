MODULE list_mod

	USE parametri

	TYPE listAtomsIndexType
		INTEGER				:: atomIndex
		TYPE(listAtomsIndexType), POINTER	:: previous => NULL()
		TYPE(listAtomsIndexType), POINTER	:: next => NULL()
	END TYPE listAtomsIndexType
	
	TYPE listAtomsIndexTypePtr
		TYPE(listAtomsIndexType), POINTER	:: head => NULL()
		INTEGER		:: counter = 0
	END TYPE listAtomsIndexTypePtr
	
	TYPE(listAtomsIndexTypePtr), DIMENSION(:,:,:) ,ALLOCATABLE	:: atomsInCellPtr
	INTEGER						:: atomsInCellCoordinates(1:3,1:nAtoms)
	
	TYPE(listAtomsIndexType), TARGET	:: listAtomsIndex(1:nAtoms)
	
	LOGICAL, PRIVATE	:: flagAllocated = .FALSE.
	
CONTAINS
	
	SUBROUTINE initList()
		
		IMPLICIT NONE
		
		IF (flagAllocated .eqv. .TRUE.) THEN
			CALL destroyList()
		END IF
		
		ALLOCATE(atomsInCellPtr(1:m,1:m,1:m))
		flagAllocated = .TRUE.
		
	END SUBROUTINE initList
	
	SUBROUTINE destroyList()
		
		IMPLICIT NONE
		
		DEALLOCATE(atomsInCellPtr)
		flagAllocated = .FALSE.
		
	END SUBROUTINE destroyList
	
	SUBROUTINE listAddHead(atomIndex,ix,iy,iz)
		
		IMPLICIT NONE
		
		INTEGER, INTENT(IN)	:: atomIndex, ix, iy, iz
		TYPE(listAtomsIndexType), POINTER	:: oldHead => NULL()
		TYPE(listAtomsIndexType), POINTER	:: newHead => NULL()
		
		newHead => listAtomsIndex(atomIndex)
		oldHead => atomsInCellPtr(ix,iy,iz)%head
		! Se esiste una lista, riarrangia gli indirizzi
		IF (ASSOCIATED(oldHead)) THEN
			newHead%next => oldHead
			NULLIFY(newHead%previous)
			oldHead%previous => newHead
		! Altrimenti inizializza gli indirizzi e poi inserisci
		ELSE
			NULLIFY(newHead%next)
			NULLIFY(newHead%previous)
		END IF
		atomsInCellPtr(ix,iy,iz)%head => newHead
		
		! Aggiorna il numero di elementi nella lista
		atomsInCellPtr(ix,iy,iz)%counter = atomsInCellPtr(ix,iy,iz)%counter + 1
		atomsInCellCoordinates(1,atomIndex) = ix
		atomsInCellCoordinates(2,atomIndex) = iy
		atomsInCellCoordinates(3,atomIndex) = iz
		
	END SUBROUTINE listAddHead
	
	SUBROUTINE listRmElement(atomIndex)
		
		IMPLICIT NONE
		
		INTEGER, INTENT(IN)	:: atomIndex
		INTEGER	:: ix,iy,iz
		
		TYPE(listAtomsIndexType), POINTER	:: previous => NULL()
		TYPE(listAtomsIndexType), POINTER	:: actual => NULL()
		TYPE(listAtomsIndexType), POINTER	:: next => NULL()
		
		! Sono ancora salvate le vecchie coordinate, le estraggo
		ix = atomsInCellCoordinates(1,atomIndex)
		iy = atomsInCellCoordinates(2,atomIndex)
		iz = atomsInCellCoordinates(3,atomIndex)
		
		actual => listAtomsIndex(atomIndex)
		previous => actual%previous
		next => actual%next
		
		IF (ASSOCIATED(previous)) THEN
			! Se è un nodo intermedio
			IF (ASSOCIATED(next)) THEN
				previous%next => actual%next
				next%previous => actual%previous
			
			! Se è l'ultimo nodo
			ELSE
				NULLIFY(previous%next)
			END IF
		! Se il nodo è il primo elemento
		ELSE
			! Se il nodo non è l'unico elemento
			IF (ASSOCIATED(next)) THEN
				NULLIFY(next%previous)
				atomsInCellPtr(ix,iy,iz)%head => actual%next
			! Se è l'unico elemento
			ELSE
				NULLIFY(atomsInCellPtr(ix,iy,iz)%head)
			END IF
		END IF
		
		NULLIFY(actual%previous)
		NULLIFY(actual%next)
		
		atomsInCellPtr(ix,iy,iz)%counter = atomsInCellPtr(ix,iy,iz)%counter - 1
		atomsInCellCoordinates(1,atomIndex) = 0
		atomsInCellCoordinates(2,atomIndex) = 0
		atomsInCellCoordinates(3,atomIndex) = 0
		
	END SUBROUTINE listRmElement

	SUBROUTINE listCountAndCollect(interactAtomsIndexFilling,interactAtomsIndex,ix,iy,iz)
		
		IMPLICIT NONE
		
		
		INTEGER, INTENT(INOUT)	:: interactAtomsIndexFilling
		INTEGER, INTENT(OUT)	:: interactAtomsIndex(1:nAtoms)
		INTEGER, INTENT(IN)		:: ix,iy,iz
				
		TYPE(listAtomsIndexType), POINTER	:: ptr => NULL()
		
		
		ptr => atomsInCellPtr(ix,iy,iz)%head
		! Esce dal DO WHILE quando ha finito di scorrere gli elementi della lista
		DO WHILE (ASSOCIATED(ptr))
			interactAtomsIndexFilling = interactAtomsIndexFilling + 1
			interactAtomsIndex(interactAtomsIndexFilling) = ptr%atomIndex
			ptr => ptr%next	
		END DO
		
		
	END SUBROUTINE listCountAndCollect
	
	SUBROUTINE fillInteractingArray(interactAtomsIndexFilling,interactAtomsIndex,ix,iy,iz)
	
		USE parametri
		
		IMPLICIT NONE
		
		INTEGER, INTENT(INOUT)	:: interactAtomsIndexFilling
		INTEGER, INTENT(OUT)	:: interactAtomsIndex(1:nAtoms)
		INTEGER, INTENT(IN)		:: ix,iy,iz
		
		INTEGER					:: tx,ty,tz
		
		
		! Pagina 32 appunti di Fisica Computazionale

		! Colleziona atomi nelle celle vicine (Sono 13 blocchi di codice quasi identici)

		! Blocco 1
		tx = ix
		ty = iy + 1
		tz = iz
		CALL modCell(tx,ty,tz)
		CALL listCountAndCollect(interactAtomsIndexFilling,interactAtomsIndex,tx,ty,tz)

		! Blocco 2
		tx = ix - 1
		ty = iy + 1
		tz = iz
		CALL modCell(tx,ty,tz)
		CALL listCountAndCollect(interactAtomsIndexFilling,interactAtomsIndex,tx,ty,tz)

		! Blocco 3
		tx = ix - 1
		ty = iy + 1
		tz = iz + 1
		CALL modCell(tx,ty,tz)
		CALL listCountAndCollect(interactAtomsIndexFilling,interactAtomsIndex,tx,ty,tz)

		! Blocco 4
		tx = ix
		ty = iy + 1
		tz = iz + 1
		CALL modCell(tx,ty,tz)
		CALL listCountAndCollect(interactAtomsIndexFilling,interactAtomsIndex,tx,ty,tz)

		! Blocco 5
		tx = ix + 1
		ty = iy + 1
		tz = iz + 1
		CALL modCell(tx,ty,tz)
		CALL listCountAndCollect(interactAtomsIndexFilling,interactAtomsIndex,tx,ty,tz)

		! Blocco 6
		tx = ix + 1
		ty = iy + 1
		tz = iz
		CALL modCell(tx,ty,tz)
		CALL listCountAndCollect(interactAtomsIndexFilling,interactAtomsIndex,tx,ty,tz)

		! Blocco 7
		tx = ix + 1
		ty = iy + 1
		tz = iz - 1
		CALL modCell(tx,ty,tz)
		CALL listCountAndCollect(interactAtomsIndexFilling,interactAtomsIndex,tx,ty,tz)

		! Blocco 8
		tx = ix
		ty = iy + 1
		tz = iz - 1
		CALL modCell(tx,ty,tz)
		CALL listCountAndCollect(interactAtomsIndexFilling,interactAtomsIndex,tx,ty,tz)

		! Blocco 9
		tx = ix - 1
		ty = iy + 1
		tz = iz - 1
		CALL modCell(tx,ty,tz)
		CALL listCountAndCollect(interactAtomsIndexFilling,interactAtomsIndex,tx,ty,tz)

		! Blocco 10
		tx = ix + 1
		ty = iy
		tz = iz
		CALL modCell(tx,ty,tz)
		CALL listCountAndCollect(interactAtomsIndexFilling,interactAtomsIndex,tx,ty,tz)

		! Blocco 11
		tx = ix + 1
		ty = iy
		tz = iz - 1
		CALL modCell(tx,ty,tz)
		CALL listCountAndCollect(interactAtomsIndexFilling,interactAtomsIndex,tx,ty,tz)

		! Blocco 12
		tx = ix
		ty = iy
		tz = iz - 1
		CALL modCell(tx,ty,tz)
		CALL listCountAndCollect(interactAtomsIndexFilling,interactAtomsIndex,tx,ty,tz)

		! Blocco 13
		tx = ix - 1
		ty = iy
		tz = iz - 1
		CALL modCell(tx,ty,tz)
		CALL listCountAndCollect(interactAtomsIndexFilling,interactAtomsIndex,tx,ty,tz)

	END SUBROUTINE fillInteractingArray
	
	SUBROUTINE modCell(x,y,z)
		
		USE parametri
		
		IMPLICIT NONE
		
		INTEGER, INTENT(INOUT)	:: x,y,z
		
		IF (x>m) THEN
			x=x-m
		ELSE IF (x<1) THEN
			x=x+m
		END IF
		IF (y>m) THEN
			y=y-m
		ELSE IF (y<1) THEN
			y=y+m
		END IF
		IF (z>m) THEN
			z=z-m
		ELSE IF (z<1) THEN
			z=z+m
		END IF
		
	END SUBROUTINE
	
END MODULE list_mod
