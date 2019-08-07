C*****************************************************************************
	REAL*8 FUNCTION B(LB,I,SGJ)
C*****************************************************************************

C	CALCULATES THE STATISTICAL TENSOR

	IMPLICIT NONE


C	DEKLARACJE FUNKCJI
	REAL*8		DCLEBG
	REAL*8 	FAZA

C	DEKLARACJE ZMIENNYCH GLOBALNYCH I BLOKOW COMMON
	REAL*8 		I
	REAL*8		LB,M,M1
	REAL*8		SGJ
	INTEGER 	IERR,IERCT
	COMMON/FGERCM/ 	IERR,IERCT

C	DEKLARACJE ZMIENNYCH LOKALNYCH
	REAL*8		SIGMA,NORM,SUM,E
	REAL*8    F,P

	REAL*8          ZERO,HALF,MHALF
	DATA            ZERO/0./,HALF/0.5/,MHALF/-0.5/


	P=2.*I
	P=FAZA(P)

	IF(SGJ.EQ.0.AND.P.EQ.1) THEN
	NORM=1.
	F=FAZA(I)
	SUM=F*DCLEBG(I,I,LB,ZERO,ZERO,ZERO)
        ENDIF

	IF(SGJ.EQ.0.AND.P.EQ.-1) THEN
	NORM=2.
	F=I+0.5
	SUM=FAZA(F)*DCLEBG(I,I,LB,MHALF,HALF,ZERO)
	F=I-0.5
	SUM=SUM+FAZA(F)*DCLEBG(I,I,LB,HALF,MHALF,ZERO)
c	WRITE(*,1000)
1000	FORMAT(' SIGMA=0 FOR HALF SPINS MEANS P(M=1/2)=P(M=-1/2)=1/2 ')
        ENDIF

	IF(SGJ.NE.0.) THEN

C	GRANICE SUMOWANIA
	SIGMA=SGJ*I

	SUM=0.
	NORM=0.
	
c Sum in equation 12.81 in Hamilton
	DO M=-I,I
c numerator
	SUM=SUM+FAZA(I+M)*DCLEBG(I,I,LB,-M,M,ZERO)*EXP(-M*M/2./SIGMA/SIGMA)
c denominator
	NORM=NORM+EXP(-M*M/2./SIGMA/SIGMA)
c	WRITE(*,*)M,I+M
c	WRITE(*,*)SUM,F,DCLEBG(I,I,LB,M1,M,ZERO),E
c	read(*,*)
	end DO

	ENDIF

c	Compute orientation parameter based on sum
c	Equation 12.81 in Hamilton
	B=SQRT(2.*I+1)*SUM/NORM

	RETURN
	END
C*****************************************************************************
	REAL*8 FUNCTION F3(LB1,LB0,LB,L,L1,IF,II)
C*****************************************************************************

	IMPLICIT NONE

C	DEKLARACJE FUNKCJI
	REAL*8 		F1
	REAL*8		DWIG3J,DWIG9J,DK
	REAL*8 	FAZA
C	DEKLARACJE ZMIENNYCH GLOBALNYCH I BLOKOW COMMON
	REAL*8 		IF,II
	REAL*8 		LB1,LB0,LB
	REAL*8 		L,L1
	INTEGER 	IERR,IERCT
	COMMON/FGERCM/ 	IERR,IERCT

C	DEKLARACJE ZMIENNYCH LOKALNYCH	
	REAL*8		NORM,FAKT
	REAL*8 	F		

	REAL*8          ZERO,ONE,MONE
	DATA            ZERO/0./,ONE/1./,MONE/-1./


	IF(LB1.EQ.0) THEN
	F3=DK(LB0,LB)*F1(LB,L,L1,IF,II)
	RETURN
	ELSE
	CONTINUE
	ENDIF
	

        F=L1+LB1+LB0+1
        F=FAZA(F)
	

        NORM=(2.*II+1.)
	NORM=NORM*(2.*IF+1.)
	NORM=NORM*(2.*L+1.)
	NORM=NORM*(2.*L1+1.)
	NORM=NORM*(2.*LB0+1.)
	NORM=NORM*(2.*LB1+1.)
 	NORM=NORM*(2.*LB+1.)
	NORM=SQRT(NORM)

	FAKT=DWIG3J(L,L1,LB,ONE,MONE,ZERO)

	FAKT=FAKT*DWIG9J(IF,L,II,IF,L1,II,LB1,LB,LB0)
	
	F3=F*NORM*FAKT
	RETURN
	END
C*****************************************************************************
	REAL*8 FUNCTION F1(LB,L,L1,IF,II)
C*****************************************************************************
	IMPLICIT NONE

C	DEKLARACJE FUNKCJI
	REAL*8       DWIG3J,DWIG6J
	REAL*8      FAZA	

C	DEKLARACJE ZMIENNYCH GLOBALNYCH
	REAL*8       IF,II
	REAL*8       LB,L,L1

C	DEKLARACJE ZMIENNYCH LOKALNYCH
	REAL*8       NORM,J3,J6
	REAL*8      F

	REAL*8          ZERO,ONE,MONE
	DATA            ZERO/0./,ONE/1./,MONE/-1./

	F=IF+II-1
	F=FAZA(F)
	NORM=(2.*LB+1.)*(2.*L+1.)*(2.*L1+1.)*(2.*II+1.)
	NORM=SQRT(NORM)
	J3=DWIG3J(L,L1,LB,ONE,MONE,ZERO)
	J6=DWIG6J(L,L1,LB,II,II,IF)

	F1=F*NORM*J3*J6


	RETURN
	END
C*****************************************************************************
	REAL*8 FUNCTION FAZA(K)
C*****************************************************************************
	
C	CALCULATES (-1)**K

	REAL*8 K
	
	FAZA=(-1.)**K
	
	END	
C*****************************************************************************
	REAL*8 FUNCTION DK(X,Y)
C*****************************************************************************

	IMPLICIT NONE

	REAL*8 X,Y
	
	IF(X.EQ.Y) THEN
	DK=1.
	ELSE
	DK=0.
	ENDIF

	RETURN
	END
