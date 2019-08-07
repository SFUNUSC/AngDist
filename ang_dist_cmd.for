	PROGRAM ang_dist_cmd

c	Program calculates directional distributions of gamma rays
c	according to eq. 12.197 of W.D Hamilton et al., 'The Electromagnetic
c	Interaction in Nuclear Spectroscopy'

c	Statistical tensor and F-coefficents are calculated using existing 
c	code by K. Starosta.

c	Declare functions
	REAL*8 A,B,DIST
c	Declare variables
	REAL*8 lambda,sigmaj
	REAL*8 l,lprime,I_final,I_init,delta
	REAL*8 angle,norm_factor,dist_val
	REAL*8 q2,q4
	REAL*8 A_val,B_val
	REAL*8 j,order
	REAL*8 ang
	
	INTEGER i,count
	CHARACTER(len=15) arg
	
	count=iargc()
	i=0

	if(count.ne.7) then
		WRITE(*,*)'./ang_dist_cmd I_final I_init l delta sigmaj q2 q4'
		stop
	endif

	WRITE(*,*)'Parsed arguments:'
	if(i.lt.count) then
		CALL getarg(i+1, arg)
		read(arg,*)I_final
		WRITE (*,*)"I_final: ",arg
		i=i+1
	endif
	if(i.lt.count) then
		CALL getarg(i+1, arg)
		read(arg,*)I_init
		WRITE (*,*)"I_init: ",arg
		i=i+1
	endif
	if(i.lt.count) then
		CALL getarg(i+1, arg)
		read(arg,*)l
		WRITE (*,*)"l: ",arg
		i=i+1
	endif
	if(i.lt.count) then
		CALL getarg(i+1, arg)
		read(arg,*)delta
		WRITE (*,*)"delta: ",arg
		i=i+1
	endif
	if(i.lt.count) then
		CALL getarg(i+1, arg)
		read(arg,*)sigmaj
		WRITE (*,*)"sigmaj: ",arg
		i=i+1
	endif
	if(i.lt.count) then
		CALL getarg(i+1, arg)
		read(arg,*)q2
		WRITE (*,*)"q2: ",arg
		i=i+1
	endif
	if(i.lt.count) then
		CALL getarg(i+1, arg)
		read(arg,*)q4
		WRITE (*,*)"q4: ",arg
		i=i+1
	endif

	if(i.ne.7) then
		WRITE(*,*)'ERROR: wrong number of arguments!'
		stop
	endif
	

	WRITE(*,*)''
	WRITE(*,*)'GAMMA RAY ANGULAR DISTRIBUTION CALCULATOR'
	WRITE(*,*)'-----------------------------------------'
	WRITE(*,*)''


c L, L' are (possible) multipolarities being considered
c eg. for M1+E2, L and L' can both be either 1 or 2
c	lambda, L, L' add according to the triangle rule*,
c	so maximum lambda is L+L', minimum lambda is |L'-L|
c but all possible combinations need to be considered
c eg. for M1+E2, L=1 L'=1, L=2,L'=1, L=2, L'=2 ...
c so the max value will always be 2*(max L)
c and the min value will always be 0
	if(delta.ne.0.0) then
c		L'=L+1 (see pg. 542, Hamilton)
		lambda=2*(l+1)
	else
		lambda=2*l
	endif
c	*from the Wigner 3j terms in eq. 12.150, Hamilton 
	
	WRITE(*,*)""
	
	WRITE(*,*)'Angular distribution coefficient and statistical tensor values:'
1	FORMAT(" A_{",F10.0,"} = ",F10.6)
2	FORMAT(" B_{",F10.0,"} = ",F10.6)
c	Print A and B factors
	do j=0,lambda,2
		A_val=A(j,l,I_final,I_init,delta)
		B_val=B(j,I_init,sigmaj)
		WRITE(*,1)j,A_val
		WRITE(*,2)j,B_val
	end do
	WRITE(*,*)""

	

c	Report angular distributions	
	WRITE(*,*)"Angular distribution function"
	WRITE(*,*)"Angle (deg), value "
	do i=0,180,5
		dist_val=0.
		do j=0,lambda,2
			order=j
			ang=i
c			Need to execute this on its own line rather than a write line 
c			to avoid recursive write statements (since there are write 
c			statements in this function)
			dist_val=dist_val+DIST(order,l,ang,I_init,I_final,sigmaj,delta,q2,q4)
c			dist_val=dist_val+LP(j,cos(i*3.14159265359/180))*B(order,I_init,sigmaj)*A(order,l,I_final,I_init,delta)
		end do
		WRITE(*,*)i,dist_val
	end do
	
c	Report angular distributions for TIGRESS angles
	WRITE(*,*)"Angular distribution at TIGRESS ring angles"
	WRITE(*,*)"Angle (deg), value "
	ang=37.524
	dist_val=0.
	do j=0,lambda,2
		order=j;
		dist_val=dist_val+DIST(order,l,ang,I_init,I_final,sigmaj,delta,q2,q4)
	end do
	WRITE(*,*)ang,dist_val
	ang=53.678
	dist_val=0.
	do j=0,lambda,2
		order=j;
		dist_val=dist_val+DIST(order,l,ang,I_init,I_final,sigmaj,delta,q2,q4)
	end do
	WRITE(*,*)ang,dist_val
	ang=81.838
	dist_val=0.
	do j=0,lambda,2
		order=j;
		dist_val=dist_val+DIST(order,l,ang,I_init,I_final,sigmaj,delta,q2,q4)
	end do
	WRITE(*,*)ang,dist_val
	ang=98.162
	dist_val=0.
	do j=0,lambda,2
		order=j;
		dist_val=dist_val+DIST(order,l,ang,I_init,I_final,sigmaj,delta,q2,q4)
	end do
	WRITE(*,*)ang,dist_val
	ang=126.322
	dist_val=0.
	do j=0,lambda,2
		order=j;
		dist_val=dist_val+DIST(order,l,ang,I_init,I_final,sigmaj,delta,q2,q4)
	end do
	WRITE(*,*)ang,dist_val
	ang=142.476
	dist_val=0.
	do j=0,lambda,2
		order=j;
		dist_val=dist_val+DIST(order,l,ang,I_init,I_final,sigmaj,delta,q2,q4)
	end do
	WRITE(*,*)ang,dist_val

	END
	
c Function calculates contribution to angular distribution
	REAL*8 FUNCTION DIST(order,l,ang,I_init,I_final,sigmaj,delta,q2,q4)

	IMPLICIT NONE

C	Declare functions
	REAL*8 A,B,LP

C	Declare global variables
	REAL*8 l,I_init,I_final,delta,sigmaj,order
	REAL*8 ang
	REAL*8 q2,q4

	if(order.eq.2) then
		DIST=q2*LP(order,cos(ang*3.14159265359/180))*B(order,I_init,sigmaj)*A(order,l,I_final,I_init,delta)
	else if(order.eq.4) then
		DIST=q4*LP(order,cos(ang*3.14159265359/180))*B(order,I_init,sigmaj)*A(order,l,I_final,I_init,delta)
	else
		DIST=LP(order,cos(ang*3.14159265359/180))*B(order,I_init,sigmaj)*A(order,l,I_final,I_init,delta)
	endif

	RETURN
	END


	
c	Function calculates angular distribution coefficients
c	according to eq. 12.185, pg. 542
	REAL*8 FUNCTION A(LB,L,IF,II,delta)

	IMPLICIT NONE

C	Declare functions
	REAL*8 F1

C	Declare global variables
	REAL*8 LB,L,IF,II,delta

	A=(F1(LB,L,L,IF,II) + 2.*delta*F1(LB,L,L+1,IF,II) + delta*delta*F1(LB,L+1,L+1,IF,II))/(1 + delta*delta)
	
	RETURN
	END

	
c	Function reports output of the Legendre polynomial of the
c	specified order at the specified value
	REAL*8 FUNCTION LP(order,val)

	IMPLICIT NONE

C	Declare global variables
	REAL*8 val
	REAL*8 order
	
	LP=0.
	
	if(order.eq.0) then
		LP=1.
		RETURN
	else if(order.eq.1) then
		LP=val
		RETURN
	else if(order.eq.2) then
		LP=0.5*(3.*val**2 - 1.)
		RETURN
	else if(order.eq.3) then
		LP=0.5*(5.*val**3 - 3.*val)
		RETURN
	else if(order.eq.4) then
		LP=(0.125)*(35.*val**4 - 30.*val**2 + 3.)
		RETURN
	else if(order.eq.5) then
		LP=(0.125)*(63.*val**5 - 70.*val**3 + 15.*val)
		RETURN
	else if(order.eq.6) then
		LP=(0.0625)*(231.*val**6 - 315.*val**4 + 105.*val**2 - 5.)
		RETURN
	else
		WRITE(*,*)"Only orders up to 6 are implemented!"
	endif
	
	
	RETURN
	END
