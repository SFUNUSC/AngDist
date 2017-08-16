	PROGRAM ang_dist

c	Program calculates directional distributions of gamma rays
c	according to eq. 12.197 of W.D Hamilton et al., 'The Electromagnetic
c	Interaction in Nuclear Spectroscopy'

c	Statistical tensor and F-coefficents are calculated using existing 
c	code by K. Starosta.
	
c	Declare functions
	REAL*8 B,F1,A,LP
c	Declare variables
	REAL*8 lambda,sigmaj
	REAL*8 l,lprime,I_final,I_init,delta
	REAL*8 angle,norm_factor,dist_val
	REAL*8 ind_val
	REAL*8 A_val,B_val
	REAL*8 j
	REAL*4 ang
	
	WRITE(*,*)''
	WRITE(*,*)'GAMMA RAY ANGULAR DISTRIBUTION CALCULATOR'
	WRITE(*,*)'-----------------------------------------'
	WRITE(*,*)''

	WRITE(*,*)'Initial and final spin'
	WRITE(*,*)'Enter [I_final, I_initial]:'
	read(*,*)I_final,I_init
	
	WRITE(*,*)'Transition multipolarity (EL, ML)'
	WRITE(*,*)'Enter [L, mixing ratio with L+1 (L+1/L ratio, 0 for no mixing)]:'
	read(*,*)l,delta

	WRITE(*,*)'Width of distribution'
	WRITE(*,*)'Enter [sigma/I]:'
	read(*,*)sigmaj
	
c	lambda, L, L' add according to the triangle rule, so maximum lambda is L+L'
	if(delta.gt.0.0) then
		lambda=2*(L+1)
	else
		lambda=2*L
	endif
	
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

	
c	Get factor to normalize angular distribution data by
	norm_factor=0.	
	do i=0,180,5
		dist_val=0.
		do j=0,lambda,2
			ind_val=j;
			dist_val=dist_val+LP(j,cos(i*3.14159265359/180))*B(ind_val,I_init,sigmaj)*A(ind_val,l,I_final,I_init,delta)
		end do
		norm_factor=norm_factor+dist_val
	end do

c	Report angular distributions	
	WRITE(*,*)"Angular distribution function"
	WRITE(*,*)"Angle (deg), value "
	do i=0,180,5
		dist_val=0.
		do j=0,lambda,2
			ind_val=j;
c			Need to execute this on its own line rather than a write line 
c			to avoid recursive write statements (since there are write 
c			statements in this function)
			dist_val=dist_val+LP(j,cos(i*3.14159265359/180))*B(ind_val,I_init,sigmaj)*A(ind_val,l,I_final,I_init,delta)
		end do
		WRITE(*,*)i,dist_val
	end do
	
	WRITE(*,*)"Angular distribution at TIGRESS ring angles"
	WRITE(*,*)"Angle (deg), value "
	ang=37.524
	dist_val=0.
	do j=0,lambda,2
		ind_val=j;
c		Need to execute this on its own line rather than a write line 
c		to avoid recursive write statements (since there are write 
c		statements in this function)
		dist_val=dist_val+LP(j,cos(ang*3.14159265359/180))*B(ind_val,I_init,sigmaj)*A(ind_val,l,I_final,I_init,delta)
	end do
	WRITE(*,*)ang,dist_val
	ang=53.678
	dist_val=0.
	do j=0,lambda,2
		ind_val=j;
c		Need to execute this on its own line rather than a write line 
c		to avoid recursive write statements (since there are write 
c		statements in this function)
		dist_val=dist_val+LP(j,cos(ang*3.14159265359/180))*B(ind_val,I_init,sigmaj)*A(ind_val,l,I_final,I_init,delta)
	end do
	WRITE(*,*)ang,dist_val
	ang=81.838
	dist_val=0.
	do j=0,lambda,2
		ind_val=j;
c		Need to execute this on its own line rather than a write line 
c		to avoid recursive write statements (since there are write 
c		statements in this function)
		dist_val=dist_val+LP(j,cos(ang*3.14159265359/180))*B(ind_val,I_init,sigmaj)*A(ind_val,l,I_final,I_init,delta)
	end do
	WRITE(*,*)ang,dist_val
	ang=98.162
	dist_val=0.
	do j=0,lambda,2
		ind_val=j;
c		Need to execute this on its own line rather than a write line 
c		to avoid recursive write statements (since there are write 
c		statements in this function)
		dist_val=dist_val+LP(j,cos(ang*3.14159265359/180))*B(ind_val,I_init,sigmaj)*A(ind_val,l,I_final,I_init,delta)
	end do
	WRITE(*,*)ang,dist_val
	ang=126.322
	dist_val=0.
	do j=0,lambda,2
		ind_val=j;
c		Need to execute this on its own line rather than a write line 
c		to avoid recursive write statements (since there are write 
c		statements in this function)
		dist_val=dist_val+LP(j,cos(ang*3.14159265359/180))*B(ind_val,I_init,sigmaj)*A(ind_val,l,I_final,I_init,delta)
	end do
	WRITE(*,*)ang,dist_val
	ang=142.476
	dist_val=0.
	do j=0,lambda,2
		ind_val=j;
c		Need to execute this on its own line rather than a write line 
c		to avoid recursive write statements (since there are write 
c		statements in this function)
		dist_val=dist_val+LP(j,cos(ang*3.14159265359/180))*B(ind_val,I_init,sigmaj)*A(ind_val,l,I_final,I_init,delta)
	end do
	WRITE(*,*)ang,dist_val

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
	REAL*4 val
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
