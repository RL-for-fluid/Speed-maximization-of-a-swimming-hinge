!*****************************************************************************!
 SUBROUTINE ibm(delta_x,delta_y,Udl,Vdl,Xl,Yl,ij,angle,avel)
!*****************************************************************************!
#ifdef  _TIMER_
USE module_timer
#endif
USE parameter, ONLY:nk,nx,ny,nl,np,dx,dy,dl,h,hi,pi,dt,rei,dxi,dyi,xo,yo
IMPLICIT NONE
REAL(nk),INTENT(OUT) :: delta_x(nx,ny,nl,np),delta_y(nx,ny,nl,np) !�f���^�֐�
REAL(nk),INTENT(OUT) :: Udl(nl,np),Vdl(nl,np)                     !Lagrange�_�̕��̑��x
REAL(nk),INTENT(OUT) :: Xl(nl,np),Yl(nl,np)                       !Lagrange�_�̍��W
INTEGER,INTENT(OUT) :: ij(nl,np,2)                                !���̕t�߂̊i�q
REAL(nk),INTENT(IN) :: angle(np),avel(np)                         !�e���̂̕�����

REAL(nk) :: absr                  !Euler�_��Lagrange�_�̋����I�Ȃ���
REAL(nk) :: phi                   !�f���^�֐��̈ꕔ
REAL(nk) :: x,y                   !Euler���W
INTEGER :: i,j                    !�i�q
INTEGER :: l                      !�\�ʗv�f
INTEGER :: m                      !���̂̐�
!=============================================================================! 
#ifdef _TIMER_
  CALL timer_start(400)
#endif

!==Lagrange�_�̕��̈ʒu�Ƒ��x==!
#ifdef _TIMER_
  CALL timer_start(410)
#endif
!$OMP PARALLEL DO COLLAPSE(2)
DO m = 1, np
DO l = 1, nl
  Xl(l,m) = xo+dl*l*cos(angle(m))
  Yl(l,m) = yo+dl*l*sin(angle(m))
  Udl(l,m) = -avel(m)*dl*l*SIN(angle(m))
  Vdl(l,m) = +avel(m)*dl*l*COS(angle(m))
ENDDO
ENDDO
!$OMP END PARALLEL DO
#ifdef _TIMER_
  CALL timer_end(410)
#endif


!==���̎���̊i�q==!
#ifdef _TIMER_
  CALL timer_start(420)
#endif
!$OMP PARALLEL DO COLLAPSE(2)
DO m = 1, np
DO l = 1, nl
  ij(l,m,1)=INT(Xl(l,m)*dxi)-5
  ij(l,m,2)=INT(Yl(l,m)*dyi)-5
  IF(ij(l,m,1)   <1 ) ij(l,m,1)=1
  IF(ij(l,m,2)   <1 ) ij(l,m,2)=1
  IF(ij(l,m,1)+10>nx) ij(l,m,1)=nx
  IF(ij(l,m,2)+10>ny) ij(l,m,2)=ny
ENDDO
ENDDO
!$OMP END PARALLEL DO
#ifdef _TIMER_
  CALL timer_end(420)
#endif

!==�f���^�֐�==!
#ifdef _TIMER_
  CALL timer_start(430)
#endif
!$OMP PARALLEL DO PRIVATE(x,y,absr,phi) COLLAPSE(2)
DO m = 1, np
DO l = 1, nl
DO j = ij(l,m,2), ij(l,m,2)+10
DO i = ij(l,m,1), ij(l,m,1)+10
  !=x�������x�̂��߂̃f���^�֐�=!
  x=i*dx ; y=(j-0.5_nk)*dy                      !���W

  absr = ABS(x-Xl(l,m))*hi                     
  phi = 0
  IF(0.5 .le. absr .and. absr .le. 1.5)THEN
    phi = (5-3*absr-sqrt(-3*(1-absr)**2+1))/6
  ENDIF
  IF(absr .le. 0.5)THEN
    phi = (1+sqrt(-3*absr**2+1))/3
  ENDIF
  
  delta_x(i,j,l,m) = phi*hi
  
  absr = ABS(y-Yl(l,m))*hi
  phi = 0
  IF(0.5 .le. absr .and. absr .le. 1.5)THEN
    phi = (5-3*absr-sqrt(-3*(1-absr)**2+1))/6
  ENDIF
  IF(absr .le. 0.5)THEN
    phi = (1+sqrt(-3*absr**2+1))/3
  ENDIF
  
  delta_x(i,j,l,m) = delta_x(i,j,l,m)*phi*hi
  
  !=y�������x�̂��߂̃f���^�֐�=!
  x=(i-0.5_nk)*dx ; y=j*dy                      !���W

  absr = ABS(x-Xl(l,m))*hi
  phi = 0
  IF(0.5 .le. absr .and. absr .le. 1.5)THEN
    phi = (5-3*absr-sqrt(-3*(1-absr)**2+1))/6
  ENDIF
  IF(absr .le. 0.5)THEN
    phi = (1+sqrt(-3*absr**2+1))/3
  ENDIF
  
  delta_y(i,j,l,m) = phi*hi
  
  absr = ABS(y-Yl(l,m))*hi
  phi = 0
  IF(0.5 .le. absr .and. absr .le. 1.5)THEN
    phi = (5-3*absr-sqrt(-3*(1-absr)**2+1))/6
  ENDIF
  IF(absr .le. 0.5)THEN
    phi = (1+sqrt(-3*absr**2+1))/3
  ENDIF
  
  delta_y(i,j,l,m) = delta_y(i,j,l,m)*phi*hi
  
ENDDO
ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO
#ifdef _TIMER_
  CALL timer_end(430)
#endif


#ifdef _TIMER_
  CALL timer_end(400)
#endif

END SUBROUTINE ibm
