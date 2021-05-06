!*****************************************************************************!
  SUBROUTINE file_make(u,v,p,pos,angle,vel,avel,drag,torque,ac,time)
!*****************************************************************************!
#ifdef  _TIMER_
USE module_timer
#endif
USE parameter, ONLY:nk,nx,ny,np,dx,dy,xo,yo,dxi,dyi,dti,max_step,pi
IMPLICIT NONE
CHARACTER(LEN=100) :: filename

REAL(nk),INTENT(IN) :: u(nx,ny),v(nx,ny)        !���x
REAL(nk),INTENT(IN) :: p(nx,ny)                 !����
REAL(nk),INTENT(IN) :: pos(2),vel(2),ac(2)      !���Ԃ̈ʒu, ���x, �����x (x,y)
REAL(nk),INTENT(IN) :: angle(np),avel(np)       !�e�_�̊p�x, �p���x (�_1,�_2)
REAL(nk),INTENT(IN) :: torque                   !�_1�ɂ�����g���N
REAL(nk),INTENT(IN) :: drag(2)                  !���̂����̂���󂯂�R�� (x,y)
INTEGER,INTENT(IN) :: time                      !������

REAL(nk) :: omg(nx,ny)                          !�Q�x

INTEGER :: i,j,im,ip,jm,jp                      !�i�q
INTEGER,SAVE :: nd=0                            !����������
!=============================================================================!
#ifdef _TIMER_
  CALL timer_start(700)
#endif

!===�Q�x===!
DO j = 1, ny
DO i = 1, nx
  ip = MOD(i,nx)+1; jp = MOD(j,ny)+1
  omg(i,j)=(-v(i,j)+v(ip,j))*dxi-(-u(i,j)+u(i,jp))*dyi
ENDDO
ENDDO

!==�e�����ʂ̗���==!
!WRITE(100) pos(1)
!WRITE(101) pos(2)
!WRITE(102) ATAN2(SIN(angle(1)),COS(angle(1)))
!WRITE(103) vel(1)
!WRITE(104) vel(2)
!WRITE(105) avel(1)
!WRITE(106) drag(1)
!WRITE(107) drag(2)
!WRITE(108) avel(1)*torque
!WRITE(109) ac(1)
!WRITE(113) torque

!==�Q�x�̓���쐬==!
!IF(time-1==nd*int(dti+1.Q-10)/2) THEN
WRITE(110) omg
WRITE(111) pos(1)
WRITE(112) ATAN2(SIN(angle(1)),COS(angle(1)))
WRITE(113) pos(2)
WRITE(114) angle(2)
!nd = nd+1
!ENDIF


#ifdef _TIMER_
  CALL timer_end(700)
#endif

END SUBROUTINE file_make
