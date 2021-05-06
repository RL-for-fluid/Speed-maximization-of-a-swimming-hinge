!*****************************************************************************!
  SUBROUTINE file_make(u,v,p,pos,angle,vel,avel,drag,torque,ac,time)
!*****************************************************************************!
#ifdef  _TIMER_
USE module_timer
#endif
USE parameter, ONLY:nk,nx,ny,np,dx,dy,xo,yo,dxi,dyi,dti,max_step,pi
IMPLICIT NONE
CHARACTER(LEN=100) :: filename

REAL(nk),INTENT(IN) :: u(nx,ny),v(nx,ny)        !速度
REAL(nk),INTENT(IN) :: p(nx,ny)                 !圧力
REAL(nk),INTENT(IN) :: pos(2),vel(2),ac(2)      !蝶番の位置, 速度, 加速度 (x,y)
REAL(nk),INTENT(IN) :: angle(np),avel(np)       !各棒の角度, 角速度 (棒1,棒2)
REAL(nk),INTENT(IN) :: torque                   !棒1にかけるトルク
REAL(nk),INTENT(IN) :: drag(2)                  !物体が流体から受ける抗力 (x,y)
INTEGER,INTENT(IN) :: time                      !反復回数

REAL(nk) :: omg(nx,ny)                          !渦度

INTEGER :: i,j,im,ip,jm,jp                      !格子
INTEGER,SAVE :: nd=0                            !無次元時間
!=============================================================================!
#ifdef _TIMER_
  CALL timer_start(700)
#endif

!===渦度===!
DO j = 1, ny
DO i = 1, nx
  ip = MOD(i,nx)+1; jp = MOD(j,ny)+1
  omg(i,j)=(-v(i,j)+v(ip,j))*dxi-(-u(i,j)+u(i,jp))*dyi
ENDDO
ENDDO

!==各物理量の履歴==!
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

!==渦度の動画作成==!
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
