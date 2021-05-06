!*****************************************************************************!
  PROGRAM main
!*****************************************************************************!
#ifdef _TIMER_
USE module_timer
#endif
USE parameter, ONLY:nk,nx,ny,np,max_step,xo,yo,dt,pi
IMPLICIT NONE
!==����==!
INTEGER :: n                                          !����step

!==���̂̕�����==!
REAL(nk) :: u(nx,ny),v(nx,ny)                         !���x
REAL(nk) :: p(nx,ny)                                  !����
REAL(nk) :: u1(nx,ny),v1(nx,ny),u2(nx,ny),v2(nx,ny)   !�����Q�N�b�^�̂��߂̑��x

!==���̂̕�����==!
REAL(nk) :: pos(2),vel(2),ac(2)     !���Ԉʒu, ���x, �����x (x����,y����)
REAL(nk) :: angle(np),avel(np)      !�e�_�̊p�x, �p���x (�_1,�_2)
REAL(nk) :: avel1(np),avel2(np)     !�����Q�N�b�^�̂��߂̊e�_�̊p���x (�_1,�_2)

REAL(nk) :: torque  !�_1���󂯂�g���N
REAL(nk) :: period  !�g���N�̎���
REAL(nk) :: drag(2) !���̂��󂯂�R�� (x,y)

!==���Ԍv��==!
INTEGER(8) :: ns,t0,t1,t2,t_rate,t_max
INTEGER :: err

!==Python==!
CHARACTER(LEN=100) :: filename
INTEGER(8) :: file_size
INTEGER :: episode
REAL(nk) :: keep_state(2,0:1000)
INTEGER :: sample,i
INTEGER :: total_step

!=============================================================================!
#ifdef _TIMER_
  CALL timer_init
  CALL timer_start(100)
  CALL timer_start(200)
#endif


!OPEN(111,FILE='position_x_movie.bin',FORM='BINARY')
!OPEN(113,FILE='position_y_movie.bin',FORM='BINARY')
!OPEN(110,FILE='omg_movie.bin',FORM='BINARY')
!OPEN(112,FILE='angle_1_movie.bin',FORM='BINARY')
!OPEN(114,FILE='angle_2_movie.bin',FORM='BINARY')


!==��������==!
!u=0;v=0                   !���̑��x
!pos(1)=xo; pos(2)=yo      !���Ԃ̈ʒu
!vel=0                     !���Ԃ̑��x
!angle(1)=0; angle(2)=pi   !�e�_�̊p�x
!avel=0                    !�e�_�̊p���x


OPEN(200,FILE='INITIAL_u.bin',FORM='BINARY',STATUS='old')
OPEN(201,FILE='INITIAL_v.bin',FORM='BINARY',STATUS='old')
OPEN(202,FILE='INITIAL_pos.bin',FORM='BINARY',STATUS='old')
OPEN(203,FILE='INITIAL_vel.bin',FORM='BINARY',STATUS='old')
OPEN(204,FILE='INITIAL_p.bin',FORM='BINARY',STATUS='old')
OPEN(205,FILE='INITIAL_angle.bin',FORM='BINARY',STATUS='old')
OPEN(206,FILE='INITIAL_avel.bin',FORM='BINARY',STATUS='old')
OPEN(207,FILE='INITIAL_keep.bin',FORM='BINARY',STATUS='old')
  READ(200) u
  READ(201) v
  READ(202) pos
  READ(203) vel
  READ(204) p
  READ(205) angle
  READ(206) avel
  READ(207) keep_state
CLOSE(200);CLOSE(201);CLOSE(202);CLOSE(203)
CLOSE(205);CLOSE(206);CLOSE(204);CLOSE(207)

!==fftwn�Ȃǂ̐ݒ�==!
CALL set_all
sample=1000

#ifdef _TIMER_
  CALL timer_end(200)
#endif


#ifdef _TIMER_
  CALL timer_start(300)
#endif
CALL system_clock(t1) 

!period=1                 !�g���N�̎���

total_step=0

DO   !�G�s�\�[�h����

  file_size=0

  READ(*,*) episode


  OPEN(1000,FILE="action.bin",FORM='BINARY',STATUS='old',ACCESS='DIRECT',RECL=8)
  OPEN(1001,FILE="observe.bin",FORM='BINARY',STATUS='replace')
  WRITE(*,*) 'INPUT0'


!==���Ԕ��W==!
DO n = 1, max_step
  !WRITE(*,*) "==========",n,"==========",pos(1),vel(1)
  
  !torque=COS(2*pi*n*dt/period) !�g���N�̐ݒ�
  
  WRITE(*,*) 'INPUT'
  
  DO WHILE(file_size<8*n)
    inquire(unit=1000, size=file_size)
    !WRITE(*,*) file_size
  ENDDO
  
  READ(1000,REC=n) torque
  total_step=total_step+1
  torque = TANH(torque)+0.1D0*SIN(total_step*pi*2/200/5) !
  
#ifdef _TIMER_
  CALL timer_start(500)
#endif
  CALL rg(u1,v1,p,u ,v ,u ,v ,pos,angle,vel,avel1,avel ,avel ,drag,torque,ac,1)    !�����Q�N�b�^1step��
  CALL rg(u2,v2,p,u1,v1,u ,v ,pos,angle,vel,avel2,avel1,avel ,drag,torque,ac,2)    !�����Q�N�b�^2step��
  CALL rg(u ,v ,p,u2,v2,u1,v1,pos,angle,vel,avel ,avel2,avel1,drag,torque,ac,3)    !�����Q�N�b�^3step��
#ifdef _TIMER_
  CALL timer_end(500)
#endif

  !===�x����W�n�̂���===!
    DO i = 0, sample-1
      keep_state(:,i) =  keep_state(:,i+1)
    ENDDO
    keep_state(1,sample)=ATAN2(SIN(angle(1)),COS(angle(1)))
    keep_state(2,sample)=avel(1)
    !keep_state(3,sample)=COS(angle(1))
    !keep_state(4,sample)=avel(1)

 ! CALL file_make(u,v,p,pos,angle,vel,avel,drag,torque,ac,n)  !�t�@�C���쐬
  
  WRITE(1001) vel(1),ATAN2(SIN(angle(1)),COS(angle(1))),avel(1),keep_state(:,800),keep_state(:,600),keep_state(:,400)
  
  WRITE(*,*) 'OUTPUT'
  
ENDDO


  WRITE(*,*) 'END'
  CLOSE(1000);CLOSE(1001)

ENDDO


CALL system_clock(t2,t_rate,t_max)
#ifdef _TIMER_
  CALL timer_end(300)
  CALL timer_end(100)
  CALL timer_finalize
#endif
WRITE(*,*) (t2-t1)/dble(t_rate)


!CLOSE(110);CLOSE(111);CLOSE(112);CLOSE(113);CLOSE(114)


!OPEN(200,FILE='INITIAL_u.bin',FORM='BINARY')
!OPEN(201,FILE='INITIAL_v.bin',FORM='BINARY')
!OPEN(202,FILE='INITIAL_pos.bin',FORM='BINARY')
!OPEN(203,FILE='INITIAL_vel.bin',FORM='BINARY')
!OPEN(204,FILE='INITIAL_p.bin',FORM='BINARY')
!OPEN(205,FILE='INITIAL_angle.bin',FORM='BINARY')
!OPEN(206,FILE='INITIAL_avel.bin',FORM='BINARY')
!OPEN(207,FILE='INITIAL_keep.bin',FORM='BINARY')
!  WRITE(200) u
!  WRITE(201) v
!  WRITE(202) pos
!  WRITE(203) vel
!  WRITE(204) p
!  WRITE(205) angle
!  WRITE(206) avel
!  WRITE(207) keep_state
!CLOSE(200);CLOSE(201);CLOSE(202);CLOSE(203)
!CLOSE(205);CLOSE(206);CLOSE(204);CLOSE(207)

END
