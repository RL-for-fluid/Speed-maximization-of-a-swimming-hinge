!*****************************************************************************!
  MODULE parameter
!*****************************************************************************!
IMPLICIT NONE
!==���x==!
INTEGER,PARAMETER :: nk = 8

!==���w�I�萔==!
REAL(16),PARAMETER :: piq = 3.1415926535897932384626433832795028_16!�~����
REAL(nk),PARAMETER :: pi = piq 

!==�n�̐ݒ�==!
INTEGER,PARAMETER :: nx = 400              !��������i�q��
INTEGER,PARAMETER :: ny = 400              !���������i�q��
REAL(16),PARAMETER :: lxq = 20             !��������̒���
REAL(16),PARAMETER :: lyq = 20             !���������̒���
REAL(nk),PARAMETER :: lx=lxq, ly=lyq
REAL(nk),PARAMETER :: xo=10, yo=10         !���W�̒����ʒu

!==���Ԋ֌W==!
REAL(16),PARAMETER :: dtq = 0.005Q0!lyq/ny/10/2           !���ԊԊu
REAL(nk),PARAMETER :: dt = dtq, dti = 1 / dtq             !���ԕ��Ƃ��̋t��
INTEGER,PARAMETER :: total_time=20                         !�v�Z����
INTEGER,PARAMETER :: max_step =total_time/dtq+1Q-8        !������
INTEGER,PARAMETER :: file_step = 1/dtq                    !�t�@�C���o�͊Ԋu

!==Reynolds���֌W==!
REAL(16),PARAMETER :: req = 100,reiq = 1/req
REAL(nk),PARAMETER :: re = req, rei=1/req

!==�v�Z��֗��Ȃ��߂�(��������}��)==!
REAL(nk),PARAMETER :: dx=lxq/nx, dy=lyq/ny, h=dx
REAL(nk),PARAMETER :: dxi=nx/lxq, dyi=ny/lyq, hi=dxi
REAL(nk),PARAMETER :: dx2i=(nx/lxq)**2, dy2i=(ny/lyq)**2
REAL(nk),PARAMETER :: hhdxi=(nx/lxq)/4, hhdyi=(ny/lyq)/4
REAL(nk),PARAMETER :: dx2rei=(nx/lxq)**2/req, dy2rei=(ny/lyq)**2/req

REAL(nk),PARAMETER :: length=1.Q0                !�_�̒��� 

!==���ߍ��݋��E�@==!
INTEGER,PARAMETER :: np = 2                     !���̂̌�
INTEGER,PARAMETER :: nl = length/h              !�\�ʗv�f�̌�
REAL(nk),PARAMETER :: dl=length/nl              !�\�ʗv�f�̊Ԋu
REAL(nk),PARAMETER :: dv = h**2                 !�\�ʗv�f�̖ʐ�

!==���̊֌W==!
REAL(nk),PARAMETER :: mass = 0.1Q0                  !�_��z�����P�ʒ���������̎���

!==�����Q�N�b�^�̌W��==!
REAL(nk) :: alpha(3),gamma(3),zeta(3)

!==FFT�ɑ΂���ݒ�==!
INTEGER(8) :: planf, plani                        !Fourier�ϊ��̕Ԃ�l
REAL(nk) :: fftp(nx,ny)                           !Fourier�ϊ��O
COMPLEX(nk) :: ffts(0:nx/2,0:ny-1)                !Fourier�ϊ���
COMPLEX(nk) :: laplacei(0:nx/2,0:ny-1,0:3)        !Laplace�̋t��
REAL(nk),PARAMETER :: nxnyi = 1.D0 / ( nx * ny )  !���K��

END MODULE parameter
