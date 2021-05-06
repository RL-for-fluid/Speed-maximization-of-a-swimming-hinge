!*****************************************************************************!
  MODULE parameter
!*****************************************************************************!
IMPLICIT NONE
!==精度==!
INTEGER,PARAMETER :: nk = 8

!==数学的定数==!
REAL(16),PARAMETER :: piq = 3.1415926535897932384626433832795028_16!円周率
REAL(nk),PARAMETER :: pi = piq 

!==系の設定==!
INTEGER,PARAMETER :: nx = 400              !流れ方向格子数
INTEGER,PARAMETER :: ny = 400              !高さ方向格子数
REAL(16),PARAMETER :: lxq = 20             !流れ方向の長さ
REAL(16),PARAMETER :: lyq = 20             !高さ方向の長さ
REAL(nk),PARAMETER :: lx=lxq, ly=lyq
REAL(nk),PARAMETER :: xo=10, yo=10         !座標の中央位置

!==時間関係==!
REAL(16),PARAMETER :: dtq = 0.005Q0!lyq/ny/10/2           !時間間隔
REAL(nk),PARAMETER :: dt = dtq, dti = 1 / dtq             !時間幅とその逆数
INTEGER,PARAMETER :: total_time=20                         !計算時間
INTEGER,PARAMETER :: max_step =total_time/dtq+1Q-8        !反復回数
INTEGER,PARAMETER :: file_step = 1/dtq                    !ファイル出力間隔

!==Reynolds数関係==!
REAL(16),PARAMETER :: req = 100,reiq = 1/req
REAL(nk),PARAMETER :: re = req, rei=1/req

!==計算上便利なために(高速化を図る)==!
REAL(nk),PARAMETER :: dx=lxq/nx, dy=lyq/ny, h=dx
REAL(nk),PARAMETER :: dxi=nx/lxq, dyi=ny/lyq, hi=dxi
REAL(nk),PARAMETER :: dx2i=(nx/lxq)**2, dy2i=(ny/lyq)**2
REAL(nk),PARAMETER :: hhdxi=(nx/lxq)/4, hhdyi=(ny/lyq)/4
REAL(nk),PARAMETER :: dx2rei=(nx/lxq)**2/req, dy2rei=(ny/lyq)**2/req

REAL(nk),PARAMETER :: length=1.Q0                !棒の長さ 

!==埋め込み境界法==!
INTEGER,PARAMETER :: np = 2                     !物体の個数
INTEGER,PARAMETER :: nl = length/h              !表面要素の個数
REAL(nk),PARAMETER :: dl=length/nl              !表面要素の間隔
REAL(nk),PARAMETER :: dv = h**2                 !表面要素の面積

!==物体関係==!
REAL(nk),PARAMETER :: mass = 0.1Q0                  !棒のz方向単位長さ当たりの質量

!==ルンゲクッタの係数==!
REAL(nk) :: alpha(3),gamma(3),zeta(3)

!==FFTに対する設定==!
INTEGER(8) :: planf, plani                        !Fourier変換の返り値
REAL(nk) :: fftp(nx,ny)                           !Fourier変換前
COMPLEX(nk) :: ffts(0:nx/2,0:ny-1)                !Fourier変換後
COMPLEX(nk) :: laplacei(0:nx/2,0:ny-1,0:3)        !Laplaceの逆数
REAL(nk),PARAMETER :: nxnyi = 1.D0 / ( nx * ny )  !正規化

END MODULE parameter
