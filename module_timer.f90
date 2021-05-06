!-----------------------------------------------------------------------
!
! 区間処理時間計測モジュール
!
! COPYRIGHT (C) 2017 FUJITSU LIMITED. ALL RIGHTS RESERVED.
!
!-----------------------------------------------------------------------
!
! 使用方法
!
!   1.本モジュールの変数T_NTIME,T_OFNAME,T_KOUTを修正する
!     T_NTIME  : 測定区間の数(実際より多くても良い)
!     T_OFNAME : 測定時間を出力するファイル名
!     T_KOUT   : 測定時間を出力するファイルのファイル番号
!                (T_KOUTが6のときは標準出力に出力される)
!
!   2.本モジュールのtimer_setnameで測定区間の名称をセットする
!     (例) T_NAME( 1) = "region1"
!     (注) T_NAMEの配列インデクスはT_NTIMEまで使用可能
!     (注) T_NAMEの配列インデクスは測定区間を示すIDになる
!
!   3.本モジュールを使用するモジュール、サブルーチンでuseする
!     (例) use module_timer
!
!   4.プログラムの最初にtimer_initサブルーチンをコールする
!     (例) call timer_init()
!     (注) MPI環境で使用する場合は、MPI_INITの後でコールする
!
!   5.プログラムの最後でtimer_finalizeサブルーチンをコールする
!     (例) call timer_finalize()
!     (注) MPI環境で使用する場合は、MPI_FINALIZEの前でコールする
!
!   6.測定したい区間をtimer_start,timer_endサブルーチンのコールで囲う
!     (例) call timer_start(1)
!               :
!             (処理)
!               :
!          call timer_end(1)
!     (注) 組み合わせとなるtimer_startとtimer_endの引数は同じ数字
!          (2のT_NAMEのインデクス)を指定する必要がある
!
!   7.MPIが無い環境で使用する場合
!     コンパイルオプションに「-D_TIMER_NOMPI_」を付与してコンパイルを行う
!
!   8.コンパイル、リンク
!     本ソースコードを含めてコンパイル、リンクを行う
!     (例) mpif90 module_timer.F90 main.f90
!
!   9.EXCELにはりつける場合
!     [データ]-[区切り位置]で[スペースによって右または左に…]を使用すると
!     綺麗に区切られる
!
!-----------------------------------------------------------------------
module module_timer
  implicit none

  ! 区間の数
  integer,parameter :: T_NTIME = 700

  ! 計測結果出力ファイル名
  character*256 :: T_OFNAME = "timer.out"

  ! 計測結果出力ファイルのファイル番号(6のとき、標準出力に出力)
  integer,parameter :: T_KOUT = 1234

  ! 区間の名前
  character*24 :: T_NAME(T_NTIME)

  ! 測定した時間
  real*8  :: T_TIMER(T_NTIME)

  ! 呼び出しカウント(timer_startの回数)
  integer :: T_SCOUNT(T_NTIME)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! 以下は内部で使用するワーク領域 !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! timer_init時の時刻
  real*8  :: T_START_TIME

  ! timer_start時の時刻格納用
  real*8  :: T_START(T_NTIME)

  ! 呼び出しカウント(timer_endの回数)
  integer :: T_ECOUNT(T_NTIME)


contains

  !---------------------------------------------------------------------
  ! 区間の名前のセット
  ! ここで測定区間の名称を記述します
  ! 配列のインデクスはそのままtimer_start、timer_endで使用するidになります
  !---------------------------------------------------------------------
  subroutine timer_setname
    T_NAME(100) = 'total-----------------'
    T_NAME(200) = '    init------------------'
    T_NAME(300) = '    main loop-------------'
    T_NAME(400) = '        ibm---------------'
    T_NAME(410) = '          l point---------'
    T_NAME(420) = '          index ----------'
    T_NAME(430) = '          delta-----------'
    T_NAME(500) = '        march-------------'
    T_NAME(510) = '          (12-a)----------'
    T_NAME(520) = '          (12-b)----------'
    T_NAME(530) = '          (12-c)----------'
    T_NAME(540) = '          (12-d)----------'
    T_NAME(550) = '          helmholtz-------'
    T_NAME(551) = '              RHS---------'
    T_NAME(552) = '              ft_f--------'
    T_NAME(553) = '              helm--------'
    T_NAME(554) = '              ft_i--------'
    T_NAME(555) = '              check-------'
    T_NAME(560) = '          poisson---------'
    T_NAME(561) = '              RHS---------'
    T_NAME(562) = '              ft_f--------'
    T_NAME(563) = '              helm--------'
    T_NAME(564) = '              ft_i--------'
    T_NAME(565) = '              check-------'
    T_NAME(570) = '          (12-gh)---------'
    T_NAME(600) = '        solution----------'
    T_NAME(700) = '        file--------------'
  end subroutine timer_setname

  !---------------------------------------------------------------------
  ! タイマーの初期化
  !---------------------------------------------------------------------
  subroutine timer_init
    integer :: i

    ! ゼロクリア
    do i=1,T_NTIME
      T_TIMER(i)  = 0.d0
      T_START(i)  = 0.d0
      T_SCOUNT(i) = 0
      T_ECOUNT(i) = 0
      T_NAME(i) = ""
    enddo

    ! 名前のセット
    call timer_setname

    ! 開始時刻のセット
    T_START_TIME = timer_gettime()

    return
  end subroutine timer_init

  !---------------------------------------------------------------------
  ! 測定開始
  !---------------------------------------------------------------------
  subroutine timer_start(id)
    integer :: id
    T_START(id) = timer_gettime()
    T_SCOUNT(id) = T_SCOUNT(id) + 1
  end subroutine timer_start

  !---------------------------------------------------------------------
  ! 測定終了
  !---------------------------------------------------------------------
  subroutine timer_end(id)
    integer :: id
    T_TIMER(id) = T_TIMER(id) + timer_spantime(T_START(id))
    T_ECOUNT(id) = T_ECOUNT(id) + 1
  end subroutine timer_end

  !---------------------------------------------------------------------
  ! 結果の出力
  !---------------------------------------------------------------------
  subroutine timer_finalize
#ifndef _TIMER_NOMPI_
    use mpi
#endif
    integer :: myrank, nprocs
#ifndef _TIMER_NOMPI_
    integer :: status(MPI_STATUS_SIZE)
    integer :: ierr
#endif
    integer :: i,j
    real*8  :: total_time

    ! 全時間
    total_time = timer_spantime(T_START_TIME)

    ! ランク番号、ランク数を取得
#ifndef _TIMER_NOMPI_
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
#else
    myrank = 0
    nprocs = 1
#endif

    ! ファイルオープン
    if (myrank.eq.0 .and. T_KOUT.ne.6) then
       open(T_KOUT,file=trim(T_OFNAME))
    endif

    ! 出力
    if( myrank.eq.0 ) then
      ! ランク0は各ランクから情報を受け取り出力
      do j=0,nprocs-1
        if( j.ne.0 ) then
          !ランク0以外から情報を受信
#ifndef _TIMER_NOMPI_
          call MPI_RECV(total_time, 1, MPI_REAL8, j, 0, MPI_COMM_WORLD, status, ierr)
          call MPI_RECV(T_TIMER, T_NTIME, MPI_REAL8, j, 0, MPI_COMM_WORLD, status, ierr)
          call MPI_RECV(T_SCOUNT, T_NTIME, MPI_INTEGER, j, 0, MPI_COMM_WORLD, status, ierr)
          call MPI_RECV(T_ECOUNT, T_NTIME, MPI_INTEGER, j, 0, MPI_COMM_WORLD, status, ierr)
#endif
        endif

        ! 出力
!       write(T_KOUT,'(a)')              "+------------------------------ ------------- ---------- ----------"
        write(T_KOUT,'(a,i0,a)')          "***** RANK", j, " *****"
        write(T_KOUT,'(a)')              "+------------------------------ ------------- ---------- ----------"
        write(T_KOUT,'(a)')              "|  ID   REGION                        TIME[s]     SCOUNT     ECOUNT"
        write(T_KOUT,'(a)')              "+------------------------------ ------------- ---------- ----------"
        write(T_KOUT,'(a,14x,1x,f12.3)') "|       TOTAL TIME", total_time
        write(T_KOUT,'(a)')              "+------------------------------ ------------- ---------- ----------"
        do i=1,T_NTIME
          if( T_SCOUNT(i).gt.0 .or. T_ECOUNT(i).gt.0 ) then
            write(T_KOUT,'(a1,i4,1x,2x,a24,1x,f12.3,1x,i10,1x,i10)') "|", i, T_NAME(i), T_TIMER(i), T_SCOUNT(i), T_ECOUNT(i)
            if( T_SCOUNT(i) .ne. T_ECOUNT(i) ) then
              write(T_KOUT,'(3a)') "Warning : SCOUNT and ECOUNT are different in region [",trim(T_NAME(i)),"]."
            endif
          endif
        enddo
!       write(T_KOUT,'(a)')              "+------------------------------ ------------- ---------- ----------"
        write(T_KOUT,'(/)')
      enddo
    else
      !ランク0以外はランク0に情報を送信
#ifndef _TIMER_NOMPI_
      call MPI_SEND(total_time, 1, MPI_REAL8, 0, 0, MPI_COMM_WORLD, ierr)
      call MPI_SEND(T_TIMER, T_NTIME, MPI_REAL8, 0, 0, MPI_COMM_WORLD, ierr)
      call MPI_SEND(T_SCOUNT, T_NTIME, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
      call MPI_SEND(T_ECOUNT, T_NTIME, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
#endif
    endif 

    ! ファイルクローズ
    if (myrank.eq.0 .and. T_KOUT.ne.6) then
       close(T_KOUT)
    endif

    return
  end subroutine timer_finalize


  !---------------------------------------------------------------------
  ! 時刻の取得
  !---------------------------------------------------------------------
  real*8 function timer_gettime()
#ifndef _TIMER_NOMPI_
    use mpi
#endif
#ifndef _TIMER_NOMPI_
    timer_gettime = MPI_WTIME()
#else
    integer :: t1
    call system_clock(t1)
    timer_gettime = t1 + 0.1
#endif
  end function timer_gettime

  !---------------------------------------------------------------------
  ! 経過時間の取得 
  !---------------------------------------------------------------------
  real*8 function timer_spantime(s)
    real*8 :: s
#ifndef _TIMER_NOMPI_
    timer_spantime = timer_gettime() - s
#else
    integer :: t1, t2, t_rate, t_max
    t1 = s
    call system_clock(t2, t_rate, t_max)
    if( t2 < t1 ) then
!      timer_spantime = (t_max - t1) + t2 + 1
      timer_spantime = 0.d0
    else
      timer_spantime = t2 - t1
    endif
    timer_spantime = timer_spantime / DBLE(t_rate)
#endif
  end function timer_spantime

end module module_timer

