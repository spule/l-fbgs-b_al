!Copyright (c) 2011, 2018, Jorge Nocedal and Jose Luis Morales, Miha Polajnar
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
!    * Redistributions of source code must retain the above copyright
!      notice, this list of conditions and the following disclaimer.
!    * Redistributions in binary form must reproduce the above copyright
!      notice, this list of conditions and the following disclaimer in the
!      documentation and/or other materials provided with the distribution.
!    * Neither the name of the <organization> nor the
!      names of its contributors may be used to endorse or promote products
!      derived from this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
!DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


!
! CHANGES
! 21.12.2018 - Miha Polajnar
! Initial version
!


module l_bfgs_b_ag
  use iso_fortran_env, only: rk => real64, & ! Only for double precision
    stdout => output_unit
  use l_bfgs_b_org, only: setulb
  implicit none
  private
!-------------------------------------------------------------------------------
  type :: optim_prob_t
    character(len=256) :: name = ''
    integer :: n = 0, neqc = 0, nieqc = 0
    real(rk), allocatable :: x(:)
    real(rk), allocatable :: lb(:)
    real(rk), allocatable :: ub(:)
    integer, allocatable :: nbd(:)
    class(*), allocatable :: dat
    procedure(func_intrf), pointer, nopass :: func => null()
    procedure(grad_intrf), pointer, nopass :: grad => null()
    procedure(eq_const_intrf), pointer, nopass :: eq_const => null()
    procedure(grad_eq_const_intrf), pointer, nopass :: grad_eq_const => null()
    procedure(ieq_const_intrf), pointer, nopass :: ieq_const => null()
    procedure(grad_ieq_const_intrf), pointer, nopass :: grad_ieq_const => null()
  end type
  type :: lbfgsb_opt_t
    real(rk) :: factr = 1.e+7_rk
    real(rk) :: pgtol = 1.e-5_rk
    integer :: iprint = 0
    integer :: m = 5
    integer :: maxeval = 99
  end type
!-------------------------------------------------------------------------------
  interface
    real(rk) pure function func_intrf(x,dat)
      import :: rk
      real(rk), intent(in) :: x(:)
      class(*), optional, intent(in) :: dat
    end function
    pure subroutine grad_intrf(x,dat,grad)
      import :: rk
      real(rk), intent(in) :: x(:)
      real(rk), intent(out) :: grad(:)
      class(*), optional, intent(in) :: dat
    end subroutine
    pure subroutine eq_const_intrf(x,dat,ce)
      import :: rk
      real(rk), intent(in) :: x(:)
      class(*), optional, intent(in) :: dat
      real(rk), intent(out) :: ce(:)
    end subroutine
   pure subroutine grad_eq_const_intrf(x,dat,grad)
      import :: rk
      real(rk), intent(in) :: x(:)
      class(*), optional, intent(in) :: dat
      real(rk), intent(out) :: grad(:,:)
    end subroutine
    pure subroutine ieq_const_intrf(x,dat,ci)
      import :: rk
      real(rk), intent(in) :: x(:)
      class(*), optional, intent(in) :: dat
      real(rk), intent(out) :: ci(:)
    end subroutine
    pure subroutine grad_ieq_const_intrf(x,dat,grad)
      import :: rk
      real(rk), intent(in) :: x(:)
      class(*), optional, intent(in) :: dat
      real(rk), intent(out) :: grad(:,:)
    end subroutine
  end interface
!-------------------------------------------------------------------------------
  interface l_bfgs_b
    module procedure :: l_bfgs_b
    !module procedure :: l_bfgs_b_al_eq_const
    !module procedure :: l_bfgs_b_al_ieq_const
  end interface
  public :: l_bfgs_b, l_bfgs_b_test, optim_prob_t, lbfgsb_opt_t
!-------------------------------------------------------------------------------
contains
!------------------------CUSTOMIZED ORIGINAL DRIVER-----------------------------
  subroutine l_bfgs_b(x,lb,ub,func,grad,factr,pgtol,iprint,nbd,m,dat,&
    maxeval)
    real(rk), intent(inout) :: x(:)
!     x is a REAL array of length n.  On initial entry
!       it must be set by the user to the values of the initial
!       estimate of the solution vector.  Upon successful exit, it
!       contains the values of the variables at the best point
!       found (usually an approximate solution).
    real(rk), intent(in) :: lb(:)
!     lb is an array of length n that must be set by
!       the user to the values of the lower bounds on the variables. If
!       the i-th variable has no lower bound, l(i) need not be defined.
    real(rk), intent(in) :: ub(:)
!     ub is a REAL array of length n that must be set by
!       the user to the values of the upper bounds on the variables. If
!       the i-th variable has no upper bound, u(i) need not be defined.
    procedure(func_intrf) :: func
!     func is a pure function to be minimized.
    procedure(grad_intrf) :: grad
!     grad is a pure subroutne to determine gradient of func.
    real(rk), optional, intent(in) :: factr
    real(rk) :: factr_local
!     factr is a REAL variable that must be set by the user.
!       It is a tolerance in the termination test for the algorithm.
!       The iteration will stop when
!
!        (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
!
!       where epsmch is the machine precision which is automatically
!       generated by the code. Typical values for factr on a computer
!       with 15 digits of accuracy in double precision are:
!       factr=1.d+12 for low accuracy;
!             1.d+7  for moderate accuracy;
!             1.d+1  for extremely high accuracy.
!       The user can suppress this termination test by setting factr=0.
    real(rk), optional, intent(in) :: pgtol
    real(rk) :: pgtol_local
!     pgtol is a REAL variable.
!       On entry pgtol >= 0 is specified by the user.  The iteration
!         will stop when
!
!                 max{|proj g_i | i = 1, ..., n} <= pgtol
!
!         where pg_i is the ith component of the projected gradient.
!       The user can suppress this termination test by setting pgtol=0.
    integer, optional, intent(in) :: iprint
    integer :: iprint_local
!     iprint is an INTEGER variable that must be set by the user.
!       It controls the frequency and type of output generated:
!        iprint<0    no output is generated;
!        iprint=0    print only one line at the last iteration;
!        0<iprint<99 print also f and |proj g| every iprint iterations;
!        iprint=99   print details of every iteration except n-vectors;
!        iprint=100  print also the changes of active set and final x;
!        iprint>100  print details of every iteration including x and g;
!       When iprint > 0, the file iterate.dat will be created to
!                        summarize the iteration.
    integer, optional, intent(in) :: nbd(:)
    integer, allocatable :: nbd_local(:)
!     nbd is an INTEGER array of dimension n that must be set by the
!       user to the type of bounds imposed on the variables:
!       nbd(i)=0 if x(i) is unbounded,
!              1 if x(i) has only a lower bound,
!              2 if x(i) has both lower and upper bounds,
!              3 if x(i) has only an upper bound.
    integer, optional, intent(in) :: m
    integer :: m_local
!     m is an INTEGER variable that must be set by the user to the
!       number of corrections used in the limited memory matrix.
!       It is not altered by the routine.  Values of m < 3  are
!       not recommended, and large values of m can result in excessive
!       computing time. The range  3 <= m <= 20 is recommended.
    class(*), optional, intent(in) :: dat
!     dat is an ULMITED POLYMORHIC VALUE for additional parameters to
!        objective function.
    integer, optional, intent(in) :: maxeval
    integer :: maxeval_local
!     maxeval is an INTEGER variable which determines maximal number of
!       function in gradient evaluations
!-----------------------------Working variables---------------------------------
    real(rk), allocatable :: wa(:)
!     wa is a REAL array of length
!       (2mmax + 5)nmax + 11mmax^2 + 8mmax used as workspace.
    integer, allocatable :: iwa(:)
!     iwa is an INTEGER  array of length 3nmax used as
!       workspace. This array must not be altered by the user.
    character(len=60) :: task
!     task is a CHARACTER string of length 60.
!       On first entry, it must be set to 'START'.
!       On a return with task(1:2)='FG', the user must evaluate the
!         function f and gradient g at the returned value of x.
!       On a return with task(1:5)='NEW_X', an iteration of the
!         algorithm has concluded, and f and g contain f(x) and g(x)
!         respectively.  The user can decide whether to continue or stop
!         the iteration.
!       When
!         task(1:4)='CONV', the termination test in L-BFGS-B has been
!           satisfied;
!         task(1:4)='ABNO', the routine has terminated abnormally
!           without being able to satisfy the termination conditions,
!           x contains the best approximation found,
!           f and g contain f(x) and g(x) respectively;
!         task(1:5)='ERROR', the routine has detected an error in the
!           input parameters;
!       On exit with task = 'CONV', 'ABNO' or 'ERROR', the variable task
!         contains additional information that the user can print.
!       This array should not be altered unless the user wants to
!          stop the run for some reason.  See driver2 or driver3
!          for a detailed explanation on how to stop the run
!          by assigning task(1:4)='STOP' in the driver.
    character(len=60) :: csave
!     csave  is a CHARACTER working array of length 60.
    logical :: lsave(4)
!     lsave is a LOGICAL working array of dimension 4.
!       On exit with task = 'NEW_X', the following information is
!         available:
!       lsave(1) = .true.  the initial x did not satisfy the bounds;
!       lsave(2) = .true.  the problem contains bounds;
!       lsave(3) = .true.  each variable has upper and lower bounds.
    integer :: isave(44)
!     isave is an INTEGER working array of dimension 44.
!       On exit with task = 'NEW_X', it contains information that
!       the user may want to access:
!         isave(30) = the current iteration number;
!         isave(34) = the total number of function and gradient
!                         evaluations;
!         isave(36) = the number of function value or gradient
!                                  evaluations in the current iteration;
!         isave(38) = the number of free variables in the current
!                         iteration;
!         isave(39) = the number of active constraints at the current
!                         iteration;
    real(rk) :: dsave(29)
!     dsave is a REAL working array of dimension 29.
!       On exit with task = 'NEW_X', it contains information that
!         the user may want to access:
!         dsave(2)  = the value of f at the previous iteration;
!         dsave(5)  = the machine precision epsmch generated by the code;
!         dsave(13) = the infinity norm of the projected gradient;
!-------------------------------------------------------------------------------
    integer :: n ! number of variables
    integer :: i
    real(rk) :: func_local
    real(rk), allocatable :: grad_local(:)
!-------------------------------------------------------------------------------
    n = size(x)
!------------------------Take care of optional arguments------------------------
    factr_local = 1.e+7_rk  ! For moderate accuracy
    pgtol_local = 1.e-5_rk
    iprint_local = 0 ! Print only one line at the last iteration
    allocate(nbd_local(n),source=2) ! Both lower and upper bounds
    m_local = 7 ! Moderate number of corrections
    maxeval_local = 99
! Change non-default values
    if (present(factr)) factr_local = factr
    if (present(pgtol)) pgtol_local = pgtol
    if (present(iprint)) iprint_local = iprint
    if (present(nbd)) nbd_local = nbd
    if (present(m)) m_local = m
!-------------------------------------------------------------------------------
    allocate(grad_local(n))
    allocate(iwa(3*n))
    allocate(wa(2*m_local*n + 5*n + 11*m_local*m_local + 8*m_local))
!-------------------------------MAIN LOOP---------------------------------------
    task = 'START'
    do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
      task.eq.'START')
      ! Call to the original code
      call setulb ( n, m_local, x, lb, ub, nbd_local, func_local, grad_local, &
        factr_local, pgtol_local, wa, iwa, task, iprint_local, &
        csave, lsave, isave, dsave )
      ! Get objective function and gradient
      if (task(1:2) .eq. 'FG') then
        if (present(dat)) then
          func_local = func(x=x,dat=dat)
          call grad(x=x,dat=dat,grad=grad_local)
        else
          func_local = func(x=x)
          call grad(x=x,grad=grad_local)
        end if
      else if (task(1:5) .eq. 'NEW_X') then
!       The minimization routine has returned with a new iterate.
!       At this point have the opportunity of stopping the iteration
!       or observing the values of certain parameters

!       1) Terminate if the total number of f and g evaluations
!            exceeds maxeval.
        if (isave(34) .ge. maxeval_local)  &
          task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'
!       2) Maybe additional stopping criteria, if needed

        ! Print final results
        write (stdout,'(2(a,i5,4x),a,1p,d12.5,4x,a,1p,d12.5)') 'Iterate', &
          isave(30),'nfg =',isave(34),'f =',func_local,'|proj g| =',dsave(13)

        if (task(1:4) .eq. 'STOP') then
          write (stdout,*) task
          write (stdout,*) 'Final X='
          write (stdout,'((1x,1p, 6(1x,d11.4)))') (x(i),i = 1,n)
        end if

      end if
    end do
!-----------------------------END OF MAIN LOOP----------------------------------
  end subroutine


!---------------------------EQUALITY CONSTRAINS---------------------------------
!  subroutine l_bfgs_b_al_eq_const(x,lb,ub,func,grad,eq_const,factr,pgtol,m,dat,&
!    iprint,nbd,maxeval)
!    real(rk), intent(inout) :: x(:)
!    real(rk), intent(in) :: lb(:)
!    real(rk), intent(in) :: ub(:)
!    procedure(func_intrf) :: func
!    procedure(grad_intrf) :: grad
!    real(rk), intent(out) :: ce(:)
!!     ce is a REAL array of equality constrains values.
!    real(rk), intent(inout) :: lagr(:)
!!     lagr is a REAL array of Lagrange multipliers values.
!    procedure(eq_const_intrf) :: eq_const
!!     eq_const is a pure subroutine to determine values of equality constrains.
!    procedure(grad_eq_const_intrf) :: grad_eq_const
!!     grad_eq_const is a pure subroutine to determine gradient values of
!!       equality constrains.
!    real(rk), optional, intent(in) :: factr
!    real(rk) :: factr_local
!    real(rk), optional, intent(in) :: pgtol
!    real(rk) :: pgtol_local
!    integer, optional, intent(in) :: iprint
!    integer :: iprint_local
!    integer, optional, intent(in) :: nbd(:)
!    integer, allocatable :: nbd_local(:)
!    integer, optional, intent(in) :: m
!    integer :: m_local
!    class(*), optional, intent(in) :: dat
!    integer, optional, intent(in) :: maxeval
!    integer :: maxeval_local
!!-----------------------------Working variables---------------------------------
!    real(rk), allocatable :: wa(:)
!    integer, allocatable :: iwa(:)
!    character(len=60) :: task
!    character(len=60) :: csave
!    logical :: lsave(4)
!    integer :: isave(44)
!    real(rk) :: dsave(29)
!!-------------------------------------------------------------------------------
!    integer :: n ! number of variables
!    integer :: i
!    real(rk) :: func_local
!    real(rk), allocatable :: grad_local(:)
!!-------------------------AUGMENTED LAGRANGIAN PARAMETERS-----------------------
!    real(rk), parameter :: omega_star = 1.e-4_rk, eta_star = 1.e-4_rk, &
!      mu_incr = 1.05_rk
!    real(rk) :: mu, omega, eta
!    mu = 1.0_rk; omega = 1._rk / mu; eta = 1._rk / mu**0.1_rk
!!-------------------------------------------------------------------------------
!    n = size(x)
!    nce = size(ce)
!!------------------------Take care of optional arguments------------------------
!    factr_local = 1.e+7_rk  ! For moderate accuracy
!    pgtol_local = 1.e-5_rk
!    iprint_local = 0 ! Print only one line at the last iteration
!    allocate(nbd_local(n),source=2) ! Both lower and upper bounds
!    m_local = 7 ! Moderate number of corrections
!    maxeval_local = 99
!! Change non-default values
!    if (present(factr)) factr_local = factr
!    if (present(pgtol)) pgtol_local = pgtol
!    if (present(iprint)) iprint_local = iprint
!    if (present(nbd)) nbd_local = nbd
!    if (present(m)) m_local = m
!    if (present(maxeval)) maxeval_local = maxeval
!!-------------------------------------------------------------------------------
!    allocate(grad_local(n))
!    allocate(iwa(3*n))
!    allocate(wa(2*m_local*n + 5*n + 11*m_local*m_local + 8*m_local))
!!-------------------------------MAIN LOOP---------------------------------------
!    task = 'START'
!    do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
!      task.eq.'START')
!      ! Call to the original code
!      call setulb ( n, m_local, x, lb, ub, nbd_local, func_local, grad_local, &
!        factr_local, pgtol_local, wa, iwa, task, iprint_local, &
!        csave, lsave, isave, dsave )
!      ! Get objective function and gradient
!      if (task(1:2) .eq. 'FG') then
!        if (present(dat)) then
!          call eq_const(x,dat,ce)
!          func_local = func(x=x,dat=dat)
!          func_local = func_local - sum(lagr*ce) + mu / 2._rk * sum(ce**2)
!          !
!          call grad(x=x,dat=dat,grad=grad_local)
!          !grad_local = grad_local - lagr * grad_ce + mu * ce * grad_ce
!        else
!          func_local = func(x=x)
!          func_local = func_local - sum(lagr*ce) + mu / 2._rk * sum(ce**2)
!          !
!          call grad(x=x,grad=grad_local)
!        end if
!      else if (task(1:5) .eq. 'NEW_X') then
!!       The minimization routine has returned with a new iterate.
!!       At this point have the opportunity of stopping the iteration
!!       or observing the values of certain parameters
!
!!       1) Terminate if the total number of f and g evaluations
!!            exceeds maxeval.
!        if (isave(34) .ge. maxeval_local)  &
!          task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'
!!       2) Maybe additional stopping criteria, if needed
!
!        ! Print final results
!        write (stdout,'(2(a,i5,4x),a,1p,d12.5,4x,a,1p,d12.5)') 'Iterate', &
!          isave(30),'nfg =',isave(34),'f =',func_local,'|proj g| =',dsave(13)
!
!        if (task(1:4) .eq. 'STOP') then
!          write (stdout,*) task
!          write (stdout,*) 'Final X='
!          write (stdout,'((1x,1p, 6(1x,d11.4)))') (x(i),i = 1,n)
!        end if
!
!      end if
!    end do
!!-----------------------------END OF MAIN LOOP----------------------------------
!  end subroutine
!
!  subroutine l_bfgs_b_al_ieq_const(m,bris1,bris2)
!    integer, optional, intent(in) :: m
!    integer, intent(in) :: bris1,bris2
!
!  end subroutine


  subroutine l_bfgs_b_test(optim_prob,lbfgsb_opt)
    type(optim_prob_t), intent(inout) :: optim_prob
    type(lbfgsb_opt_t), intent(in) :: lbfgsb_opt
!-----------------------------Working variables---------------------------------
    real(rk), allocatable :: wa(:)
    integer, allocatable :: iwa(:)
    character(len=60) :: task
    character(len=60) :: csave
    logical :: lsave(4)
    integer :: isave(44)
    real(rk) :: dsave(29)
!-------------------------------------------------------------------------------
    integer :: i
    real(rk) :: f
    real(rk), allocatable :: g(:)
!-------------------------AUGMENTED LAGRANGIAN PARAMETERS-----------------------
    real(rk), parameter :: omega_star = 1.e-4_rk, eta_star = 1.e-4_rk, &
      mu_incr = 1.05_rk
    real(rk) :: mu, omega, eta
!---------------------------------CHECKS----------------------------------------
    if (.not.associated(optim_prob%func)) error stop 'Function not associated &
      &in optimization problem.'
    if (.not.associated(optim_prob%grad)) error stop 'Gradient not associated &
      &in optimization problem.'
!-------------------------------------------------------------------------------
    associate(n => optim_prob%n, m => lbfgsb_opt%m)
    allocate(g(n))
    allocate(iwa(3*n))
    allocate(wa(2*m*n + 5*n + 11*m*m + 8*m))
!-------------------------------------------------------------------------------
    if (lbfgsb_opt%iprint >= 0) then
      write(stdout,'(a,a)') 'Running the optimization problem ', optim_prob%name
    end if
!-------------------------------MAIN LOOP---------------------------------------
    task = 'START'
    do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
      task.eq.'START')
      ! Call to the original code
      call setulb (n, m, optim_prob%x, &
        optim_prob%lb, optim_prob%ub, optim_prob%nbd, &
        f, g, &
        lbfgsb_opt%factr, lbfgsb_opt%pgtol, &
        wa, iwa, task, &
        lbfgsb_opt%iprint, &
        csave, lsave, isave, dsave )
      if (task(1:2) .eq. 'FG') then
        if (allocated(optim_prob%dat)) then
          f = optim_prob%func(x=optim_prob%x,dat=optim_prob%dat)
          call optim_prob%grad(x=optim_prob%x,dat=optim_prob%dat,grad=g)
        else
          f = optim_prob%func(x=optim_prob%x)
          call optim_prob%grad(x=optim_prob%x,grad=g)
        end if





      else if (task(1:5) .eq. 'NEW_X') then
        if (isave(34) .ge. lbfgsb_opt%maxeval)  &
          task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'
        write (stdout,'(2(a,i5,4x),a,1p,d12.5,4x,a,1p,d12.5)') 'Iterate', &
          isave(30),'nfg =',isave(34),'f =',f,'|proj g| =',dsave(13)
        if (task(1:4) .eq. 'STOP') then
          write (stdout,*) task
          write (stdout,*) 'Final X='
          write (stdout,'((1x,1p, 6(1x,d11.4)))') (optim_prob%x(i),i = 1,n)
        end if
      end if
    end do
    end associate
!-----------------------------END OF MAIN LOOP----------------------------------
  end subroutine

end module
