! Test program for using customized calling subroutine, use
! the Rosenbrock function from driver 1, added constants

!
! CHANGES
! 21.12.2018 - Miha Polajnar
! Initial version
!
program custom1
  use iso_fortran_env, only: rk => real64
  use l_bfgs_b_ag, only: l_bfgs_b, optim_prob_t, lbfgsb_opt_t
  implicit none
  integer, parameter :: n = 25
  integer :: i
  type :: rconst_t
    real(rk) :: a, b
  end type
  type(optim_prob_t) :: optim_prob
  type(lbfgsb_opt_t) :: lbfgsb_opt
  type(rconst_t) :: rconst
  !
  lbfgsb_opt = lbfgsb_opt_t(iprint=1)
  rconst = rconst_t(a=4._rk,b=8._rk)
  !
  optim_prob%name = 'Rosenbrock'
  optim_prob%n = n
  optim_prob%x = [(3._rk,i=1,n)]
  optim_prob%nbd=[(2,i=1,n)]
  allocate(optim_prob%lb(n),optim_prob%ub(n))
  do i = 1, n, 2
    optim_prob%lb(i)   = 1.0_rk
    optim_prob%ub(i)   = 1.0e2_rk
  end do
  do i =2, n, 2
    optim_prob%lb(i)   = -1.0e2_rk
    optim_prob%ub(i)   =  1.0e2_rk
  end do
  optim_prob%func => rosenbrock_func
  optim_prob%grad => grad_rosenbrock_func
  allocate(optim_prob%dat,source=rconst)
  ! Run the algorithm
  call l_bfgs_b(optim_prob,lbfgsb_opt)
contains
  !
  pure function rosenbrock_func(x,dat) result(res)
    real(rk), intent(in) :: x(:)
    class(*), intent(in) :: dat
    real(rk) :: res
    integer :: i
    res = .25_rk * (x(1) - 1._rk)**2
    do i = 2, size(x)
      res = res + (x(i) - x(i - 1)**2)**2
    end do
    select type(dat)
    type is (rconst_t)
      res = dat%a * res
    class default
      error stop 'Rosenbrock function, incorrect type.'
    end select
  end function
  !
  pure subroutine grad_rosenbrock_func(x,dat,grad)
    real(rk), intent(in) :: x(:)
    class(*), intent(in) :: dat
    real(rk), intent(out) :: grad(:)
    integer :: i, n
    real(rk) :: t1, t2
    n = size(x)
    t1 = x(2) - x(1)**2
    grad(1) = 2._rk * (x(1) - 1._rk) - 1.6e1_rk * x(1) * t1
    do i = 2, n - 1
       t2   = t1
       t1   = x(i+1) - x(i)**2
       grad(i) = 8._rk * t2 - 1.6e1_rk * x(i) * t1
    end do
    select type (dat)
    type is (rconst_t)
      grad(n) = dat%b * t1
    class default
      error stop 'Rosenbrock function gradient, incorrect type.'
    end select
  end subroutine
  !
end program
