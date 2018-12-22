! Test program for using customized calling subroutine, use
! the Rosenbrock function from driver 1, added constants

!
! CHANGES
! 21.12.2018 - Miha Polajnar
! Initial version
!

program custom1
  use iso_fortran_env, only: rk => real64
  use l_bfgs_b_ag, only: l_bfgs_b
  implicit none
  integer,  parameter :: n = 25, m = 5, iprint = 1
  integer :: i
  real(rk), parameter :: factr  = 1.0e+7_rk, pgtol  = 1.0e-5_rk
  real(rk), allocatable  :: x(:), lb(:), ub(:)

  type :: rconst_t
    real(rk) :: a, b
  end type
  type(rconst_t) :: rconst

  rconst = rconst_t(a=4._rk,b=8._rk)

  allocate(x(n), lb(n), ub(n))

  do i = 1, n, 2
    lb(i)   = 1.0_rk
    ub(i)   = 1.0e2_rk
  end do

  do i =2, n, 2
    lb(i)   = -1.0e2_rk
    ub(i)   =  1.0e2_rk
  end do

  do i = 1, n
    x(i) = 3.0_rk
  end do

  call l_bfgs_b(x,lb,ub,rosenbrock_func,grad_rosenbrock_func,&
    iprint=iprint,m=m,dat=rconst)


contains

  pure function rosenbrock_func(x,dat) result(res)
    real(rk), intent(in) :: x(:)
    class(*), optional, intent(in) :: dat
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

  pure subroutine grad_rosenbrock_func(x,dat,grad)
    real(rk), intent(in) :: x(:)
    class(*), optional, intent(in) :: dat
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
    select type(dat)
    type is (rconst_t)
      grad(n) = dat%b * t1
    class default
      error stop 'Rosenbrock function gradient, incorrect type.'
    end select
  end subroutine

end program
