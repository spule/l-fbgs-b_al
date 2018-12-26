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
  integer, parameter :: n = 2 + 2
  type :: rconst_t
    real(rk) :: a, b
  end type
  type(optim_prob_t) :: optim_prob
  type(lbfgsb_opt_t) :: lbfgsb_opt
  !
  lbfgsb_opt = lbfgsb_opt_t(iprint=200,factr=1.e7_rk,pgtol=0._rk,&
    maxeval=500,m=5)
  !
  optim_prob%name = 'ROSENBROCK - CONSTRANED WITH CUBIC AND LINE'
  optim_prob%n = n
  optim_prob%x = [-1._rk, 2._rk, 1._rk, 1._rk]
  optim_prob%nbd=[ 2, 2, 1, 1 ]
  optim_prob%lb = [ -1.5_rk, -0.5_rk, 0._rk, 0._rk ]
  optim_prob%ub = [  1.5_rk,  2.5_rk, 0._rk, 0._rk ]
  optim_prob%func => rosenbrock
  optim_prob%grad => grad_rosenbrock
  optim_prob%neqc = 2
  optim_prob%eq_const => rosenbrock_eq_const
  optim_prob%grad_eq_const => grad_rosenbrock_eq_const
  ! Run the algorithm
  call l_bfgs_b(optim_prob,lbfgsb_opt)
contains

  real(rk) pure function rosenbrock(x,dat)
    real(rk), intent(in) :: x(:)
    class(*), intent(in) :: dat
    rosenbrock = (1._rk - x(1))**2 + 100._rk *(x(2) - x(1))**2
  end function

  pure subroutine grad_rosenbrock(x,dat,grad)
    real(rk), intent(in) :: x(:)
    class(*), intent(in) :: dat
    real(rk), intent(out) :: grad(:)
    grad(1) = -2._rk * (1._rk - x(1)) - 200._rk * (x(2) - x(1))
    grad(2) = 200._rk * (x(2) - x(1))
    grad(3) = 0._rk
    grad(4) = 0._rk
  end subroutine

  pure subroutine rosenbrock_eq_const(x,dat,ce)
    real(rk), intent(in) :: x(:)
    class(*), intent(in) :: dat
    real(rk), intent(out) :: ce(:)
    ce(1) = -((x(1) - 1._rk)**3 -x(2) + 1._rk) - x(3)
    ce(2) = -(x(1) + x(2) - 2._rk) - x(4)
  end subroutine

  pure subroutine grad_rosenbrock_eq_const(x,dat,grad)
    real(rk), intent(in) :: x(:)
    class(*), intent(in) :: dat
    real(rk), intent(out) :: grad(:,:)
    grad(1,1) = -3._rk * (x(1) - 1._rk)**2
    grad(1,2) = 1._rk
    grad(1,3) = -1._rk
    grad(1,4) = 0._rk
    grad(2,1) = -1._rk
    grad(2,2) = -1._rk
    grad(2,3) = 0._rk
    grad(2,4) = -1._rk
  end subroutine
end program
