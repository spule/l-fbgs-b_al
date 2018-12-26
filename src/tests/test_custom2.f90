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
  integer, parameter :: n = 2
  integer :: i
  type :: rconst_t
    real(rk) :: a, b
  end type
  type(optim_prob_t) :: optim_prob
  type(lbfgsb_opt_t) :: lbfgsb_opt
  !
  lbfgsb_opt = lbfgsb_opt_t(iprint=-200,factr=1.e7_rk,pgtol=1.e-5_rk,&
    maxeval=500)
  !
  optim_prob%name = 'NOCEDAL EX. 17.1'
  optim_prob%n = n
  optim_prob%x = [(0._rk,i=1,n)]
  optim_prob%nbd=[(2,i=1,n)]
  allocate(optim_prob%lb(n),optim_prob%ub(n))
  optim_prob%lb = [(-5.12_rk,i=1,n)]
  optim_prob%ub = [(5.12_rk,i=1,n)]
  optim_prob%func => nocedal_ex_17_1
  optim_prob%grad => grad_nocedal_ex_17_1
  optim_prob%neqc = 1
  optim_prob%eq_const => nocedal_ex_17_1_eq_const
  optim_prob%grad_eq_const => grad_nocedal_ex_17_1_eq_const
  allocate(optim_prob%lagr(optim_prob%neqc),source=10._rk)
  allocate(optim_prob%dat,source=2._rk)
  ! Run the algorithm
  call l_bfgs_b(optim_prob,lbfgsb_opt)
contains

  real(rk) pure function nocedal_ex_17_1(x,dat)
    real(rk), intent(in) :: x(:)
    class(*), intent(in) :: dat
    nocedal_ex_17_1 = x(1) + x(2)
  end function

  pure subroutine grad_nocedal_ex_17_1(x,dat,grad)
    real(rk), intent(in) :: x(:)
    class(*), intent(in) :: dat
    real(rk), intent(out) :: grad(:)
    grad(1) = 1._rk
    grad(2) = 1._rk
  end subroutine

  pure subroutine nocedal_ex_17_1_eq_const(x,dat,ce)
    real(rk), intent(in) :: x(:)
    class(*), intent(in) :: dat
    real(rk), intent(out) :: ce(:)
    select type(dat)
    type is (real(rk))
      ce(1) = x(1)**2 + x(2)**2 - dat
    class default
      error stop 'Nocedal EX. 17.1 wrong type.'
    end select
  end subroutine

  pure subroutine grad_nocedal_ex_17_1_eq_const(x,dat,grad)
    real(rk), intent(in) :: x(:)
    class(*), intent(in) :: dat
    real(rk), intent(out) :: grad(:,:)
    grad(1,1) = 2._rk * x(1)
    grad(1,2) = 2._rk * x(2)
  end subroutine
end program
