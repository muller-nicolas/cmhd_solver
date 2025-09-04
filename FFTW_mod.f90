module fftw_wrapper
    use, intrinsic :: iso_c_binding
    implicit none

    ! FFTW constants
    integer, parameter :: FFTW_MEASURE = 0
    integer, parameter :: FFTW_ESTIMATE = 64

    interface
        ! Plan creation (double precision, R2C)
        function fftw_plan_dft_r2c_2d(n0, n1, in, out, flags) bind(C, name="fftw_plan_dft_r2c_2d")
            import :: c_int, c_ptr
            integer(c_int), value :: n0, n1, flags
            type(c_ptr), value :: in, out
            type(c_ptr) :: fftw_plan_dft_r2c_2d
        end function

        ! Plan creation (double precision, C2R)
        function fftw_plan_dft_c2r_2d(n0, n1, in, out, flags) bind(C, name="fftw_plan_dft_c2r_2d")
            import :: c_int, c_ptr
            integer(c_int), value :: n0, n1, flags
            type(c_ptr), value :: in, out
            type(c_ptr) :: fftw_plan_dft_c2r_2d
        end function

        ! Execution routines (generic, takes plan + arrays)
        subroutine fftw_execute_dft_r2c(plan, in, out) bind(C, name="fftw_execute_dft_r2c")
            import :: c_ptr
            type(c_ptr), value :: plan, in, out
        end subroutine

        subroutine fftw_execute_dft_c2r(plan, in, out) bind(C, name="fftw_execute_dft_c2r")
            import :: c_ptr
            type(c_ptr), value :: plan, in, out
        end subroutine

        ! Destroy plan
        subroutine fftw_destroy_plan(plan) bind(C, name="fftw_destroy_plan")
            import :: c_ptr
            type(c_ptr), value :: plan
        end subroutine
    end interface
end module fftw_wrapper

