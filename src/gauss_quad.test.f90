!!**************************************************************************************************
!! Created on 2024-06-04 at 18:15:13 CEST by David Gaspard <gaspard.dave@gmail.com>
!! This program is distributed under the Creative Commons (CC) BY-NC-SA license.
!! Fortran program to test the Gaussian quadrature utilities.
!!**************************************************************************************************
program gauss_quad_test
    use gauss_quad
	implicit none
    
	!! Useful mathematical constants:
    real(wp), parameter :: pi = acos(-1.0_wp)    !! Pi constant at working precision.
    
    call test_gauss_jacobi_quad_1()
    call test_gauss_jacobi_quad_2a()
    call test_gauss_jacobi_quad_2b()
	call test_gauss_jacobi_quad_3()
	
    !! Maybe also test the integration of functions such as: f(x) = sin(9*x)/(x*sqrt(1-x))
    
	contains
    
    !!**********************************************************************************************
    !! Test the first few moments of the Gauss-Jacobi quadrature for various values of the parameters.
    !!**********************************************************************************************
    subroutine test_gauss_jacobi_quad_1()
        call test_gauss_jacobi_quad_mom(12,  0.0_wp, -0.5_wp,  0.0_wp, 1.0_wp)
        call test_gauss_jacobi_quad_mom(12, -0.6_wp, -0.8_wp, -0.2_wp, 0.9_wp)
        call test_gauss_jacobi_quad_mom(12, -0.9_wp,  0.4_wp, -0.7_wp, 1.6_wp)
    end subroutine
	
    !!**********************************************************************************************
    !! Test the first few moments of the Gauss-Jacobi quadrature for the given parameters.
    !!**********************************************************************************************
    subroutine test_gauss_jacobi_quad_mom(np, a, b, xmin, xmax)
        integer, intent(in) :: np
        real(wp), intent(in) :: a, b, xmin, xmax
        real(wp), allocatable :: xpoints(:), weights(:)
        real(wp) :: sumw, sumw_ex, sumwx, sumwx_ex, sumwx2, sumwx2_ex
        
        print '(a)', "====== TEST Gauss-Jacobi quadrature: First moments ======"
        print '(a,i0,4(a,g0.3))', "[TEST] Parameters: np=", np, ", a=", a, ", b=", b, ", xmin=", xmin, ", xmax=", xmax
        
        allocate(xpoints(np))
        allocate(weights(np))
        
        call gauss_jacobi_quad(np, a, b, xmin, xmax, xpoints, weights)
        
        sumw = sum(weights)
        sumw_ex = ((xmax-xmin)**(a+b+1)) * (gamma(a+1)*gamma(b+1)/gamma(a+b+2))
        print '(3(a,g13.6))', "[TEST] sum_k w_k       =", sumw, &
            "expected =", sumw_ex, "relative error =", abs((sumw - sumw_ex)/sumw_ex)
        
        sumwx = sum(weights*xpoints)
        sumwx_ex = sumw_ex * ((a+1)*xmax + (b+1)*xmin)/(a+b+2)
        print '(3(a,g13.6))', "[TEST] sum_k w_k x_k   =", sumwx, &
            "expected =", sumwx_ex, "relative error =", abs((sumwx - sumwx_ex)/sumwx_ex)
        
        sumwx2 = sum(weights*xpoints*xpoints)
        sumwx2_ex = sumw_ex * ((a+1)*(a+2)*xmax**2 + 2*(a+1)*(b+1)*xmin*xmax + (b+1)*(b+2)*xmin**2)/((a+b+2)*(a+b+3))
        print '(3(a,g13.6))', "[TEST] sum_k w_k x_k^2 =", sumwx2, &
            "expected =", sumwx2_ex, "relative error =", abs((sumwx2 - sumwx2_ex)/sumwx2_ex)
        
        deallocate(xpoints)
        deallocate(weights)
        
    end subroutine
    
	!!**********************************************************************************************
	!! Test the Gauss-Jacobi quadrature in the Chebyshev case (first kind) using exactly known values.
	!!**********************************************************************************************
	subroutine test_gauss_jacobi_quad_2a()
        real(wp), parameter :: xmin = -1.0_wp, xmax = 1.0_wp, a = -0.5_wp, b = -0.5_wp
        integer, parameter :: np = 10
        real(wp) :: xpoints(np), weights(np), xpoints_ex(np), weights_ex(np), xpoints_err, weights_err
        integer :: i
        
        !! Expected results (analytically known in this case):
        do i = 1, np
            xpoints_ex(i) = -cos(pi*(i - 0.5_wp)/np)
            weights_ex(i) = pi/np
        end do
        
        print '(a)', "====== TEST Gauss-Jacobi quadrature: Chebyshev 1st kind ======"
        
        call gauss_jacobi_quad(np, a, b, xmin, xmax, xpoints, weights)
        xpoints_err = norm2(xpoints - xpoints_ex)
        weights_err = norm2(weights - weights_ex)
        
        !! Print the errors:
        print '(a,g0.3)', "[TEST] xpoints_err = ", xpoints_err
        print '(a,g0.3)', "[TEST] weights_err = ", weights_err
        
	end subroutine
    
	!!**********************************************************************************************
	!! Test the Gauss-Jacobi quadrature in the Chebyshev case (second kind) using exactly known values.
	!!**********************************************************************************************
	subroutine test_gauss_jacobi_quad_2b()
        real(wp), parameter :: xmin = -1.0_wp, xmax = 1.0_wp, a = 0.5_wp, b = 0.5_wp
        integer, parameter :: np = 10
        real(wp) :: xpoints(np), weights(np), xpoints_ex(np), weights_ex(np), xpoints_err, weights_err
        integer :: i
        
        !! Expected results (analytically known in this case):
        do i = 1, np
            xpoints_ex(i) = -cos(i*pi/(np + 1.0_wp))
            weights_ex(i) = (pi/(np + 1.0_wp)) * sin(i*pi/(np + 1.0_wp))**2
        end do
        
        print '(a)', "====== TEST Gauss-Jacobi quadrature: Chebyshev 2nd kind ======"
        
        call gauss_jacobi_quad(np, a, b, xmin, xmax, xpoints, weights)
        xpoints_err = norm2(xpoints - xpoints_ex)
        weights_err = norm2(weights - weights_ex)
        
        !! Print the errors:
        print '(a,g0.3)', "[TEST] xpoints_err = ", xpoints_err
        print '(a,g0.3)', "[TEST] weights_err = ", weights_err
        
	end subroutine
    
	!!**********************************************************************************************
	!! Test the Gauss-Jacobi quadrature in the Legendre case using tabulated values.
	!!**********************************************************************************************
	subroutine test_gauss_jacobi_quad_3()
        real(wp), parameter :: xmin = -1.0_wp, xmax = 1.0_wp, a = 0.0_wp, b = 0.0_wp
        integer, parameter :: np = 10
        real(wp) :: xpoints(np), weights(np), xpoints_err, weights_err
        
        !! Expected results (see: https://pomax.github.io/bezierinfo/legendre-gauss.html#n10):
        real(wp), parameter :: xpoints_ex(np) = [-0.9739065285171717_wp, -0.8650633666889845_wp,     &
            -0.6794095682990244_wp, -0.4333953941292472_wp, -0.1488743389816312_wp, 0.1488743389816312_wp, &
            0.4333953941292472_wp, 0.6794095682990244_wp, 0.8650633666889845_wp, 0.9739065285171717_wp]
        real(wp), parameter :: weights_ex(np) = [0.0666713443086881_wp, 0.1494513491505806_wp,    &
            0.2190863625159820_wp, 0.2692667193099963_wp, 0.2955242247147529_wp, 0.2955242247147529_wp, &
            0.2692667193099963_wp, 0.2190863625159820_wp, 0.1494513491505806_wp, 0.0666713443086881_wp]
        
        print '(a)', "====== TEST Gauss-Jacobi quadrature: Legendre case ======"
        
        call gauss_jacobi_quad(np, a, b, xmin, xmax, xpoints, weights)
        xpoints_err = norm2(xpoints - xpoints_ex)
        weights_err = norm2(weights - weights_ex)
        
        !! Print the errors:
        print '(a,g0.3)', "[TEST] xpoints_err = ", xpoints_err
        print '(a,g0.3)', "[TEST] weights_err = ", weights_err
        
	end subroutine
	
end program
