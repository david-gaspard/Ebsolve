!!**************************************************************************************************
!! Created on 2024-06-04 at 16:18:03 CEST by David Gaspard <david.gaspard@espci.fr>
!! This program is distributed under the MIT License.
!! Fortran module to provide Gaussian quadrature utilities. Note that this module depends on LAPACK.
!!**************************************************************************************************
module gauss_quad
    use base_utils
	implicit none
	
    !! Constant strings:
    character(len=*), parameter :: module_name_gauss_quad = "GaussQuad"   !! Name of the present module.
    character(len=*), parameter :: errmsg_gauss_quad = module_name_gauss_quad // ": Invalid arguments, aborting..."
	
	contains
	
	!!**********************************************************************************************
	!! Subroutine to compute the nodes and weights of the Gauss-Jacobi quadrature.
    !! The purpose of this routine is to prepare for the integration of the function f(x) weighted by W(x):
	!! Pseudocode: integrate(f(x)*W(x), for x in [xmin, xmax]) ≃ sum(w_k*f(x_k), for k in {0, ..., npoints-1})
	!! The weight function of the Gauss-Jacobi quadrature reads: W(x) = (x - xmin)^a * (xmax - x)^b.
	!! Therefore, be aware that the exponent at "xmin" is "a" and the exponent at "xmax" is "b".
	!! Other conventions exist in the literature because the Jacobi polynomials are defined with the reverse exponents!
    !! 
    !! This subroutine is based on the Golub-Welsch matrix algorithm which consists in solving an eigenvalue problem
    !! such that the eigenvalues are the nodes of the Gaussian quadrature and the first component of the corresponding
    !! eigenvector gives the weight.
    !! For details about the Golub-Welsch algorithm, see: W. H. Press et al., Numerical Recipes, 3rd ed. (2007), page 188.
    !! 
    !! Arguments:
	!! np       = (IN) Desired number of points of the Gaussian quadrature.
	!! a        = (IN) Exponent of the weight at the point "xmin".
	!! b        = (IN) Exponent of the weight at the point "xmax".
	!! xmin     = (IN) Minimum abscissa of the integral.
	!! xmax     = (IN) Maximum abscissa of the integral.
    !! xpoints  = (OUT) Real array containing the nodes of the Gauss-Jacobi quadrature.
    !! weights  = (OUT) Real array containing the weights in the same order as "xpoints".
    !! 
    !! Licensing: This subroutine is distributed under the MIT License.
    !! Author: David Gaspard <david.gaspard@espci.fr>
    !! Created: June 4-5, 2024, based on a previous Python function written on September 19, 2023.
	!!**********************************************************************************************
	subroutine gauss_jacobi_quad(np, a, b, xmin, xmax, xpoints, weights)
        integer, intent(in) :: np
        real(wp), intent(in) :: a, b, xmin, xmax
        real(wp), intent(out) :: xpoints(:), weights(:)
        character, parameter :: compz = "I"   !! Tells LAPACK that, on exit, the matrix "Z" should contain the eigenvectors of the tridiagonal matrix.
        real(wp), allocatable :: diag(:), sub(:), z(:, :), work(:)  !! Diagonal, subdiagonal, orthogonal matrix, and LAPACK workspace.
        integer :: n, info   !! LAPACK info flag.
        
        external dsteqr  !! LAPACK subroutine to compute the eigenvalues and eigenvectors of a symmetric tridiagonal matrix.
        
        !! 1. First check for possible errors:
        if (np <= 0) then
            write (stderr, '(a,i0,a)') tag_error // module_name_gauss_quad // ": Invalid number of points (np=", &
                np, "), aborting..."
            stop errmsg_gauss_quad
        else if (xmin > xmax) then
            write (stderr, '(2(a,g0.3),a)') tag_error // module_name_gauss_quad // ": Invalid interval [xmin=", &
                xmin, ", xmax=", xmax, "], aborting..."
            stop errmsg_gauss_quad
        else if (a <= -1.) then
            write (stderr, '(a,g0.3,a)') tag_error // module_name_gauss_quad // ": Invalid exponent a = ", &
                a, ", the integral will not converge..."
            stop errmsg_gauss_quad
        else if (b <= -1.) then
            write (stderr, '(a,g0.3,a)') tag_error // module_name_gauss_quad // ": Invalid exponent b = ", &
                b, ", the integral will not converge..."
            stop errmsg_gauss_quad
        else if (size(xpoints) /= np .or. size(weights) /= np) then
            write (stderr, '(3(a,i0),a)') tag_error // module_name_gauss_quad // ": Invalid size of arrays xpoints(", &
                size(xpoints), ") and weights(", size(weights), "), expected size is ", np, ", aborting..."
            stop errmsg_gauss_quad
        end if
        
        !! 2. Builds the diagonal and subdiagonal of the symmetric tridiagonal matrix (aka the Jacobi matrix):
        allocate(diag(np))   !! Allocate space for the diagonal.
        allocate(sub(np-1))  !! Allocate space for the subdiagonal.
        allocate(z(np, np))  !! Allocate space for eigenvectors.
        allocate(work(max(1, 2*np - 2)))  !! Allocate space for the "work" array.
        
        do n = 1, np   !! Loop along the diagonal elements to fill "diag" and "sub".
            if (n == 1) then   !! An simplification occurs for the first elements. This exception is necessary to avoid division by zero in special cases.
                diag(n) = (a-b)/(a+b+2)
                sub(n) = 2.0_wp*sqrt(((a+1)*(b+1))/((a+b+2)*(a+b+2)*(a+b+3)))
            else
                diag(n) = ((a-b)*(a+b))/((a+b+2*n-2)*(a+b+2*n))  !! Diagonal elements.
                if (n < np) then
                    sub(n) = 2.0_wp*sqrt((n*(a+n)*(b+n)*(a+b+n))/((a+b+2*n-1)*(a+b+2*n)*(a+b+2*n)*(a+b+2*n+1)))  !! Subdiagonal elements.
                end if
            end if
        end do
        
        !! 3. Diagonalizes the Jacobi matrix and deduces the nodes and weights:
        if (np == 1) then
            z = 1.0_wp  !! In the case np=1, everything is trivial and LAPACK must not be called.
        else
            call dsteqr(compz, np, diag, sub, z, np, work, info)  !! Diagonalizes the Jacobi tridiagonal matrix (equivalent to scipy's eigh_tridiagonal).
            if (info < 0) then  !! Check the LAPACK status flag:
                write (stderr, '(a,i0,a,i0,a)') tag_error // module_name_gauss_quad // ": dsteqr(): The ", -info, "-th &
                    &argument had an illegal value (LAPACK info=", info, "), aborting..."
                stop errmsg_gauss_quad
            else if (info > 0) then
                write (stderr, '(2(a,i0),a)') tag_error // module_name_gauss_quad // ": dsteqr(): Failed to find &
                    &all the eigenvalues in a total of ", 30*np, " iterations (LAPACK info=", info, "), aborting..."
                stop errmsg_gauss_quad
            end if
        end if
        
        xpoints = 0.5_wp*((xmax-xmin)*diag + (xmax+xmin))
        weights = ((xmax-xmin)**(a+b+1)) * (gamma(a+1)*gamma(b+1)/gamma(a+b+2)) * z(1, :) * z(1, :)
        
        !! 4. Frees the allocated memory:
        deallocate(diag)
        deallocate(sub)
        deallocate(z)
        deallocate(work)
        
    end subroutine
    
end module
