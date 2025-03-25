!!**************************************************************************************************
!! Created on 2024-07-21 at 12:51:13 CEST by David Gaspard <david.gaspard@espci.fr>
!! This program is distributed under the MIT License.
!! Fortran program to test the Eilenberger integrator.
!!**************************************************************************************************
program integrator_test
    use integrator
    use testing
    use base_utils
	implicit none
	
    real(wp), parameter :: tol = 1.0e-13_wp    !! Global tolerance per real (or complex) number.
    
	call test_apply_mexp()
    call test_ebintegrate("test/ebintegrate_sol_1.nml")
    call test_ebintegrate("test/ebintegrate_sol_2.nml")
    call test_ebintegrate("test/ebintegrate_sol_3.nml")
    call test_ebintegrate("test/ebintegrate_sol_4.nml")
    call test_ebintegrate("test/ebintegrate_sol_5.nml")
    call test_ebintegrate("test/ebintegrate_sol_6.nml")
	!!call test_ebintegrate_error()
    
	contains
	
	!!**********************************************************************************************
	!! Subroutine to test the matrix exponential transformation (hyperbolic rotation).
    !! This subroutine tests the apply_mexp() subroutine with various values of the parameters.
	!!**********************************************************************************************
	subroutine test_apply_mexp()
        complex(wp) :: argls(3, 7)
        integer :: i
        
        call print_test_title("apply_mexp()")
        
        argls = reshape([ &
            (0.0_wp, -2._wp), (3._wp, 1._wp), (0.0_wp, 0.0_wp), (1._wp, -4._wp), (2._wp, 0.5_wp), &
                (-1.5527686200792445_wp, 2.3465604163078308_wp), (-1.2869423865071256_wp, 1.6105214353777922_wp), &
            (1._wp, -3._wp), (0._wp, 2._wp), (-2._wp, 1._wp), (-0.5_wp, 1.5_wp), (2.5_wp, -3.5_wp), &
                (-0.8894938643351842_wp, -10.7115121264736_wp), (-2.852841963244258_wp, 3.3061366751200874_wp), &
            (0.0_wp, 0.0_wp), (0.0_wp, 2.), (0.0_wp, 0.0_wp), (-1._wp, 2._wp), (1._wp, 0.5_wp), (-2._wp, 4._wp), (1._wp, 0.5_wp) &
        ], shape(argls), order=[2,1])
        
        do i = 1, size(argls, 1)
            
            call apply_mexp(argls(i,1), argls(i,2), argls(i,3), argls(i,4), argls(i,5))
            
            if (assert_equal_complex(argls(i,4), argls(i,6), tol) .and. assert_equal_complex(argls(i,5), argls(i,7), tol)) then
                print '(a,i0,a,g0.2,a)', tag_pass // "#", i, " | Expected result within tolerance limits (tol=", tol, ")."
            else 
                print '(a,i0,a)', tag_fail // "#", i, " | Significant differences detected."
            end if
            
        end do
        
	end subroutine
    
	!!**********************************************************************************************
    !! Test the solution of the Eilenberger equation given by the ebintegrate() subroutine with a
    !! reference solution found in the "test/" folder.
	!!**********************************************************************************************
    subroutine test_ebintegrate(filename)
        character(len=*), intent(in) :: filename
        integer :: sgn, fu, openstat, readstat
        integer, parameter :: nx = 7
        complex(wp) :: mu
        real(wp) :: xmesh(nx), xa, xb, normtol, normerr
        complex(wp) :: pmesh(3, nx-1), ca(3), cb(3), amesh(2, nx), bmesh(2, nx), gmesh(3, nx)
        complex(wp) :: amesh_expc(2, nx), bmesh_expc(2, nx), gmesh_expc(3, nx)
        
        !! 1. Read the benchmark solution from a namelist (external file):
		namelist /ebintegrate_sol/ sgn, mu, xmesh, pmesh, xa, ca, xb, cb, amesh_expc, bmesh_expc, gmesh_expc
		
        call print_test_title("ebintegrate() [" // trim(filename) // "]")
        
        open(newunit=fu, file=filename, action='read', iostat=openstat)
        read(unit=fu, nml=ebintegrate_sol, iostat=readstat)
        close(fu)
        
        if (openstat /= 0) then 
			print '(a)', tag_error // "File not found: '" // trim(filename) // "'"
			stop
		else if (readstat /= 0) then
			print '(a)', tag_error // "Invalid namelist format in file: '" // trim(filename) // "'"
			stop
		end if
		
        !! 2. Compute the solution predicted by ebintegrate():
        call ebintegrate(sgn, mu, xmesh, pmesh, xa, ca, xb, cb, amesh, bmesh, gmesh)
        
        !! 3. Compare with the benchmark solution:
		normtol = tol*nx  !! Total tolerance for the norm.
        normerr = norm_complex_2darray(amesh - amesh_expc)
		if (normerr > normtol) then
			print '(2(a,g0.2),a)', tag_fail // "Errors found in a(x). Error=", normerr, " (> ", normtol, ")"
            call print_complex_2darray("amesh", "g0.6", amesh)
            call print_complex_2darray("amesh_expc", "g0.6", amesh_expc)
        else
            print '(2(a,g0.2),a)', tag_pass // "a(x) is correct. Error=", normerr, " (< ", normtol, ")"
		end if
        
        normerr = norm_complex_2darray(bmesh - bmesh_expc)
        if (normerr > normtol) then
			print '(2(a,g0.2),a)', tag_fail // "Errors found in b(x). Error=", normerr, " (> ", normtol, ")"
            call print_complex_2darray("bmesh", "g0.6", bmesh)
            call print_complex_2darray("bmesh_expc", "g0.6", bmesh_expc)
        else
            print '(2(a,g0.2),a)', tag_pass // "b(x) is correct. Error=", normerr, " (< ", normtol, ")"
		end if
        
        normerr = norm_complex_2darray(gmesh - gmesh_expc)
        if (normerr > normtol) then
			print '(2(a,g0.2),a)', tag_fail // "Errors found in g(x). Error=", normerr, " (> ", normtol, ")"
            call print_complex_2darray("gmesh", "g0.6", gmesh)
            call print_complex_2darray("gmesh_expc", "g0.6", gmesh_expc)
        else
            print '(2(a,g0.2),a)', tag_pass // "g(x) is correct. Error=", normerr, " (< ", normtol, ")"
		end if
        
    end subroutine
    
	!!**********************************************************************************************
    !! Test the error detection of the ebintegrate_test() subroutine in a situation it is expected to trigger.
	!!**********************************************************************************************
    subroutine test_ebintegrate_error()
        integer, parameter :: sgn = +1, nx = 7
        complex(wp), parameter :: mu = (1.0_wp, 0.0_wp)
        real(wp), parameter :: xmesh(nx) = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
        real(wp), parameter :: xa = 1.5_wp, xb = 6.5_wp
        complex(wp), parameter :: ca(3) = [0.0, 1.0, 0.0], cb(3) = [0.0, 0.0, 1.0]
        complex(wp), allocatable :: pmesh(:, :) 
        complex(wp) :: amesh(2, nx), bmesh(2, nx), gmesh(3, nx)
        integer :: ix
        
        allocate(pmesh(3, nx))  !! <- Invalid length. Expected dimensions: pmesh(3, nx-1).
        do ix = 1, nx
            pmesh(:, ix) = [1.0, 2.0, 3.0]
        end do
        
        call print_test_title("ebintegrate() - Error detection")
        
        call print_real_1darray("xmesh", "g0.2", xmesh)
        call print_complex_2darray("pmesh", "g0.2", pmesh)
        
        call ebintegrate_test(sgn, mu, xmesh, pmesh, xa, ca, xb, cb, amesh, bmesh, gmesh)
        
        deallocate(pmesh)
        
    end subroutine
    
end program
