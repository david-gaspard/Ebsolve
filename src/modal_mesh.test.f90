!!**************************************************************************************************
!! Created on 2024-07-26 at 17:33:00 CEST by David Gaspard <david.gaspard@espci.fr>
!! Based on the same module written on 2024-01-23 at 15:20:17 CET by David Gaspard for the Ebsolve program (v1).
!! This program is distributed under the MIT License.
!! Fortran program to test the modal_mesh module which is supposed to generate the modal lattices "mumesh" and "cmesh".
!!**************************************************************************************************
program modal_mesh_test
    use modal_mesh
    use testing
	implicit none
	real(wp), parameter :: tol = 1.0e-14_wp  !! General tolerance on all floating point operations.
    
    !! === TEST PERIODIC BOUNDARY CONDITIONS ===
    
    !! Test the internal lattice distance subroutine:
    call test_distance_square_lattice_periodic_1d()
    call test_distance_square_lattice_periodic_2d()
    
    !! Test the mumesh and cmesh for a waveguide with periodic boundary conditions:
	call test_fill_mumesh_periodic_2d()
	call test_fill_mumesh_periodic_3d()
    
    !! === TEST INFINITE SLAB CASE ===
    
    call test_fill_mumesh_infinite()
    
	contains
    
    !! ========================== TEST PERIODIC BOUNDARY CONDITIONS ================================
    
	!!**********************************************************************************************
	!! Subroutine to test the generation of distances on a square lattice in 1D.
	!!**********************************************************************************************
	subroutine test_distance_square_lattice_periodic_1d()
		integer, parameter :: dlat = 1
		real(wp), parameter :: rmax = 3.01_wp  !! Close to an integer, but just above.
		integer, parameter :: ndist_expc = 7   !! Expected number of distances.
		real(wp), parameter :: dists_expc(ndist_expc) = [3._wp, 2._wp, 1._wp, 0._wp, 1._wp, 2._wp, 3._wp]
		real(wp), allocatable :: dists(:)
		integer :: ndist
		
		write (stdout, '(a)') tag_test // "====== TEST 1D LATTICE (PERIODIC) ======"
		
		call distance_square_lattice_periodic(dlat, rmax, ndist, dists)
		
		if (assert_equal_real_1darray(dists(1:ndist_expc), dists_expc, tol)) then
			print '(a)', tag_pass // "OK, dists is correct within tolerance limits."
		else 
			print '(a)', tag_fail // "Errors detected in dists."
            call print_real_1darray("Distance array", "g0.2", dists)
            call print_real_1darray("Expected array", "g0.2", dists_expc)
		end if
		
	end subroutine
    
	!!**********************************************************************************************
	!! Subroutine to test the generation of distances on a square lattice in 2D.
	!!**********************************************************************************************
	subroutine test_distance_square_lattice_periodic_2d()
		integer, parameter :: dlat = 2
		real(wp), parameter :: rmax = 4.48_wp  !! Close to sqrt(4^2 + 2^2) = 2*sqrt(5) = 4.47214..., but just above.
		integer, parameter :: ndist_expc = 69   !! Expected number of distances.
		real(wp), parameter :: dists_expc(ndist_expc) = [ &
				4.472135954999579_wp, 4.123105625617661_wp, 4.000000000000000_wp, 4.123105625617661_wp,  &
				4.472135954999579_wp, 4.242640687119285_wp, 3.605551275463989_wp, 3.162277660168379_wp,  &
				3.000000000000000_wp, 3.162277660168379_wp, 3.605551275463989_wp, 4.242640687119285_wp,  & 
				4.472135954999579_wp, 3.605551275463989_wp, 2.828427124746190_wp, 2.236067977499790_wp,  &
				2.000000000000000_wp, 2.236067977499790_wp, 2.828427124746190_wp, 3.605551275463989_wp,  &
				4.472135954999579_wp, 4.123105625617661_wp, 3.162277660168379_wp, 2.236067977499790_wp,  &
				1.414213562373095_wp, 1.000000000000000_wp, 1.414213562373095_wp, 2.236067977499790_wp,  &
				3.162277660168379_wp, 4.123105625617661_wp, 4.000000000000000_wp, 3.000000000000000_wp,  &
				2.000000000000000_wp, 1.000000000000000_wp, 0.000000000000000_wp,                        &
				1.000000000000000_wp, 2.000000000000000_wp,                                              &
				3.000000000000000_wp, 4.000000000000000_wp, 4.123105625617661_wp, 3.162277660168379_wp,  &
				2.236067977499790_wp, 1.414213562373095_wp, 1.000000000000000_wp, 1.414213562373095_wp,  &
				2.236067977499790_wp, 3.162277660168379_wp, 4.123105625617661_wp, 4.472135954999579_wp,  &
				3.605551275463989_wp, 2.828427124746190_wp, 2.236067977499790_wp, 2.000000000000000_wp,  &
				2.236067977499790_wp, 2.828427124746190_wp, 3.605551275463989_wp, 4.472135954999579_wp,  &
				4.242640687119285_wp, 3.605551275463989_wp, 3.162277660168379_wp, 3.000000000000000_wp,  &
				3.162277660168379_wp, 3.605551275463989_wp, 4.242640687119285_wp, 4.472135954999579_wp,  &
				4.123105625617661_wp, 4.000000000000000_wp, 4.123105625617661_wp, 4.472135954999579_wp ]
		
		real(wp), allocatable :: dists(:)
		integer :: ndist
		
		write (stdout, '(a)') tag_test // "====== TEST 2D LATTICE (PERIODIC) ======"
		
		call distance_square_lattice_periodic(dlat, rmax, ndist, dists)
		
		!!print *, "Distance array = ", dists
		
		if (assert_equal_real_1darray(dists(1:ndist_expc), dists_expc, tol)) then
			print '(a)', tag_pass // "OK, dists is correct within tolerance limits."
		else 
			print '(a)', tag_fail // "Errors detected in dists."
            call print_real_1darray("Distance array", "g0.2", dists)
            call print_real_1darray("Expected array", "g0.2", dists_expc)
		end if
		
	end subroutine
    
	!!**********************************************************************************************
	!! Subroutine to test the "muspace" object of "periodic" type in 2D.
	!!**********************************************************************************************
	subroutine test_fill_mumesh_periodic_2d()
        integer, parameter :: d = 2
		real(wp), parameter :: wol = 5.90_wp
		integer, parameter :: nmu_expc = 6
		complex(wp), parameter :: mumesh_expc(nmu_expc) = [ &
				(1._wp,0._wp),(0.9855316447529919_wp,0._wp),(0.9407924804324012_wp,0._wp), &
				(0.8610770031105449_wp,0._wp),(0.7350931675322523_wp,0._wp),(0.5308630428259602_wp,0._wp) ]
		complex(wp), parameter :: cmesh_expc(nmu_expc) = [ &
				(1.,0.),(2.,0.),(2.,0.),(2.,0.),(2.,0.),(2.,0.) ]
        
        complex(wp), allocatable :: mumesh(:), cmesh(:)  !! Gets the "mumesh" and "cmesh" from "muspace".
        
        !! 1. Initialize the muspace object and extract the arrays:
        print '(a)', tag_test // "====== MUSPACE PERIODIC 2D ======"
        call fill_mumesh_periodic(d, wol, mumesh, cmesh)
        
        !! 2. Compare with the expected arrays:
		if (assert_equal_complex_1darray(mumesh, mumesh_expc, tol)) then
			print '(a)', tag_pass // "OK, mumesh is correct within tolerance limits."
		else 
			print '(a)', tag_fail // "Errors detected in mumesh."
            call print_complex_1darray("Mumesh array",   "g0.2", mumesh)
            call print_complex_1darray("Expected array", "g0.2", mumesh_expc)
		end if
		if (assert_equal_complex_1darray(cmesh, cmesh_expc, tol)) then
			print '(a)', tag_pass // "OK, cmesh is correct within tolerance limits."
		else 
			print '(a)', tag_fail // "Errors detected in cmesh."
            call print_complex_1darray("Cmesh array",    "g0.2", cmesh)
            call print_complex_1darray("Expected array", "g0.2", cmesh_expc)
		end if
        
	end subroutine
    
	!!**********************************************************************************************
	!! Subroutine to test the "muspace" object of "periodic" type in 3D.
	!!**********************************************************************************************
	subroutine test_fill_mumesh_periodic_3d()
        integer, parameter :: d = 3
		real(wp), parameter :: wol = 5.10_wp
		integer, parameter :: nmu_expc = 15
		complex(wp), parameter :: mumesh_expc(nmu_expc) = [ &
				(1._wp,0.),(0.980588215690195_wp,0.),(0.9607843137254901_wp,0.),(0.919898361234502_wp,0.), &
				(0.898758167558105_wp,0.),(0.8321213793695271_wp,0.),(0.8086898285216187_wp,0.),(0.7845587852448059_wp,0.), &
				(0.7072426979165267_wp,0.),(0.6203643929237792_wp,0.),(0.588562000776613_wp,0.),(0.5549400665915645_wp,0.), &
				(0.48069218322083346_wp,0.),(0.1970563847278602_wp,0.),(0.019607843137252536_wp,0.) ]
		complex(wp), parameter :: cmesh_expc(nmu_expc) = [ &
				(1.,0.),(4.,0.),(4.,0.),(4.,0.),(8.,0.),(4.,0.),(4.,0.),(8.,0.),(8.,0.),(4.,0.), &
				(8.,0.),(4.,0.),(8.,0.),(12.,0.),(8.,0.) ]
		
		complex(wp), allocatable :: mumesh(:), cmesh(:)
        
        !! 1. Initialize the muspace object and extract the arrays:
        print '(a)', tag_test // "====== MUSPACE PERIODIC 3D ======"
        call fill_mumesh_periodic(d, wol, mumesh, cmesh)
        
        !! 2. Compare with the expected arrays:
		if (assert_equal_complex_1darray(mumesh, mumesh_expc, tol)) then
			print '(a)', tag_pass // "OK, mumesh is correct within tolerance limits."
		else 
			print '(a)', tag_fail // "Errors detected in mumesh."
            call print_complex_1darray("Mumesh array",   "g0.2", mumesh)
            call print_complex_1darray("Expected array", "g0.2", mumesh_expc)
		end if
		if (assert_equal_complex_1darray(cmesh, cmesh_expc, tol)) then
			print '(a)', tag_pass // "OK, cmesh is correct within tolerance limits."
		else 
			print '(a)', tag_fail // "Errors detected in cmesh."
            call print_complex_1darray("Cmesh array",    "g0.2", cmesh)
            call print_complex_1darray("Expected array", "g0.2", cmesh_expc)
		end if
        
    end subroutine
    
    !! ========================== TEST INFINITE SLAB CASE ================================
    
    !!**********************************************************************************************
    !! Test the muspace object of "infinite" type in arbitrary dimension "d" using the first moments
    !! which are known in closed form for the continuum limit (nmu -> infinity).
    !!**********************************************************************************************
    subroutine test_fill_mumesh_infinite()
        call test_fill_mumesh_infinite_mom(2, 20, 1.0_wp)
        call test_fill_mumesh_infinite_mom(3, 12, 1.0_wp)
    end subroutine
    
    !!**********************************************************************************************
    !! Test the mumesh/cmesh arrays of "infinite" type in arbitrary dimension "d" using the first moments
    !! which are known in closed form for the continuum limit (nmu -> infinity).
    !!**********************************************************************************************
    subroutine test_fill_mumesh_infinite_mom(d, nmu, ash)
        integer, intent(in) :: d, nmu
        real(wp), intent(in) :: ash
        complex(wp), allocatable :: mumesh(:), cmesh(:)
        complex(wp) :: sumc, sumc_ex, sumcdmu, sumcdmu_ex, sumcmu, sumcmu_ex, sumcmu2, sumcmu2_ex
        
        !! 1. Compute the mumesh/cmesh arrays:
        print '(a,i0,a)', tag_test // "====== MUSPACE INFINITE ", d, "D ======"
        print '(2(a,i0),a,g0.3)', tag_test // "Parameters: d=", d, ", nmu=", nmu, ", ash=", ash
        
        call fill_mumesh_infinite(d, nmu, ash, mumesh, cmesh)
        
        !! 2. Test the first moments of the coefficients:
        sumc = sum(cmesh)
        sumc_ex = 2.0_wp/(d - 1.0_wp)
        print '(3(a,g0.6,sp,g0.6,ss,"i"))', tag_test // "Sum_i c_i        = ", sumc, &
            ", expected = ", sumc_ex, ", relative error = ", abs((sumc - sumc_ex)/sumc_ex)
        
        sumcdmu = sum(cmesh/mumesh)
        sumcdmu_ex = sqrt(pi)*gamma((d-1.0_wp)/2.0_wp)/gamma(d/2.0_wp)
        print '(3(a,g0.6,sp,g0.6,ss,"i"))', tag_test // "Sum_i c_i/mu_i   = ", sumcdmu, &
            ", expected = ", sumcdmu_ex, ", relative error = ", abs((sumcdmu - sumcdmu_ex)/sumcdmu_ex)
        
        sumcmu = sum(cmesh*mumesh)
        sumcmu_ex = sqrt(pi)*gamma((d-1.0_wp)/2.0_wp)/(2.0_wp*gamma(d/2.0_wp + 1.0_wp))
        print '(3(a,g0.6,sp,g0.6,ss,"i"))', tag_test // "Sum_i c_i*mu_i   = ", sumcmu, &
            ", expected = ", sumcmu_ex, ", relative error = ", abs((sumcmu - sumcmu_ex)/sumcmu_ex)
        
        sumcmu2 = sum(cmesh*mumesh*mumesh)
        sumcmu2_ex = 4.0_wp/((d - 1)*(d + 1))
        print '(3(a,g0.6,sp,g0.6,ss,"i"))', tag_test // "Sum_i c_i*mu_i^2 = ", sumcmu2, &
            ", expected = ", sumcmu2_ex, ", relative error = ", abs((sumcmu2 - sumcmu2_ex)/sumcmu2_ex)
        
    end subroutine
    
end program
