!!**************************************************************************************************
!! Created on 2024-07-26 at 16:42:18 CEST by David Gaspard <david.gaspard@espci.fr>
!! This program is distributed under the MIT License.
!! Fortran module to generate the modal lattice of the Ebsolve program (v2).
!!**************************************************************************************************
module modal_mesh
    use gauss_quad
    implicit none
    
    contains
    
    !!**********************************************************************************************
    !! Subroutine to parse the modal meshes from the given file "funit".
    !!**********************************************************************************************
    subroutine parse_modal_mesh(funit, modaltype, mumesh, cmesh, modestring)
        integer, intent(in) :: funit
        character(len=*), intent(in) :: modaltype
        complex(wp), allocatable, intent(out) :: mumesh(:), cmesh(:)
        character(len=*), intent(out) :: modestring
        
        select case (modaltype)
            case ("periodic")
                call parse_periodic(funit, mumesh, cmesh, modestring)
            case ("infinite")
                call parse_infinite(funit, mumesh, cmesh, modestring)
            case default
                write (stderr, '(a)') tag_error // "Unknown modal type: '" // trim(modaltype) // "'"
                stop errmsg_invalid_arg
        end select
        
    end subroutine
    
    !!=================================== PERIODIC TYPE ============================================
    
    !!**********************************************************************************************
    !! Subroutine to parse the modal meshes from the given file "funit" assuming the selected
    !! transverse boundary conditions are "periodic".
    !!**********************************************************************************************
    subroutine parse_periodic(funit, mumesh, cmesh, modestring)
        integer, intent(in) :: funit
        complex(wp), allocatable, intent(out) :: mumesh(:), cmesh(:)
        character(len=*), intent(out) :: modestring
        integer :: d, istat
        real(wp) :: wol
        
        namelist /periodic/ d, wol
        
        !! 1. Read the namelist and check for errors:
        read(unit=funit, nml=periodic, iostat=istat)
		call check_and_rewind_nml(funit, istat)  !! Check for reading error and rewind the namelist file (for further reading).
		
        !! 2. Fill the mumesh/cmesh for the 'periodic' boundary conditions:
        call fill_mumesh_periodic(d, wol, mumesh, cmesh)
        
        !! 3. Write the modestring:
        write (modestring, '(a,i0,a,f0.6,a,i0,a)') &
            "Mode space (periodic, d=", d, ", W/lambda=", wol, ", Nmu=", size(mumesh), ")"
        
    end subroutine
	
	!!**********************************************************************************************
	!! Compute the arrays 'mumesh' and 'cmesh' for a square waveguide with periodic boundary conditions
	!! and according to the given parameters.
	!! This subroutine assumes that the arrays 'mumesh' and 'cmesh' are not allocated.
	!! d   = Total number of spatial dimensions of the waveguide (usually d=2 or d=3).
	!! wol = Width-to-wavelength ratio, W/lambda.
	!!**********************************************************************************************
	subroutine fill_mumesh_periodic(d, wol, mumesh, cmesh)
        integer, intent(in) :: d
        real(wp), intent(in) :: wol
        complex(wp), allocatable, intent(out) :: mumesh(:), cmesh(:)
		real(wp), allocatable :: dists(:), dists_uniq(:)
		integer, allocatable :: dists_cnts(:)
		integer :: ndist
		
		!! Generate the list of all values of 'mu':
		call distance_square_lattice_periodic(d - 1, wol, ndist, dists)  !! Compute the distances from the origin to lattice points.
		call ssort_real(dists(1:ndist))  !! Sort the distances in ascending order to detect and count duplicates.
		call uniq_real(dists(1:ndist), dists_uniq, dists_cnts)  !! Stores the distinct values of distances and count the multiplicities.
		
        mumesh = sqrt(1.0_wp + tiny(1.0_wp)*iu - (dists_uniq/wol)**2)   !! Implicit allocation here.
		cmesh = dists_cnts
		
	end subroutine
	
	!!**********************************************************************************************
	!! Compute the list of the distances between the origin (0, 0, ...) and the points on an integer lattice.
    !! In other words, this subroutine gives the list of the transverse wavenumber, k_(perp,n)*W/(2*pi),
    !! of the eigenmodes of a square of side "W" with periodic boundary conditions.
	!! dlat  = Number of spatial dimensions of the integer lattice (dlat = d - 1).
	!! rmax  = Maximum radius of the lattice in the lattice units.
	!! ndist = Returned number of distances.
	!! dists = Returned list of distances, stored in a 1D array.
	!!**********************************************************************************************
	subroutine distance_square_lattice_periodic(dlat, rmax, ndist, dists)
		integer, intent(in) :: dlat
		real(wp), intent(in) :: rmax
		integer, intent(out) :: ndist
		real(wp), intent(out), allocatable :: dists(:)
		integer :: nmax, vcenter, vside, i, j, k, kdiv
		real(wp) :: dist2, rmax2
		
		!! Prepare the array:
		vcenter = floor(rmax)     !! Position of the center since the lattice is computed on the positive integers.
		vside = 2*vcenter + 1     !! Number of modes along a side of the lattice.
		nmax = vside**dlat       !! Estimated size of the array (strict upper bound).
		if (allocated(dists)) deallocate(dists)
		allocate(dists(nmax))
		
		rmax2 = rmax*rmax
		ndist = 0  !! At the beginning, there is no point in the array.
		do i = 1, nmax
			!! Compute the square norm of the lattice point:
			dist2 = 0._wp  !! Square norm of the lattice vector.
			k = i - 1
			do j = 1, dlat
				kdiv = k/vside  !! Integer division (same algo as for base conversion, but base is 'vside').
				dist2 = dist2 + (k - vside*kdiv - vcenter)**2
				k = kdiv
			end do
			!! Check if the point is inside the sphere of radius "rmax":
			if (dist2 <= rmax2) then
				ndist = ndist + 1
				dists(ndist) = sqrt(dist2)
			end if
		end do
	end subroutine
    
    !!=================================== INFINITE TYPE ============================================
    
    !!**********************************************************************************************
    !! Subroutine to parse the modal meshes from the given file "funit" assuming the selected
    !! transverse boundary conditions are "infinite", i.e., an infinite slab.
    !!**********************************************************************************************
    subroutine parse_infinite(funit, mumesh, cmesh, modestring)
        integer, intent(in) :: funit
        complex(wp), allocatable, intent(out) :: mumesh(:), cmesh(:)
        character(len=*), intent(out) :: modestring
        integer :: d, nmu, istat
        real(wp) :: ash
        
        namelist /infinite/ d, nmu, ash
        
        !! 1. Read the namelist and check for errors:
        read(unit=funit, nml=infinite, iostat=istat)
		call check_and_rewind_nml(funit, istat)  !! Check for reading error and rewind the namelist file (for further reading).
		
        !! 2. Fill the mumesh/cmesh for the 'periodic' boundary conditions:
        call fill_mumesh_infinite(d, nmu, ash, mumesh, cmesh)
        
        !! 3. Write the modestring:
        write (modestring, '(2(a,i0),a,f0.6)') &
            "Mode space (infinite, d=", d, ", Nmu=", nmu, ", ash=", ash, ")"
        
    end subroutine
    
    !!**********************************************************************************************
    !! Function defining the integration path mu(t) used to integrate over the direction subspace
    !! in the subroutine "fill_mumesh_infinite()".
    !! Note that this path must go through the points: mu(-1)=-1, mu(0)=0, and mu(1)=1.
    !!**********************************************************************************************
    elemental function mu_path(ash, t) result(mu)
        real(wp), intent(in) :: ash, t
        complex(wp) :: mu
        mu = t + iu*ash*t*(1.0_wp - t)*(1.0_wp + t)
    end function
    
    !!**********************************************************************************************
    !! Function defining the derivative of the integration path mu(t) with respect to the parameter "t".
    !! This function must be the exact analytical derivative of mu(t) given by "mu_path()".
    !!**********************************************************************************************
    elemental function mu_der(ash, t) result(dmu)
        real(wp), intent(in) :: ash, t
        complex(wp) :: dmu
        dmu = 1.0_wp + iu*ash*(1.0_wp - 3*t*t)
    end function
    
	!!**********************************************************************************************
	!! Compute the arrays 'mumesh' and 'cmesh' for an infinite slab.
	!! Note that the integration path mu(t), initially covering the interval mu=[-1, 1], is shifted in the complex plane of "mu"
    !! in order to avoid the singularities close to the real axis of "mu".
	!! d   = Total number of spatial dimensions of the waveguide (usually d=2 or d=3).
	!! nmu = Number of simulated pseudo-modes. Number of points of the Gaussian quadrature in the directional subspace.
	!! ash = Shift factor of the integration path in the complex 'mu' plane. The integration path is: mu(t) = t + i*ash*t*(1-t)*(1+t) for t in [0, 1].
	!!**********************************************************************************************
	subroutine fill_mumesh_infinite(d, nmu, ash, mumesh, cmesh)
        integer, intent(in) :: d, nmu
        real(wp), intent(in) :: ash
        complex(wp), allocatable, intent(out) :: mumesh(:), cmesh(:)
        real(wp), allocatable :: tpoints(:), weights(:)
        real(wp) :: a, b, xmin, xmax
        
        allocate(tpoints(nmu))
        allocate(weights(nmu))
        
        !! Parameters of the Gauss-Jacobi quadrature:
        a = 0.0_wp
        b = (d - 3.0_wp)/2
        xmin = 0.0_wp
        xmax = 1.0_wp
        
        call gauss_jacobi_quad(nmu, a, b, xmin, xmax, tpoints, weights)
        mumesh = mu_path(ash, tpoints)  !! Implicit allocation.
        cmesh = 2.0_wp * mumesh * mu_der(ash, tpoints) &
            * (((1.0_wp - mumesh*mumesh)/(1.0_wp - tpoints))**b) * weights
        
        deallocate(tpoints)
        deallocate(weights)
        
    end subroutine
    
end module
