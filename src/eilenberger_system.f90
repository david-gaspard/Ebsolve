!!**************************************************************************************************
!! Created on 2024-07-24 at 17:24:46 CEST by David Gaspard <david.gaspard@espci.fr>
!! This program is distributed under the MIT License.
!! Fortran module providing the Eilenberger_System object for the Ebsolve (v2) program.
!!**************************************************************************************************
module eilenberger_system
    use modal_mesh   !! Import the modal mesh generator.
    use integrator   !! Import the Eilenberger integrator.
	implicit none
	
    real(wp), parameter :: x0 = 0.0_wp  !! Coordinate of the leftmost edge of the disordered region on the mesh (do not change).
    real(wp), parameter :: x1 = 1.0_wp  !! Coordinate of the rightmost edge of the disordered region on the mesh (do not change).
    real(wp), parameter :: xshift = 10*epsilon(1._wp)  !! Small shift of the 'xmesh' in the special case 'xa' or 'xb' coincides with it.
    
	!!**********************************************************************************************
	!! Type defining the Eilenberger system. This type contains all the variables describing
	!! the current state of the Eilenberger equation, in particular, the matrix field Qn(x) and the matrix radiance g(mu,x).
	!!**********************************************************************************************
	type eilenberger_system_t
		private   !! All the type members should be private. This ensures encapsulation, i.e., 
		          !! members should not be modified from outside to maintain the object consistency.
            
            !! Directional space (muspace):
            character(len=:), allocatable :: modetype     !! Type of modal space. "periodic" (periodic boundary conditions) or "infinite" (infinite slab).
            character(len=:), allocatable :: modestring   !! Description of the modal space. Typically: "Mode space (periodic, d=2, W/lambda=3.5, Nmu=4)".
            integer :: nmu                                !! Number of "mu" points. The true number of modes is n_mode = 2*nmu - 1.
            complex(wp), allocatable :: mumesh(:)         !! List of direction cosines mu=cos(theta) of the modes. Dimensions: mumesh(nmu).
            complex(wp), allocatable :: cmesh(:)          !! Weight coefficients "c" associated to each "mu" in the same order. Dimensions: cmesh(nmu).
            
            !! Positional space (xspace):
            integer :: nxdiso                             !! Number of points of the disordered region in the longitudinal "x" direction.
            integer :: nxfree                             !! Number of points of the free region on one side (the other side also has "nxfree" too). Should be at least nxfree=1.
            integer :: nx                                 !! Total number of points in the "x" direction used in the discretization (nx = nxdiso + 2*nxfree).
            real(wp), allocatable :: xmesh(:)             !! List of "x" points on the mesh in increasing order. Dimensions: xmesh(nx).
            real(wp), allocatable :: sigma_scat(:)        !! Inverse scattering mean free path, 1/l_scat(x), in intervals of xmesh. Dimensions: (nx-1).
            real(wp), allocatable :: sigma_abso(:)        !! Inverse absorption mean free path, 1/l_abso(x), in intervals of xmesh. Dimensions: (nx-1).
            real(wp) :: dscat                             !! Total scattering depth, integral dx/l_scat(x) for x over xmesh.
            real(wp) :: dabso                             !! Total absorption depth, integral dx/l_abso(x) for x over xmesh.
            
            !! Contact parameters:
            real(wp) :: xa, xb                            !! Positions of the contact interactions (xa, xb).
            real(wp) :: naper_a, naper_b                  !! Numerical apertures of contact interactions. naper = sin(theta_max), 0 < naper <= 1.
            real(wp) :: mumin_a, mumin_b                  !! Minimum real part of mu=cos(theta) imposed by the numericall aperture: mumin=sqrt(1 - naper^2).
            complex(wp) :: gamma_a, gamma_b               !! Corresponding coupling parameters, i.e., variables of the generating function.
            real(wp) :: tm, geps                          !! Transmission value, 'tm', and small shift over the gamma parameter, 'geps'.
            
            !! Generating function and density (output):
            complex(wp) :: fa, fb                   !! Generating functions at positions "xa" and "xb", respectively.
            complex(wp) :: fm                       !! Sought generating function, F(gammaz) = fm = (fa + fb)/2
            real(wp) :: rhom                        !! Eigenvalue density computed from "fm".
            
            !! Fields (output):
            complex(wp), allocatable :: qn(:,:)     !! Field Qn(x). Dim: qn(3,nx). Triplets order: (Qn11, Qn12, Qn21).
            complex(wp), allocatable :: jn(:,:)     !! Matrix current Jn(x). Dim: jn(3,nx). Triplets  order: (Jn11, Jn12, Jn21).
            complex(wp), allocatable :: gp(:,:,:)   !! Radiance g_+(mu,x). Dim: gp(3,nmu,nx). Triplets order: (gp11, gp12, gp21).
            complex(wp), allocatable :: gm(:,:,:)   !! Radiance g_-(mu,x). Dim: gm(3,nmu,nx). Triplets order: (gm11, gm12, gm21).
            
		contains
			!! Initialization:
            procedure :: parse => parse_ebsys
			procedure :: init => init_ebsys
			procedure :: del => del_ebsys
            !! Setters:
            procedure :: set_transmission => set_transmission_ebsys
			procedure :: set_qn => set_qn_ebsys
            !! Getters:
            procedure :: get_nx => get_nx_ebsys
            procedure :: get_xmesh => get_xmesh_ebsys
			procedure :: get_qn => get_qn_ebsys
			procedure :: get_jn => get_jn_ebsys
            procedure :: get_rho => get_rho_ebsys
			procedure :: get_genfun => get_genfun_ebsys
            procedure :: to_string => to_string_ebsys
            procedure :: output_directory => output_directory_ebsys
            !! Computation:
            procedure :: next_state => next_state_ebsys
			
	end type
    
	contains
    
    !!====================================== PARSERS ===============================================
    
    !!**********************************************************************************************
	!! Subroutine to parse the Eilenberger_System object from the given settings (namelist) file.
	!!**********************************************************************************************
    subroutine parse_ebsys(self, funit)
        class(eilenberger_system_t), intent(inout) :: self
        integer, intent(in) :: funit
        character(len=50) :: modetype
        integer :: nxdiso, nxfree, istat
		real(wp) :: dscat, dabso, xa, xb, naper_a, naper_b
        complex(wp), allocatable :: mumesh(:), cmesh(:)
        character(len=150) :: modestring
        
        namelist /eilenberger_system/ modetype, nxdiso, nxfree, dscat, dabso, xa, xb, naper_a, naper_b
        
        !! 1. Read the main Eilenberger_System settings and check for errors:
        read(unit=funit, nml=eilenberger_system, iostat=istat)
		call check_and_rewind_nml(funit, istat)  !! Check for reading error and rewind the namelist file (for further reading).
		
        !! 2. Parse and fill the mumesh/cmesh:
        call parse_modal_mesh(funit, modetype, mumesh, cmesh, modestring)
        
        !! 3. Initialize the Eilenberger_System object:
        call init_ebsys(self, modetype, modestring, mumesh, cmesh, nxdiso, nxfree, dscat, dabso, xa, naper_a, xb, naper_b)
        
        deallocate(mumesh)  !! Free the memory allcaoted by parse_modal_mesh().
        deallocate(cmesh)
        
    end subroutine
    
	!!==================================== INITIALIZERS ============================================
    
	!!**********************************************************************************************
	!! Subroutine to initialize the Eilenberger_System object.
	!!**********************************************************************************************
	subroutine init_ebsys(self, modetype, modestring, mumesh, cmesh, nxdiso, nxfree, dscat, dabso, xa, naper_a, xb, naper_b)
		class(eilenberger_system_t), intent(inout) :: self
        complex(wp), intent(in) :: mumesh(:), cmesh(:)
        character(len=*), intent(in) :: modetype, modestring
        integer, intent(in) :: nxdiso, nxfree
        real(wp), intent(in) :: dscat, dabso, xa, xb, naper_a, naper_b
        integer :: ix, nx
        
        !! 1. Check for possible errors in the arguments:
        call check_ebsys(mumesh, cmesh, nxdiso, nxfree, dscat, dabso, naper_a, naper_b)
        
        !! 2. Copy the state variables in the current Eilenberger_System object:
        self%nmu = size(mumesh)  !! Number of "mu" points (modes or pseudo-modes).
        self%mumesh = mumesh     !! Implicit allocations and copy.
        self%cmesh = cmesh
        self%modetype = trim(modetype)
        self%modestring = trim(modestring)
        self%nxdiso = nxdiso
        self%nxfree = nxfree
        self%nx = nxdiso + 2*nxfree
        self%dscat = dscat
        self%dabso = dabso
        self%xa = xa   !! Contact interactions.
        self%xb = xb
        self%naper_a = naper_a
        self%naper_b = naper_b
        self%mumin_a = sqrt(1.0_wp - naper_a**2)
        self%mumin_b = sqrt(1.0_wp - naper_b**2)
        
        !! 3. Build the positional lattice:
        nx = self%nx
        allocate(self%xmesh(nx))
		do ix = 1, nx  !! Fill the xmesh with a uniform mesh:
			
            !!self%xmesh(ix) = x0 + (ix - 2)*(x1 - x0)/(nx - 3)
            
            self%xmesh(ix) = x0 + (ix - nxfree - 1)*(x1 - x0)/(nxdiso - 1) 
            
            if (self%xmesh(ix) == xa .or. self%xmesh(ix) == xb) then  !! Slightly shift the xmesh if 'xa' or 'xb' belongs to it.
				self%xmesh(ix) = self%xmesh(ix) - xshift
			end if
		end do
		if (xa <= self%xmesh(1))  self%xmesh(1)  = 2*xa - x0  !! Shifts the boundaries of xmesh to include the points 'xa' and 'xb'.
		if (xb <= self%xmesh(1))  self%xmesh(1)  = 2*xb - x0
		if (xa >= self%xmesh(nx)) self%xmesh(nx) = 2*xa - x1
		if (xb >= self%xmesh(nx)) self%xmesh(nx) = 2*xb - x1
        
        !! 4. Build the inverse mean-free-path lattices:
        allocate(self%sigma_scat(nx-1))
        allocate(self%sigma_abso(nx-1))
        self%sigma_scat = 0.0_wp
        self%sigma_abso = 0.0_wp
        do ix = nxfree+1, nxdiso+nxfree-1
            self%sigma_scat(ix) = dscat
            self%sigma_abso(ix) = dabso
        end do
        
        !! 5. Allocate space for the fields:
        allocate(self%qn(3, nx))
        allocate(self%jn(3, nx))
        allocate(self%gp(3, self%nmu, nx))
        allocate(self%gm(3, self%nmu, nx))
        
        !! 6. Zero-fill the uninitialized variables to avoid memory garbage:
        self%qn = 0.0_wp
		self%jn = 0.0_wp
		self%gp = 0.0_wp
		self%gm = 0.0_wp
        self%gamma_a = 0.0_wp
        self%gamma_b = 0.0_wp
        self%tm = 0.0_wp
        self%geps = 0.0_wp
        self%fa = 0.0_wp
        self%fb = 0.0_wp
        self%fm = 0.0_wp
        self%rhom = 0.0_wp
        
	end subroutine
    
    !!**********************************************************************************************
    !! Check the validity of the parameters of Eilenberger_System object. Stops the program on error.
    !!**********************************************************************************************
    subroutine check_ebsys(mumesh, cmesh, nxdiso, nxfree, dscat, dabso, naper_a, naper_b)
        complex(wp), intent(in) :: mumesh(:), cmesh(:)
        integer, intent(in) :: nxdiso, nxfree
        real(wp), intent(in) :: dscat, dabso, naper_a, naper_b
        integer, parameter :: nxdiso_min = 2, nxfree_min = 1
        logical :: error
        
        error = .false.
        if (size(mumesh) /= size(cmesh)) then
            write (stderr, '(2(a,i0),a)') tag_error // "Inconsistent lengths of 'mumesh' (", &
                size(mumesh), ") and 'cmesh' (", size(cmesh), ")..."
            error = .true.
        else if (nxdiso < nxdiso_min) then
            write (stderr, '(2(a,i0),a)') tag_error // "Too few x points in the bulk, received nxdiso=", &
                nxdiso, ", expected at least ", nxdiso_min, "..."
            error = .true.
        else if (nxfree < nxfree_min) then
            write (stderr, '(2(a,i0),a)') tag_error // "Too few x points at the edge, received nxfree=", &
                nxfree, ", expected at least ", nxfree_min, "..."
            error = .true.
        else if (dscat < 0.0_wp .or. dabso < 0.0_wp) then
            write (stderr, '(a)') tag_error // "Invalid optical depth, cannot be negative..."
            error = .true.
        else if (naper_a < 0.0_wp .or. naper_a > 1.0_wp .or. naper_b < 0.0_wp .or. naper_b > 1.0_wp) then
            write (stderr, '(a)') tag_error // "Invalid numerical aperture, must be within [0, 1]..."
            error = .true.
        end if
        
        if (error) stop errmsg_invalid_arg
        
    end subroutine
    
	!!**********************************************************************************************
	!! Frees the memory allocated by the arrays of the type Eilenberger_System().
	!! Calling this subroutine is only useful when memory consumption is critical.
	!! However, it should be reminded that Fortran automatically deallocates
	!! allocated arrays when out of scope.
	!!**********************************************************************************************
	subroutine del_ebsys(self)
		class(eilenberger_system_t), intent(inout) :: self
        if (allocated(self%mumesh)) deallocate(self%mumesh)
        if (allocated(self%cmesh)) deallocate(self%cmesh)
		if (allocated(self%xmesh)) deallocate(self%xmesh)
		if (allocated(self%sigma_scat)) deallocate(self%sigma_scat)
		if (allocated(self%sigma_abso)) deallocate(self%sigma_abso)
		if (allocated(self%qn)) deallocate(self%qn)
		if (allocated(self%jn)) deallocate(self%jn)
		if (allocated(self%gp)) deallocate(self%gp)
		if (allocated(self%gm)) deallocate(self%gm)
	end subroutine
    
    !!**********************************************************************************************
    !! Assigns the transmission value "tm" (=T) and the small shift "geps" over gamma = gamma_a*gamma_b.
    !! Typically, "geps" is of the order of the machine precision, but can be even smaller.
    !! In fact, the parameter "geps" roughly controls the zero level of the transmission-eigenvalue density "rho(T)".
    !!**********************************************************************************************
    subroutine set_transmission_ebsys(self, tm, geps)
        class(eilenberger_system_t), intent(inout) :: self
        real(wp), intent(in) :: tm, geps
        
        !! 1. Check for possible errors:
        if (tm == 0.) then
            write (stderr, '(a)') tag_error // "Transmission value 'tm' cannot be zero."
            stop errmsg_invalid_arg
        end if
        
        !! 2. Assigns the values:
        self%tm = tm
        self%geps = geps
        self%gamma_a = sqrt(1./tm + geps*iu)  !! The product gamma_a*gamma_b must be equal to 1./tm + geps*iu :
        self%gamma_b = self%gamma_a
        
    end subroutine
    
	!!**********************************************************************************************
	!! Assigns the value of the Qn(x) field. This is useful when using the relaxation algorithm to
    !! accelerate the convergence of the fixed-point iterative for Qn(x).
	!!**********************************************************************************************
	subroutine set_qn_ebsys(self, qn)
		class(eilenberger_system_t), intent(inout) :: self
		complex(wp), intent(in) :: qn(:, :)
		if (size(qn, 1) /= 3 .or. size(qn, 2) /= self%nx) then
			write (stderr, '(a)') tag_error // "Invalid dimension of the Qn field."
			stop errmsg_invalid_arg
		end if
		self%qn = qn  !! In principle, no reallocation is needed here.
	end subroutine
    
    !!====================================== GETTERS ===============================================
    
    !!**********************************************************************************************
    !! Gets "nx", i.e., the total number of points of "xmesh" used by the Eilenberger_System.
    !!**********************************************************************************************
    function get_nx_ebsys(self) result(nx)
        class(eilenberger_system_t), intent(in) :: self
        integer :: nx
        nx = self%nx
    end function
    
    !!**********************************************************************************************
    !! Returns the "xmesh", i.e., the array of points in the direction longitudinal to the waveguide ("x" direction).
    !!**********************************************************************************************
    subroutine get_xmesh_ebsys(self, xmesh)
        class(eilenberger_system_t), intent(in) :: self
		real(wp), intent(out) :: xmesh(:)
		xmesh = self%xmesh
    end subroutine
    
	!!**********************************************************************************************
	!! Gets the values of the generating functions: "fa" is the value at "xa", "fb" at "xb", and
	!! "fm" is the sought mean generating function F(gammaz) = fm = (fa + fb)/2.
	!! In principle, this subroutine should only be called after the Eilenberger_System has been
	!! solved and convergence has been reached.
	!!**********************************************************************************************
	function get_rho_ebsys(self) result(rho)
		class(eilenberger_system_t), intent(in) :: self
        real(wp) :: rho
		rho = self%rhom
	end function
    
	!!**********************************************************************************************
	!! Gets the values of the generating functions: "fa" is the value at "xa", "fb" at "xb", and
	!! "fm" is the sought mean generating function F(gammaz) = fm = (fa + fb)/2.
	!! In principle, this subroutine should only be called after the Eilenberger_System has been
	!! solved and convergence has been reached.
	!!**********************************************************************************************
	subroutine get_genfun_ebsys(self, fa, fb, fm)
		class(eilenberger_system_t), intent(in) :: self
		complex(wp), intent(out) :: fa, fb, fm
		fa = self%fa
		fb = self%fb
		fm = self%fm
	end subroutine
	
	!!**********************************************************************************************
	!! Gets the matrix field Qn(x) from outside the class. The "qn" array is supposed to be already
    !! allocated to the appropriate dimensions: (3,nx).
	!!**********************************************************************************************
	subroutine get_qn_ebsys(self, qn)
		class(eilenberger_system_t), intent(in) :: self
		complex(wp), intent(out) :: qn(:,:)
		qn = self%qn
	end subroutine
	
	!!**********************************************************************************************
	!! Gets the matrix current Jn(x) from outside the class. The "jn" array is supposed to be already
    !! allocated to the appropriate dimensions: (3,nx).
	!!**********************************************************************************************
	subroutine get_jn_ebsys(self, jn)
		class(eilenberger_system_t), intent(in) :: self
		complex(wp), intent(out) :: jn(:,:)
        call compute_jn_from_g_ebsys(self)  !! First compute the matrix current Jn(x) [not computed by default].
		jn = self%jn
	end subroutine
    
	!!**********************************************************************************************
	!! Returns a string representation of the current Eilenberger_System "self"
    !! using a given prefix string "prefix" (typically a comment character).
	!!**********************************************************************************************
    function to_string_ebsys(self, prefix, newline) result(str)
        class(eilenberger_system_t), intent(in) :: self
        character(len=*), intent(in) :: prefix, newline
        character(len=400+3*len(prefix)+3*len(newline)) :: str
        
        write (str, '(a)') &
            trim(prefix) // " " // trim(self%modestring) // newline
        write (str, '(2(a,f0.6),2(a,i0),a)') &
            trim(str) // trim(prefix) // " Medium (L/lscat=", &
            self%dscat, ", L/labso=", self%dabso, ", Nxdiso=", self%nxdiso, ", Nxfree=", self%nxfree, ")" // newline
        write (str, '(4(a,f0.6),a)') &
            trim(str) // trim(prefix) // " Contacts (xa=", &
                self%xa, ", xb=", self%xb, ", NApera=", self%naper_a, ", NAperb=", self%naper_b, ")" // newline
        
    end function
    
    !!**********************************************************************************************
    !! Generate a path for output files derived from the current Eilenberger_System.
    !! pathsep = Path separator used by the filesystem (typically "/" on UNIX).
    !!**********************************************************************************************
    subroutine output_directory_ebsys(self, outputdir, pathsep)
        class(eilenberger_system_t), intent(in) :: self
        character(len=*), intent(out) :: outputdir
        character(len=*), intent(in) :: pathsep
        character(len=30) :: string_dscat, string_dabso, folder_dscat, folder_dabso
        
        string_dscat = "dscat_"
        string_dabso = "dabso_"
        if (self%dscat < 1.) string_dscat = trim(string_dscat) // "0"   !! Append a zero just before the decimal point.
        if (self%dabso < 1.) string_dabso = trim(string_dabso) // "0"
        
        write (folder_dscat, '(a,f0.6)') trim(string_dscat), self%dscat
        write (folder_dabso, '(a,f0.6)') trim(string_dabso), self%dabso
        
        !! Prepare the path for the output directory:
        outputdir = "out" // pathsep // &
                trim_zero(folder_dscat) // pathsep // &
                trim_zero(folder_dabso) // pathsep // &
                self%modetype // pathsep
        
    end subroutine
    
    !!==================================== COMPUTATIONS ============================================
	
	!!**********************************************************************************************
	!! Compute the Qn(x) field based on the current values of the g_+(mu,x) and g_-(mu,x) radiances.
	!! This subroutine assumes that the radiances g_+(mu,x) and g_-(mu,x) have been computed already.
	!! The resulting matrix field is stored in the class member: "qn".
	!!**********************************************************************************************
	subroutine compute_qn_from_g_ebsys(self)
		class(eilenberger_system_t) :: self
		complex(wp) :: denom, cdmu
		integer :: i
		!! First compute the normalization factor (denominator):
		denom = (0.0_wp, 0.0_wp) 
		self%qn = 0.0_wp  !! Reset the matrix field Qn(x) to zero before summing.
		do i = 1, self%nmu
			cdmu = self%cmesh(i)/self%mumesh(i)  !! Compute the quotient c_i/mu_i.
			denom = denom + cdmu  !! denom = sum(c_n/mu_n, over n)
			self%qn(:,:) = self%qn(:,:) + cdmu*(self%gp(:,i,:) + self%gm(:,i,:))/2.
		end do
		!! Then normalize the Qn(x) field by 'denom':
		self%qn = self%qn/denom
	end subroutine
	
	!!**********************************************************************************************
	!! Compute the matrix current Jn(x) based on the current values of the g_+(mu,x) and g_-(mu,x)
	!! radiances. This subroutine assumes that the radiances g_+(mu,x) and g_-(mu,x) have been
	!! computed already. The resulting matrix current is stored in the class member: "jn".
	!!**********************************************************************************************
	subroutine compute_jn_from_g_ebsys(self)
		class(eilenberger_system_t) :: self
		complex(wp) :: denom
		integer :: i
		!! First compute the normalization factor (denominator):
		denom = (0.0_wp, 0.0_wp)
		self%jn = 0.0_wp  !! Reset the matrix current Jn(x) to zero before summing.
		do i = 1, self%nmu
			denom = denom + self%cmesh(i)  !! denom = sum(c_n, over n)
			self%jn(:,:) = self%jn(:,:) + self%cmesh(i)*(self%gp(:,i,:) - self%gm(:,i,:))/2.
		end do
		!! Then normalize the Jn(x) current by 'denom':
		self%jn = self%jn/denom
	end subroutine
	
    !!**********************************************************************************************
    !! Integrates the Eilenberger equation along the 'x' axis assuming hte Qn(x) field is fixed.
    !! This subroutine finds the new matrix radiances g_+(mu,x) and g_-(mu,x) for all values of "mu" using
    !! the subroutine: "ebintegrate()".
    !!**********************************************************************************************
    subroutine integrate_x_ebsys(self)
        class(eilenberger_system_t) :: self
        complex(wp), allocatable :: amesh(:, :), bmesh(:, :), pmesh(:, :)
        complex(wp) :: mu, ca(3), cb(3)
		integer :: ix, imu
        
        !! Allocate buffer arrays (Schopohl spinors):
        allocate(amesh(2, self%nx))
        allocate(bmesh(2, self%nx))
        
        !! Compute the field P(x) = -Qn(x)/(2*lscat(x)) - Lambda_3/(2*labso(x)):
        allocate(pmesh(3, self%nx-1))
        do ix = 1, self%nx-1  !! Loop over the positions ("x" axis).
            pmesh(:, ix) = - self%sigma_scat(ix) * (self%qn(:, ix) + self%qn(:, ix+1))/4  &
                           - (self%sigma_abso(ix)/2) * [1.0, 0.0, 0.0]
        end do
        
        do imu = 1, self%nmu  !! Loop over the directions, i.e., values of "mu".
            
			mu = self%mumesh(imu)  !! Extract the value of "mu" .
            
			if (mu%re > self%mumin_a) then  !! Account for the limited numerical aperture at x=xa.
                ca = iu * mu * self%gamma_a * [0.0, 1.0, 0.0]   !! Contact matrix: Ca = i * mu * gamma_a * Lambda_+
            else
                ca = 0.0   !! No interaction for extreme rays beyond the numerical aperture.
            end if
            if (mu%re > self%mumin_b) then  !! Account for the limited numerical aperture at x=xb.
                cb = iu * mu * self%gamma_b * [0.0, 0.0, 1.0]   !! Contact matrix: Cb = i * mu * gamma_b * Lambda_-
            else
                cb = 0.0   !! No interaction for extreme rays beyond the numerical aperture.
            end if
            
            !! Integrates the Eilenberger equation along the rays (in the position space):
            call ebintegrate(+1, mu, self%xmesh, pmesh, self%xa, +ca, self%xb, +cb, amesh, bmesh, self%gp(:,imu,:))  !! Integrates g_+(mu,x).
            call ebintegrate(-1, mu, self%xmesh, pmesh, self%xa, -ca, self%xb, -cb, amesh, bmesh, self%gm(:,imu,:))  !! Integrates g_-(mu,x).
		end do
        
        deallocate(amesh)
        deallocate(bmesh)
        deallocate(pmesh)
        
    end subroutine
    
	!!**********************************************************************************************
	!! Compute the next state of the Eilenberger_System. The computations proceed in the following order:
    !! 1. From a given matrix field Qn(x), compute the matrix radiance g_+(mu,x) and g_-(mu,x).
    !! 2. From g_+(mu,x) and g_-(mu,x), deduce the new matrix field Qn(x). The latter is ready for the next iteration.
    !! 3. From g_+(mu,x) and g_-(mu,x), compute the generating function F(gamma) and the eigenvalue density rho(T).
    !! Note that this subroutine assumes the Qn(x) field is already set to the best possible guess (typically zero).
	!!**********************************************************************************************
    subroutine next_state_ebsys(self)
        class(eilenberger_system_t) :: self
        call integrate_x_ebsys(self)
        call compute_qn_from_g_ebsys(self)
		call compute_rho_from_g_ebsys(self)
    end subroutine
    
	!!**********************************************************************************************
	!! Computes the generating function F(gammaz) = (1/N_modes) * Tr[ t⁺t / (1 - gammaz * t⁺t) ]
	!! from the current value of the matrix radiances g_+(mu,x) and g_-(mu,x) and accounting for the
    !! limited numerical aperture of the observables at x=xa and x=xb.
    !! Therefore, this subroutine assumes that the matrix radiances has been computed already.
    !! 
    !! The important base formulas are:
    !! 
    !!   1|  rho(T) = Im F(1/T + 0*i)/(pi*T^2)
    !! 
    !!   2|  F(gamma) = (i/2*min(N_a, N_b)) * ( J_21(xa)/gamma_b + J_12(xb)/gamma_a )
    !!       N_a, N_b = Number of modes involved at points xa and xb, respectively.
    !!       N_a, N_b = sum_imu c_imu, where the sum is limited by the corresponding numerical aperture.
    !! 
    !!   3|  J(x) = sum_imu c_imu * (g_+(imu, x) - g_-(imu, x))/2
    !!       Note that the sum over the modes (imu) is limited by the corresponding numerical aperture. 
    !! 
	!! The results are stored in the class members: "fa", "fb", "fm", "rhom".
	!!**********************************************************************************************
	subroutine compute_rho_from_g_ebsys(self)
		class(eilenberger_system_t) :: self
		complex(wp) :: mu, ja, jb, delta, denom
        integer :: imu
        
        !! 1. Compute the matrix currents at x=xa and x=xb using a loop over the modes.
        denom = (0.0_wp, 0.0_wp)
        !!denomb = (0.0_wp, 0.0_wp)
		ja = (0.0_wp, 0.0_wp) !! Reset the matrix current J(x) to zero before summing.
        jb = (0.0_wp, 0.0_wp)
        
		do imu = 1, self%nmu   !! Loop over the modes.
            
            mu = self%mumesh(imu)  !! Extract the direction cosine of the mode, mu=cos(theta).
            
            if (mu%re > self%mumin_a) then
                !!denoma = denoma + self%cmesh(imu)  !! denom = sum(c_n, over n)
                call interp1(self%xmesh, self%gp(3, imu, :) - self%gm(3, imu, :), self%xa, delta)
                ja = ja + self%cmesh(imu)*delta/2
            end if
            if (mu%re > self%mumin_b) then
                !!denomb = denomb + self%cmesh(imu)  !! denom = sum(c_n, over n)
                call interp1(self%xmesh, self%gp(2, imu, :) - self%gm(2, imu, :), self%xb, delta)
                jb = jb + self%cmesh(imu)*delta/2
            end if
            if (mu%re > max(self%mumin_a, self%mumin_b)) then
                denom = denom + self%cmesh(imu)  !! denom = sum(c_n, over n)
            end if
            
		end do
        
        !! 2. Normalize the matrix currents:
		!jna = jna/denoma
		!jnb = jnb/denomb
        
		!! 2. Computes the generating functions:
		self%fa = iu*ja/(denom*self%gamma_b)
		self%fb = iu*jb/(denom*self%gamma_a)
		self%fm = (self%fa + self%fb)/2.
        self%rhom = (self%fm%im)/(pi*self%tm*self%tm)  !! Compute the eigenvalue density.
        
    end subroutine
    
end module
