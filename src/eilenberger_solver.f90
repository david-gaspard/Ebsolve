!!**************************************************************************************************
!! Created on 2024-07-26 at 12:48:37 CEST by David Gaspard <gaspard.dave@gmail.com>
!! Based on the same module written on 2024-01-24 at 11:10:04 CET by David Gaspard for the Ebsolve program (v1).
!! This program is distributed under the Creative Commons (CC) BY-NC-SA license.
!! Fortran module providing the Eilenberger solver type for the Ebsolve program (v2).
!! This object provides the essential subroutines to solve the Eilenberger_System object iteratively.
!! This object also contains all the parameters of the iterative solver.
!!**************************************************************************************************
module eilenberger_solver
    use eilenberger_system
	implicit none
	
	!!**********************************************************************************************
    !! Type defining the Eilenberger solver. This type contains all the parameters and subroutine
    !! needed to solver the Eilenberger_System using the iterative relaxation method.
	!!**********************************************************************************************
    type :: eilenberger_solver_t
        private  !! All type members should be private to ensures encapsulation.
            !! Constant parameters:
            character(len=:), allocatable :: method  !! Iterative method ("fpi", "relax", ...).
            integer :: maxit    !! Maximum number of iterations used to solver the Eilenberger system.
                                !! Typically, "maxit" is much larger than 100 but not as large as 10^6.
            real(wp) :: qtol    !! Tolerance on the relative variation of the Qn(x) field. Typically "qtol" is between 1e-6 and 1e-15. 
	                            !! Too large tolerance, such as 1e-3, are not recommended because it does not directly represent
	                            !! the tolerance on the results, in particular the generating function and rho(T).
            real(wp) :: frelax  !! Relaxation factor. Main parameter of the relaxation method.
                                !! Empirically, the optimal factor is between 0.5 and 0.9 in the quasiballistic regime (L/l_scat < 10).
                                !! In the diffusive regime, over-relaxation with frelax > 2 can also be used to speed up convergence.
                                !! The special case frelax=1 corresponds to the standard fixed-point iteration.
                                !! This parameter is ignored when using the solve_fpi() subroutine.
            integer :: verbose  !! Verbosity level of the solver (0=Quiet, 1=Verbose).
            
            !! Variable parameters:
            real(wp), allocatable :: qerr_history(:)  !! History of the relative errors on the Qn(x) field. Total length=maxit, useful length=niter.
            real(wp), allocatable :: rho_history(:)   !! History of the values of rho(T) through the iterations. Total length=maxit, useful length=niter.
            integer :: niter    !! Final number of iteration needed to meet the convergence criterion. Also the useful length of "qerr_history".
            logical :: conv     !! Convergence flag (True=Converged, False=Not converged).
            
        contains
            !! Initializers:
            procedure :: parse => parse_ebsolve
            procedure :: init => init_ebsolve
            procedure :: del => del_ebsolve
            !! Getters:
            procedure :: converged => converged_ebsolve
            procedure :: get_niter => get_niter_ebsolve
            procedure :: get_qerr => get_qerr_ebsolve
            procedure :: get_rho => get_rho_ebsolve
            procedure :: to_string => to_string_ebsolve
            !! Solvers:
            procedure :: solve => solve_ebsolve
            
    end type
    
	contains
    
	!!**********************************************************************************************
	!! Parse the solver from a file "funit".
	!!**********************************************************************************************
	subroutine parse_ebsolve(self, funit)
        class(eilenberger_solver_t), intent(inout) :: self
		integer, intent(in) :: funit
        character(len=100) :: method
		integer :: istat, maxit, verbose
		real(wp) :: qtol, frelax
		
		namelist /solver/ method, maxit, qtol, frelax, verbose
		
		read(unit=funit, nml=solver, iostat=istat)  !! Read the namelist.
		call check_and_rewind_nml(funit, istat)  !! Check for input errors in the namelist, and rewind.
		call init_ebsolve(self, method, maxit, qtol, frelax, verbose)  !! Calls the initializer.
		
    end subroutine
    
	!!**********************************************************************************************
	!! Initialize the solver with the needed parameters.
	!!**********************************************************************************************
	subroutine init_ebsolve(self, method, maxit, qtol, frelax, verbose)
		class(eilenberger_solver_t), intent(inout) :: self
        character(len=*), intent(in) :: method
        integer, intent(in) :: maxit, verbose
        real(wp), intent(in) :: qtol, frelax
        logical :: error
        
        !! 1. Check for possible errors:
        error = .false.
        if (maxit <= 0) then
			write (stderr, '(a)') tag_error // "Invalid number of iterations."
			error = .true.
		else if (qtol > 1. .or. qtol < epsilon(1._wp)) then
			write (stderr, '(a)') tag_error // "Invalid tolerance, it must be in [meps, 1]."
			error = .true.
        else if (method /= "fpi" .and. method /= "relax") then
            write (stderr, '(a)') tag_error // "Unknown iterative method: " // trim(method)
			error = .true.
		end if
        if (error) stop errmsg_invalid_arg  !! Stop on error.
        
        !! 2. Assigns the parameters:
        self%method = trim(method)  !! Implicit allocation.
        self%maxit = maxit
        self%qtol = qtol
        self%frelax = frelax
        self%verbose = verbose
        self%niter = 0
        self%conv = .false.
        
        !! 3. Allocates the arrays:
        if (allocated(self%qerr_history)) deallocate(self%qerr_history)
        if (allocated(self%rho_history)) deallocate(self%rho_history)
        allocate(self%qerr_history(maxit))
        allocate(self%rho_history(maxit))
        self%qerr_history = 0.  !! Sets all the relative errors to a default value, a virtually unreachable value.
        self%rho_history = 0.   !! Sets all the densities to zero, just to avoid garbage.
        
	end subroutine
    
	!!**********************************************************************************************
	!! Frees the memory occupied by the arrays in the present solver "self".
	!!**********************************************************************************************
	subroutine del_ebsolve(self)
        class(eilenberger_solver_t), intent(inout) :: self
        if (allocated(self%method)) deallocate(self%method)
        if (allocated(self%qerr_history)) deallocate(self%qerr_history)
        if (allocated(self%rho_history)) deallocate(self%rho_history)
    end subroutine
	
    !!====================================== GETTERS ===============================================
    
	!!**********************************************************************************************
	!! Returns the convergence status of the solver "self".
	!!**********************************************************************************************
    function converged_ebsolve(self) result(conv)
        class(eilenberger_solver_t), intent(in) :: self
        logical :: conv
        conv = self%conv
    end function
    
	!!**********************************************************************************************
	!! Returns the final number of iterations used to solver the Eilenberger_System.
	!!**********************************************************************************************
    function get_niter_ebsolve(self) result(niter)
        class(eilenberger_solver_t), intent(in) :: self
        integer :: niter
        niter = self%niter
    end function
    
	!!**********************************************************************************************
	!! Checks if the given iteration index "iter" is in the correct interval of existing history.
    !! If not, this subroutine stops the program.
	!!**********************************************************************************************
    subroutine check_historical_index(iter, itermin, itermax)
        integer, intent(in) :: iter, itermin, itermax
        if (iter < itermin .or. itermax < iter) then
            write (stderr, '(a,i0,a,i0,a,i0,a)') tag_error // "Invalid iteration ", iter, &
                ", expected within interval ", itermin, "..", itermax, "."
            stop errmsg_invalid_arg
        end if
    end subroutine
    
	!!**********************************************************************************************
	!! Returns the relative error on the Qn(x) field at a given iteration index "iter".
	!!**********************************************************************************************
    function get_qerr_ebsolve(self, iter) result(qerr)
        class(eilenberger_solver_t), intent(in) :: self
        integer, intent(in) :: iter
        real(wp) :: qerr
        
        call check_historical_index(iter, itermin=1, itermax=self%niter)
        qerr = self%qerr_history(iter)
        
    end function
    
	!!**********************************************************************************************
	!! Returns the value of the density rho(T) at a given iteration index "iter".
	!!**********************************************************************************************
    function get_rho_ebsolve(self, iter) result(rho)
        class(eilenberger_solver_t), intent(in) :: self
        integer, intent(in) :: iter
        real(wp) :: rho
        
        call check_historical_index(iter, itermin=1, itermax=self%niter)
        rho = self%rho_history(iter)
        
    end function
    
	!!**********************************************************************************************
	!! Returns a string summary of the parameters of the present solver "self" to the output stream "output"
    !! using a given prefix string "prefix" (typically a comment character).
	!!**********************************************************************************************
    function to_string_ebsolve(self, prefix) result(str)
        class(eilenberger_solver_t), intent(in) :: self
        character(len=*), intent(in) :: prefix
        character(len=400+len(prefix)) :: str
        
        select case (self%method)
            case ("fpi")
                write (str, '(a,i0,a,g0.3,a)') trim(prefix) // " Solver (" // trim(self%method) // &
                    ", maxit=", self%maxit, ", qtol=", self%qtol, ")"
            case ("relax")
                write (str, '(a,i0,a,g0.3,a,g0.3,a)') trim(prefix) // " Solver (" // trim(self%method) // &
                    ", maxit=", self%maxit, ", qtol=", self%qtol, ", frelax=", self%frelax, ")"
        end select
        
    end function
    
	!!**********************************************************************************************
	!! Prints the given parameters representing the current iteration status.
    !! ebsys = Current Eilenberger_System.
	!!**********************************************************************************************
	subroutine print_iteration_ebsolve(self, ebsys)
        class(eilenberger_solver_t), intent(inout) :: self
        class(eilenberger_system_t), intent(in) :: ebsys
		
		print '(a,i0,a,g10.3,a,g0)', tag_info // "#", self%niter, & 
				char_tab // "| qerr=", self%qerr_history(self%niter), char_tab // "rho=", ebsys%get_rho()
        
	end subroutine
    
	!!**********************************************************************************************
	!! Solve the Eilenberger_System "ebsys" using the prescribed method ("fpi", "relax").
	!! The Eilenberger_System must have been initialized with the init() subroutine before calling this subroutine.
    !! After calling this subroutine, the solution is stored in the final "ebsys" object, as far as convergence
    !! has been reached.
	!!**********************************************************************************************
    subroutine solve_ebsolve(self, ebsys)
        class(eilenberger_solver_t), intent(inout) :: self
        class(eilenberger_system_t), intent(inout) :: ebsys
        
        !! Print some information at the beginning:
        if (self%verbose >= 1) then
			print '(a,i0,a)', tag_info // "Running the method '" // self%method // "' with maxit=", self%maxit, "..."
		end if
        
        !! Run the appropriate method:
        select case (self%method)
        case ("fpi")
            call solve_fpi_ebsolve(self, ebsys)
        case ("relax")
            call solve_relax_ebsolve(self, ebsys)
        end select
        
        !! Print some information at the end:
		if (self%verbose >= 1) then
			call print_iteration_ebsolve(self, ebsys)
			if (self%conv) then
				print '(a,i0,a)', tag_info // "Successfully converged after ", self%niter, " iterations."
			else
				print '(a,i0,a)', tag_warn // "Failed to reach convergence after ", self%niter, " iterations."
			end if
		end if
        
    end subroutine
    
	!!**********************************************************************************************
	!! Solve the Eilenberger_System "ebsys" using an elementary fixed-point iteration (FPI).
	!! The Eilenberger_System must have been initialized with the init() subroutine before calling this subroutine.
    !! After calling this subroutine, the solution is stored in the final "ebsys" object, as far as convergence
    !! has been reached.
	!!**********************************************************************************************
    subroutine solve_fpi_ebsolve(self, ebsys)
        class(eilenberger_solver_t), intent(inout) :: self
        class(eilenberger_system_t), intent(inout) :: ebsys
        complex(wp), allocatable :: qn_old(:,:), qn_new(:,:)
		real(wp) :: qerr
		integer :: s, nx
        
        !! A. Allocate some space for the arrays:
        nx = ebsys%get_nx()   !! Extract the number of x points.
        allocate(qn_old(3, nx))
        allocate(qn_new(3, nx))
        
		!! B. Prepare for the iteration:
		self%conv = .false.  !! Resets the convergence status to False at the beginning.
        self%niter = 0       !! Reset the number of iteration to zero.
		qn_old = 0.          !! Sets the initial Qn(x) field to the best empirical guess, which is simply Qn=0.
        
		!! C. Loop on the iterations (cannot be parallelized):
		do s = 1, self%maxit
            
            !! 1. Perform a single iteration: From the current Qn(x) field, compute the next Qn(x) field.
            call ebsys%next_state()
            call ebsys%get_qn(qn_new)
            
            !! 2. Compute the relative variation of the Qn(x) field:
			qerr = norm_complex_2darray(qn_new - qn_old)/norm_complex_2darray(qn_new)
            self%niter = self%niter + 1  !! Increments the number of iterations.
            self%qerr_history(self%niter) = qerr  !! Saves the relative variation of Qn(x) in the history.
            self%rho_history(self%niter) = ebsys%get_rho()  !! Saves the current value of rho(T) in the history.
            
            !! 3. Print the current iteration:
			if (self%verbose >= 1 .and. ismultipower(s, 10)) then
				call print_iteration_ebsolve(self, ebsys)
			end if
            
            !! 4. Stopping criterion:
			if (qerr < self%qtol) then
				self%conv = .true.
				exit
			end if
            
            !! 5. Save the current field just before the next iteration:
            qn_old = qn_new
            
        end do
        
        !! D. Frees allocated memory:
        deallocate(qn_old)
        deallocate(qn_new)
		
    end subroutine
    
	!!**********************************************************************************************
	!! Solve the Eilenberger_System "ebsys" using the relaxed fixed-point iteration method.
	!! The Eilenberger_System must have been initialized with the init() subroutine before calling this subroutine.
    !! After calling this subroutine, the solution is stored in the final "ebsys" object, as far as convergence
    !! has been reached.
	!!**********************************************************************************************
    subroutine solve_relax_ebsolve(self, ebsys)
        class(eilenberger_solver_t), intent(inout) :: self
        class(eilenberger_system_t), intent(inout) :: ebsys
        complex(wp), allocatable :: qn_old(:,:), qn_new(:,:), delta_qn(:,:)
		real(wp) :: qerr
		integer :: s, nx
        
        !! A. Allocate some space for the arrays:
        nx = ebsys%get_nx()   !! Extract the number of x points.
        allocate(qn_old(3, nx))
        allocate(qn_new(3, nx))
        allocate(delta_qn(3, nx))
        
		!! B. Prepare for the iteration:
		self%conv = .false.  !! Resets the convergence status to False at the beginning.
        self%niter = 0       !! Reset the number of iteration to zero.
		qn_old = 0.          !! Sets the initial Qn(x) field to the best empirical guess, which is simply Qn=0.
        
		!! C. Loop on the iterations (cannot be parallelized):
		do s = 1, self%maxit
            
            !! 1. Perform a single iteration: From the current Qn(x) field, compute the next Qn(x) field.
            call ebsys%set_qn(qn_old)
            call ebsys%next_state()
            call ebsys%get_qn(qn_new)
            delta_qn = qn_new - qn_old  !! Compute the residual of the Qn(x) field.
            
            !! 2. Compute the relative variation of the Qn(x) field:
			qerr = norm_complex_2darray(delta_qn)/norm_complex_2darray(qn_new)
            self%niter = self%niter + 1  !! Increments the number of iterations.
            self%qerr_history(self%niter) = qerr  !! Saves the relative variation of Qn(x) in the history.
            self%rho_history(self%niter) = ebsys%get_rho()  !! Saves the current value of rho(T) in the history.
            
            !! 3. Print the current iteration:
			if (self%verbose >= 1 .and. ismultipower(s, 10)) then
				call print_iteration_ebsolve(self, ebsys)
			end if
            
            !! 4. Stopping criterion:
			if (qerr < self%qtol) then
				self%conv = .true.
				exit
			end if
            
            !! 5. Relaxation (do not accept the full step, only a fraction of it):
            qn_old = qn_old + self%frelax * delta_qn
            
        end do
        
        !! D. Frees allocated memory:
        deallocate(qn_old)
        deallocate(qn_new)
        deallocate(delta_qn)
		
    end subroutine
    
end module
