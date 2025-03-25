!!**************************************************************************************************
!! Created on 2024-07-26 at 12:54:33 CEST by David Gaspard <gaspard.dave@gmail.com>
!! This program is distributed under the Creative Commons (CC) BY-NC-SA license.
!! Fortran program Ebsolve (v2). This program solves the Eilenberger equation, a self-consistent matrix
!! transport equation, associated to the full counting statistics problem, i.e., the determination of
!! the eigenvalue distribution of the transmission matrix t⁺t in a disordered waveguide or a slab.
!!**************************************************************************************************
program main
    use omp_lib       !! OpenMP library for parallelization (compile with -fopenmp).
    use eilenberger_solver
	implicit none
	character(len=150) :: arg  !! This assumes that filenames are not too long.
	integer :: i, funit, rstat
	logical :: fexist
	
	write (stdout, '(a)') "====== This is " // trim(program_copyright) // " ======"
	if (command_argument_count() == 0) then
		write (stderr, '(a)') tag_error // "No settings file provided, doing nothing..."
	end if
	do i = 1, command_argument_count()
		call get_command_argument(i, arg)
		inquire(file=arg, exist=fexist)
		if (fexist) then
			open(newunit=funit, file=arg, action='read', iostat=rstat)
			if (rstat == 0) then
				write (stdout, '(a)') tag_info // "Reading file: " // trim(arg)
				call execute_file(funit)
				close(funit)
			else
				write (stderr, '(a)') tag_error // "Cannot open file: " // trim(arg)
			end if
		else
			write (stderr, '(a)') tag_error // "File not found: " // trim(arg)
		end if
	end do
	
	contains
	
	!!**********************************************************************************************
	!! Executes the instructions given in the given namelist file, assuming this file exists and can be opened.
	!!**********************************************************************************************
	subroutine execute_file(funit)
		integer, intent(in) :: funit
		type(eilenberger_system_t) :: ebsys
        type(eilenberger_solver_t) :: solver
		integer :: istat
		character(len=150) :: task
		namelist /main_settings/ task
		
		!! 1. Read main settings and check for errors:
		read(unit=funit, nml=main_settings, iostat=istat)
		call check_and_rewind_nml(funit, istat)  !! Check for reading error and rewind the namelist file (for further reading).
		
        !! 2. Read the settings from the file and print them:
        call ebsys%parse(funit)
        call solver%parse(funit)
        write (stdout, '(a)') trim(ebsys%to_string("|", char_nl)) // trim(solver%to_string("|"))
        
		!! 3. Execute the appropriate task:
        print '(a)', tag_info // "Running task '" // trim(task) // "'..."
		select case (task)
            case ("distrib")
                call task_distrib(funit, ebsys, solver)  !! Compute the distribution rho(T).
            case ("fields")
                call task_fields(funit, ebsys, solver)  !! Compute only the fields for given "T".
            case default
                write (stderr, '(a)') tag_error // "Unknown task: " // trim(task)
                stop errmsg_invalid_arg
		end select
        
        !! 4. Finally, deallocates memory (see parse(funit) methods):
        call ebsys%del()
        call solver%del()
		
	end subroutine
	
    !!**********************************************************************************************
    !! Parse the parameters of the "distrib" task from the file "funit".
    !! Arguments:
    !!   funit    = (in)  Input file where the namelist is read.
    !!   rhodata  = (out) 2D array "rhodata" containing the points [T, rho(T), Niter(T)]. On output, the column "T" is initialized.
    !!   geps     = (out) Imaginary part of "gamma", the parameter of the generating function of rho(T). gamma = 1/T + geps*i.
    !!   nthreads = (out) Number of threads used for multithreading.
    !!**********************************************************************************************
    subroutine parse_task_distrib(funit, rhodata, geps, nthreads)
        integer, intent(in) :: funit
        real(wp), allocatable, intent(out) :: rhodata(:, :)
        real(wp), intent(out) :: geps
        integer, intent(out) :: nthreads
        integer, parameter :: ntm_min = 1, nthreads_min = 1
        integer :: ntm, istat
        real(wp) :: tmin, tmax
        namelist /distrib/ ntm, tmin, tmax, geps, nthreads
        
        !! 1. Parse the "distrib" namelist:
        read(unit=funit, nml=distrib, iostat=istat)
		call check_and_rewind_nml(funit, istat)
        
        !! 2. Check for possible invalid parameters:
        if (ntm < ntm_min) then
            write (stderr, '(2(a,i0),a)') tag_error // &
                "Invalid number of samples of T, received ", ntm, ", expected at least ", ntm_min, "."
            stop errmsg_invalid_arg
        else if (tmax <= tmin) then
            write (stderr, '(a)') tag_error // "Invalid interval of T, expected tmin < tmax."
            stop errmsg_invalid_arg
        else if (nthreads < nthreads_min) then
            write (stderr, '(2(a,i0),a)') tag_error // &
                "Invalid number of threads, received ", nthreads, ", expected at least ", nthreads_min, "."
            stop errmsg_invalid_arg
        end if
        
        !! 3. Prepare the "rhodata" array for thread columns [T, rho(T), Niter(T)]:
        allocate(rhodata(ntm, 3))
        call fill_chebyshev_1(tmin, tmax, rhodata(:, 1))
        
    end subroutine
    
	!!**********************************************************************************************
	!! Compute the transmission-eigenvalue density rho(T) for various values of the transmission eigenvalue T
    !! and writes the results to files.
	!!**********************************************************************************************
    subroutine task_distrib(funit, ebsys, solver)
        integer, intent(in) :: funit
		class(eilenberger_system_t), intent(inout) :: ebsys
		class(eilenberger_solver_t), intent(inout) :: solver
        real(wp), allocatable :: rhodata(:, :)
        real(wp) :: geps, time, speed
        integer :: nthreads, time_start(8), time_end(8), delta(8)
        
        !! 1. Parse the "distrib" task parameters and prepare the arrays:
        call parse_task_distrib(funit, rhodata, geps, nthreads)
        
        !! 2. Main computational loop to compute the distribution rho(T):
        call date_and_time(values=time_start)
        if (nthreads == 1) then
            call compute_distrib_serial(ebsys, solver, geps, rhodata)
        else
            call compute_distrib_omp(ebsys, solver, geps, nthreads, rhodata)
        end if
        
        !! 3. Determine the computation time (in seconds):
        call date_and_time(values=time_end)
        delta = time_end - time_start
        time = (( delta(3)*24 + delta(5) )*60 + delta(6) )*60 + delta(7) + 0.001_wp*delta(8)
        speed = size(rhodata, 1)/time
        print '(2(a,f0.3),a)', tag_info // "Done, computation time is ", time, " s, avg speed is ", speed, " job/s."
        
        !! 3. Output the distribution:
        call output_task_distrib(ebsys, solver, geps, nthreads, time, rhodata)
        
        deallocate(rhodata)  !! Deallocates the memory.
        
    end subroutine
    
    !!**********************************************************************************************
    !! Compute the transmission-eigenvalue distribution using a serial loop.
    !!**********************************************************************************************
    subroutine compute_distrib_serial(ebsys, solver, geps, rhodata)
        class(eilenberger_system_t), intent(inout) :: ebsys
		class(eilenberger_solver_t), intent(inout) :: solver
        real(wp), intent(in) :: geps
        real(wp), intent(inout) :: rhodata(:, :)
        integer :: i, ntm, niter
        real(wp) :: tm, rho, qerr
        character(len=30) :: info
        
        ntm = size(rhodata, 1)  !! Number of points to compute.
        
        do i = 1, ntm
            tm = rhodata(i, 1)                !! Extract the desired transmission eigenvalue from the data array.
            call ebsys%set_transmission(tm, geps)  !! Assigns the desired transmission T to the Eilenberger_System.
            call solver%solve(ebsys)               !! Solve the Eilenberger equation iteratively.
            rho = ebsys%get_rho()                  !! Extract the resulting density rho(T).
            niter = solver%get_niter()             !! Extract the number of iterations needed to solve the Eilenberger equation.
            qerr = solver%get_qerr(niter)          !! Extract the final relative variation of the Qn(x) field.
            rhodata(i, 2) = rho               !! Saves the resulting density rho(T) to the array distrib_data.
            rhodata(i, 3) = niter             !! Saves the number of iterations to the array distrib_data.
            
            !! Use color to help visually detect convergence problems:
            if (solver%converged()) then
                info = tcolor_lightgreen // "Found solution" // tcolor_nc
            else
                info = tcolor_bold // tcolor_lightred // "Not converged!"// tcolor_nc
            end if
            
            print '(a,i0,a,g13.6,3x,a,3x,a,g19.12,3x,a,g10.3,3x,a,i0,a)', tag_info // "#", i, &
                char_tab // "| T=", tm, trim(info), "rho=", rho, "qerr=", qerr, "niter=", niter
            
        end do
        
    end subroutine
    
	!!**********************************************************************************************
	!! Compute the transmission-eigenvalue density rho(T) for various values of the transmission eigenvalue T
    !! and writes the results to files.
	!!**********************************************************************************************
    subroutine compute_distrib_omp(ebsys_shared, solver_shared, geps, nthreads, rhodata)
		class(eilenberger_system_t), intent(in) :: ebsys_shared
		class(eilenberger_solver_t), intent(in) :: solver_shared
        real(wp), intent(in) :: geps
        integer, intent(in) :: nthreads
        real(wp), intent(inout) :: rhodata(:, :)
        type(eilenberger_solver_t) :: solver  !! Declare local (OMP-Private) objects of type Eilenberger_Solver.
        type(eilenberger_system_t) :: ebsys   !! Declare local (OMP-Private) objects of type Eilenberger_System.
        real(wp) :: tm, rho, qerr  !! Various OMP-Private variables.
        integer :: i, niter        !! Various OMP-Private variables.
        character(len=30) :: msg   !! Short information message (OMP-Shared).
        integer :: ntm, cjob, start_time(8)  !! Various OMP-Shared variables.
        
        ntm = size(rhodata, 1)  !! Number of points to compute.
        
        write (msg, '(a,i0,a)') "distrib_omp, ", nthreads, " thr"
        call date_and_time(values=start_time)   !! Start the time measurement. date_and_time() is the true wall time even with parallelization.
        cjob = 0    !! Initialize the number of completed jobs (OMP-Shared).
        
        !! 3. Main computational loop to compute the distribution rho(T):
        !$omp parallel num_threads(nthreads) default(shared) private(ebsys, solver, tm, rho, niter, qerr)
            
            solver = solver_shared  !! Create a local copy of the solver, one per thread (implicit allocation).
                                    !! Needed because this object stores the history of the iterations, and this is OMP-Private.
            ebsys = ebsys_shared    !! Initialize the Eilenberger_System (memory allocation). This is also an OMP-Private step.
            
            !$omp do schedule(dynamic, 1)
            do i = 1, ntm
                tm = rhodata(i, 1)                     !! Extract the desired transmission eigenvalue from the data array.
                call ebsys%set_transmission(tm, geps)  !! Assigns the desired transmission T to the Eilenberger_System.
                call solver%solve(ebsys)               !! Solve the Eilenberger equation iteratively.
                rho = ebsys%get_rho()                  !! Extract the resulting density rho(T).
                niter = solver%get_niter()             !! Extract the number of iterations needed to solve the Eilenberger equation.
                qerr = solver%get_qerr(niter)          !! Extract the final relative variation of the Qn(x) field.
                rhodata(i, 2) = rho                    !! Saves the resulting density rho(T) to the array distrib_data.
                rhodata(i, 3) = niter                  !! Saves the number of iterations to the array distrib_data.
                
                !! Critical section to deal with the progress bar:
                !$omp critical
                cjob = cjob + 1      !! Use a critical section to print a progress bar.
                call print_progress_bar(cjob, ntm, start_time, msg)  !! Print the progress bar with the expected time of arrival.
                !$omp end critical
                
            end do
            !$omp end do
            
            call ebsys%del()  !! Deallocate memory allocated by all the Eilenberger_System (OMP-Private step).
            call solver%del() !! Deallocate memory allocated by all the Eilenberger_Solver (OMP-Private step).
            
        !$omp end parallel
        
        print '()'  !! Print a new line after the progress bar.
        
    end subroutine
    
    !!**********************************************************************************************
    !! Save the output of the "distrib" task to a CSV file, generate a TikZ file, and compile the figure.
    !!**********************************************************************************************
    subroutine output_task_distrib(ebsys, solver, geps, nthreads, time, rhodata)
        class(eilenberger_system_t), intent(in) :: ebsys
		class(eilenberger_solver_t), intent(in) :: solver
        real(wp), intent(in) :: geps, time
        integer, intent(in) :: nthreads
        real(wp), intent(in) :: rhodata(:, :)
        character(len=*), parameter :: prefix = "%%", delimiter = ", ", fmtreal = "g0.16"
        character(len=100) :: outputdir, filename
        character(len=500) :: cmd, title_tikz
        integer :: outfp
        
        !! 1. Create a unique filename and open it:
        call ebsys%output_directory(outputdir, pathsep)  !! Extract the suggested output directory name.
        outputdir = trim(outputdir) // "distrib" // pathsep
        call create_directory(outputdir)  !! Ensure that the directory exists.
        call unique_filename(trim(outputdir) // "result_", ".csv", filename)
        open(newunit=outfp, file=filename, status='new', action='write')
        
        !! 2. Write the parameters in the header:
        call print_timestamp(outfp, prefix)
        write (outfp, '(a)') trim(ebsys%to_string(prefix, char_nl)) // trim(solver%to_string(prefix))
        write (outfp, '(a,i0,a,g0.6,a,i0,a,f0.3,a)') &
            prefix // " Task (ntm=", size(rhodata, 1), ", geps=", geps, ", nthreads=", nthreads, ", time=", time, "s)"
        write (outfp, '(a)') "tm, rho, niter"
        
        !! 3. Write the data and close the file:
        call write_real_2darray(outfp, fmtreal, delimiter, rhodata)
        close(outfp)
        
        print '(a)', tag_info // "Data written to file: '" // trim(filename) // "'..."
        
        !! 4. Write a TikZ file and compile it:
        write (title_tikz, '(a,i0,a,g0.6,a,i0,a,f0.3,a)') &
            trim(ebsys%to_string("", "\\ ")) // trim(solver%to_string("")) // &
            "\\ Task (ntm=", size(rhodata, 1), ", geps=", geps, ", nthreads=", nthreads, ", time=", time, "s)"
        
        cmd = "./tikz/task_distrib.py " // trim(filename) // " '" // trim(title_tikz) // "' "
        
        !!print '(a)', tag_info // "TITLE: '" // trim(title_tikz) // "'."
        !!print '(a)', tag_info // "COMMAND: '" // trim(cmd) // "'."
        
        call execute_command_line(cmd)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Execute the task "fields", i.e., compute the fields Qn(x) and Jn(x) for a given transmission eigenvalue "T".
    !!**********************************************************************************************
    subroutine task_fields(funit, ebsys, solver)
        integer, intent(in) :: funit
		class(eilenberger_system_t), intent(inout) :: ebsys
		class(eilenberger_solver_t), intent(inout) :: solver
        real(wp), allocatable :: xmesh(:), fieldsdata(:, :)
        complex(wp), allocatable :: qn(:, :), jn(:, :)
        real(wp) :: tm, geps, rho, qerr
        integer :: i, nx, niter, outfp, istat
        character(len=*), parameter :: prefix = "%%", delimiter = ", ", fmtreal = "g0.16"
        character(len=*), parameter :: script_name = "./tikz/task_fields.py"  !! Location of script to be executed (relative to the root directory of the program).
        character(len=150) :: outputdir, filename
        character(len=800) :: title_tikz
        character(len=len(script_name)+len(filename)+len(title_tikz)+20) :: cmd
        
        namelist /fields/ tm, geps
        
        !! 1. Parse the "distrib" namelist:
        read(unit=funit, nml=fields, iostat=istat)
		call check_and_rewind_nml(funit, istat)
        
        !! 2. Check for possible invalid parameters:
        if (tm < 0.0_wp .or. 1.0_wp < tm) then
            write (stderr, '(2(a,i0),a)') tag_error // &
                "Invalid transmission eigenvalue, received ", tm, ", expected in [0, 1], aborting..."
            stop errmsg_invalid_arg
        end if
        
        !! 3. Solves the Eilenberger equation and extract the fields Qn(x) and Jn(x):
        call ebsys%set_transmission(tm, geps)  !! Assigns the desired transmission T to the Eilenberger_System.
        call solver%solve(ebsys)               !! Solve the Eilenberger equation iteratively.
        nx = ebsys%get_nx()                    !! Extract the number of "x" points for allocation.
        allocate(xmesh(nx))                    !! Allocate space for the xmesh.
        allocate(qn(3, nx))                    !! Allocate space for Qn(x).
        allocate(jn(3, nx))                    !! Allocate space for Jn(x).
        call ebsys%get_xmesh(xmesh)            !! Extract the xmesh.
        call ebsys%get_qn(qn)                  !! Extract the Qn(x) field.
        call ebsys%get_jn(jn)                  !! Extract the Jn(x) matrix current.
        rho = ebsys%get_rho()                  !! Extract the resulting density rho(T).
        niter = solver%get_niter()             !! Extract the number of iterations needed to solve the Eilenberger equation.
        qerr = solver%get_qerr(niter)          !! Extract the final relative variation of the Qn(x) field.
        
        !! 4. Put the data in an exportable array:
        allocate(fieldsdata(nx, 13))  !! x, ReQn11, ImQn11, ReQn12, ImQn12, ReQn21, ImQn21, ReJn11, ImJn11, ReJn12, ImJn12, ReJn21, ImJn21.
        do i = 1, nx
            fieldsdata(i, :) = [xmesh(i), &
                qn(1, i)%re, qn(1, i)%im, qn(2, i)%re, qn(2, i)%im, qn(3, i)%re, qn(3, i)%im, &
                jn(1, i)%re, jn(1, i)%im, jn(2, i)%re, jn(2, i)%im, jn(3, i)%re, jn(3, i)%im  ]
        end do
        
        !! 5. Create a unique filename and open it:
        call ebsys%output_directory(outputdir, pathsep)  !! Extract the suggested output directory name.
        outputdir = trim(outputdir) // "fields" // pathsep
        call create_directory(outputdir)  !! Ensure that the directory exists.
        call unique_filename(trim(outputdir) // "result_", ".csv", filename)
        open(newunit=outfp, file=filename, status='new', action='write')
        
        !! 6. Write the parameters in the header:
        call print_timestamp(outfp, prefix)
        write (outfp, '(a)') trim(ebsys%to_string(prefix, char_nl)) // trim(solver%to_string(prefix))
        write (outfp, '(3(a,g0.6),a,i0,a,g0.6,a)') &
            prefix // " Task (Fields, tm=", tm, ", geps=", geps, ")  Results (rho=", rho, ", niter=", niter, ", qerr=", qerr, ")"
        write (outfp, '(a)') "x, reqn11, imqn11, reqn12, imqn12, reqn21, imqn21, rejn11, imjn11, rejn12, imjn12, rejn21, imjn21"
        
        !! 7. Write the data and close the file:
        call write_real_2darray(outfp, fmtreal, delimiter, fieldsdata)
        close(outfp)
        
        print '(a)', tag_info // "Data written to file: '" // trim(filename) // "'..."
        
        !! 8. Write a TikZ file and compile it:
        write (title_tikz, '(3(a,g0.6),a,i0,a,g0.6,a)') &
            trim(ebsys%to_string("", "\\ ")) // trim(solver%to_string("")) // &
            "\\ Task (Fields, tm=", tm, ", geps=", geps, ")  Results (rho=", rho, ", niter=", niter, ", qerr=", qerr, ")"
        
        cmd = script_name // " " // trim(filename) // " '" // trim(title_tikz) // "' "
        
        call execute_command_line(cmd)
        
        deallocate(xmesh)
        deallocate(qn)
        deallocate(jn)
        deallocate(fieldsdata)
        
    end subroutine
    
end program
