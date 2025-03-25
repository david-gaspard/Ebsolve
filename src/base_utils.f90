!!**************************************************************************************************
!! Created on 2024-07-22 at 14:02:02 CEST by David Gaspard <gaspard.dave@gmail.com>
!! This program is distributed under the Creative Commons (CC) BY-NC-SA license.
!! Fortran module providing general and useful base routines.
!!**************************************************************************************************
module base_utils
    use constants
	implicit none
	
	contains
	
	!!**********************************************************************************************
	!! Computes the Frobenius norm of the given complex 1D array.
	!!**********************************************************************************************
	function norm_complex_1darray(array) result(norm)
		complex(wp), intent(in) :: array(:)
		real(wp) :: n2, norm
		integer :: i
		n2 = 0._wp
		do i = 1, size(array)
			n2 = n2 + (array(i)%re)**2 + (array(i)%im)**2
		end do
		norm = sqrt(n2)
	end function
    
	!!**********************************************************************************************
	!! Computes the Frobenius norm of the given complex 2D array.
	!!**********************************************************************************************
	function norm_complex_2darray(array) result(norm)
		complex(wp), intent(in) :: array(:, :)
		real(wp) :: n2, norm
		integer :: i, j
		n2 = 0._wp
		do i = 1, size(array, 1)
			do j = 1, size(array, 2)
				n2 = n2 + (array(i, j)%re)**2 + (array(i, j)%im)**2
			end do
		end do
		norm = sqrt(n2)
	end function
	
    !!**********************************************************************************************
	!! Subroutine to print a 1D real array to standard output.
	!!**********************************************************************************************
    subroutine print_real_1darray(strname, fmtreal, array)
        real(wp), intent(in) :: array(:)
		character(len=*), intent(in) :: strname, fmtreal
        character(len=len(fmtreal)+40) :: fmtstr
        write (fmtstr, '(a,i0,3a)') "(a,1x,'= [',", size(array), "(1x,", fmtreal, ",1x),']')"
        write (stdout, fmtstr) trim(strname), array
    end subroutine
    
    !!**********************************************************************************************
	!! Subroutine to print a 1D complex array to standard output.
	!!**********************************************************************************************
    subroutine print_complex_1darray(strname, fmtreal, array)
        complex(wp), intent(in) :: array(:)
		character(len=*), intent(in) :: strname, fmtreal
        character(len=len(fmtreal)+40) :: fmtstr
        write (fmtstr, '(a,i0,a)') "(a,1x,'= [',", size(array), &
            "(1x," // trim(fmtreal) // ",sp," // trim(fmtreal) // ",ss,'i',1x),']')"
        write (stdout, fmtstr) trim(strname), array
    end subroutine
    
    !!**********************************************************************************************
	!! Subroutine to print a 2D real array to standard output.
	!!**********************************************************************************************
    subroutine print_real_2darray(strname, fmtreal, matrix)
        real(wp), intent(in) :: matrix(:, :)
		character(len=*), intent(in) :: strname, fmtreal
        character(len=len(fmtreal)+40) :: fmtstr
		integer :: nrow, ncol
		nrow = size(matrix, 1)
		ncol = size(matrix, 2)
        write (stdout, '(2a)') trim(strname), " = "
        
        write (fmtstr, '(a,i0,a,i0,3a)') &
            "(", nrow, "('|',", ncol, "(1x,", fmtreal, ",1x),'|',/))"
        
        !!print '(a)', "[INFO] fmtstr = " // fmtstr
        
        write (stdout, fmtstr) transpose(matrix)
        
    end subroutine
    
    !!**********************************************************************************************
	!! Subroutine to print a 2D complex array to standard output.
	!!**********************************************************************************************
	subroutine print_complex_2darray(strname, fmtreal, a)
		complex(wp), intent(in) :: a(:, :)
		character(len=*), intent(in) :: strname, fmtreal
		integer :: i, j, nrow, ncol
		nrow = size(a, 1)
		ncol = size(a, 2)
		write (*, '(2a)') trim(strname), " = "
		do i = 1, nrow
			write (*, '(a)', advance='no') "|"
			do j = 1, ncol
				write (*, "(a,"//fmtreal//",sp,"//fmtreal//",ss,'i')", advance='no') char_tab, a(i, j)
			end do
			write (*, '(a)') " |"
		end do
	end subroutine
    
	!!**********************************************************************************************
	!! Write a 2D real array "matrix" with format "fmtreal" to the output unit "ounit".
    !! Arguments:
    !!   ounit     = (in) Output stream. Typically a file or "stdout".
    !!   fmtreal   = (in) Format used for all real numbers. Typically: "g0.16".
    !!   delimiter = (in) Delimiter used to separate number in the file. Typically: ",".
    !!   matrix    = (in) Real 2D array to be written to the output stream.
	!!**********************************************************************************************
    subroutine write_real_2darray(ounit, fmtreal, delimiter, matrix)
        integer, intent(in) :: ounit
        character(len=*), intent(in) :: fmtreal, delimiter
        real(wp), intent(in) :: matrix(:, :)
        character(len=2*len(fmtreal)+len(delimiter)+35) :: fmtstr
        integer :: nrow, ncol
        
        !! Prepare the formatting string:
        nrow = size(matrix, 1)
        ncol = size(matrix, 2)
        
        write (fmtstr, '(2(a,i0),a)') &
            "(", nrow, "(", ncol-1, "(" // fmtreal // ",'" // delimiter // "')," // fmtreal // ",/))"
        
        write (ounit, fmtstr) transpose(matrix)
        
    end subroutine
    
	!!**********************************************************************************************
	!! Interpolates linearly the abscissa/ordinate data in the arrays "xmesh" and "ymesh"
	!! at position "x". The result is written in the variable "y".
	!! Inspired by: https://scicomp.stackexchange.com/questions/20960/linear-interpolation-in-fortran
	!! xmesh = Array of abscissas, not necessarily in increasing order (real numbers, input).
	!! ymesh = Array of ordinates (complex numbers, input).
	!! x     = Abscissa at which the inerpolated value is sought (real number, input).
	!!         This coordinate must belong to the "x" interval covered by xmesh (no extrapolation allowed).
	!! y     = Interpolated value of the ordinate (complex number, output).
	!!**********************************************************************************************
	subroutine interp1(xmesh, ymesh, x, y)
		real(wp), intent(in) :: x, xmesh(:)
		complex(wp), intent(in) :: ymesh(:)
		complex(wp), intent(out) :: y
		real(wp) :: x0, x1
		complex(wp) :: y0, y1
		integer :: i, ix
		!! Find the index of the point just below "x":
		do i = 1, size(xmesh) - 1
			x0 = min(xmesh(i), xmesh(i+1))
			x1 = max(xmesh(i), xmesh(i+1))
			if (x0 <= x .and. x < x1) then
				ix = i
				!!print '(a,i0)', "Found at index ix=", ix
				exit
			end if
		end do
		x0 = xmesh(ix)
		x1 = xmesh(ix + 1)
		y0 = ymesh(ix)
		y1 = ymesh(ix + 1)
		y = (y1 - y0)*(x - x0)/(x1 - x0) + y0
    end subroutine
    
    !!**********************************************************************************************
    !! Swap two real variables. This kind of subroutine can be used by all possible sorting subroutines.
    !!**********************************************************************************************
    subroutine swap_real(x, y)
		real(wp), intent(inout) :: x, y
		real(wp) :: tmp
		tmp = x
		x = y
		y = tmp
    end subroutine
    
    !!**********************************************************************************************
    !! Sort a real array using the selection sort algorithm ("ssort" stands for selection sort).
    !!**********************************************************************************************
    subroutine ssort_real(array)
		real(wp), intent(inout) :: array(:)
		real(wp) :: amin
		integer :: i, j, k, n
		n = size(array)
		do i = 1, n-1
			k = i
			amin = array(k)
			do j = i+1, n   !! Search for the minimum value of the sub-array.
				if (array(j) < amin) then
					k = j
					amin = array(j)
				end if
			end do
			if (k /= i) then
				!!print *, "Array before: ", array
				!!print '(a,i0,a,i0)', "Swap indices: ", i, " <-> ", k
				call swap_real(array(i), array(k))
			end if
		end do
    end subroutine
    
    !!**********************************************************************************************
    !! Finds the distinct values of the given real array and stores them in "array_uniq".
    !! Also counts their multiplicities and stores them into "array_cnts".
    !! The output arrays will be allocated to the appropriate length.
    !! WARNING: This subroutine assumes that "array" is sorted so that duplicate values are
    !! next to each other ! Otherwise, the behavior of this subroutine is undefined.
    !!**********************************************************************************************
	subroutine uniq_real(array, array_uniq, array_cnts)
		real(wp), intent(in) :: array(:)
		real(wp), allocatable, intent(out) :: array_uniq(:)
		integer, allocatable, intent(out) :: array_cnts(:)
		integer :: i, n, nuniq
		if (allocated(array_uniq)) deallocate(array_uniq)
		if (allocated(array_cnts)) deallocate(array_cnts)
		n = size(array)
		!! 1. First compute the number of unique elements:
		nuniq = 0
		if (n /= 0) nuniq = 1
		do i = 2, n
			if (array(i) /= array(i-1)) nuniq = nuniq + 1
		end do
		allocate(array_uniq(nuniq))
		allocate(array_cnts(nuniq))
		!! 2. Then stores the unique elements:
		nuniq = 0
		if (n /= 0) then
			array_uniq(1) = array(1)
			array_cnts = 0  !! Initialize the counts to zero.
			array_cnts(1) = 1
			nuniq = 1
		end if
		do i = 2, n
			if (array(i) /= array(i-1)) then
				nuniq = nuniq + 1
				array_uniq(nuniq) = array(i)
			end if
			array_cnts(nuniq) = array_cnts(nuniq) + 1
		end do
    end subroutine
    
	!!**********************************************************************************************
	!! Test if the given integer "i" is a multiple of a power of "base".
	!! Returns True if it is so, False otherwise.
	!! This function is mainly used to reduce the amount of printing in an iteration sequence.
	!!**********************************************************************************************
	function ismultipower(i, base) result(bool)
		integer, intent(in) :: i, base
		logical :: bool
		integer :: k, kdiv, nnz
		nnz = 0  !! Count the number of non-zero digit.
		k = i
		do while(k /= 0)
			kdiv = k/base  !! Integer division.
			if (k /= kdiv*base) nnz = nnz + 1
			k = kdiv
		end do
		bool = (nnz == 1)
	end function
    
    !!**********************************************************************************************
    !! Fill the array "xmesh" of Chebyshev nodes of the first kind in the interval [xmin, xmax].
    !! By construction, the nodes never reach exactly the bounds "xmin" and "xmax".
    !! The array "xmesh" is assumed to be already allocated.
    !!**********************************************************************************************
    subroutine fill_chebyshev_1(xmin, xmax, xmesh)
        real(wp), intent(in) :: xmin, xmax
        real(wp), intent(out) :: xmesh(:)
        real(wp) :: center, extent 
        integer :: i, n
        
        center = (xmax + xmin)/2
        extent = (xmax - xmin)/2.
        n = size(xmesh)
        do i = 1, n
            xmesh(n - i + 1) = extent*cos(pi*(i - 0.5)/n) + center
        end do
        
    end subroutine
    
    !!**********************************************************************************************
    !! Trim the spaces and zeros at the end of the string, and the spaces at the beginning of the string.
    !! This function does not need to be passed to trim().
    !!**********************************************************************************************
    function trim_zero(string) result(res)
        character(len=*), intent(in) :: string
        character(len=:), allocatable :: res
        integer :: imin, imax
        imin = 1
        do while (string(imin:imin) == " ")
            imin = imin + 1
        end do
        imax = len(string)
        do while (string(imax:imax) == "0" .or. string(imax:imax) == " ")
            imax = imax - 1
        end do
        res = string(imin:imax)
    end function
    
	!!**********************************************************************************************
	!! Returns a string representation of the current time.
	!!**********************************************************************************************
    function datetime_string() result(str)
        character(len=30) :: str
        character(len=6) :: zone
        integer :: dt(8)
        
        call date_and_time(values=dt, zone=zone)
        write (str, '(i4,a,5(i0.2,a))') &
            dt(1), "-", dt(2), "-", dt(3), " at ", &
            dt(5), ":", dt(6), ":", dt(7), " " // trim(zone)
        
    end function
    
    !!**********************************************************************************************
    !! Print a timestamp in the stream "ounit", typically the file containing the simulation results.
    !! The string "prefix" is the comment character.
    !!**********************************************************************************************
    subroutine print_timestamp(ounit, prefix)
        integer, intent(in) :: ounit
        character(len=*), intent(in) :: prefix
        
        write (ounit, '(a)') &
            trim(prefix) // " Computed on " // trim(datetime_string()) // " by " // trim(program_copyright) 
        
    end subroutine
    
    !!**********************************************************************************************
	!! Returns the given time duration "time" (in seconds) in a string of the format "DDDd HH:MM:SS".
	!!**********************************************************************************************
    function format_time(time) result(str)
        real(wp), intent(in) :: time   !! Time duration in seconds.
        character(len=30) :: str
        real(wp) :: t
        integer :: day, hour, minute, second
        
        day = floor(time/86400)
        t = time - 86400*day
        hour = floor(t/3600)
        t = t - 3600*hour
        minute = floor(t/60)
        second = floor(t - 60*minute)
        
        if (day == 0) then  !! If day=0, then print a short version:
            write (str, '(i0.2,":",i0.2,":",i0.2)') hour, minute, second
        else  !! If day is not zero, then print a longer version:
            write (str, '(i0,"d ",i0,"h ",i0,"m ",i0,"s")') day, hour, minute, second
        end if
        
    end function
	
    !!**********************************************************************************************
    !! Prints the progress bar to standard output with the expected time of arrival.
    !! cjob  = Index of the current completed job [1, ..., njob].
    !! njob  = Total number of jobs.
    !! start = Date and time array given by the subroutine: date_and_time(values=start).
    !! msg   = Short message printed before the progress bar.
    !!**********************************************************************************************
    subroutine print_progress_bar(cjob, njob, start, msg)
        integer, intent(in) :: cjob, njob, start(8)
        character(len=*), intent(in) :: msg
        integer, parameter :: progress_bar_length = 50
        integer :: nchar, now(8)
        real(wp) :: delta(8), time, speed, eta
        
        !! 1. Compute the expected time of arrival:
        call date_and_time(values=now)
        delta = real(now - start, kind=wp)
        time = 86400*delta(3) + 3600*delta(5) + 60*delta(6) + delta(7) + delta(8)/1000  !! Elapsed time in seconds.
        speed = cjob/time           !! Compute the mean speed in job per second.
        eta = (njob - cjob)/speed   !! Compute the expected time of arrival (in seconds).
        
        !! 2. Print the progress bar:
        nchar = (progress_bar_length*cjob)/njob   !! Integer division.
        write(stdout, '(a,f0.2,a,f0.1,a)', advance='no') &
            tag_info // trim(msg) // " [" // repeat("#", nchar) // &
            repeat(" ", progress_bar_length - nchar) // "] ", 100.0*cjob/njob, "% ETA " // &
            trim(format_time(eta)) // " (", speed, " job/s)" // char_cr
        
    end subroutine
    
    !!**********************************************************************************************
    !! Ends the progress bar by printing the total computation time to std output.
    !! This subroutine also compute the total computation time "time" (in seconds).
    !!**********************************************************************************************
    subroutine end_progress_bar(njob, start, time)
        integer, intent(in) :: njob, start(8)
        real(wp), intent(out) :: time
        integer :: end_time(8)
        real(wp) :: delta(8), speed
        
        call date_and_time(values=end_time)      !! End the time measurement.
        delta = real(end_time - start, kind=wp)  !! Time difference.
        time = 86400*delta(3) + 3600*delta(5) + 60*delta(6) + delta(7) + delta(8)/1000  !! Convert the computation time in seconds.
        speed = njob/time   !! Average computation speed.
        print '(/,a,2(f0.2,a))', tag_info // "Done, computation time is ", time, " s, avg speed is ", speed, " job/s"
        
    end subroutine
    
	!!**********************************************************************************************
	!! Subroutine to generate a new filename by checking the availability of the filename.
	!! The template is "base_x.ext" where "x" is a number from 1 to MAX_INT.
	!! If the file already exists, then increments the index "x" until the filename does not exist.
	!! The purpose is of course to avoid overwriting existing data.
	!! base     = Basename of the file, including the folder.
	!! ext      = File extension, including the dot prefix.
	!! filename = Resulting output filename, guaranteed to exist if no error is sent.
	!!**********************************************************************************************
	subroutine unique_filename(base, ext, filename)
		character(len=*), intent(in) :: base, ext
		character(len=*), intent(out) :: filename
		integer, parameter :: nmax = 1000000000  !! Large value still representable on int32.
		integer :: i
		logical :: fexist
		do i = 1, nmax
			write (filename, '(a,i0,a)') trim(base), i, trim(ext)
			!!print '(a)', tag_info // "Testing filename '" // filename // "'..."
			inquire(file=filename, exist=fexist)
			if (.not. fexist) return
		end do
		!! The loop should never reach this point:
		write (stderr, '(a,i0,a)') tag_error // "Too many files (nmax=", nmax, "), no valid filename !"
		stop errmsg_invalid_arg
	end subroutine
	
	!!**********************************************************************************************
	!! Check if the given output directory "dirname" exists. If not, it will be created.
    !! WARNING : This subroutine requires the "mkdir" command !
	!!**********************************************************************************************
	subroutine create_directory(dirname)
		character(len=*), intent(in) :: dirname
		call execute_command_line("mkdir -p " // trim(dirname))
	end subroutine
    
	!!**********************************************************************************************
	!! If the IO status flag 'istat' is zero, then tries to rewind the file.
	!! Otherwise, displays the invalid entry in the namelist.
	!! In principle, this routine should be called after each namelist reading.
	!!**********************************************************************************************
    subroutine check_and_rewind_nml(funit, istat)
		integer, intent(in) :: funit
		integer, intent(in) :: istat
		integer :: rewind_stat
		character(len=150) :: line
		if (istat == 0) then
			!!print '(a)', tag_info // "Rewinding the namelist file..."
			rewind (unit=funit, iostat=rewind_stat)
			if (rewind_stat /= 0) then
				print '(a)', tag_error // "Error in rewinding."
				stop errmsg_invalid_arg
			end if
		else
			backspace(funit)  !! Go back one line.
			read (funit, '(a)') line
			write (stderr, '(a)') tag_error // "Invalid entry in namelist: " // trim(line) // "..."
			if (line == "/") then
				write (stderr, '(a)') tag_warn // "Try checking the spelling of the namelist groups..."
			end if
			stop errmsg_invalid_arg
		end if
	end subroutine
    
end module
