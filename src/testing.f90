!!**************************************************************************************************
!! Created on 2024-07-21 at 13:02:08 CEST by David Gaspard <david.gaspard@espci.fr>
!! This program is distributed under the MIT License.
!! Fortran module providing subroutine useful for unit testing.
!!**************************************************************************************************
module testing
    use constants
	implicit none
    
	!! Constant tags reserved for tests:
    character(len=*), parameter :: tag_test  = "["//tcolor_bold//tcolor_cyan//"TEST"//tcolor_nc//"] "
    character(len=*), parameter :: tag_fail  = "["//tcolor_bold//tcolor_lightred//"FAIL"//tcolor_nc//"] "
	character(len=*), parameter :: tag_pass  = "["//tcolor_bold//tcolor_lightgreen//" OK "//tcolor_nc//"] "
    
	contains
	
    !!**********************************************************************************************
    !! Subroutint to print the title of a test, hence avoiding code duplication in tests.
    !!**********************************************************************************************
    subroutine print_test_title(title)
        character(len=*), intent(in) :: title
        
        write(stdout, '(a)') tag_test // "====== TEST : " // trim(title) // " ======"
        
    end subroutine
    
	!!**********************************************************************************************
	!! Safely compare two real numbers and check if they are equal to within the given tolerance 'tol'.
	!! If not, prints them.
	!! This function returns True if the numbers are equal, False otherwise.
	!!**********************************************************************************************
	function assert_equal_real(x1, x2, tol) result(equal)
		real(wp), intent(in) :: x1, x2, tol
		logical :: equal
		equal = .true.
		if (abs(x1 - x2) > tol*(abs(x1) + abs(x2))/2) then
			print '(a,g0,a,g0,a)', tag_fail // "Numbers are too different: x1=", x1, ", and x2=", x2, "."
			equal = .false.
		end if
	end function
	
	!!**********************************************************************************************
	!! Safely compare two complex numbers and check if they are equal to within the given tolerance 'tol'.
	!! If not, prints them.
	!! This function returns True if the numbers are equal, False otherwise.
	!!**********************************************************************************************
	function assert_equal_complex(z1, z2, tol) result(equal)
		complex(wp), intent(in) :: z1, z2
		real(wp), intent(in) :: tol
		logical :: equal
		equal = .true.
		if (abs(z1 - z2) > tol*(abs(z1) + abs(z2))/2) then
			print '(a,g0,sp,g0,ss,"i",a,g0,sp,g0,ss,"i",a)', &
                tag_fail // "Numbers are too different: z1=", z1, ", and z2=", z2, "."
            !!print '(a,g0,a,g0,a)', tag_fail // "Relative error=", abs(z1 - z2)/((abs(z1) + abs(z2))/2), ", tolerance=", tol, "."
			equal = .false.
		end if
	end function
	
	!!**********************************************************************************************
	!! Safely compare two real 1D arrays and check if they are equal to within the given tolerance 'tol'.
	!! If not prints a summary of the differences.
	!! This function returns True if the arrays are equal, False otherwise.
	!!**********************************************************************************************
	function assert_equal_real_1darray(array1, array2, tol) result(equal)
		real(wp), intent(in) :: array1(:), array2(:), tol
		integer :: i, len1, len2
		logical :: equal
		equal = .true.
		len1 = size(array1)
		len2 = size(array2)
		if (len1 /= len2) then
			print '(a,i0,a,i0,a)', tag_fail // "Arrays are of different size: len(array1)=", len1, &
				" and len(array2)=", len2, "."
			equal = .false.
			return
		end if
		do i = 1, len1
			!! Safe test avoiding possible division by zero:
			if (abs(array1(i) - array2(i)) > tol*(abs(array1(i)) + abs(array2(i)))/2) then
				print '(a,i0,a,g0,a,i0,a,g0,a)', &
					tag_fail // "Elements are too different: array1(", i, ")=", array1(i), &
					", and array2(", i, ")=", array2(i), "."
				equal = .false.
			end if
		end do
	end function
	
	!!**********************************************************************************************
	!! Safely compare two complex 1D arrays and check if they are equal to within the given tolerance 'tol'.
	!! If not prints a summary of the differences.
	!! This function returns True if the arrays are equal, False otherwise.
	!!**********************************************************************************************
	function assert_equal_complex_1darray(array1, array2, tol) result(equal)
		complex(wp), intent(in) :: array1(:), array2(:)
		real(wp), intent(in) :: tol
		integer :: i, len1, len2
		logical :: equal
		equal = .true.
		len1 = size(array1)
		len2 = size(array2)
		if (len1 /= len2) then
			print '(a,i0,a,i0,a)', tag_fail // "Arrays are of different size: len(array1)=", len1, &
				" and len(array2)=", len2, "."
			equal = .false.
			return
		end if
		do i = 1, len1
			!! Safe test avoiding possible division by zero:
			if (abs(array1(i) - array2(i)) > tol*(abs(array1(i)) + abs(array2(i)))/2) then
				print '(a,i0,a,g0,sp,g0,"i",a,i0,a,g0,sp,g0,"i",a)', &
					tag_fail // "Elements are too different: array1(", i, ")=", array1(i), &
					", and array2(", i, ")=", array2(i), "."
				equal = .false.
			end if
		end do
	end function
	
end module
