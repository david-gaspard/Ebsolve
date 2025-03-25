!!**************************************************************************************************
!! Created on 2024-07-21 at 13:07:32 CEST by David Gaspard <david.gaspard@espci.fr>
!! This program is distributed under the MIT License.
!! Fortran module providing the global constants of the Ebsolve program (v2).
!!**************************************************************************************************
module constants
    use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64, input_unit, output_unit, error_unit
    implicit none
    
    !! Intrinsic constants:
    integer, parameter :: wp     = real64       !! Current working precision (wp) is double precision.
    integer, parameter :: stdin  = input_unit   !! Standard input.
    integer, parameter :: stdout = output_unit  !! Standard output display.
    integer, parameter :: stderr = error_unit   !! Standard error display.
    
    !! Useful mathematical constants:
    real(wp), parameter :: pi = acos(-1.0_wp)      !! Pi constant at working precision.
    complex(wp), parameter :: iu = (0._wp, 1._wp)  !! Imaginary unit (iu).
    
    !! Global strings regarding the program:
    character(len=*), parameter :: program_shortname = "Ebsolve v1.0"
    character(len=*), parameter :: program_fullname  = program_shortname // " - Program to solve the Eilenberger equation"
    character(len=*), parameter :: program_copyright = program_shortname // " (c) 2024 David GASPARD <david.gaspard@espci.fr>"
    
    !! Special characters:
    character, parameter :: char_tab = char(9)       !! Tab character.
    character, parameter :: char_nl = new_line('a')  !! New line character.
    character, parameter :: char_cr = achar(13)      !! Carriage return character (used for progress bar).
    character, parameter :: pathsep = '/'            !! Path separator used on UNIX filesystems.
    
    !! Define colors (for TTY only):
    character(len=*), parameter :: tcolor_nc         = achar(27)//"[0m"
    character(len=*), parameter :: tcolor_bold       = achar(27)//"[1m"
    character(len=*), parameter :: tcolor_red        = achar(27)//"[31m"
    character(len=*), parameter :: tcolor_lightred   = achar(27)//"[91m"
	character(len=*), parameter :: tcolor_lightgreen = achar(27)//"[92m"
	character(len=*), parameter :: tcolor_yellow     = achar(27)//"[93m"
	character(len=*), parameter :: tcolor_blue       = achar(27)//"[94m"
	character(len=*), parameter :: tcolor_purple     = achar(27)//"[95m"
	character(len=*), parameter :: tcolor_cyan       = achar(27)//"[96m"
	
    !! Constant tags:
    character(len=*), parameter :: tag_info  = "[INFO] "
    character(len=*), parameter :: tag_warn  = "["//tcolor_bold//tcolor_yellow//"WARN"//tcolor_nc//"] "
    character(len=*), parameter :: tag_error = "["//tcolor_bold//tcolor_red//"ERROR"//tcolor_nc//"] "
    
	!! Generic error messages:
	character(len=*), parameter :: errmsg_invalid_arg   = "Ebsolve: Invalid arguments, aborting..."
	character(len=*), parameter :: errmsg_unbound_index = "Ebsolve: Index out of bounds, aborting..."
	
end module
