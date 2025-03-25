!!**************************************************************************************************
!! Created on 2024-07-19 at 11:58:15 CEST by David Gaspard <david.gaspard@espci.fr>
!! This program is distributed under the MIT License.
!! Fortran module to integrate the Eilenberger equation, a matrix transport equation, in a single direction along a single ray.
!!**************************************************************************************************
module integrator
    use constants
	implicit none
	
    !! Define some global constants:
    logical, parameter :: run_check_ebintegrate = .false.   !! Execute some test before perfomring computations (recommended for debugging).
    
	contains
    
	!!**********************************************************************************************
	!! Subroutine to integrate the Eilenberger equation (a matrix transport equation) in a single direction along a single ray.
    !! More precisely, this subroutine solves in a given direction "mu" the linear Eilenberger equation of the form:
    !! 
    !!   ± mu * dg(±,mu,x)/dx = [ P(x), g(±,mu,x) ] + delta(x - x_a) * [ C_a(±,mu), g(±,mu,x) ]
    !!                                              + delta(x - x_b) * [ C_b(±,mu), g(±,mu,x) ] ,
    !! 
    !! with the boundary conditions (at infinity): g(in) = Lambda_3 + g_21 * Lambda_-, and g(out) = Lambda_3 + g_12 * Lambda_+ ,
    !! g_12 and g_21 being unknown, and Lambda_±, Lambda_3 the usual Pauli matrices. It should be noted that
    !! these boundary conditions only fix 3 out of 4 matrix element of g(mu,x) and thus require
    !! the specialized ODE integrator which follows.
    !! 
    !! The symbols in the above Eilenberger equation are:
    !!   ±         = Sign of the propagation direction.
    !!               Plus is the increasing-x direction, and minus is the decreasing-x direction.
    !!   mu        = Direction cosine of the current ray with respect to the direction of the positional mesh ("xmesh").
    !!   x         = Position along the ray.
    !!   g(±,mu,x) = Matrix radiance. Unknown function to be integrated by the present subroutine.
    !!   P(x)      = Matrix defining the distributed interactions. Typically: P(x) = - Qn(x)/(2*lscat(x)) - Lambda_3/(2*labso(x)).
    !!               In principle, Qn(x) depends on g(±,mu,x), but this dependence is neglected in this subroutine.
    !!   delta(x)  = Usual Dirac delta.
    !!   x_a, x_b  = Positions of the contact interactions.
    !!   C_a, C_b  = Matrices defining the contact interactions.
    !!               Note that they may depend on the direction "mu" (e.g. in case of limited numerical aperture).
    !! 
    !! Arguments of this subroutine:
    !!   mu           = (in)  Direction cosine of the current ray with respect to the direction of the positional mesh ("xmesh").
    !!   sgn          = (in)  Sign of the propagation direction (plus or minus).
    !!                        For sgn=+1, compute g_+(mu,x), for sgn=-1, compute g_-(mu,x).
    !!   xmesh        = (in)  Array of the positional lattice ("x" lattice) associated to g(±,mu,x). Dimensions: xmesh(nx).
    !!   pmesh        = (in)  Array containing the distributed interaction P(x) at the midpoints.
    !!                        Dimensions: pmesh(3, nx-1). Triplet order: (P_11, P_12, P_21).
    !!   xa, xb       = (in)  Positions of the contact interactions relative to the "xmesh".
    !!   ca, cb       = (in)  Matrices of the contact interactions. Dimensions: ca(3), cb(3). Triplet order: (C_11, C_12, C_21).
    !!   amesh, bmesh = (out) Array of the computed vector radiances a(±,mu,x) and b(±,mu,x).
    !!                        Dimensions: amesh(2, nx), bmesh(2, nx). Triplet order: (a_1, a_2).
    !!   gmesh        = (out) Array of the computed matrix radiance g(±,mu,x).
    !!                        Dimensions: gmesh(3, nx). Triplet order: (g_11, g_12, g_21).
	!!**********************************************************************************************
	subroutine ebintegrate(sgn, mu, xmesh, pmesh, xa, ca, xb, cb, amesh, bmesh, gmesh)
		integer, intent(in) :: sgn
        complex(wp), intent(in) :: mu
        real(wp), intent(in) :: xmesh(:), xa, xb
        complex(wp), intent(in) :: pmesh(:, :), ca(:), cb(:)
        complex(wp), intent(out) :: amesh(:, :), bmesh(:, :), gmesh(:, :)
        complex(wp) :: denom
        integer :: ix, nx
        
        !! 1. Check for possible input errors:
        if (run_check_ebintegrate) then
            call check_ebintegrate(sgn, mu, xmesh, pmesh, xa, ca, xb, cb, amesh, bmesh, gmesh)
        end if
        
        !! 2. Solve the vector form of the Eilenberger equation (integrates in two directions):
        call ebintegrate_vec('a', sgn, mu, xmesh, pmesh, xa, ca, xb, cb, amesh)  !! Integrate in the direction given by 'sgn'.
        call ebintegrate_vec('b', sgn, mu, xmesh, pmesh, xa, ca, xb, cb, bmesh)  !! Integrate in the direction opposite to 'sgn'.
        
        !! 3. Compute g(±,mu,x) from the vectors:
        nx = size(xmesh)
        do ix = 1, nx
            denom = amesh(2, ix)*bmesh(1, ix) - amesh(1, ix)*bmesh(2, ix)
            gmesh(1, ix) = (amesh(2, ix)*bmesh(1, ix) + amesh(1, ix)*bmesh(2, ix))/denom  !! (a2 b1+a1 b2)/(a2 b1-a1 b2)
            gmesh(2, ix) = -2*amesh(1, ix)*bmesh(1, ix)/denom  !! (-2 a1 b1)/(a2 b1-a1 b2)
            gmesh(3, ix) =  2*amesh(2, ix)*bmesh(2, ix)/denom  !! (2 a2 b2)/(a2 b1-a1 b2)
        end do
        
	end subroutine
	
    !!**********************************************************************************************
    !! Checks the validity of the arguments of ebintegrate(). Stops the program on error.
    !! This subroutine is called by ebintegrate() when the global flag "run_test_ebintegrate" is enabled.
    !! Nevertheless its strongly recommended in case of debugging.
    !!**********************************************************************************************
    subroutine check_ebintegrate(sgn, mu, xmesh, pmesh, xa, ca, xb, cb, amesh, bmesh, gmesh)
        integer, intent(in) :: sgn
        real(wp), intent(in) :: xmesh(:), xa, xb
        complex(wp), intent(in) :: mu, pmesh(:, :), ca(:), cb(:), amesh(:, :), bmesh(:, :), gmesh(:, :)
        integer :: nx
        logical :: error
        
        error = .false.  !! Becomes true on error.
        nx = size(xmesh)
        if (sgn /= 1 .and. sgn /= -1) then
            write (stderr, '(a,i0,a)') tag_error // "Invalid direction sign, received ", &
                sgn, ", must be +1 or -1..."
            error = .true.
        else if (mu == 0.) then
            write (stderr, '(a)') tag_error // "Invalid direction cosine, cannot be zero..."
            error = .true.
        else if (.not. isvalidxmesh(xmesh)) then
            write (stderr, '(a)') tag_error // "Invalid 'xmesh', must be sorted in increasing order with nonzero intervals..."
            error = .true.
        else if (.not. (islocatedinxmesh(xa, xmesh) .and. islocatedinxmesh(xb, xmesh))) then
            write (stderr, '(a)') tag_error // "Contact points 'xa' and 'xb' must be within the interval covered by 'xmesh'..."
            error = .true.
        else if (any(xa == xmesh) .or. any(xb == xmesh)) then
            write (stderr, '(a)') tag_error // "Contact points 'xa' and 'xb' cannot coincide with the points &
                &in 'xmesh' due to the discontinuity..."
            error = .true.
        else if (size(pmesh, 1) /= 3 .or. size(pmesh, 2) /= nx-1) then
            write (stderr, '(3(a,i0),a)') tag_error // "Invalid dimensions of 'pmesh', received (", &
                size(pmesh, 1), ", ", size(pmesh, 2), "), expected (3, ", nx-1, ")..."
            write (stderr, '(a)') tag_info // "Note that 'pmesh' contains the values of P(x) in the intervals &
                &between of the points of 'xmesh', not at the points of 'xmesh' themselves..."
            error = .true.
        else if (size(ca) /= 3 .or. size(cb) /= 3) then
            write (stderr, '(2(a,i0),a)') tag_error // "Invalid dimensions of 'ca' and 'cb', received ", &
                size(ca), " and ", size(cb), ", expected 3..."
            error = .true.
        else if (size(gmesh, 1) /= 3 .or. size(gmesh, 2) /= nx) then
            write (stderr, '(3(a,i0),a)') tag_error // "Invalid dimensions of 'gmesh', received (", &
                size(gmesh, 1), ", ", size(gmesh, 2), "), expected (3, ", nx, ")..."
            error = .true.
        else if (size(amesh, 1) /= 2 .or. size(amesh, 2) /= nx) then
            write (stderr, '(3(a,i0),a)') tag_error // "Invalid dimensions of 'amesh', received (", &
                size(amesh, 1), ", ", size(amesh, 2), "), expected (2, ", nx, ")..."
            error = .true.
        else if (size(bmesh, 1) /= 2 .or. size(bmesh, 2) /= nx) then
            write (stderr, '(3(a,i0),a)') tag_error // "Invalid dimensions of 'bmesh', received (", &
                size(bmesh, 1), ", ", size(bmesh, 2), "), expected (2, ", nx, ")..."
            error = .true.
        end if
        
        if (error) stop errmsg_invalid_arg
        
    end subroutine
    
    !!**********************************************************************************************
    !! Solves the vector form of the Eilenberger equation:
    !! 
    !!   ± mu * dv(±,mu,x)/dx = P(x) * v(±,mu,x) + delta(x - x_a) * C_a(±,mu) * v(±,mu,x)
    !!                                           + delta(x - x_b) * C_b(±,mu) * v(±,mu,x) ,
    !! 
    !! where the symbols are the same as in the previous subroutine ebintegrate().
    !! This subroutine does not perform error checks.
    !! 
    !! Arguments:
    !!   aob    = (in)  Character 'a' or 'b' deciding the type of vector and the boundary conditions.
    !!   sgn    = (in)  Sign of the propagation direction (plus or minus).
    !!                  For sgn=+1, compute g_+(mu,x), for sgn=-1, compute g_-(mu,x).
    !!   xmesh  = (in)  Array of the positional lattice ("x" lattice) associated to g(±,mu,x). Dimensions: xmesh(nx).
    !!   pmesh  = (in)  Array containing the distributed interaction P(x) at the midpoints.
    !!                  Dimensions: pmesh(3, nx-1). Triplet order: (P_11, P_12, P_21).
    !!   xa, xb = (in)  Positions of the contact interactions relative to the "xmesh".
    !!   ca, cb = (in)  Matrices of the contact interactions. Dimensions: ca(3), cb(3). Triplet order: (C_11, C_12, C_21).
    !!   vmesh  = (out) Array of the computed vector radiances. Dimensions: vmesh(2, nx). Doublet order: (v_1, v_2).
    !!**********************************************************************************************
    subroutine ebintegrate_vec(aob, sgn, mu, xmesh, pmesh, xa, ca, xb, cb, vmesh)
        character, intent(in) :: aob
        integer, intent(in) :: sgn
        complex(wp), intent(in) :: mu
        real(wp), intent(in) :: xmesh(:), xa, xb
        complex(wp), intent(in) :: pmesh(:, :), ca(:), cb(:)
        complex(wp), intent(out) :: vmesh(:, :)
        complex(wp) :: u, v, sgdx
        real(wp) :: x0, x1, xn
        integer :: nx, i0, i1, di, imin, coll
        
        !! A. Initialize the integration over x (boundary conditions):
		if (aob == 'a') then
			u = 0.  !! Boundary condition.
			v = 1.
			di = sgn  !! Integrate in the same direction as 'sgn'.
		else if (aob == 'b') then
			u = 1.  !! Boundary condition.
			v = 0.
			di = -sgn !! Integrate in the direction opposite to 'sgn'.
		else
			write (stderr, '(3a)') tag_error // "Unknown option: aob='", aob, "', expected 'a' or 'b'."
			stop errmsg_invalid_arg
		end if
        
        nx = size(xmesh)
		if (di == +1) then
			i0 = 1   !! Index of the initial position.
		else if (di == -1) then
			i0 = nx   !! Index of the initial position.
		else 
			write (stderr, '(a,i0,a)') tag_error // "Invalid parameter: sgn=", sgn, ", expected +1 or -1."
			stop errmsg_invalid_arg
		end if
        vmesh(1, i0) = u  !! Initialize the array.
        vmesh(2, i0) = v
		x0 = xmesh(i0)  !! Current position.
		i1 = i0 + di  !! Index of the next point planned in xmesh.
		
		!! B. Integration over x:
		do while (1 <= i1 .and. i1 <= nx)
			
			!! 1. Determine the actual next point (either the next xmesh point, x1, or xa, or xb):
			x1 = xmesh(i1) !! Next point planned in xmesh.
			if (between_excl(xa, x0, x1) .and. .not. between_excl(xb, x0, xa)) then
				xn = xa   !! If xa in [x0, x1] and xb not in [x0, xa], then next collision point is xa.
				coll = 1  !! Collision status (0=No collision, 1=Collision at xa, 2=Collision at xb).
			else if (between_excl(xb, x0, x1) .and. .not. between_excl(xa, x0, xb)) then 
				xn = xb   !! If xb in [x0, x1] and xa not in [x0, xb], then next collision point is xb.
				coll = 2  !! Collision status (0=No collision, 1=Collision at xa, 2=Collision at xb).
			else
				xn = x1   !! Next point is x1.
				coll = 0  !! Collision status (0=No collision, 1=Collision at xa, 2=Collision at xb).
			end if
			
			!! 2. Integrate over the interval [x0, xn]:
			imin = min(i0, i1)    !! The P(x) matrix in the current interval is pmesh(:, imin).
            sgdx = sgn*(xn - x0)/mu    !! Factor of the P(x) matrix before exponentiation.
            call apply_mexp(sgdx*pmesh(1, imin), sgdx*pmesh(2, imin), sgdx*pmesh(3, imin), u, v)
            
			!! 3. Deals with the contact interactions (xa or xb) at the end of the interval.
			if (coll == 1) then  !! Collision with xa.
                sgdx = sgn*di/mu
                call apply_mexp(sgdx*ca(1), sgdx*ca(2), sgdx*ca(3), u, v)
			else if (coll == 2) then  !! Collision with xb.
                sgdx = sgn*di/mu
				call apply_mexp(sgdx*cb(1), sgdx*cb(2), sgdx*cb(3), u, v)
			else  !! No collision detected.
				vmesh(1, i1) = u   !! Store the components (u,v) in the vmesh.
				vmesh(2, i1) = v
				i0 = i1       !! Update the current position index.
				i1 = i1 + di  !! Go to the next mesh point.
			end if
			x0 = xn  !! Update the current position.
            
        end do
        
    end subroutine
    
    !!**********************************************************************************************
    !! Check if the given value 'x' is contained in the interval [a, b], regardless of the order
    !! of 'a' and 'b'. In other words, 'b' can be smaller than 'a'.
    !! Returns True if 'x' is in the given interval (excluding the bounds), False otherwise.
    !!**********************************************************************************************
    function between_excl(x, a, b) result(bool)
		logical :: bool
		real(wp), intent(in) :: x, a, b
		real(wp) :: xmin, xmax
		if (a < b) then
			xmin = a
			xmax = b
		else
			xmin = b
			xmax = a
		end if
		bool = (xmin < x .and. x < xmax)
    end function
	
    !!**********************************************************************************************
    !! Computes efficiently the transformation exp(A) * vec, where:
    !!   A   = Traceless 2x2 matrix. A = (a11, a12, a21, -a11).
    !!   vec = Complex 2-component vector (spinor). vec = (u, v).
    !! 
    !! This subroutine uses the exact analytical expression of the matrix exponential of a traceless 2x2 matrix.
    !! This routine is very stable and takes into account the special case of a singular matrix [ typically det(A) = 0 ].
    !! 
    !! Arguments:
    !!   a11, a12, a21 = Matrix elements of the matrix A to be exponentiated.
    !!                   Note that, since A is traceless, we have necessarily: a22 = -a11.
    !!   u, v          = On input, the initial components of the vector "vec".
    !!                   On output, the components after the transformation exp(A) * vec.
    !!**********************************************************************************************
    subroutine apply_mexp(a11, a12, a21, u, v)
        complex(wp), intent(in) :: a11, a12, a21
        complex(wp), intent(inout) :: u, v
        complex(wp) :: u0, v0, sdet, ch, sh
        
        sdet = sqrt(a11 * a11 + a12 * a21)
        if (sdet /= 0.) then
            ch = cosh(sdet)
            sh = sinh(sdet)/sdet
        else
            ch = 1._wp
            sh = 1._wp
        end if
        
        u0 = u   !! Buffer values.
        v0 = v
        
        u = ch * u0 + sh * (a11 * u0 + a12 * v0)    !! Apply the transformation.
        v = ch * v0 + sh * (a21 * u0 - a11 * v0)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Check if the given "xmesh" (positional lattice) is a valid for the current program.
    !! Returns True if the "xmesh" is valid, False otherwise.
    !!**********************************************************************************************
    function isvalidxmesh(xmesh) result(bool)
        real(wp), intent(in) :: xmesh(:)
        integer :: ix
        logical :: bool
        bool = .true.
        do ix = 2, size(xmesh)
            if (xmesh(ix-1) >= xmesh(ix)) then
                bool = .false.
                exit
            end if
        end do
    end function
    
    !!**********************************************************************************************
    !! Check if the given point "x" is contained in the interval covered by "xmesh" (positional lattice).
    !! Note that this function assumes the "xmesh" is sorted in increasing order (see other tests).
    !! Returns True if "x" is contained in "xmesh" is valid, False otherwise.
    !!**********************************************************************************************
    function islocatedinxmesh(x, xmesh) result(bool)
        real(wp), intent(in) :: x, xmesh(:)
        logical :: bool
        bool = (xmesh(1) < x .and. x < xmesh(size(xmesh)))
    end function
    
end module
