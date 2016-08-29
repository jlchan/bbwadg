function test_burgers

x = linspace(-1,1,100);
v = burgers_solution(1e-3,100,x,1,.32);
plot(x,v)
end

function vu = burgers_solution ( nu, vxn, vx, vtn, vt )

%*****************************************************************************80
%
%% BURGERS_SOLUTION evaluates a solution to the Burgers equation.
%
%  Discussion:
%
%    The form of the Burgers equation considered here is
%
%      du       du        d^2 u
%      -- + u * -- = nu * -----
%      dt       dx        dx^2
%
%    for -1.0 < x < +1.0, and 0 < t.
%
%    Initial conditions are u(x,0) = - sin(pi*x).  Boundary conditions
%    are u(-1,t) = u(+1,t) = 0.  The viscosity parameter nu is taken
%    to be 0.01 / pi, although this is not essential.
%
%    The authors note an integral representation for the solution u(x,t),
%    and present a better version of the formula that is amenable to
%    approximation using Hermite quadrature.
%
%    This program library does little more than evaluate the exact solution
%    at a user-specified set of points, using the quadrature rule.
%    Internally, the order of this quadrature rule is set to 8, but the
%    user can easily modify this value if greater accuracy is desired.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    17 November 2011
%
%  Author:
%
%    John Burkardt.
%
%  Reference:
%
%    Claude Basdevant, Michel Deville, Pierre Haldenwang, J Lacroix,
%    J Ouazzani, Roger Peyret, Paolo Orlandi, Anthony Patera,
%    Spectral and finite difference solutions of the Burgers equation,
%    Computers and Fluids,
%    Volume 14, Number 1, 1986, pages 23-41.
%
%  Parameters:
%
%    Input, real NU, the viscosity.
%
%    Input, integer VXN, the number of spatial grid points.
%
%    Input, real VX(VXN), the spatial grid points.
%
%    Input, integer VTN, the number of time grid points.
%
%    Input, real VT(VTN), the time grid points.
%
%    Output, real VU(VXN,VTN), the solution of the Burgers
%    equation at each space and time grid point.
%
qn = 8;
%
%  Compute the rule.
%
[ qx, qw ] = hermite_ek_compute ( qn );
%
%  Evaluate U(X,T) for later times.
%
vu = zeros ( vxn, vtn );

for vti = 1 : vtn
    
    if ( vt(vti) == 0.0 )
        
        vu(1:vxn,vti) = - sin ( pi * vx(1:vxn) );
        
    else
        
        for vxi = 1 : vxn
            
            top = 0.0;
            bot = 0.0;
            
            for qi = 1 : qn
                
                c = 2.0 * sqrt ( nu * vt(vti) );
                
                top = top - qw(qi) * c * sin ( pi * ( vx(vxi) - c * qx(qi) ) ) ...
                    * exp ( - cos ( pi * ( vx(vxi) - c * qx(qi)  ) ) ...
                    / ( 2.0 * pi * nu ) );
                
                bot = bot + qw(qi) * c ...
                    * exp ( - cos ( pi * ( vx(vxi) - c * qx(qi)  ) ) ...
                    / ( 2.0 * pi * nu ) );
                
                vu(vxi,vti) = top / bot;
                
            end
            
        end
        
    end
    
end

return
end


function [ x, w ] = hermite_ek_compute ( n )

%*****************************************************************************80
%
%% HERMITE_EK_COMPUTE computes a Gauss-Hermite quadrature rule.
%
%  Discussion:
%
%    The code uses an algorithm by Elhay and Kautsky.
%
%    The abscissas are the zeros of the N-th order Hermite polynomial.
%
%    The integral:
%
%      integral ( -oo < x < +oo ) exp ( - x * x ) * f(x) dx
%
%    The quadrature rule:
%
%      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    19 April 2011
%
%  Author:
%
%    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Sylvan Elhay, Jaroslav Kautsky,
%    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
%    Interpolatory Quadrature,
%    ACM Transactions on Mathematical Software,
%    Volume 13, Number 4, December 1987, pages 399-415.
%
%  Parameters:
%
%    Input, integer N, the number of abscissas.
%
%    Output, real X(N), the abscissas.
%
%    Output, real W(N), the weights.
%

%
%  Define the zero-th moment.
%
zemu = gamma ( 0.5 );
%
%  Define the Jacobi matrix.
%
bj = zeros ( n, 1 );
for i = 1 : n
    bj(i) = i / 2.0;
end
bj(1:n) = sqrt ( bj(1:n) );

x = zeros ( n, 1 );

w = zeros ( n, 1 );
w(1) = sqrt ( zemu );
%
%  Diagonalize the Jacobi matrix.
%
[ x, w ] = imtqlx ( n, x, bj, w );

w(1:n) = w(1:n).^2;

return
end


function [ d, z ] = imtqlx ( n, d, e, z )

%*****************************************************************************80
%
%% IMTQLX diagonalizes a symmetric tridiagonal matrix.
%
%  Discussion:
%
%    This routine is a slightly modified version of the EISPACK routine to
%    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
%
%    The authors thank the authors of EISPACK for permission to use this
%    routine.
%
%    It has been modified to produce the product Q' * Z, where Z is an input
%    vector and Q is the orthogonal matrix diagonalizing the input matrix.
%    The changes consist (essentially) of applying the orthogonal transformations
%    directly to Z as they are generated.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    04 January 2010
%
%  Author:
%
%    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Sylvan Elhay, Jaroslav Kautsky,
%    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
%    Interpolatory Quadrature,
%    ACM Transactions on Mathematical Software,
%    Volume 13, Number 4, December 1987, pages 399-415.
%
%    Roger Martin, James Wilkinson,
%    The Implicit QL Algorithm,
%    Numerische Mathematik,
%    Volume 12, Number 5, December 1968, pages 377-383.
%
%  Parameters:
%
%    Input, integer N, the order of the matrix.
%
%    Input, real D(N), the diagonal entries of the matrix.
%
%    Input, real E(N), the subdiagonal entries of the
%    matrix, in entries E(1) through E(N-1). 
%
%    Input, real Z(N), a vector to be operated on.
%
%    Output, real D(N), the diagonal entries of the diagonalized matrix.
%
%    Output, real Z(N), the value of Q' * Z, where Q is the matrix that 
%    diagonalizes the input symmetric tridiagonal matrix.
%
  itn = 30;

  prec = eps;

  if ( n == 1 )
    return
  end

  e(n) = 0.0;

  for l = 1 : n

    j = 0;

    while ( 1 )

      for m = l : n

        if ( m == n )
          break
        end

        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) )
          break
        end

      end

      p = d(l);

      if ( m == l )
        break
      end

      if ( j == itn )
        fprintf ( 1, '\n' );
        fprintf ( 1, 'IMTQLX - Fatal error!\n' );
        fprintf ( 1, '  Iteration limit exceeded.\n' );
        error ( 'IMTQLX - Fatal error!' );
      end

      j = j + 1;
      g = ( d(l+1) - p ) / ( 2.0 * e(l) );
      r =  sqrt ( g * g + 1.0 );
      g = d(m) - p + e(l) / ( g + r8_sign ( g ) * abs ( r ) );
      s = 1.0;
      c = 1.0;
      p = 0.0;
      mml = m - l;

      for ii = 1 : mml

        i = m - ii;
        f = s * e(i);
        b = c * e(i);

        if ( abs ( f ) >= abs ( g ) )
          c = g / f;
          r =  sqrt ( c * c + 1.0 );
          e(i+1) = f * r;
          s = 1.0 / r;
          c = c * s;
        else
          s = f / g;
          r =  sqrt ( s * s + 1.0 );
          e(i+1) = g * r;
          c = 1.0 / r;
          s = s * c;
        end

        g = d(i+1) - p;
        r = ( d(i) - g ) * s + 2.0 * c * b;
        p = s * r;
        d(i+1) = g + p;
        g = c * r - b;
        f = z(i+1);
        z(i+1) = s * z(i) + c * f;
        z(i) = c * z(i) - s * f;

      end

      d(l) = d(l) - p;
      e(l) = g;
      e(m) = 0.0;

    end

  end

  for ii = 2 : n

     i = ii - 1;
     k = i;
     p = d(i);

     for j = ii : n
       if ( d(j) < p )
         k = j;
         p = d(j);
       end
     end

     if ( k ~= i )
       d(k) = d(i);
       d(i) = p;
       p = z(i);
       z(i) = z(k);
       z(k) = p;
     end

  end

  return
end


function value = r8_sign ( x )

%*****************************************************************************80
%
%% R8_SIGN returns the sign of an R8.
%
%  Discussion:
%
%    The value is +1 if the number is positive or zero, and it is -1 otherwise.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    21 March 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the number whose sign is desired.
%
%    Output, real VALUE, the sign of X.
%
  if ( 0 <= x )
    value = +1.0;
  else
    value = -1.0;
  end

  return
end