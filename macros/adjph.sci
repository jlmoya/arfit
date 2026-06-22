function ox=adjph(x)
//  Normalization of columns of a complex matrix.
//  Calling Sequence
//   OX=adjph(X)
//   Parameters
//   X: complex matrix
//   OX: complex matrix
//  Description
//
//  Given a complex matrix X, OX=adjph(X) returns the complex matrix OX
//  that is obtained from X by multiplying column vectors of X with
//  phase factors exp(%i*phi) such that the real part and the imaginary
//  part of each column vector of OX are orthogonal and the norm of the
//  real part is greater than or equal to the norm of the imaginary
//  part.
//
//  adjph is called by armode.
//
//  See also 
//  armode
// Authors
//  Modified 16-Dec-99
//  Tapio Schneider   tapio@gps.caltech.edu

  for j = 1:size(x,2)				
    a       = real(x(:,j));                     // real part of jth column of x
    b       = imag(x(:,j));                     // imag part of jth column of x
    phi     = .5*atan( 2*mtlb_sum(a.*b)/(b'*b-a'*a) );
    bnorm   = norm(sin(phi).*a+cos(phi).*b);    // norm of new imaginary part
    anorm   = norm(cos(phi).*a-sin(phi).*b);    // norm of new real part
    if bnorm > anorm 
      if phi < 0
	phi = phi-%pi/2;
      else
	phi = phi+%pi/2;
      end
    end
    ox(:,j) = x(:,j).*exp(%i*phi);
  end






endfunction