function [R, rscale]=arqr(v, p, mcor)
// QR factorization for least squares estimation of AR model.
// Calling Sequence
// [R, rscale]=arqr(v, p, mcor)
// Description
//  [R, SCALE]=ARQR(v,p,mcor) computes the QR factorization needed in
//  the least squares estimation of parameters of an AR(p) model. If
//  the input flag mcor equals one, a vector of intercept terms is
//  being fitted. If mcor equals zero, the process v is assumed to have
//  mean zero. The output argument R is the upper triangular matrix
//  appearing in the QR factorization of the AR model, and SCALE is a
//  vector of scaling factors used to regularize the QR factorization.
//
//  ARQR is called by ARFIT. 
//
//  See also 
//  arfit
// Authors
//  Modified 29-Dec-99
//           24-Oct-10 Tim Mullen (added support for multiple realizatons)
//
//  Author: Tapio Schneider     tapio@gps.caltech.edu

  // n:   number of time steps (per realization)
  // m:   number of variables (dimension of state vectors) 
  // ntr: number of realizations (trials)
  if length(size(v))>2 then
     [n,m,ntr]   = size(v);     
  else
     ntr=1; [n,m]=size(v);
  end
  
   

  ne    = ntr*(n-p);                 // number of block equations of size m
  np    = m*p+mcor;             // number of parameter vectors of size m

  // If the intercept vector w is to be fitted, least squares (LS)
  // estimation proceeds by solving the normal equations for the linear
  // regression model
  //
  //                  v(k,:)' = Aaug*u(k,:)' + noise(C)
  //
  // with Aaug=[w A] and `predictors' 
  //
  //              u(k,:) = [1 v(k-1,:) ...  v(k-p,:)]. 
  //
  // If the process mean is taken to be zero, the augmented coefficient
  // matrix is Aaug=A, and the regression model
  //
  //                u(k,:) = [v(k-1,:) ...  v(k-p,:)]
  //
  // is fitted. 
  // The number np is the dimension of the `predictors' u(k). 

  // Assemble the data matrix K (of which a QR factorization will be computed)
  K = zeros(ne,np+m);                 // initialize K
  if (mcor == 1)
    // first column of K consists of ones for estimation of intercept vector w
    K(:,1) = ones(ne,1);
  end
  
  
//   for j=1:p
//     K(:, mcor+m*(j-1)+1:mcor+m*j) = [v(p-j+1:n-j, :)];
//   end
//   // Add `observations' v (left hand side of regression model) to K
//   K(:,np+1:np+m) = [v(p+1:n, :)];
  

  // Assemble `predictors' u in K 
  for itr=1:ntr
    for j=1:p
      K((n-p)*(itr-1) + 1 : (n-p)*itr, mcor+m*(j-1)+1 : mcor+m*j) = squeeze(v(p-j+1:n-j, :, itr));
    end
    // Add `observations' v (left hand side of regression model) to K
    K((n-p)*(itr-1) + 1 : (n-p)*itr, np+1 : np+m) = squeeze(v(p+1:n, :, itr));
  end


  // Compute regularized QR factorization of K: The regularization
  // parameter delta is chosen according to Higham's (1996) Theorem
  // 10.7 on the stability of a Cholesky factorization. Replace the
  // regularization parameter delta below by a parameter that depends
  // on the observational error if the observational error dominates
  // the rounding error (cf. Neumaier, A. and T. Schneider, 2001:
  // "Estimation of parameters and eigenmodes of multivariate
  // autoregressive models", ACM Trans. Math. Softw., 27, 27--57.).
  q     = np + m;             // number of columns of K
  delta = (q^2 + q + 1)*%eps;  // Higham's choice for a Cholesky factorization
  rscale = sqrt(delta)*sqrt(mtlb_sum(K.^2));   
  [tmp,tmpA]=qr([K; diag(rscale)],"e");
  R     = mtlb_triu(tmpA);
endfunction

function B = squeeze(A)
   //  a squeeze function likes the Matlab 's one ?
   if argn(2) ~= 1 then
      error("squeeze: bad number of input arg")
   end
   B_dims = size(A)
   if length(B_dims) <= 2 then
      B = A
   else  // length(B_dims)>=3  => A is an hypermat
      B_dims(B_dims == 1) = []
      if length(B_dims) >= 2 then
         B = matrix(A.entries, B_dims)
      else
         B = A.entries
      end
   end
endfunction

