function [siglev,res, lmp, dof_lmp]=arres(w,A,v,k)
// Test of residuals of fitted AR model.
// Calling Sequence
// [siglev,res, lmp, dof_lmp]=arres(w,A,v,k)
// Description
//  [siglev,res]=ARRES(w,A,v) computes the time series of residuals
//
//        res(k,:)' = v(k+p,:)'- w - A1*v(k+p-1,:)' - ... - Ap*v(k,:)'
//
//  of an AR(p) model with A=[A1 ... Ap]. If v has three dimensions,
//  the 3rd dimension is treated as indicating multiple realizations
//  (trials) of the time series, that is, v(:,:,itr) is taken as the
//  itr-th trial. 
//
//
//  Also returned is the significance level siglev of the modified
//  Li-McLeod portmanteau (LMP) statistic.
//
//  Correlation matrices for the LMP statistic are computed up to lag
//  k=20, which can be changed to lag k by using
//  [siglev,res]=ARRES(w,A,v,k).
// 
//  Bibliography
//    Li, W. K., and A. I. McLeod, 1981: Distribution of the
//        Residual Autocorrelations in Multivariate ARMA Time Series
//        Models, J. Roy. Stat. Soc. B, 43, 231--239.
// Authors
//  Modified 11-Jun-07 German Gomez-Herrero, german.gomezherrero@ieee.org
//  - Input data can be a ddataset object (to save memory resources)
//
//  Modified 17-Dec-99
//       24-Oct-10 Tim Mullen (added support for multiple realizations)  
//
//  Author: Tapio Schneider   tapio@gps.caltech.edu
//


  [nargout,nargin]=argn();

//if isa(v,'ddataset'),
//    v = load(v);
//end
[n,m] = size(v);
if n<m,
    v = v';
    [n,m] = size(v);
end

n     = size(v,1);                    // number of observations
m     = size(v,2);                    // dimension of state vectors
  if length(size(v))>2 then
     ntr   = size(v,3);                    // number of realizations (trials)
  else
     ntr=1;
  end


p     = size(A,2)/m;                  // order of model

nres  = n-p;                          // number of residuals

// Default value for k
if (nargin < 4)
    //k   = 60;
    k   = min(20, nres-1);

end
if (k <= p) & nargout >1                          // check if k is in valid range
    error('Maximum lag of residual correlation matrices too small.');
end
if (k >= nres)
    error('Maximum lag of residual correlation matrices too large.');
end

w     = w(:)';                        // force w to be row vector

// Get time series of residuals
 res = zeros(nres,m,ntr);

  l = 1:nres;                           // vectorized loop l=1,...,nres  
  res(l,:,:) = v(l+p,:,:) - mtlb_repmat(w, [nres, 1, ntr]);
  for itr=1:ntr
    for j=1:p
      res(l,:,itr) = res(l,:,itr) - v(l-j+p,:,itr)*A(:, (j-1)*m+1:j*m)';
    end
  end
  // end of loop over l
  
  // For computation of correlation matrices, center residuals by
  // subtraction of the mean  
  resc  = res - mtlb_repmat(mtlb_mean(res), [nres, 1, 1]);
  
  // Compute correlation matrix of the residuals
  // compute lag zero correlation matrix
  c0    = zeros(m, m);
  for itr=1:ntr
    c0  = c0 + resc(:,:,itr)'*resc(:,:,itr);
  end
  d     = diag(c0);
  dd    = sqrt(d*d');
  c0    = c0./dd;

  // compute lag l correlation matrix
  cl    = zeros(m, m, k);
  for l=1:k
    for itr=1:ntr
      cl(:,:,l) = cl(:,:,l)    + resc(1:nres-l, :, itr)'*resc(l+1:nres, :, itr);
    end
    cl(:,:,l) = cl(:,:,l)./dd;  
  end
  
  // Get "covariance matrix" in LMP statistic
  c0_inv= inv(c0);                      // inverse of lag 0 correlation matrix
  rr    = kron(c0_inv, c0_inv);         // "covariance matrix" in LMP statistic

  // Compute modified Li-McLeod portmanteau statistic
  lmp   = 0;                            // LMP statistic initialization
  x     = zeros(m*m,1);                 // correlation matrix arranged as vector
  for l=1:k
    x   = matrix(cl(:,:,l), m^2, 1);    // arrange cl as vector by stacking columns
    lmp = lmp + x'*rr*x;                // sum up LMP statistic
  end
  ntot  = n*ntr;                        // total number of observations
  lmp   = ntot*lmp + m^2*k*(k+1)/2/ntot;// add remaining term and scale
  dof_lmp = m^2*(k-p);                  // degrees of freedom for LMP statistic
      
  // Significance level with which hypothesis of uncorrelatedness is rejected
  //siglev = 1 - gammainc(lmp/2, dof_lmp/2);
   

  g_x=lmp/2;
  g_a=dof_lmp/2;
   [ma,na] = size(g_a);
   [mx,nx] = size(g_x);
   if ( prod([mx,nx])==1 & prod([ma,na])<>1 ) then
	    g_x=g_x*ones(g_a);
    elseif ( prod([mx,nx])<>1 & prod([ma,na])==1 ) then
	    g_a=g_a*ones(g_x);
   end;
   Rate = ones(g_a);
  // Make a special case for entries with a=0 
  kzeroa = find(g_a==0);
  g_a(kzeroa) = 1;
  [g_p,g_q] = cdfgam("PQ",g_x,g_a,Rate);
   g_p(kzeroa)=1;
   

   siglev = 1 - g_p;
  
endfunction








// l = 1:nres; 	// vectorized loop l=1,...,nres
// //res=zeros(nres,2);			
// res(l,:) = v(l+p,:) - ones(nres,1)*w;
// for j=1:p
//     res(l,:) = res(l,:) - v(l-j+p,:)*A(:, (j-1)*m+1:j*m)';
// end
// // end of loop over l
// 
// // Center residuals by subtraction of the mean
// res   = res - ones(nres,1)*mtlb_mean(res);
// 
// if nargout >1,
//     // Compute lag zero correlation matrix of the residuals
//     c0    = res'*res;
//     d     = diag(c0);
//     dd    = sqrt(d*d');
//     c0    = c0./dd;
// 
//     // Get "covariance matrix" in LMP statistic
//     c0_inv= inv(c0);                      // inverse of lag 0 correlation matrix
//     rr    = (c0_inv .*. c0_inv);         // "covariance matrix" in LMP statistic
// 
//     // Initialize LMP statistic and correlation matrices
//     lmp   = 0;                            // LMP statistic
//     cl    = zeros(m,m);                   // correlation matrix
//     x     = zeros(m*m,1);                 // correlation matrix arranged as vector
// 
//     // Compute modified Li-McLeod portmanteau statistic
//     for l=1:k
//         cl  = (res(1:nres-l, :)'*res(l+1:nres,:))./dd;  // lag l correlation matrix
//         x   = matrix(cl,m*m, 1);           // arrange cl as vector by stacking columns
//         lmp = lmp + x'*rr*x;                // sum up LMP statistic
//     end
//     lmp   = n*lmp + m^2*k*(k+1)/2/n;      // add remaining term and scale
//     dof_lmp = m^2*(k-p);                  // degrees of freedom for LMP statistic
// 
//     // Significance level with which hypothesis of uncorrelatedness is rejected
//     siglev = 1 - gammainc(lmp/2, dof_lmp/2);
// 
// end
// 
// 
// 
// endfunction

