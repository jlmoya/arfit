function t=tquant(n, p)
//  Quantiles of Student's t distribution
//  Calling Sequence
//  t=tquant(n, p)
//  Description
//  tquant(n, p) is the p-quantile of a t distributed random variable
//  with n degrees of freedom; that is, tquant(n, p) is the value below
//  which 100p percent of the t distribution with n degrees of freedom
//  lies.
//  Bibliography  
//  L. Devroye, 1986: "Non-Uniform Random Variate Generation", Springer
//  M. Abramowitz and I. A. Stegun, 1964: "Handbook of Mathematical  Functions" 
//  
//    Authors
//   H. Nahrstaedt - Jan 2011
//  Tapio Schneider      tapio@gps.caltech.edu


  if (n ~= round(n) | n < 1)
    error('Usage: TQUANT(n,p) - Degrees of freedom n must be positive integer.')
  end

  if (p<0 | p>1)
    error('Usage: TQUANT(n,p) - Probability p must be in [0,1].')  
  elseif p == 1
    t   = %inf;
    return
  elseif p == 0
    t   = -%inf;
    return
  end




   [mn,nn] = size(n);
   [mp,np] = size(p);
   if ( prod([mn,nn])==1 & prod([mp,np])<>1 ) then
	    n=n*ones(p);
    elseif ( prod([mn,nn])<>1 & prod([mp,np])==1 ) then
	    p=p*ones(n);
   end;
     
   

  // Make a special case for entries with a=0 
  I0 = p==0;
  I1 = p==1;

  ind =  ~(I0 | I1);
  if length(p(ind))>0 then 	
   t=cdft("T",n(ind),p(ind),1-p(ind));
  else
   t=zeros(p);
   t(p==0)=-%inf;
   t(p==1)=%inf;
  end;
  

// 
//   
//   if n == 1
//     // Cauchy distribution (cf. Devroye [1986, pp. 29 and 450])
//     t   = tan(%pi*(p-.5));
//   elseif p >= 0.5 
//     // positive t-values (cf. M. Abramowitz and I. A. Stegun [1964,
//     // Chapter 26])
//     b0  = [0, 1];
//     //f   = inline('1 - betainc(b, n/2, .5)/2 - p', 'b', 'n', 'p'); 
//     //f   = inline('1 - betainc(b, n/2, .5)/2 - p', 'b', 'n', 'p'); 
// 
// // function F=f(b,n,p)
// //  F=1 - betainc(b, n/2, .5)/2 - p;
// // endfunction
//   execstr("function F=f(b),F=1 - betainc(b, "+string(n)+"/2, .5)/2 - "+string(p)+";,endfunction");
// 
// 
//     opt = optimset('Display', 'off'); 
//     b   = fzero('f', b0, opt); 
//     // old calling sequence (on Windows/Mac with older Matlab versions):
//     //b   = fzero(f, b0, eps, 0, n, p); 
//     
//     t   = sqrt(n/b-n);
//   else
//     // negative t-values
//     b0  = [0, 1];
//     //f   = inline('betainc(b, n/2, .5)/2 - p', 'b', 'n', 'p'); 
// 
// //     function F=f(b,n,p)
// //     F= betainc(b, n/2, .5)/2 - p;
// //     endfunction
//     execstr("function F=f(b),F=betainc(b, "+string(n)+"/2, .5)/2 - "+string(p)+";,endfunction");
// 
//     opt = optimset('Display', 'off'); // does not work on Windows/Mac
//     b   = fzero('f', b0, opt); 
//     // old calling sequence (for Windows/Mac compatibility):
//     //b   = fzero(f, b0, eps, 0, n, p);  
//     t   = -sqrt(n/b-n);
//   end

endfunction

// 
// 
// function [F]=betainc(x,a,b)
// //  The beta cumulative distribution function
// //  Calling Sequence
// //         F = pbeta(x,a,b)
// //	Input	x	matrix (elements between 0 and 1)
// //		a,b	positive reals, parameters of the beta distribution
// //		(x,a,b  can be scalar or matrix with common size)
// //
// //	Output	F	for each element of x, F=Prob(X<x) where X is a random 
// //			variable with beta density :
// //			x--> x.^(a-1) .* (1-x).^(b-1) ./ beta(a,b)1_{0<x<1}
// //       Anders Holtsberg, 18-11-93
// //       Copyright (c) Anders Holtsberg
// //       last update: dec 2001 (jpc)
// //       completely changed in order to directly use cdfbet 
// if length(a)==1 then a = a*ones(x);end 
// if length(b)==1 then b = b*ones(x);end 
// 
// Ii = find(x>0 & x<1);
// F = 0*ones(x);
// Iu = find(x>=1);
// if Iu<>[] then F(Iu)=1 ;end 
// if Ii<>[] then
//   F1=cdfbet("PQ",x(Ii),1-x(Ii),a(Ii),b(Ii))
//   F(Ii)=F1;
// end
// endfunction
