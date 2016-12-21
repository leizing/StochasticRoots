function prj_x = spgm_proj(x)
%SPGM_PROJ  Project A = reshape(x) to the set of stochastic matrices
%   PRJ_X = SPGM_PROJ(X) compute the nearest stochastic matrix to given matrix
%   A = reshape(X), where the distance is measured by matrix Frobenius norm.

  x = real(x);  % discard the imaginary part;
  n = length(x);
%   m = int32(sqrt(double(n)));  % not necessary in MATLAB2010a
  m = sqrt(n);
  A = reshape(x,m,m);
  prj_A = zeros(m);
  for i = 1:m   % regularize on a row by row basis;
      a = A(i,:);
      a = a - 1/m * (sum(a)-1);  % projection of a on the hyperplane;
      
      while any(a<0)  % iterative algorithm published by Kreinin et al.
          a(a<0) = 0;
          ind = find(a>0);
          b = a(ind);          
          k = length(b);
          b = b - 1/k *(sum(b)-1);
          a(ind) = b;
      end
      
      prj_A(i,:) = a;
  end
   prj_x = prj_A(:);
  