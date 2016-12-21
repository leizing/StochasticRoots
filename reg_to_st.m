function rA = reg_to_st(A)
%REG_TO_ST  Regularize A to its nearest stochastic matrix.
%   rA = REG_TO_ST(A) compute the nearest stochastic matrix to given matrix
%   A, where the distance is measured by matrix Frobenius norm.

  A = real(A);  % discard the imaginary part;
  [m n] = size(A);
  rA =[];
  for i = 1:m   % regularize on a row by row basis;
      a = A(i,:);
      a = a - 1/n * (sum(a)-1);  % projection of a on the hyperplane;
      
      while any(a<0)  % iterative algorithm published by Kreinin et al.
          a(a<0) = 0;
          ind = find(a>0);
          b = a(ind);          
          k = length(b);
          b = b - 1/k *(sum(b)-1);
          a(ind) = b;
      end
      
      rA(i,:) = a;
   end
  