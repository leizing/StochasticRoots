function sA = near_sto_fro(A)
%NEAR_GEN  Regularize L to its nearest intensity matrix in Frobenius norm.
%   G = NEAR_GEN(L) compute the nearest intensity matrix to given matrix
%   L, where the distance is measured by matrix Frobenius norm.

  sA = real(A);  % discard the imaginary part;
  [m n] = size(sA);

  
  for i = 1:m   % regularize on a row by row basis;
      a = sA(i,:);
      a = a - 1/n * (sum(a)-1);  % projection of a on the hyperplane;
      if a >= 0 
          sA(i,:) = a;
          continue
      end      
      [a ind] = sort(a,'descend');
      for k = 1:n
          c(k) = sum(a(1:k)) - k*a(k);
      end
      l = find(c <= 1, 1, 'last' );
      a = a + (1 - sum(a(1:l)))/l;
      a(l+1:n) = 0;
      [temp rev_ind] = sort(ind,'ascend');
      a = a(rev_ind);
         
      sA(i,:) = a;
   end
  