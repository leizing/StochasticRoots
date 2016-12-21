function G = near_gen_fro(L)
%NEAR_GEN  Regularize L to its nearest intensity matrix in Frobenius norm.
%   G = NEAR_GEN(L) compute the nearest intensity matrix to given matrix
%   L, where the distance is measured by matrix Frobenius norm.

  L = real(L);  % discard the imaginary part;
  [m n] = size(L);
  G = zeros(m,n);
  
  for i = 1:m   % regularize on a row by row basis;
      a = L(i,:);
      a = a - (sum(a))/n;  % projection of a on the hyperplane;
      [a ind] = sort(a,'descend');
      flag = 0;
      for l = 2:n-1
          if a(l+1) >= sum(a(l:n))/(n-l+1)
              flag = 1;
              break
          end
      end
      if flag 
          g = a - (sum(a) - sum(a(2:l)))/(n-l+1);
          g(2:l) = 0;
      else
          g = a;
      end
      [temp rev_ind] = sort(ind,'ascend');
      g = g(rev_ind);
      G(i,:) = g;
   end
  