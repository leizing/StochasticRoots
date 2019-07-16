function G = near_gen_WA(L)
%NEAR_GEN_WA  Regularize L to an intensity matrix.
%   G = NEAR_GEN_WA(L) compute an intensity matrix by weighted adjustment.

  G = real(L);  % discard the imaginary part;
  [m n] = size(G);
  dG = diag(G);
  G(G<0) = 0;
  G = G + diag(dG); % keep the diagonal of G even though it is negative.
% the line below will make error for zeros rows
%   G = G - abs(G) .* repmat(sum(G,2)./sum(abs(G),2),1,n);  
  for i = 1:m
      g = G(i,:);
      n1g = norm(g,1);
      if  n1g > eps*n  % nonzero row
          lbda = sum(g)/n1g;
          for j = 1:n
              if G(i,j) ~= 0
                  G(i,j) = G(i,j) - abs(G(i,j))*lbda;
              end
          end
      end
  end
          
  


  
  