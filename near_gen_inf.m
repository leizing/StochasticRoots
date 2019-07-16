function G = near_gen_inf(L)
%NEAR_GEN_INF  Regularize L to its nearest intensity matrix in inf norm.
%   G = NEAR_GEN_INF(L) compute the nearest intensity matrix to given 
%   matrix L, where the distance is measured by matrix inf norm.

  G = real(L);  % discard the imaginary part;
  [m n] = size(G);
  G(G<0) = 0;
  G = G - diag(diag(G));
  G = G - diag(sum(G,2));
  