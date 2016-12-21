function [f g] = objfun_MB(x,A,p)
%OBJFUN_MB calculate the objective function and (optionally) its gradient
%for MATLAB function fmincon. 
% 
% global A p;
n = length(x);
m = int32(sqrt(double(n)));
M = reshape(x,m,m);
Diff = M^p - A;

f = trace( Diff' * Diff);

if nargout > 1 % gradient required
    Grad = 0;
    for i=1:p
       Grad = Grad + M'^(i-1)* Diff * M'^(p-i);
    end
    g = 2 * Grad(:);
end




