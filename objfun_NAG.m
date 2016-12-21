function [mode, objf, objgrd, user] = objfun_NAG(mode, n, x, objgrd,...
                                      nstate, user)
%OBJFUN_NAG calculate the objective function and (optionally) its gradient
%for e04uc. 
%   For details of e04uc, see:
%   http://www.nag.co.uk/numeric/MB//manual_21_1/pdf/E04/e04uc.pdf

A = user{1};
p = user{2};
tol = user{3};
m = int32(sqrt(double(n)));
M = reshape(x,m,m);
Diff = M^p - A;

if (mode == 0 || mode == 2)
objf = trace( Diff' * Diff);
else
objf = 0;
end

if (mode == 1 || mode == 2)  % calculate the gradient of the object funtion
    Grad = 0;
    for i=1:p
       Grad = Grad + M'^(i-1)* Diff * M'^(p-i);
    end
objgrd = 2 * Grad(:);

if norm(spgm_proj(x-objgrd) - x) < tol
	mode = int32(-1);
end
end