function [X0, nlnpA, iter, time, objf, objgrd, r, x, istate, c, cjac, clamda,  ...
          user, lwsav, iwsav, rwsav, ifail] = rootm_nlnp_NAG(A,p,start)
%ROOTM_NLNP_NAG compute an approximate pth root of matrix A
%   [X1 nlnpA] = ROOTM_NLNP_NAG(A,p) compute an approximate pth root of
%   a stochastic matrix A by solving nonlinear programming problem via
%   e04uc in NAG Toolbox for MATLAB. X1 returns an approximate root obtained
%   from adjusting the principal root to a nearest stochastic matrix if the
%   principal root is not stochastic. nlnpA returns the final optimal
%   solution by solving NLP using X1 as an initial estimate.
%   For details of e04uc, see:
%   http://www.nag.co.uk/numeric/MB//manual_21_1/pdf/E04/e04uc.pdf
global tol_global
n = int32(length(A));
m = int32(n^2);  % n, the number of variables in nonlinear programming;
switch start
    case 'PrincRoot'  % rA from regularization as starting point
        rA = rootm(A,p);  %  the principal root of A
        rA = reg_to_st(rA);   % adjusting to a nearest stochastic matrix
        X0 = rA;      
    case 'Ident' % identity matrices5 as starting point
        X0 = eye(n);    
    case 'UTri'
        D = diag(diag(A)); % initial estimate: Case 2
        X0 = rootm(D,p);
        for i = 1:n-1
            X0(i,n) = 1-X0(i,i);
        end
        X0(n,n) = 1;
    case 'FullRow'
        D = diag(diag(A)); % initial estimate: Case 4
        X0 = rootm(D,p);
        for i = 1:n
            for j = 1:n
                 if i ~= j 
                     X0(i,j) = (1-X0(i,i))/(n-1);
                 end
            end 
        end
    case 'StoRand'
        X0 = rand(n);
        X0 = reg_to_st(X0);
    case 'GenFro'
        L = logm(A);
        G = near_gen_fro(L);
        X0 = expm(1/p*G);
    case 'GenInf'
        L = logm(A);
        G = near_gen_inf(L);
        X0 = expm(1/p*G);
    case 'GenWA'
        L = logm(A);
        G = near_gen_WA(L);
        X0 = expm(1/p*G);       
end
x = X0(:);
a = kron(ones(1,n),eye(n));
bl = [zeros(m,1);ones(n,1)];
bu = ones(n+m,1);

istate = zeros(n+m,1,'int32');   % initialization for e04uc
cjac = 0;
clamda = zeros(n+m);

r = zeros(m,m);

%I = eye(n); 
%x= I(:);   % use identity matrix as initial value
%[cwsav,lwsav,iwsav,rwsav] = e04wb('e04uc');
[temp,lwsav,iwsav,rwsav] = e04wb('e04uc');
tic;
[iter, istate, c, cjac, clamda, objf, objgrd, r, x, user, lwsav, iwsav, rwsav, ifail] = ...
    e04uc(a, bl, bu, 'd04udm', 'objfun_NAG', ...
     istate, cjac, clamda, r, x, lwsav, iwsav, rwsav, 'n', m, 'nclin', n, ...
     'ncnln', int32(0), 'user',{A p tol_global});
time = toc;

nlnpA = reshape(x,n,n);