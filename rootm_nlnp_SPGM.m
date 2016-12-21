function [X0 nlnpA niter time resultsf,resultsg,resultsa,resultsg0, nfval,ngval,xs] ...
            = rootm_nlnp_SPGM(A,p,start)
%ROOTM_NLNP_SPGM compute an approximate pth root of matrix A
%   [rA nlnpA] = ROOTM_NLNP_MB(A,p) compute an approximate stochastic pth 
%   root of a stochastic matrix A by solving nonlinear programming problem 
%   via MATLAB function fmincon in Optimization Toolbox. rA returns an
%   approximate root obtained from adjusting the principal root to a 
%   nearest stochastic matrix if the principal root is not stochastic. 
%   nlnpA returns the final optimal solution by solving NLP via fmincon 
%   using rA as an initial estimate.
% 

m = int32(length(A));

% n = int32(m^2);  % n, the number of variables in nonlinear programming;
switch start
    case 'PrincRoot'  % rA from regularization as starting point
        rA = rootm(A,p);  %  the principal root of A
        rA = reg_to_st(rA);   % adjusting to a nearest stochastic matrix
        X0 = rA;      
    case 'Ident' % identity matrices5 as starting point
        X0 = eye(m);    
    case 'UTri'
        D = diag(diag(A)); % initial estimate: Case 2
        X0 = rootm(D,p);
        for i = 1:m-1
            X0(i,m) = 1-X0(i,i);
        end
        X0(m,m) = 1;
    case 'FullRow'
        D = diag(diag(A)); % initial estimate: Case 4
        X0 = rootm(D,p);
        for i = 1:m
            for j = 1:m
                 if i ~= j 
                     X0(i,j) = (1-X0(i,i))/(m-1);
                 end
            end 
        end
    case 'StoRand'
        X0 = rand(m);
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

x0 = X0(:);

tol = 1e-6;
sn = 1e-3;
sx = 1e3;
Mspgm = 10;

%     resultsa = func_grad2(x0,A,p,B);
     [resultsa g0] = spgm_obj(x0,A,p);
     resultsg0 = norm(spgm_proj(x0-g0)-x0);   % initial q(X)
    %one iteration of the 'altk' method, to ensure a decrease
    
    tic;
%         [resultsf,resultsg,niter,nfval,xs] = ...
%             SPGM(@(x)fmincon_grad2(x,ln,kn),@(x,g)compute_proj_grad(x,ln,kn,'nl',g),...
%             @(x)projq(x,ln,kn),x0,accuracy,ln*10*kn,0.1/(kn*ln),Mspgm);
                     
        [resultsf,resultsg,niter,nfval,xs] = ...
            SPGM_2( @(x)spgm_obj(x,A,p), @(x)spgm_proj(x),...
                    x0,sn,sx,tol,Mspgm);
    time = toc;
    ngval = nfval;
nlnpA = reshape(xs,m,m);


% From Rued's code
% function [resultsa,resultsf,resultsg,nfval,ngval,niter,time,xs]=spgmspec(x0,ln,kn,Mspgm)
%     global accuracy;
%     resultsa = fmincon_grad2(x0,ln,kn);
%     %one iteration of the 'altk' method, to ensure a decrease
%      
%     x0 = reshape(x0,[],1);
%     tic;
%         [resultsf,resultsg,niter,nfval,xs] = ...
%             SPGM(@(x)fmincon_grad2(x,ln,kn),@(x,g)compute_proj_grad(x,ln,kn,'nl',g),...
%             @(x)projq(x,ln,kn),x0,accuracy,ln*10*kn,0.1/(kn*ln),Mspgm);
%     time = toc;
%     ngval = nfval;