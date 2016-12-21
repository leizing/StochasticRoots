function [X0, nlnpA, niter,time, fval,grad, exitflag,output,lambda]...
                = rootm_nlnp_MB(A,p,alg,start)
%ROOTM_NLNP_MB compute an approximate pth root of matrix A
%   [rA nlnpA] = ROOTM_NLNP_MB(A,p) compute an approximate stochastic pth 
%   root of a stochastic matrix A by solving nonlinear programming problem 
%   via MATLAB function fmincon in Optimization Toolbox. rA returns an
%   approximate root obtained from adjusting the principal root to a 
%   nearest stochastic matrix if the principal root is not stochastic. 
%   nlnpA returns the final optimal solution by solving NLP via fmincon 
%   using rA as an initial estimate.
global tol_global
m = int32(length(A));
n = int32(m^2);  % n, the number of variables in nonlinear programming;
Aineq = [];      % no linear inequality constraint
bineq = [];
Aeq = kron(ones(1,m),eye(m));   %  linear equality constraints
beq = ones(m,1);
lb = zeros(n,1);  %  lower bound of the variables
ub = ones(n,1);   %  upper bound of the variables 
nonlcon = [];     %  no nonlinear constraint
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

options = optimset('Algorithm',alg,'GradObj','on','MaxIter',40000,...
                   'DerivativeCheck','off','TolFun',1e-8,'TolCon',1e-8,...
                   'OutputFcn',@outfun);

tic;
[x,fval,exitflag,output,lambda,grad] = fmincon(@(x)objfun_MB(x,A,p),x0,...
                            Aineq,bineq,Aeq,beq,lb,ub,nonlcon,options);
time = toc;
niter = output.iterations;
    % nested output function, to set my own stopping criterion
    function stop = outfun(x,optimValues,state)
        stop = false;
        % Check if directional derivative is less than tol.
        [fk gk] = objfun_MB(x,A,p);
        if norm(spgm_proj(x-gk) - x) < tol_global
            stop = true;
        end 
    end

nlnpA = reshape(x,m,m);
end