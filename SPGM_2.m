function [fk,resultsg,k,nf,xk]=SPGM_2(func,proj,x0,sign,sigx,tol,M)
% function [fk,resultsg,k,nf,xk,temp]=SPGM_2(func,proj,x0,ln,lx,tol,M)
% [f,grad] = func(x);
% d = proj(x)

% Algorithm based on the paper:
% Algorithm 813: {SPG}---{Software} for Convex-Constrained Optimization
% by Birgin, Ernesto G. and Mart\'{\i}nez, Jos\'{e} Mario and Raydan, Marcos
% 2001

k = 1;
if (nargin < 7) || (isempty(M))
    M = 4;
end
gamma =  1e-2;%0.5;
tau1 = 0.1;
tau2 = 0.9;
farray = zeros(M,1);
xk = x0; [fk,gk] = func(xk);
sigma = min(sigx,max(sign,norm(proj(xk-gk)-xk,inf)));
% lbda = min(lx,max(ln,1/norm(proj(xk-gk)-xk)));
nf = 1;

% temp(k) = lbda;
while (k < 40000)
    farray = [fk;farray(1:M-1)];
    fmax = max(farray(1:min(k,M)));
    dk = proj(xk-1/sigma*gk)-xk;
%     dk = proj(xk-lbda*gk)-xk;
    delta = gk'*dk;
    alpha = 1;
   
    %Test nonmonotone Armijo-like criterion
    while 1
        xp = xk + alpha * dk;
        [fp,gp] = func(xp); 
%         fp = func(xp); 
        nf = nf + 1;
        if fp <= fmax + gamma * alpha * delta
            break
        end
        %compute a safeguarded new trial steplength
        alphatmp = -0.5*alpha^2*delta/(fp-fk-alpha*delta);
        if (alphatmp >= tau1) && (alphatmp <= alpha * tau2)
            alpha = alphatmp;
        else
            alpha = alpha/2;
        end
    end
    
    sk = xp - xk; yk = gp - gk;%sk'*yk/(sk'*sk)
    sigma = min(sigx,max(sign,sk'*yk/(sk'*sk)));
%     lbda = min(lx,max(ln,(sk'*sk)/(sk'*yk)));

    fk = fp; gk = gp; xk = xp;
%     resultsg = projg(xk,gk);
%     if (resultsg < tol)
%         break
%     end
    resultsg = norm(proj(xk-gk)-xk);
    if (resultsg < tol)
        break
    end

    k = k + 1;
end
    





