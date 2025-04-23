function [x, prob] = normal_discrete_1(n, mu, sigma)
% Creates n points and probabilities for a normal distribution .
%
% PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
%     CompEcon toolbox of Mario Miranda and Paul Fackler available under
%     http://www4.ncsu.edu/~pfackler/compecon/toolbox.html
%
% REFERENCE: Miranda, M. & Fackler, P. (2002). Applied Computational Economics
%            and Finance. Cambridge: MIT Press.
% INPUTS
%     n:     Num. of points for discretization
%     mu:    Mean of x
%     sigma: Standard deviation of x
% OUTPUTS
%     x:     Grid (a vector)
%     prob:  Probability mass function (a vector)

maxit = 200;

x = zeros(n,1);
prob = zeros(n,1);

% calculate 1/pi^0.25
pim4 = 1/pi^0.25;

% get number of points (best if n is odd)
m = round((n+1)/2);

% start iteration
for i = 1:m

    % set reasonable starting values
    if(i == 1)
        z = sqrt(2*n+1)-1.85575*((2*n+1)^(-1/6));
    elseif(i == 2)
        z = z - 1.14d0*(n^0.426)/z;
    elseif(i == 3)
        z = 1.86*z+0.86*x(1);
    elseif(i == 4)
        z = 1.91*z+0.91*x(2);
    else
        z = 2*z+x(i-2);
    end

    % root finding iterations
    its = 0;
    while (its < maxit)
        its = its+1;
        p1 = pim4;
        p2 = 0;
        for j = 1:n
            p3 = p2;
            p2 = p1;
            p1 = z*sqrt(2/j)*p2-sqrt((j-1)/j)*p3;
        end
        pp = sqrt(2*(n))*p2;
        z1 = z;
        z  = z1-p1/pp;
        if(abs(z-z1) < 1e-14)
            break
        end
        if(its >= maxit)
            error('normal_discrete: Could not discretize normal distribution')
        end
        x(n+1-i) = z;
        x(i) = -z;
        prob(i) = 2/pp^2;
        prob(n+1-i) = prob(i);
    end
end

% set output data
prob = prob/sqrt(pi);
x = x*sqrt(2)*sigma + mu;

end %end function