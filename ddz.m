function [d,dsq]=ddz(z)
% derivative matricies for independent variable z.
% Three point stencil.
% One-sided derivatives at boundaries.
% no assumptions about z spacing.
N=length(z);
d=zeros(N,N);
dsq=zeros(N,N);
for n=2:N-1
    % a is the point before
    % b is the point after
    dela = z(n-1)-z(n);
    delb = z(n+1)-z(n);
    denom = delb*dela^2-dela*delb^2;
    d(n,n) = -(dela^2-delb^2)/denom;
    d(n,n-1) = -(delb^2)/denom;
    d(n,n+1)= (dela^2)/denom;
    dsq(n,n)   = (2*dela-2*delb)/denom;
    dsq(n,n-1) = 2*delb/denom;
    dsq(n,n+1) = -2*dela/denom;
end
    % a is the first point after
    % b is the second point after
    dela = z(2)-z(1);
    delb = z(3)-z(1);
    denom = delb*dela^2-dela*delb^2;
    d(1,1) = -(dela^2-delb^2)/denom;
    d(1,2) = -(delb^2)/denom;
    d(1,3)= (dela^2)/denom;
    dsq(1,1)   = (2*dela-2*delb)/denom;
    dsq(1,2) = 2*delb/denom;
    dsq(1,3) = -2*dela/denom;
    % a is the first point before
    % b is the second point before
    dela = z(N-1)-z(N);
    delb = z(N-2)-z(N);
    denom = delb*dela^2-dela*delb^2;
    d(N,N) = -(dela^2-delb^2)/denom;
    d(N,N-1) = -(delb^2)/denom;
    d(N,N-2)= (dela^2)/denom;
    dsq(N,N)   = (2*dela-2*delb)/denom;
    dsq(N,N-1) = 2*delb/denom;
    dsq(N,N-2) = -2*dela/denom;
return
end
