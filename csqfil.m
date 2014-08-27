function [y] = csqfil(x,t,s,T)
% [y,fil] = csqfil(x,t,T)
% Filter the time series x with a weighting function
%
% weight(t,s) = (1+cos(2*pi*(s-t)/T))^2 for |s-t| < T/2
% This filter has a 1/2 power point of f = 1/T
% first zero power at f~3/T
% Roll off of (T*f)^(-10)
% if T is not given T = 3*mean(diff(s))
%
 if(nargin==0)
 % plot the filter frequency response
 % compare a moving average of width T/3
  N = 2^16;
  f = (1:N)-1-ceil(N/2);
  t = f/N; t = t-mean(t);
  f = f/sqrt(N);
  t = t*sqrt(N);
  T = 1.0;
  wa = (1+sign(1-(6*t/T).^2));
  wa = wa./sum(wa);
  Fwa = fftshift(fft(wa));
  Pwa = abs(Fwa).^2;
  w = csqwt(t,T);
  w = w/sum(w);
  Fw = fftshift(fft(w));
  Pw = abs(Fw).^2;
%
% Pass band
  subplot(2,2,1);
  plot(f,Pw,"b;cosq  ;",f,Pwa,"r;T/3 Moving avg. ;");
  grid on;
  axis([0,1/T,.4,1]);
%  title(["cosq filter 1/T=" num2str(1/T)]);
  xlabel("freq.");
  ylabel("Filter");
% Filtered band
  subplot(2,2,2);
  loglog(f,Pw,"b;cosq ;",f,(T*f).^-10,"k;fT^{-10};",f,Pwa,"r;T/3 moving avg. ;",f,(T*f).^-2,"k;fT^{-2};");
  axis([1/T,128/T,128^(-6),2]);
%  title(["cosq filter 1/T=" num2str(1/T)]);
  xlabel("freq.");
  ylabel("Filter");
% weights
  subplot(2,1,2);
  plot(t,w,"b;cosq ;",t,wa,"r;T/3 moving avg. ;");
  axis([-T/2,T/2,0,1.05*max([w,wa])]);
  title(["cosq filter, T=" num2str(T)]);
  xlabel("Time diff.");
  ylabel("Weight")
  print("csqfilldemo.png","-dpng","-S1280,1024","-F:8")
 else
  if(nargin<4)
   T = 3*mean(abs(diff(s)));
  end%if
  N = length(t);
  M = length(s);
% ensure there is nothing to alias into the new sampling
  y = zeros(1,M);
  for i=1:M
   w = zeros(1,N);
   for j=1:N
    if((abs(s(i)-t(j))<T/2)*(1-isnan(x(j))))
     w(j) = csqwt(s(i)-t(j),T);
     y(i) = y(i)+x(j)*w(j);
    end%if
   end%for	
    wsum = sum(w);
    if(wsum==0)
     y(i) = NaN;
    else
     y(i) = y(i)/wsum;
    end%if
  end%for
 end%if
end%function

function w = csqwt(dt,T)
 w=(1+sign(1-(2*dt/T).^2)).*(.5+.5.*cos(2.*pi.*(dt./T))).^2;
end%function
