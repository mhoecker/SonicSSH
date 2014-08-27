t = 1:2^10;
N = length(t);
dt = mean(diff(t));
f = (cumsum(ones(size(t)))-1-floor(N/2))/dt;
ftest =256.25;
Nbin = 16;
x = sin(2*pi*ftest*t/max(t));
xb = floor(x*Nbin)./Nbin;
x2b = floor(x*2*Nbin)./(2*Nbin);
x4b = floor(x*4*Nbin)./(4*Nbin);
Fx = fftshift(fft(x))/N;
Fxb = fftshift(fft(xb))/N;
Fx2b = fftshift(fft(x2b))/N;
Fx4b = fftshift(fft(x4b))/N;
loglog(abs(f),abs(Fx4b),abs(f),abs(Fx2b),abs(f),abs(Fxb))
