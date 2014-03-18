t = linspace(0,2,8192);
x = sin(2*pi*t)+(rand(size(t))-.5+rand(size(t))-.5+rand(size(t))-.5+rand(size(t))-.5+rand(size(t))-.5);
y = x;
xm = x;
ym = y;
R = (2)^2;
for i=2:length(x)
 y(i) = (x(i)+y(i-1)/R)*(R-1)/R;
end%for
 M = 32;
 ysub = y(1:M:end);
 tsub = t(1:M:end);
 ym = cumsum(ysub)*mean(diff(tsub));
 ym = cumsum(ym)*mean(diff(tsub));
 xm = cumsum(x)*mean(diff(t));
 xm = cumsum(xm)*mean(diff(t));
 xmsub = xm(1:M:end);
subplot(4,1,1)
plot(t,x,";x;")
subplot(4,1,2)
plot(tsub,ysub,";y;")
subplot(4,1,3)
plot(tsub,xmsub,";x_m;",tsub,ym,";y_m;")
subplot(4,1,4)
semilogy(t,abs(x-y)./(abs(x)+abs(y)),tsub,abs(xmsub-ym)./(abs(xmsub)+abs(ym)))
