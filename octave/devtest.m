N = 10;
#clear x,y,dxy,ordx,ordy;
if(N>0)
 x = [0];
 y = [0];
end%if
for i=2:N
 for k=1:i
  x = [x,k-.5-i/2.0];
  y = [y,-(i-1)*cos(pi/6)];
 end%for
end%for
# Remove the mean so that x and dx are the same
y = y-mean(y);
# Scale so that the aray size is 1
x = x/N;
y = y/N;
#
dxy = [];
ordx = [];
ordy = [];
for i=0:N-1
 for j=i:N-1
  ordx = [ordx,j-i];
  ordy = [ordy,i];
  dxy = [dxy;(x.^(j-i).*y.^(i))/(factorial(j-i)*factorial(i))];
 end%for
end%for
ixy = inv(dxy');
z = ones(size(x));
z = z+x+y;
z = z+x.*x/2+x.*y+y.*y/2;
z = z+x.*x.*x/6+x.*x.*y/2+x.*y.*y/2+y.*y.*y/6;
z = exp(x)+exp(y);
f = [ixy*z'];
F = zeros([N,N]);
for i=1:length(ordx)
  F(ordx(i)+1,ordy(i)+1)=f(i);
end%for
for i=2:N
   F(i,N-i+2:end)=NaN;
end%for
F
f = ixy*x';
f = ixy*y';
f = ixy*(x'.*x');
f = ixy*(x'.*y');
f = ixy*(y'.*y');
#
subplot(2,2,1)
plot(x,y,"+")
axis([-1,1,-1,1])
subplot(2,2,2)
semilogy(abs(eig(dxy)),"-+")

