# Array spacing (m)
R = 1;
# Array size
# order gives the highest order derivative estimated
Order = 3;
# Record Length (sec)
T = 240;
# Wavelengths (m)
L = R*(55+45*rand(1,4));
# Amplitude
A = rand(size(L));
# Orientation (degrees)
#theta = pi*rand(1,1);
theta = pi*rand(size(L));
# create an array
x = [];
for i=1:Order
 for j=1:i
   x = [x;j-(i+1)/2,(i-1)*sin(pi/3)];
 end#for
end#for
x = R*x;
Nx = length(x(:,1));
x0 = mean(x,1);
# Assume 10Hz sample rate
t = [0:.1:T];
Nt = length(t);
# Add noise
h = rand(Nx,Nt)/50;
#Add all the waves
for i=1:length(L)
# Wave charachteristics
 k = [cos(theta(i)),sin(theta(i))]*2*pi/L(i);
 omega = sqrt(10*sqrt(sum(k.^2)));
# Wave height
 h = h+A(i).*sin(x*k'.-omega*t);
end#for
#
# matrix of diferences
# A =
#     1   ,1   ,1   ,...
#     dx  ,dx  ,dx  ,...
#     dy  ,dy  ,dy  ,...
#     dx^2,dx^2,dx^2,...
#     dy^2,dy^2,dy^2,...
#     dxdy,dxdy,dxdy,...
#     .   ,.   ,.   ,.
#     .   ,.   ,.   ,.
#     .   ,.   ,.   ,.
#
Taylor = ones(Nx,Nx);
Taylor(1,:) = 1;
# dx,dy
if(Nx>2)
 Taylor(2,:) = x(:,1)-x0(1);
 Taylor(3,:) = x(:,2)-x0(2);
end#if
# dx^2, dx*dy, dy^2
if(Nx>5)
 Taylor(4,:) = .5*(x(:,1)-x0(1)).^2;
 Taylor(5,:) = .5*(x(:,2)-x0(2)).^2;
 Taylor(6,:) = (x(:,2)-x0(2)).*(x(:,1)-x0(1));
end#if
#
# oops that should be transpossed
Taylor = Taylor';
# Taylor expansion looks like
#
#         mean(h)
#         dh/dx
#         dh/dy
#Taylor * d^2h/dx^2 = h
#         d^2h/dy^2
#         d^2h/dxdy
#         .
#         .
#         .
#
#
# Invert the matrix
Taylorinv = inv(Taylor);
#
# Calculate derrivatives
hders = Taylorinv*h;
#
#
figure(1)
subplot(2,1,1)
plot(t,atan2(hders(3,:),hders(2,:))*180/pi)
axis([0,T,-180,180])
xlabel("time (sec)")
ylabel("Wave face direction")
subplot(2,2,3)
plot(hders(2,:),hders(3,:),".;;");
axis(length(L)*pi*[-1,1,0,1]/min(L));
axis equal;
xlabel("A*k_x");
ylabel("A*k_y");
drawnow
subplot(2,2,4)
plot(+A.*2.*pi.*cos(theta)./L,+A.*2.*pi.*sin(theta)./L,"+;+wave vectors;",-A.*2.*pi.*cos(theta)./L,-A.*2.*pi.*sin(theta)./L,".;-wave vectors;");
axis(2*pi*[-1,1,0,1]/min(L));
axis equal;
xlabel("A*k_x");
ylabel("A*k_y");
drawnow
