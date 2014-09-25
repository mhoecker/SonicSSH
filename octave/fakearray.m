# Array spacing (m)
R = 1;
# Array size
Order = 3;
# Record Length (sec)
T = 60;
# Wavelength (m)
L = 100*R;
# Orientation (degrees)
#theta = pi*rand(1,1);
theta = pi/2;
# create an array
x = [];
for i=1:Order
 for j=1:i
   x = [x;j-(i+1)/2,(i-1)*sin(pi/3)];
 end%for
end%for
x = R*x;
Nx = length(x(:,1));
x0 = mean(x,1);
# Assume 10Hz sample rate
t = [0:.1:T];
# Wave charachteristics
k = [cos(theta),sin(theta)]*2*pi/L;
omega = sqrt(10*sqrt(sum(k.^2)));
#
# Wave height
h = sin(x*k'.-omega*t);
# Add noise
h = h+rand(size(h))/20;
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
A = ones(Nx,Nx);
A(1,:) = 1;
# dx,dy
if(Nx>2)
 A(2,:) = x(:,1)-x0(1);
 A(3,:) = x(:,2)-x0(2);
end%if
# dx^2, dx*dy, dy^2
if(Nx>5)
 A(4,:) = .5*(x(:,1)-x0(1)).^2;
 A(5,:) = .5*(x(:,2)-x0(2)).^2;
 A(6,:) = (x(:,2)-x0(2)).*(x(:,1)-x0(1));
end%if
#
# oops that should be transpossed
A = A';
# Taylor expansion looks like
#
#       mean(h)
#       dh/dx
#       dh/dy
# A  *  d^2h/dx^2 = h
#       d^2h/dy^2
#       d^2h/dxdy
#         .
#         .
#         .
#
#
# Invert the matrix
Ainv = inv(A);
#
# Calculate derrivatives
hders = Ainv*h;
#
#
figure(1)
plot(t,h-hders(1,:))
figure(2)
subplot(2,1,1)
plot(t,atan2(hders(3,:),hders(2,:))*180/pi)
axis([0,T,-180,180])
xlabel("time (sec)")
ylabel("Wave face direction")
subplot(2,1,2)
plot(hders(2,:),hders(3,:),"o;;");
axis(4*pi*[-1,1,-1,1]/L);
axis equal;
xlabel("A*k_x");
ylabel("A*k_y");
drawnow
