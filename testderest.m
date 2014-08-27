#test the derivative fitter
load /Users/mhoecker/Documents/SonicSSH/scripts/AvgVel75kHz150to200m.mat
who
N = 300;
Ni = 500;
Nf = Ni+N;
#t = linspace(0,1,N);
#course = -pi/4;
#deviation = course+pi/2;
#devamp = .01;
#x = 45+sin(course)*t+devamp*sin(2*pi*t)*sin(deviation);
#y = 45+cos(course)*t+devamp*sin(2*pi*t)*cos(deviation);
#u = -y./sqrt(x.^2+y.^2);
#v = x./sqrt(x.^2+y.^2);
t = t(Ni:Nf);
x = lon(Ni:Nf)';
y = lat(Ni:Nf)';
u = uavg(Ni:Nf)';
v = vavg(Ni:Nf)';
[dudx,dudy,dvdx,dvdy,vor,div,vbyf] = strainrateest(u,v,x,y,15000,5000,32);
subplot(2,2,1)
quiver(x,y,u,v)
hold on
plot(x,y)
hold off
subplot(2,2,2)
plot(vor,y,div,y)
subplot(2,2,3)
plot(x,vor,x,div)
subplot(2,2,4)
plot(t,vor,t,div)
