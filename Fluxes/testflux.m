addpath ~/Documents/cor3_5/
u = 8; #wind speed m/s
zu = 16; #height of wind speed measurement
t = 30.2; #air temp
zt = 16; #height of air temp
rh = 78.5; #Relative humidity %
zq = 16; #height of relative humidity
P = 1003.3; #Pressure 
ts = 30.1; #Water Temp
Rs = 0; #Downward Short Wave
Rl = 439.2; #Downward Long Wave Rad
lat = 22.14; #Lattitude
zi = NaN; #PBL height?
rain = 0; #rain rate
cp = 10; #wave speed
sigH = 1 ; #wave height
A=coare35vn(u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi,rain,cp,sigH)
