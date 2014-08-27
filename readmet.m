function [A,B,X,PT,T] = readmet(metname)
addpath ~/Documents/cor3_5/
metdir = "/Volumes/current_cruise/met/data/";
fluxdir = "/Volumes/scienceparty_share/FluxData/";
nowutc = gmtime(time());
if(nargin<1)
 metname = [num2str(mod(nowutc.year,100))];
 if(nowutc.mon<9) 
  metname = [metname "0" num2str(nowutc.mon+1)];
 else
  metname = [metname num2str(nowutc.mon+1)];
 end%if
 if(nowutc.mday<10) 
  metname = [metname "0" num2str(nowutc.mday)];
 else
  metname = [metname num2str(nowutc.mday)];
 end%if
end%if
zu = 17.06;
zt = 17.06;
zq = 17.06;
MET = load([metdir metname ".MET"]);
secs = mod(MET(:,1),100);
mins = mod(floor(MET(:,1)/100),100);
hous = floor(MET(:,1)/10000);
#Time
T = secs+60*(mins+60*hous);
#Air Temp
AT = MET(:,2);
# Barometric Pressure (Pa)
BP = MET(:,3);
# Barometric Pressure Temp (C)
BC = MET(:,4);
#
BS = MET(:,5);
# Accumulated Precipitation (mm)
PR = MET(:,6);
# Precipitation Rate
PT = diff([PR;PR(end)])./mean(diff(T));
# Relative Humidity
RH = MET(:,7);
# Relative humidity (temp)
RT = MET(:,8);
# Dew Point
DP = MET(:,9);
# LW body Temp (K)
LB = MET(:,10);
# LWR Thermopile (V)
LT = MET(:,11);
# Long Wave Radiation (W/m^2)
LW = MET(:,12);
# Short Wave Radiation (W/m^w)
SW = MET(:,13);
# Photosynthetically available Radiation (uE/s/m^2)
PA = MET(:,14);
# Relative Wind Speed
WS = MET(:,15);
# Relative Wind Directions
WD = MET(:,16);
# True Wind Speed
TW = MET(:,17);
# True Wind Direction
TI = MET(:,18);
# Relative Wind speed 2
WS2 = MET(:,19);
# Relative wind direction 2
WD2 = MET(:,20);
# True wind speed 2
TW2 = MET(:,21);
# True wind direction 2
TI2 = MET(:,22);
# Thermosalinograph Temp (C)
TT = MET(:,23);
# Thermosalinograph conductivity (mS/cm)
TC = MET(:,24);
# Salinity (psu)
SA = MET(:,25);
# sigma_t
SD = MET(:,26);
# Sound Velocity (m/s)
SV = MET(:,27);
# Thermosalinograph uncorrected conductivity (mS/cm)
TG = MET(:,28);
# USW Flow Meter (L/min)
FI = MET(:,29);
# Pressure (psi)
PS = MET(:,30);
# Thermosalinograph Temp 2 (C)
TT2 = MET(:,31);
# Thermosalinigraph Conductivity (mS/cm)
TC2 = MET(:,32);
# Salinity-2 (psu)
SA2 = MET(:,33);
# sigma_t 2
SD2 = MET(:,34);
# Sound Velocity 2 (m/s)
SV2 = MET(:,35);
# Thermosalinograph uncorrected conductivity (mS/cm)
TG2 = MET(:,36);
# Oxygen Current (uA)
OC = MET(:,37);
# Oxygen Temperature (C)
OT = MET(:,38);
# Oxygen (mL/L)
OX = MET(:,39);
# Oxygen Saturation Value (mL/L)
OS = MET(:,40);
# Flourometer
FL = MET(:,41);
# USW Flow meter 2
FI2 = MET(:,42);
# Pressure 2
PS_2 = MET(:,43);
# VRU Pitch (deg)
VP = MET(:,44);
# VRU Roll (deg)
VR = MET(:,45);
# VRU heave (m)
VH = MET(:,46);
# Ship List
VX = MET(:,47);
# Ship Trim
VY = MET(:,48);
# Ship Gyrocompass Heading (deg)
GY = MET(:,49);
#
MB = MET(:,50);
# Bottom Depth (m)
BT = MET(:,51);
#
LF = MET(:,52);
#
HF = MET(:,53);
# Lattitude
LA = MET(:,54);
# Longitude
LO = MET(:,55);
#GPS time
GT = MET(:,56);
# Ships Course over ground (deg)
CR = MET(:,57);
#Ship Speed over ground (kts)
SP = MET(:,58);
#GPS Datetime
ZD = MET(:,59);
#GPS altitude (m)
GA = MET(:,60);
# GPS Status and # of satellites
GS = MET(:,61);
#Ashtech Heading
SH = MET(:,62);
#Ashtech Ptch
SM = MET(:,63);
#Ashtech roll
SR = MET(:,64);
#Winch wire out
ZO = MET(:,65);
#Winch Speed
ZS = MET(:,66);
#Winch Tension
ZT = MET(:,67);
#
ZI = MET(:,68);
# Winch wire out 2
ZO2 = MET(:,69);
#winch Speed 2
ZS2 = MET(:,70);
# Winch tension 2
ZT2 = MET(:,71);
#Input into core3.5
Flux1=coare35vn(TW,zu,AT,zt,RH,zq,BP,(TT+TT2)/2,SW,LW,LA,NaN,PT,NaN,NaN);
Flux2=coare35vn(TW2,zu,AT,zt,RH,zq,BP,(TT+TT2)/2,SW,LW,LA,NaN,PT,NaN,NaN);
save("-v7",[fluxdir metname ".mat"],"Flux1","Flux2","MET","PT","T")