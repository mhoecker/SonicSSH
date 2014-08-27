function ADCPstrainrate(filenamein,filenameout)
#load /Users/mhoecker/Documents/SonicSSH/scripts/AvgVel75kHz150to200m.mat
load(filenamein)
x = min(lon):.2:max(lon);
y = min(lat):.2:max(lat);
[xx,yy] = meshgrid(x,y);
lonf = xx(:);
latf = yy(:);
tf = spatialfilter(lat,lon,t,lonf,latf,20000); 
u = spatialfilter(lat,lon,uavg,lonf,latf,20000);
v = spatialfilter(lat,lon,vavg,lonf,latf,20000);
[dudx,dudy,dvdx,dvdy,vor,div,vbyf] = strainrateest(u,v,x,y,200000,20000,128);
save("-v7",filenameout,"latf","lonf","tf","u","v","dudx","dudy","dvdx","dvdy","vor","div","vbyf")
end%function