function	RevelleADCPShear(outloc,daterange)
# Load ADCP data
nc150 = netcdf("/Volumes/current_cruise/adcp/RR1308/proc/nb150/contour/nb150.nc","r");
t150 = nc150{'time'}(:);
nc75 = netcdf("/Volumes/current_cruise/adcp/RR1308/proc/os75nb/contour/os75nb.nc","r");
t75 = nc75{'time'}(:);
if(nargin()<2)
daterange = [max([min(t150),min(t75)]),min([max(t150),max(t75)])];
end%if
t150idx = find((t150>min(daterange)).*(t150<max(daterange)));
t150 = nc150{'time'}(t150idx)+1;
lat150 = nc150{'lat'}(t150idx);
lon150 = nc150{'lon'}(t150idx);
u150 = nc150{'u'}(t150idx,:);
v150 = nc150{'v'}(t150idx,:);
d150 = nc150{'depth'}(t150idx,:);
ncclose(nc150);
bad150  = find((abs(u150)+abs(v150)>=1e38));
u150(bad150) = NaN;
v150(bad150) = NaN;
NaNcount = sum(isnan(u150));
somegood150 = find(NaNcount<length(t150));
u150 = u150(:,somegood150);
v150 = v150(:,somegood150);
d150 = d150(:,somegood150);
Su150 = zeros(size(u150));
Sv150 = zeros(size(v150));
tt150 = meshgrid(t150,d150(1,:));
#
for i=1:length(t150)
 utmp = u150(i,:)';
 vtmp = v150(i,:)';
 ztmp = d150(i,:)';
 tmpbad = find(isnan(utmp));
 utmpclean = utmp;
 vtmpclean = vtmp;
 utmpclean(tmpbad) = 0;
 vtmpclean(tmpbad) = 0;
 ddzbad = nanstencil(utmp,3,1);
 ddzmat = ddz(ztmp);
 Su150(i,:) = ddzbad+ddzmat*utmpclean;
 Sv150(i,:) = ddzbad+ddzmat*vtmpclean;
end%for
S150 = sqrt(Su150.^2+Sv150.^2);
#
t75idx = find((t75>min(daterange)).*(t75<max(daterange)));
t75 = nc75{'time'}(t75idx)+1;
lat75 = nc75{'lat'}(t75idx);
lon75 = nc75{'lon'}(t75idx);
u75 = nc75{'u'}(t75idx,:);
v75 = nc75{'v'}(t75idx,:);
d75 = nc75{'depth'}(t75idx,:);
ncclose(nc75);
bad75  = find(abs(u75)+abs(v75)>=1e38);
u75(bad75) = NaN;
v75(bad75) = NaN;
NaNcount75 = sum(isnan(u75));
somegood75 = find(NaNcount75<length(t75));
u75 = u75(:,somegood75);
v75 = v75(:,somegood75);
d75 = d75(:,somegood75);
Su75 = zeros(size(u75));
Sv75 = zeros(size(v75));
tt75 = meshgrid(t75,d75(1,:));
#
for i=1:length(t75)
 utmp = u75(i,:)';
 vtmp = v75(i,:)';
 ztmp = d75(i,:)';
 tmpbad = find(isnan(utmp));
 utmpclean = utmp;
 vtmpclean = vtmp;
 utmpclean(tmpbad) = 0;
 vtmpclean(tmpbad) = 0;
 ddzbad = nanstencil(utmp,3,1);
 ddzmat = ddz(ztmp);
 Su75(i,:) = ddzmat*utmpclean;
 Sv75(i,:) = ddzmat*vtmpclean;
 dbad = find(isnan(ddzbad));
 Su75(i,dbad) = 0;
 Sv75(i,dbad) = 0;
end%for
S75 = sqrt(Su75.^2+Sv75.^2);
#
save("-v7",[outloc "ADCPshear.mat"],"d150","d75","t150","t75","lat150","lon150","lat75","lon75","u150","v150","S150","Su150","Sv150","u75","v75","S75","Su75","Sv75")
dbad = find(isnan(S150));
S150(dbad) = 0;
Su150(dbad) = 0;
Sv150(dbad) = 0;
dbad = find(isnan(S75));
S75(dbad) = 0;
Su75(dbad) = 0;
Sv75(dbad) = 0;
z150 = -nanmean(d150);
z75 = -nanmean(d75);
#TZ plots
binmatrix(t150',(z150),S150',[outloc "S150kHzTZ.dat"])
binmatrix(t150',(z150),Su150',[outloc "Su150kHzTZ.dat"])
binmatrix(t150',(z150),Sv150',[outloc "Sv150kHzTZ.dat"])
binmatrix(t75',(z75),S75',[outloc "S75kHzTZ.dat"])
binmatrix(t75',(z75),Su75',[outloc "Su75kHzTZ.dat"])
binmatrix(t75',(z75),Sv75',[outloc "Sv75kHzTZ.dat"])
binmatrix(t150',(z150),u150',[outloc "u150kHzTZ.dat"])
binmatrix(t150',(z150),v150',[outloc "v150kHzTZ.dat"])
binmatrix(t75',(z75),u75',[outloc "u75kHzTZ.dat"])
binmatrix(t75',(z75),v75',[outloc "v75kHzTZ.dat"])
#LatZ plots
binmatrix(lat150',(z150),u150',[outloc "u150kHzLatZ.dat"])
binmatrix(lat150',(z150),v150',[outloc "v150kHzLatZ.dat"])
binmatrix(lat75',(z75),u75',[outloc "u75kHzLatZ.dat"])
binmatrix(lat75',(z75),v75',[outloc "v75kHzLatZ.dat"])
binmatrix(lat150',(z150),S150',[outloc "S150kHzLatZ.dat"])
binmatrix(lat150',(z150),Su150',[outloc "Su150kHzLatZ.dat"])
binmatrix(lat150',(z150),Sv150',[outloc "Sv150kHzLatZ.dat"])
binmatrix(lat75',(z75),S75',[outloc "S75kHzLatZ.dat"])
binmatrix(lat75',(z75),Su75',[outloc "Su75kHzLatZ.dat"])
binmatrix(lat75',(z75),Sv75',[outloc "Sv75kHzLatZ.dat"])
#LonZ
binmatrix(lon150',(z150),u150',[outloc "u150kHzLonZ.dat"])
binmatrix(lon150',(z150),v150',[outloc "v150kHzLonZ.dat"])
binmatrix(lon75',(z75),u75',[outloc "u75kHzLonZ.dat"])
binmatrix(lon75',(z75),v75',[outloc "v75kHzLonZ.dat"])
binmatrix(lon150',(z150),S150',[outloc "S150kHzLonZ.dat"])
binmatrix(lon150',(z150),Su150',[outloc "Su150kHzLonZ.dat"])
binmatrix(lon150',(z150),Sv150',[outloc "Sv150kHzLonZ.dat"])
binmatrix(lon75',(z75),S75',[outloc "S75kHzLonZ.dat"])
binmatrix(lon75',(z75),Su75',[outloc "Su75kHzLonZ.dat"])
binmatrix(lon75',(z75),Sv75',[outloc "Sv75kHzLonZ.dat"])
#
pltid = fopen([outloc "ADCP.plt"],'w');
[termtype,termtext]=termselect('pdfarticle');
# Set comon color palette
fprintf(pltid,'set palette defined (0 "dark-red", .25 "red", .4 "yellow", .5 "light-grey" , .6 "cyan", .75 "blue",1 "magenta")\n');
fprintf(pltid,'set view map\n')
fprintf(pltid,'unset surface\n')
fprintf(pltid,'unset key\n')
fprintf(pltid,'set pm3d\n')
fprintf(pltid,'set term %s\n',termtype)
# Set common velocity scale
vscale = sqrt(max([max(max(u150.^2+v150.^2)),max(max(u75.^2+v75.^2))]));
fprintf(pltid,'set cbrange[-%f:%f]\n',vscale,vscale)
# plot u150, v150, u75, v75
fprintf(pltid,'set ylabel "depth (m)"\n')
# TZ
fprintf(pltid,'set xlabel "date (2013 yearday)"\n')
fprintf(pltid,'set title" 150 kHz u"\n')
fprintf(pltid,'set output "%s"\n',[outloc "u150kHzTZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "u150kHzTZ.dat"])
fprintf(pltid,'set title" 150 kHz v"\n')
fprintf(pltid,'set output "%s"\n',[outloc "v150kHzTZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "v150kHzTZ.dat"])
fprintf(pltid,'set title" 75 kHz u"\n')
fprintf(pltid,'set output "%s"\n',[outloc "u75kHzTZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "u75kHzTZ.dat"])
fprintf(pltid,'set title" 75 kHz v"\n')
fprintf(pltid,'set output "%s"\n',[outloc "v75kHzTZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "v75kHzTZ.dat"])
# LatZ
fprintf(pltid,'set xlabel "Latitude (deg)"\n')
fprintf(pltid,'set title" 150 kHz u"\n')
fprintf(pltid,'set output "%s"\n',[outloc "u150kHzLatZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "u150kHzLatZ.dat"])
fprintf(pltid,'set title" 150 kHz v"\n')
fprintf(pltid,'set output "%s"\n',[outloc "v150kHzLatZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "v150kHzLatZ.dat"])
fprintf(pltid,'set title" 75 kHz u"\n')
fprintf(pltid,'set output "%s"\n',[outloc "u75kHzLatZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "u75kHzLatZ.dat"])
fprintf(pltid,'set title" 75 kHz v"\n')
fprintf(pltid,'set output "%s"\n',[outloc "v75kHzLatZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "v75kHzLatZ.dat"])
# LonZ
fprintf(pltid,'set xlabel "Longitude (deg)"\n')
fprintf(pltid,'set title" 150 kHz u"\n')
fprintf(pltid,'set output "%s"\n',[outloc "u150kHzLonZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "u150kHzLonZ.dat"])
fprintf(pltid,'set title" 150 kHz v"\n')
fprintf(pltid,'set output "%s"\n',[outloc "v150kHzLonZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "v150kHzLonZ.dat"])
fprintf(pltid,'set title" 75 kHz u"\n')
fprintf(pltid,'set output "%s"\n',[outloc "u75kHzLonZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "u75kHzLonZ.dat"])
fprintf(pltid,'set title" 75 kHz v"\n')
fprintf(pltid,'set output "%s"\n',[outloc "v75kHzLonZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "v75kHzLonZ.dat"])
# Set common Shear scale
Sscale = 2*sqrt(max([nanmean(nanmean(S150.^2,2)),nanmean(nanmean(S75.^2,2))]));
fprintf(pltid,'set cbrange[-%f:%f]\n',Sscale,Sscale)
# plot Su150, Sv150, Su75, Sv75
#TZ plots
fprintf(pltid,'set xlabel "date (2013 yearday)"\n')
fprintf(pltid,'set title" 150 kHz du/dz"\n')
fprintf(pltid,'set output "%s"\n',[outloc "Su150kHzTZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "Su150kHzTZ.dat"])
fprintf(pltid,'set title" 150 kHz dv/dz"\n')
fprintf(pltid,'set output "%s"\n',[outloc "Sv150kHzTZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "Sv150kHzTZ.dat"])
fprintf(pltid,'set title" 75 kHz du/dz"\n')
fprintf(pltid,'set output "%s"\n',[outloc "Su75kHzTZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "Su75kHzTZ.dat"])
fprintf(pltid,'set title" 75 kHz dv/dz"\n')
fprintf(pltid,'set output "%s"\n',[outloc "Sv75kHzTZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "Sv75kHzTZ.dat"])
#LatZ plots
fprintf(pltid,'set xlabel "Latitude (deg)"\n')
fprintf(pltid,'set title" 150 kHz du/dz"\n')
fprintf(pltid,'set output "%s"\n',[outloc "Su150kHzLatZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "Su150kHzLatZ.dat"])
fprintf(pltid,'set title" 150 kHz dv/dz"\n')
fprintf(pltid,'set output "%s"\n',[outloc "Sv150kHzLatZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "Sv150kHzLatZ.dat"])
fprintf(pltid,'set title" 75 kHz du/dz"\n')
fprintf(pltid,'set output "%s"\n',[outloc "Su75kHzLatZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "Su75kHzLatZ.dat"])
fprintf(pltid,'set title" 75 kHz dv/dz"\n')
fprintf(pltid,'set output "%s"\n',[outloc "Sv75kHzLatZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "Sv75kHzLatZ.dat"])
#LonZ plots
fprintf(pltid,'set xlabel "Longitude (deg)"\n')
fprintf(pltid,'set title" 150 kHz du/dz"\n')
fprintf(pltid,'set output "%s"\n',[outloc "Su150kHzLonZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "Su150kHzLonZ.dat"])
fprintf(pltid,'set title" 150 kHz dv/dz"\n')
fprintf(pltid,'set output "%s"\n',[outloc "Sv150kHzLonZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "Sv150kHzLonZ.dat"])
fprintf(pltid,'set title" 75 kHz du/dz"\n')
fprintf(pltid,'set output "%s"\n',[outloc "Su75kHzLonZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "Su75kHzLonZ.dat"])
fprintf(pltid,'set title" 75 kHz dv/dz"\n')
fprintf(pltid,'set output "%s"\n',[outloc "Sv75kHzLonZ" termtext])
fprintf(pltid,'splot "%s" binary matrix \n',[outloc "Sv75kHzLonZ.dat"])
# Set common positive definite Shear scale
# plot S150 and S75
# close plt file
fclose(pltid);
unix(["gnuplot " outloc "ADCP.plt"]);