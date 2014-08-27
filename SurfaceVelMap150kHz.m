	addpath /Users/mhoecker/Documents/SonicSSH/scripts
while(1)
# Plot upper 170 m (minus 50m surface layer) mean currents from ADCP
for i=1:12
middepth = 50*i-25;
depthrange = 25;
depthbound = middepth+[-depthrange,depthrange];
latrange = [18.4,23.0];
lonrange = [127.0,132.2];
ar = 1/cos(mean(latrange)*pi/180);;
velscale = .4;
#ar = 1/cos(80*pi/180);
#
#
# Load ADCP data
if(middepth<175)
nc = netcdf("/Volumes/current_cruise/adcp/RR1308/proc/nb150/contour/nb150.nc","r"); adcpfq = "150kHz";
fileroot = ["/Volumes/scienceparty_share/ADCP/AvgVel150kHz" num2str(middepth-depthrange) "to" num2str(middepth+depthrange) "m"];
else
nc = netcdf("/Volumes/current_cruise/adcp/RR1308/proc/os75nb/contour/os75nb.nc","r"); adcpfq = "75kHz";
fileroot = ["/Volumes/scienceparty_share/ADCP/AvgVel75kHz" num2str(middepth-depthrange) "to" num2str(middepth+depthrange) "m"];
end#if
u = nc{'u'}(:);
idxubad = find(abs(u)>10);
u(idxubad) = NaN;
v = nc{'v'}(:);
idxvbad = find(abs(v)>10);
v(idxvbad) = NaN;
lat = nc{'lat'}(:);
lon = nc{'lon'}(:);
depth = nc{'depth'}(:);
t = nc{'time'}(:);
ncclose(nc);
nc = netcdf("/Volumes/current_cruise/ssha_files/ssha_global_wrt_mean_2013171.nc",'r');
sshdate = "Jun 20 2013";
sshlat = nc{'lat'}(:);
idxlat = find((sshlat>min(latrange)-1).*(sshlat<max(latrange)+1));
sshlat = sshlat(idxlat);
sshlon = nc{'lon'}(:);
idxlon = find((sshlon>min(lonrange)-1).*(sshlon<max(lonrange)+1));
sshlon = sshlon(idxlon);
ssh = nc{'ssh'}(idxlat,idxlon);
ncclose(nc)
sshdat =[fileroot "ssh.dat"];
binmatrix(sshlon',sshlat,ssh,sshdat);
sshpcntr = [fileroot "sshp.cntr"];
ssh0cntr = [fileroot "ssh0.cntr"];
sshmcntr = [fileroot "sshm.cntr"];
uavg = zeros(size(t));
vavg = zeros(size(t));
for i=1:length(t)
	binrange = find(depthrange<abs(depth(i,:)-middepth));
	uavg(i) = nanmean(u(i,binrange)');
	vavg(i) = nanmean(v(i,binrange)');
end%for
ts = min(t):1/48:max(t);
#[vor,div,vbyf] = strainrateest(uavg',vavg',lon',lat',55500,555,400);
					  vor = zeros(size(uavg))';
					  div = zeros(size(uavg))';
					  vbyf = zeros(size(uavg))';
matfile = [fileroot ".mat"];
save("-v7",matfile,"lat","lon","uavg","vavg","t","vor","div","vbyf")
#
binfile = [fileroot ".dat"];
form = binarray(t',[lat';lon';uavg';vavg';vor;div;vbyf],binfile);
# plot in gnuplot using U, V, heading, and speed as pc
[termtxt,termsfx] = termselect('pdfarticle');
pltfile = [fileroot ".plt"];
pltidx = fopen(pltfile,'w');
fprintf(pltidx,'set contour base\n')
fprintf(pltidx,'unset surface\n')
fprintf(pltidx,'set table "%s"\n',sshpcntr)
fprintf(pltidx,'set cntrparam levels incremental 4,4,28 \n')
fprintf(pltidx,'splot "%s" binary matrix\n',sshdat)	
fprintf(pltidx,'set table "%s"\n',sshmcntr)
fprintf(pltidx,'set cntrparam levels incremental -28,4,-4 \n')
fprintf(pltidx,'splot "%s" binary matrix\n',sshdat)
fprintf(pltidx,'set table "%s"\n',ssh0cntr)
fprintf(pltidx,'set cntrparam levels discrete 0 \n')
fprintf(pltidx,'splot "%s" binary matrix\n',sshdat)	
fprintf(pltidx,'unset table\n')
fprintf(pltidx,'reset\n')
fprintf(pltidx,'set term %s\n',termtxt)
fprintf(pltidx,'set output "%s%s"\n',fileroot,termsfx)
fprintf(pltidx,'set title "SSH from %s, Avg Velocity %3.0fm-%3.0fm %s ADCP"\n',sshdate,min(depthbound),max(depthbound),adcpfq)
fprintf(pltidx,'set grid\n')
fprintf(pltidx,'set xrange [%f:%f]\n',min(lonrange),max(lonrange))
fprintf(pltidx,'set xtics %f,1\n',floor(min(lonrange)))
fprintf(pltidx,'set mxtics \n')
fprintf(pltidx,'set xlabel "Longitude"\n')
fprintf(pltidx,'set cbtics 0,30\n')
fprintf(pltidx,'set cblabel "Heading"\n')
fprintf(pltidx,'set yrange [%f:%f]\n',min(latrange),max(latrange))
fprintf(pltidx,'set ytics %f,1\n',floor(min(latrange)))
fprintf(pltidx,'set mytics \n')
fprintf(pltidx,'set ylabel "Latitude"\n')
fprintf(pltidx,'set palette mode HSV\n')
fprintf(pltidx,'set palette function gray,1,1\n')
fprintf(pltidx,'set cbrange [0:360]\n')
fprintf(pltidx,'unset colorbox\n')
fprintf(pltidx,'set multiplot\n')
fprintf(pltidx,'set size 1,1\n')
fprintf(pltidx,'set size ratio %f\n',-ar)
fprintf(pltidx,'set key\n')
fprintf(pltidx,'set arrow from %f,%f to %f,%f nohead lw 2 front\n',min(lonrange)+.5*velscale,max(latrange)-1.5*velscale,min(lonrange)+.5*velscale,max(latrange)-.5*velscale)
fprintf(pltidx,'set label "1 m/s" at %f,%f front\n',min(lonrange)+.5*velscale,max(latrange)-velscale)
fprintf(pltidx,'set lmargin 7\n')
fprintf(pltidx,'set rmargin 2\n')
fprintf(pltidx,'set tmargin 2\n')
fprintf(pltidx,'set bmargin 3\n')
fprintf(pltidx,'plot ')
fprintf(pltidx,' "%s" with lines lc -1 lw .5 notitle\\\n',sshpcntr)
fprintf(pltidx,', "%s" with lines lc "grey" lw .5 notitle\\\n',sshmcntr)
fprintf(pltidx,', "%s" with lines lt -1 lw 1 notitle\\\n',ssh0cntr)
fprintf(pltidx,', "%s" %s u 3:2:($4*%f):($5*%f):((360+5*floor(36*atan2($4,$5)/pi))%s) w vectors nohead lw .5 lc pal not \\\n',binfile,form,ar*velscale,velscale,'%360')
fprintf(pltidx,', "%s" %s u 3:2 w lines lw .5 lc rgbcolor "#404080" title "Ship Track"\n',binfile,form)
fprintf(pltidx,'unset key\n')
fprintf(pltidx,'set size .2,.2\n')
fprintf(pltidx,'set origin .6,.65\n')
fprintf(pltidx,'set parametric\n')
fprintf(pltidx,'set pm3d\n')
fprintf(pltidx,'set pm3d corners2color c1 \n')
fprintf(pltidx,'set view map\n')
fprintf(pltidx,'set size ratio -1\n')
fprintf(pltidx,'unset xlabel\n')
fprintf(pltidx,'unset title\n')
fprintf(pltidx,'unset ylabel\n')
fprintf(pltidx,'unset tics\n')
fprintf(pltidx,'unset surface\n')
fprintf(pltidx,'set border 0\n')
fprintf(pltidx,'set lmargin 0\n')
fprintf(pltidx,'set rmargin 0\n')
fprintf(pltidx,'set tmargin 0\n')
fprintf(pltidx,'set bmargin 0\n')
fprintf(pltidx,'set isosamples 180,180\n')
fprintf(pltidx,'set xrange [-1:1]\n')
fprintf(pltidx,'set yrange [-1:1]\n')
fprintf(pltidx,'set urange [.5:1]\n')
fprintf(pltidx,'set vrange [0:2*pi]\n')
fprintf(pltidx,'set label "E" at 1,0 center offset character 1,0 front\n')
fprintf(pltidx,'set label "N" at 0,1 center offset character 0,.5 front\n')
fprintf(pltidx,'set label "W" at -1,0 center offset character -1,0 front\n')
fprintf(pltidx,'set label "S" at 0,-1 center offset character 0,-.5 front\n')
#fprintf(pltidx,'set label "NE" at sqrt(2)/2,sqrt(2)/2 center offset character 1,0 front\n')
#fprintf(pltidx,'set label "NW" at -sqrt(2)/2,sqrt(2)/2 center offset character -1,0 front\n')
#fprintf(pltidx,'set label "SW" at -sqrt(2)/2,-sqrt(2)/2 center offset character -1,-1 front\n')
#fprintf(pltidx,'set label "SE" at sqrt(2)/2,-sqrt(2)/2 center offset character 1,-1 front\n')
fprintf(pltidx,'splot u*sin(v),u*cos(v),floor(v*180/pi)%s lc pal\n',"%360")
fprintf(pltidx,'unset multiplot\n')
#fprintf(pltidx,'set output "%svorticity%s"\n',fileroot,termsfx)
#fprintf(pltidx,'set title "SSH, Avg Vorticity %3.0fm-%3.0fm %s ADCP"\n',min(depthbound),max(depthbound),adcpfq)
#fprintf(pltidx,'set grid\n')
#fprintf(pltidx,'set xrange [%f:%f]\n',min(lonrange),max(lonrange))
#fprintf(pltidx,'set xtics %f,1\n',floor(min(lonrange)))
#fprintf(pltidx,'set mxtics \n')
#fprintf(pltidx,'set xlabel "Longitude"\n')
#fprintf(pltidx,'set cbtics auto\n')
#fprintf(pltidx,'set cblabel "Heading"\n')
#fprintf(pltidx,'set yrange [%f:%f]\n',min(latrange),max(latrange))
#fprintf(pltidx,'set ytics %f,1\n',floor(min(latrange)))
#fprintf(pltidx,'set mytics \n')
#fprintf(pltidx,'set ylabel "Latitude"\n')
#fprintf(pltidx,'set palette mode RGB\n')
#fprintf(pltidx,'set palette defined (0 "dark-red", .25 "red", .4 "yellow", .5 "light-grey" , .6 "cyan", .75 "blue",1 "magenta")\n');
#fprintf(pltidx,'set cbrange [-2:2]\n')
#fprintf(pltidx,'set colorbox\n')
#fprintf(pltidx,'set cblabel "Vorticity/f"\n')
#fprintf(pltidx,'set multiplot\n')
#fprintf(pltidx,'set origin 0,0\n')
#fprintf(pltidx,'set size 1,1\n')
#fprintf(pltidx,'set size ratio %f\n',-ar)
#fprintf(pltidx,'set key\n')
#fprintf(pltidx,'set border 31\n')
#fprintf(pltidx,'unset arrow\n')
#fprintf(pltidx,'unset label\n')
#fprintf(pltidx,'set lmargin 7\n')
#fprintf(pltidx,'set rmargin 2\n')
#fprintf(pltidx,'set tmargin 2\n')
#fprintf(pltidx,'set bmargin 3\n')
#fprintf(pltidx,'plot ')
#fprintf(pltidx,' "%s" with lines lc -1 lw 1 notitle\\\n',sshpcntr)
#fprintf(pltidx,', "%s" with lines lc "grey" lw 1 notitle\\\n',sshmcntr)
#fprintf(pltidx,', "%s" with lines lt -1 lw 2 notitle\\\n',ssh0cntr)
#fprintf(pltidx,', "%s" %s u 3:2:8 w lines lw 4 lc pal title "Vorticity/f"\\\n',binfile,form)
fclose(pltidx);
unix(["gnuplot " pltfile]);
				ncclose(nc);
end%for
RevelleADCPShear("/Volumes/scienceparty_share/ADCP/");
pause(900);
end%while