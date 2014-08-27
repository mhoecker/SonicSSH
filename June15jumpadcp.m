# Plot the Jump on June 14 2013 ~20:30 UTC
# Also plot the Lat/Lon
#
adcpfile = "/Volumes/current_cruise/hdss/Matfiles/50k/50k_2013-06-15_60s.mat";
Fluxdir = "/Volumes/scienceparty_share/FluxData/";
load([ Fluxdir "130615.mat"]);
load(adcpfile);
jmpidx = 1:length(T);
Thours = T(jmpidx)/3600;
xrange = [min(Thours),max(Thours)];
crange = [-0.75,.75];
zrange = [-800,-50];
subplot(3,1,2);
plot(Thours,MET(jmpidx,23),";SST_1;",Thours,MET(jmpidx,31),";SST_2;");
title("SST Jump on June 15 2013 UTC")
ylabel("T (C)")
axis(xrange)
#subplot(6,1,2);
#Calculate the net heat flux using both wind speed measurements
#Net1 = MET(jmpidx,13)+MET(jmpidx,12) + Flux1(jmpidx,3) + Flux1(jmpidx,4);
#Net2 = MET(jmpidx,13)+MET(jmpidx,12) + Flux2(jmpidx,3) + Flux2(jmpidx,4);
#plot(Thours,Net1,";Flux_1;",Thours,Net2,";Flux_2;");
#axis([xrange,-100,1650])
#ylabel("Net Heat Flux (W/m^2)")
#subplot(6,1,3); 
#plot(Thours,MET(jmpidx,54));axis(xrange);ylabel("Latitude (deg)")
#subplot(6,1,4);
#plot(Thours,MET(jmpidx,55)); axis(xrange); ylabel("Longitude (deg)"); xlabel("Hours UTC")
subplot(3,1,1)
[UTChour,z] = meshgrid(24*(sonar.datenum-datenum([2013,06,15])),-sonar.depths);
pcolor(UTChour,z,real(sonar.U));shading flat;
axis([xrange,zrange])
caxis(crange)
subplot(3,1,3)
pcolor(UTChour	,z,imag(sonar.U));shading flat;
axis([xrange,zrange])
caxis(crange)
print([ Fluxdir "Plots/June15jumpadcp.png"],"-dpng")