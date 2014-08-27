function SST_HeatFlux_Lat_Lon_plot(yymmdd)
# Plot the Jump on year = 20yy month mm day dd UTC
# Also plot the Lat/Lon
#
Fluxdir = "/Volumes/scienceparty_share/FluxData/";
load([ Fluxdir num2str(yymmdd) ".mat"])
jmpidx = 1:length(T);
Thours = T(jmpidx)/3600;
xrange = [min(Thours),max(Thours)]
subplot(4,1,1);
plot(Thours,MET(jmpidx,23),";SST_1;",Thours,MET(jmpidx,31),";SST_2;");
title(["SST Jump on " num2str(yymmdd) " UTC"])
ylabel("T (C)")
axis(xrange)
subplot(4,1,2);
#Calculate the net heat flux using both wind speed measurements
Net1 = MET(jmpidx,13) + Flux1(jmpidx,3) + Flux1(jmpidx,4) - Flux1(jmpidx,24);
Net2 = MET(jmpidx,13) + Flux2(jmpidx,3) + Flux2(jmpidx,4) - Flux2(jmpidx,24);
plot(Thours,Net1,";Flux_1;",Thours,Net2,";Flux_2;");
axis([xrange,-200,1600])
plot(Thours,Net1,";Flux_1;",Thours,Net2,";Flux_2;");
axis([xrange,-200,1600])
ylabel("Net Heat Flux (W/m^2)")
subplot(4,1,3);
plot(Thours,MET(jmpidx,54));
axis(xrange)
ylabel("Latitude (deg)")
subplot(4,1,4);
plot(Thours,MET(jmpidx,55))
axis(xrange)
ylabel("Longitude (deg)")
xlabel("Hours UTC")

print([ Fluxdir "Plots/SST_Heat_LatLon" num2str(yymmdd) ".png"],"-dpng")