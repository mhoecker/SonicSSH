#feed kpp model
CTDfile = ;
ADCPfile = ;
Fluxfile = ;
# Parameters
Cp = 4000; #J/K kg Heat Capacity of water (fix me)
F0 = 4*pi/(24*3600); #Corriolis parameter at North Pole
rho0 = 1000;
hbl = 0;
r1 = .6854;
amu1 = 1./1.116;
r2 = 1-r1;
amu2 = 1/22.59;
imod = 0;
# Get T-S Profile
#T = ;
#S = ;
#ZTS = ;
#TStime = ;
#Cor = F0*sin(Lat);
# Get Flux data close to the profile
#fluxidx  = find((fluxtime-TStime)<dtime)
#WUsurf = mean(-tau(fluxidx).*sin(winder(fluxidx)));
#WVsurf = mean(-tau(fluxidx).*cos(winder(fluxidx)));
#WTsurf = (Sensible+Latent+Longwave)/Cp;
#WSsurf = (Latent/L)*SurfaceS;
# Get ADCP data close to the profile
#adcpidx = find((adcptime-TStime)<dtime)
# Fill velocity to surface
# interpolate onto a grid
# invoke KPP
#kpp(U,V,T,S,ZZ,ZP,WUSURF,WVSURF,WTSURF,WSSURF,SWRD,COR,hbl,r1,amu1,r2,amu2,imod)