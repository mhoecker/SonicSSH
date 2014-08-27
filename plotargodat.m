function plotargodat(dirloc)
	addpath /Users/mhoecker/Documents/TEOS-10
	addpath /Users/mhoecker/Documents/TEOS-10/library
	load([dirloc "ingestedArgo.mat"])
	T = Argodata{1};
	S = Argodata{2};
	P = Argodata{3};
	Lat = Argodata{4};
	Lon = Argodata{5};
	Jday = Argodata{6};
	clear Argodata;
	Z = gsw_z_from_p(P,mean(Lat));
	load("/Volumes/scienceparty_share/CTDdata/EddyMix_CTDstations.mat","TEMP","LAT","LON","SALT","DEPTH");
	Latctd = LAT;
	Lonctd = LON;
	Tctd = [];
	Sctd = [];
	for i=1:length(TEMP(1,:))
		Tctd = [Tctd;csqfil(TEMP(:,i),-DEPTH(:,i),Z)];
		Sctd = [Sctd;csqfil(SALT(:,i),-DEPTH(:,i),Z)];
	end%for
	clear TEMP LAT LON SALT DEPTH
	Tmean = nanmean(T);
	Smean = nanmean(S);
	save("-v7",[dirloc "CTD_Argo.mat"],"Tmean","Smean","Z","Sctd","Tctd","Latctd","Lonctd")
	[Latp,Plat] = meshgrid(Lat,P);
	[Lonp,Plon] = meshgrid(Lon,P);
	[Jdayp,Pjday] = meshgrid(Jday,P);
	figure(1)
	subplot(1,2,1)
	plot(Tmean,-P)
	xlabel("Temperatire C")
	ylabel("Pressure dBar")
	subplot(1,2,2)
	plot(Smean,-P)
	xlabel("Salinity PSU")
	print([dirloc "Mean_T_S.png"],"-dpng")
	figure(2)
	subplot(1,2,1)
	pcolor(Latp,-Plat,T'); shading flat
	title("Temperature")
	xlabel("Longitude deg")
	ylabel("Pressure dBar")
	subplot(1,2,2)
	pcolor(Latp,-Plat,S'); shading flat
	title("Salinity")
	xlabel("Latitude deg")
	figure(3)
	subplot(1,1,1)
	plot(Lon,Lat,"o")
	ylabel("Latitude deg")
	xlabel("Longitude deg")
	title("Argo Float Locations")
	print([dirloc "Argo_Locations.png"],"-dpng")
	figure(4)
	subplot(2,1,1)
	plot(Latctd(:)*ones(size(Z))+(Tctd-(Tmean(:)*ones(1,6))')/5,Z)
	axis([19,24])
	ylabel("depth m")
	xlabel("Temperature (C)")
	title("Anomaly")
	subplot(2,1,2)
	plot(Latctd(:)*ones(size(Z))+(Sctd-(Smean(:)*ones(1,6))')*2,Z)
	xlabel("Salinity (PSU)")
	axis([19,24])
	print([dirloc "CTD-Argo_Lat.png"],"-dpng")
	figure(5)
	subplot(2,1,1)
	plot(Lonctd(:)*ones(size(Z))+(Tctd-(Tmean(:)*ones(1,6))')/5,Z)
	#axis([19,24])
	ylabel("depth m")
	xlabel("Temperature (C)")
	title("Anomaly")
	subplot(2,1,2)
	plot(Lonctd(:)*ones(size(Z))+(Sctd-(Smean(:)*ones(1,6))')*2,Z)
	xlabel("Salinity (PSU)")
	#axis([19,24])
	print([dirloc "CTD-Argo_Lon.png"],"-dpng")
end%function