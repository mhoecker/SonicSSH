#Isopycnal Plots
load /Volumes/scienceparty_share/CTDdata/EddyMix_CTDstations.mat 
Z = -512:2:0;
Sigmasmooth = zeros(length(LAT),length(Z));
Tempsmooth = zeros(length(LAT),length(Z));
potTempsmooth = zeros(length(LAT),length(Z));
Saltsmooth = zeros(length(LAT),length(Z));
for i=1:length(LAT)
	Sigmasmooth(i,:) = csqfil(SIGMA(:,i),-DEPTH(:,i),Z);
	Tempsmooth(i,:) = csqfil(TEMP(:,i),-DEPTH(:,i),Z);
	potTempsmooth(i,:) = csqfil(potTEMP(:,i),-DEPTH(:,i),Z);
	Saltsmooth(i,:) = csqfil(SALT(:,i),-DEPTH(:,i),Z);
#	Sigmasmooth(i,:) = interp1(SIGMA(:,i),-DEPTH(:,i),Z);
#	Tempsmooth(i,:) = interp1(TEMP(:,i),-DEPTH(:,i),Z);
#	potTempsmooth(i,:) = interp1(potTEMP(:,i),-DEPTH(:,i),Z);
#	Saltsmooth(i,:) = interp1(SALT(:,i),-DEPTH(:,i),Z);
end%for
T = nanmean(TIME);
save("-v7","SmoothCTD.mat","Sigmasmooth","Tempsmooth","Saltsmooth","T","Z","LAT","LON");
binmatrix(T-ones(size(T)).*(datenum([2013,1,1])+1),Z',Sigmasmooth',"SmoothCTDvT.dat");
binmatrix(LAT',Z',Sigmasmooth',"SmoothCTDvLat.dat");
binmatrix(LON',Z',Sigmasmooth',"SmoothCTDvLon.dat");