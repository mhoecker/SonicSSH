function [Argodata] = ingestArgoDir(dirloc,P)
	T = [];
	S = [];
	Lat = [];
	Lon = [];
	Jday = [];
	srcfile = {};
	if(nargin<2)
		P = 0:10:1000;
	end%if
	Argodirectory = dir(dirloc);
	for i=1:length(Argodirectory)
#	for i=1:min([length(Argodirectory),500])
		thefile = Argodirectory(i).name;
		if(strfind(thefile,".nc")==length(thefile)-2)
			thefile = [dirloc Argodirectory(i).name];
			nc = netcdf(thefile,'r');
			dy = abs(nc{'LATITUDE'}(:)-20.4);
			dx = abs(nc{'LONGITUDE'}(:)-129.5);
			drsq = dy^2+(dx*cos(pi*nc{'LATITUDE'}(:)/180))^2;
			if(drsq<9)
 			srcfile = {srcfile{:},thefile};
				try
					T = [T;csqfil(nc{'TEMP_ADJUSTED'}(:),nc{'PRES_ADJUSTED'}(:),P)];
				catch
					T = [T;csqfil(nc{'TEMP'}(:),nc{'PRES'}(:),P)];
				end%try_catch
				try
					S = [S;csqfil(nc{'PSAL_ADJUSTED'}(:),nc{'PRES_ADJUSTED'}(:),P)];
				catch
					S = [S;csqfil(nc{'PSAL'}(:),nc{'PRES'}(:),P)];
				end%try_catch
				Lat = [Lat;nc{'LATITUDE'}(:)];
				Lon = [Lon;nc{'LONGITUDE'}(:)];
					Jday = [Jday;nc{'JULD'}(:)];
			end%if
			ncclose(nc)
		end%if
	end%for
	T(find(T>100))= NaN;
	T(find(T<-4))= NaN;
	S(find(S<0))= NaN;
	S(find(S>100))= NaN;
	Argodata = {T,S,P,Lat,Lon,Jday,srcfile};
	save("-v7",[dirloc "ingestedArgo.mat"],"Argodata")
end%function