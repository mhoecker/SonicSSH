function [Jasondata] = ingestJasonDir(dirloc,latrange,lonrange)
Lat = [];
Lon = [];
SSHa = [];
# Loop over ssh files
	Jasondirectory = dir(dirloc);
	for i=1:length(Jasondirectory)
#	for i=1:min([800,length(Jasondirectory)])
		thefile = Jasondirectory(i).name;
		if(strfind(thefile,".nc")==length(thefile)-2)
			thefile = [dirloc Jasondirectory(i).name];
			nc = netcdf(thefile,'r');
			latI = nc{'lat'}(:)*1e-6; 
			lonI = nc{'lon'}(:)*1e-6;
			sshI = nc{'ssha'}(:);
# Find points in region of interest
			locidx = 	find((latI>min(latrange)).*(latI<max(latrange)).*(lonI>min(lonrange)).*(lonI<max(lonrange)).*(sshI!=32767));
			Lat = [Lat,latI(locidx)']; 
			Lon = [Lon,lonI(locidx)'];
			SSHa = [SSHa,sshI(locidx)'];
			clear latI lonI sshI
			ncclose(nc)
		end%if
	end%for
Jasondata = [Lat;Lon;SSHa];
end%function