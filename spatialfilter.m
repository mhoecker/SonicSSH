function fil = spatialfilter(Lain,Loin,val,Laout,Loout,R)
# Filter using a gaussian filter of width R
# uses the Gibbs sea water package
 addpath /Users/mhoecker/Documents/TEOS-10
 addpath /Users/mhoecker/Documents/TEOS-10/library
#
 for i=1:length(Laout)
  dist = zeros(size(val));
  w = dist;
  for j=1:length(Lain)
   dist(j) = gsw_distance([Loout(i),Loin(j)],[Laout(i),Lain(j)])
  end%for
  w = exp(-(dist.^2)/(2*R^2));
  if(sum(w)>0)
   fil(i) = sum(w.*val)./sum(w);
  else
   fil(i) = NaN;
  end%if
 end%for
end%function