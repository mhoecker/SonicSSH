function themean = nanmean(hasnan,dim)
 if(nargin<2)
  dim=1;
 end%if
nonan = hasnan;
Ntot = size(hasnan)(dim);
Nnan = sum(isnan(hasnan),dim);
nanidx = find(isnan(hasnan));
nonan(nanidx) = 0;
 if(Ntot-Nnan>0)
themean = sum(nonan,dim)./(Ntot-Nnan);
 else
  themean = NaN;
end%function
