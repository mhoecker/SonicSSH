function bad = nanstencil(data,rw,cw)
% Which points are poluted by NaNs
% using a stencil of size rows by cols
% assumed the stencil
% *is not periodic
% *is rectagular
% *allways uses the same number of points
 sizedat = size(data);
 bad = zeros(sizedat);
 for i=1:sizedat(1)
  for j=1:sizedat(2)
   if(isnan(data(i,j)))
    bad(i,j)=NaN;
    if(i<=ceil(rw/2))
     istart = 0;
    elseif(i>=sizedat(1)-floor(rw/2))
      istart = sizedat(1)-rw;
    else
      istart = i-ceil(rw/2);
    end%if
    if(j<=ceil(cw/2))
      jstart = 0;
    elseif(j>=sizedat(2)-floor(cw/2))
      jstart = sizedat(2)-cw;
    else
      jstart = j-ceil(cw/2);
    end%if
    for k=1:rw
     for l=1:cw
      bad(istart+k,jstart+l)=NaN;
     end%for
    end%for
   end%if
  end%for
 end%for
end%function
