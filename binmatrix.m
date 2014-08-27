function binmatrix(x,y,zz,filename)
% write a binary matrix readble by gnuplot
	Nx = length(x);
	Ny = length(y);
	Sx = size(x);
	Sy = size(y);
	Sz = size(zz);
	if((length(Sx)<3)*(length(Sy)<3)*(length(Sz)<3))
		fid = fopen(filename,"w");
		plottable = zeros(Ny+1,Nx+1);
		plottable(1,:) = [Nx,x];
		plottable(2:Ny+1,1) = y';
		plottable(2:Ny+1,2:Nx+1) = zz;
		plottable = plottable';
		fwrite(fid,plottable,'float');
		fclose(fid);
	else
		printf("binmatrix(x,y,zz,filename) was given a mishapen argument\n")
		printf("Function assumes x is a [1,Nx] matrix\n")
		printf("Function assumes y is a [1,Ny] matrix\n")
		printf("Function assumes zz is a [Nx,Ny] matrix\n")
	endif
endfunction
