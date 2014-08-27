function [form,S] = binarray(x,y,filename)
% write a binary array readble by gnuplot
% returns the format which needs to be given to gnuplot to read the file
	xy = [x;y];
	S = size(xy);
	if((length(S)==2)*(nargin()>1))
		fid = fopen(filename,"w");
		fwrite(fid,xy,'float');
		fclose(fid);
		form = [' binary format="'];
		for i=1:S(1)
			form = [form '%float'];
		end
		form = [form '" '];
	else
		printf("binarray(x,y,filename) was given a mishapen arguments\n")
		printf("Function assumes [x;y] is a [M,N] matrix\n")
	endif
endfunction
