function [termtxt,termsfx] = termselect(termtype)
# function [temtype,termsfx] = termselect(termtype)
# choose a consistent terminal type for gnuplot
# the options are
#
#		'epsposter' 12"x12" 48pt text
#		'epsarticle' 9"x6" 12pt text
#		'epsarticlesmall' 4.5"x3" 10pt font
#		'pngposter' is intended to be blown up to 12"x12"
#		'canvas' cerates an interactive web page
#
#		if the term type is unrecognized the funtion returns values for a png
#
if(nargin<1)
	termtype = "";
endif
	switch(termtype)
		case 'epsposter'
			termtxt = 'postscript eps enhanced color size 12in,12in font "Helvetica,48" blacktext linewidth 2';
			termsfx = '.eps';
		case 'pdfposter'
			termtxt = 'pdf enhanced color size 12in,12in font "Helvetica,48" linewidth 2';
			termsfx = '.pdf';
		case 'pngposter'
			# if blown up to 12inx12in creates a 300dpi image
			termtxt = 'png enhanced truecolor size 3600,3600 nocrop font "Helvetica,96" linewidth 6';
			termsfx = '.png';
		case 'epsarticle'
			termtxt = 'postscript eps enhanced color size 9in,6in font "Helvetica,12" blacktext linewidth 2';
			termsfx = '.eps';
		case 'pdfarticle'
			termtxt = 'pdf enhanced color size 9in,6in font "Helvetica,12" linewidth 2';
			termsfx = '.pdf';
		case 'epsarticlesmall'
			termtxt = 'postscript eps enhanced color size 4.5in,3in font "Helvetica,10" blacktext linewidth 2';
			termsfx = '.eps';
		case 'pdfarticlesmall'
			termtxt = 'ppdf enhanced color size 4.5in,3in font "Helvetica,10" linewidth 2';
			termsfx = '.eps';
		case 'canvas'
			termtxt = ' canvas size 1536,1024 jsdir "js" mousing enhanced linewidth 1';
			termsfx = '.html';
		otherwise
			termtxt = "png enhanced size 1536,1024 truecolor nocrop linewidth 2";
			termsfx = '.png';
	end
end
