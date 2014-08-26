module sonicsensor(origin=[0,0,0],orientation=[0,0,0],width=.0625){
	facets = 64;
	translate(v=origin){
		rotate(a=orientation){
		 union(){
		  translate([2*width,0,-1.25*width]){
		   cylinder(r=.125*width,h=2.5*width,$fn=facets);
		  }
		  translate([0,0,.6*width]){
					rotate([0,90,0]){
						cylinder(r=.5*width,h=4*width,$fn=facets);
						cylinder(r=.6*width,h=width,$fn=facets);
						translate([0,0,3*width]){
 						cylinder(r=.6*width,h=width,$fn=facets);
						}
					}
					sphere(r=.6*width,$fn=4*facets);
					cylinder(r=.6*width,h=width,$fn=facets);
					cylinder(r=.5*width,h=2*width,$fn=facets);
					translate([0,0,2*width]){
					 difference(){
							cylinder(r1=.5*width,r2=width,h=2*width,$fn=facets);
							translate([0,0,.5*width]){
								cylinder(r1=.5*width,r2=width,h=2*width,$fn=facets);
							}
						}
					}
				}
			}
		}
	}
}
module trapezoid(width=.0625,base=2,side=1,angle=60.0,orientation=[0,0,0],origin=[0,0,0],sensors=[0,1,2,3,4]){
	cang = cos(angle);
	sang = sin(angle);
	thang = tan(angle/2);
	tang = tan(angle);
	yb = -thang*base/2;
	yt = sang*side+yb;
	xb = base/2;
	xt = xb-cang*side;
	dws = sang*width;
	dwc = cang*width;
	dlb = width/thang;
	dlt = width/tang;
	bot = [
	[-xb+dlt,yb],
	[-xb+dlt,yb+width],
	[+xb-dlt,yb+width],
	[+xb-dlt,yb]
	];
	top = [
	[-xt,yt],
	[-xt,yt-width],
	[+xt,yt-width],
	[+xt,yt]
	];
	lsi = [
	[-xb,yb+dlt],
	[-xb+dws,yb],
	[-xt+dws,yt-dlt],
	[-xt,yt]
	];
	rsi = [
	[+xb,yb+dlt],
	[+xb-dws,yb],
	[+xt-dws,yt-dlt],
	[+xt,yt]
	];
	mid = [
	[0-width/2,yb],
	[0+width/2,yb],
	[0+width/2,yt],
	[0-width/2,yt],
	];
	translate(v=origin){
		rotate(orientation){
			union(){
				linear_extrude(height=width){
					polygon(points=bot,convexity=2);
				}
				linear_extrude(height=width){
					polygon(points=top,convexity=2);
				}
				translate([0,0,-width]){
					linear_extrude(height=width){
						polygon(points=lsi,convexity=2);
					}
					linear_extrude(height=width){
						polygon(points=rsi,convexity=2);
					}
					linear_extrude(height=width){
						polygon(points=mid,convexity=2);
					}
				}
				color("Gray",a=.5){
				 for(i=sensors){
				  if(i==0){
	 				 sonicsensor(origin=[xb,yb+width/2,width],orientation=[0,0,180],width=width);
	 				}
				  if(i==1){
 				  sonicsensor(origin=[0,yb+width/2,width],orientation=[0,0,0],width=width);
	 				}
				  if(i==2){
	 				 sonicsensor(origin=[-xb,yb+width/2,width],orientation=[0,0,00],width=width);
	 				}
				  if(i==3){
 				 	sonicsensor(origin=[-xt,yt-width/2,width],orientation=[0,0,0],width=width);
	 				}
				  if(i==4){
 				 	sonicsensor(origin=[xt,yt-width/2,width],orientation=[0,0,180],width=width);
 				 }
 				}
				}
			}
		}
	}
}
c60 = cos(60);
s60 = sin(60);
width = 1;
base = 200/3;
side = 100/3;
union(){
 trapezoid(width=width,base=base,side=side,angle=60.0,orientation=[0,0,-120],origin=[.5*side,s60*side,0],sensors=[0,1,2]);
 trapezoid(width=width,base=base,side=side,angle=60.0,orientation=[0,0,0],origin=[0,0,0],sensors=[0,1,2,4]);
 trapezoid(width=width,base=base,side=side,angle=60.0,orientation=[0,0,120],origin=[1*side,0,0],sensors=[0,1,2]);
}
