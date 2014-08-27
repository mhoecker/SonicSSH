module fakechanel(origin=[0,0,0],orientation=[0,0,0],length=96){
 spacing=48;
 r=7;
 holes = ceil(length/spacing);
 w = 41;
 h = 41;
 t = 2.7;
 dw = 10;
 dh = 7.1;
 profile = [
  [-w/2+dw,h-dh],
  [-w/2+dw,h],
  [-w/2,h],
  [-w/2,0],
  [+w/2,0],
  [+w/2,h],
  [+w/2-dw,h],
  [+w/2-dw,h-dh],
  [+w/2-dw+t,h-dh],
  [+w/2-dw+t,h-t],
  [+w/2-t,h-t],
  [+w/2-t,+t],
  [-w/2+t,+t],
  [-w/2+t,h-t],
  [-w/2+dw-t,h-t],
  [-w/2+dw-t,h-dh]
 ];
 translate(origin){
  rotate(orientation){
   difference(){
    rotate([90,0,90]){
     linear_extrude(height=length){
      polygon(points=profile,convexity=8);
     }
    }
    union(){
     for(i=[1:holes]){
      translate([spacing*(i-.5),0,-t]){
       cylinder(r=r,h=3*t,$fn=16);
      }
     }
    }
   }
  }
 }
}
module sonicsensor(origin=[0,0,0],orientation=[0,0,0]){
	facets = 8;
	width = 20;
	r = 7;
	translate(origin){
		rotate(orientation){
		 union(){
		  translate([72,0,-1*width]){
		   cylinder(r=r,h=2.5*width,$fn=facets);
		  }
		  translate([0,0,.6*width]){
					rotate([0,90,0]){
						cylinder(r=.5*width,h=120,$fn=facets);
						cylinder(r=.6*width,h=width,$fn=facets);
						translate([0,0,120]){
 						cylinder(r=.6*width,h=width,$fn=facets);
						}
					}
					sphere(r=.6*width,$fn=facets);
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
module trapezoid(side=480,orientation=[0,0,0],origin=[0,0,0],sensors=[0,1,2,3,4]){
 xb = side;
 yb = -side*tan(30);
 xt = side/2;
 yt = side*sin(60)+yb;
	translate(v=origin){
		rotate(orientation){
			union(){
			 fakechanel(origin=[-xb+48,yb,-41],orientation=[0,0,0],length=2*side-96);
			 fakechanel(origin=[-xt+48,yt,-41],orientation=[0,0,0],length=side-72);
			 fakechanel(origin=[-120,yb+21,-41],orientation=[0,0,90],length=side*sin(60)-42);
			 fakechanel(origin=[xb-side/2,yb,-82],orientation=[0,0,60],length=side/2);
			 fakechanel(origin=[xb-side/3,yb,-82],orientation=[0,0,60],length=side/3);
				color("Gray",a=.5){
				 for(i=sensors){
				  if(i==0){
	 				 sonicsensor(origin=[xb,yb,0],orientation=[0,0,180],width=width);
	 				}
				  if(i==1){
 				  sonicsensor(origin=[0,yb,0],orientation=[0,0,0],width=width);
	 				}
				  if(i==2){
	 				 sonicsensor(origin=[-xb,yb,0],orientation=[0,0,00],width=width);
	 				}
				  if(i==3){
 				 	sonicsensor(origin=[-xt,yt,0],orientation=[0,0,0],width=width);
	 				}
				  if(i==4){
 				 	sonicsensor(origin=[xt,yt,0],orientation=[0,0,180],width=width);
 				 }
 				}
				}
			}
		}
	}
}
c60 = cos(60);
s60 = sin(60);
side = 48*3*8;
union(){
 //trapezoid();
 trapezoid(side=side,orientation=[0,0,-120],origin=[.5*side,s60*side,0],sensors=[0,1,2,3]);
 trapezoid(side=side,orientation=[0,0,0],origin=[0,0,0],sensors=[1,2,4]);
 trapezoid(side=side,orientation=[0,0,120],origin=[1*side,0,0],sensors=[1,2,3]);
}
