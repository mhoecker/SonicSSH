function [x,y] = sensorarray(a,b,iter)
 x = [];
 y= [];
 for i=1:iter
  x = [x,a(1)*(-.5)*i+b(1)*(-.5)*i];
  x = [x,a(1)*(+.5)*i+b(1)*(-.5)*i];
  x = [x,a(1)*(+.5)*i+b(1)*(+.5)*i];
  y = [y,a(2)*(i)+b(2)*(i)];
  y = [y,a(2)*(i-1)+b(2)*(i)];
  y = [y,a(2)*(i)+b(2)*(i-1)];
 end%for
end%function
