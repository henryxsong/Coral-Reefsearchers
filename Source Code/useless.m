syms x y z
% v = -5:0.01:5;  % plotting range from -5 to 5
% [x y] = meshgrid(v);  % get 2-D mesh for x and y
% cond1 = x+y+x.^2 < 3;  % check conditions for these values
% cond2 = y+x+y.^2 < 3;
% cond1 = double(cond1);  % convert to double for plotting
% cond2 = double(cond2);
% cond1(cond1 == 0) = NaN;  % set the 0s to NaN so they are not plotted
% cond2(cond2 == 0) = NaN;
% cond = cond1.*cond2;  % multiply the two condaces to keep only the common points
% surf(x,y,cond)
% view(0,90) 

ineqplot('x.^2+y.^2-z.^3<10',[-5 5], 'r');
%h = ineqplot('y<x+3',[0 10 -5 5]);
set(h,'MarkerFaceColor','r'); % Change color