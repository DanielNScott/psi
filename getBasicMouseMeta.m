function mouse = getBasicMouseMeta()
%BASICMOUSE Summary of this function goes here
%   Detailed explanation goes here

   % Fixed parameters
   mouse.fixed.lapse.min = 0.05;
   mouse.fixed.lapse.max = 0.05;
   mouse.fixed.lapse.tau = 2;

   mouse.fixed.guess.min = 0.0;
   mouse.fixed.guess.max = 0.0;
   mouse.fixed.guess.tau = 2;
   
   % Parameters govorning by-day learning
   % Threshold
   mouse.meta.thresh.maxbeg = 0.2;
   mouse.meta.thresh.maxfin = 0.08;
   
   mouse.meta.thresh.minbeg = 0.09;
   mouse.meta.thresh.minfin = 0.06;
   
   mouse.meta.thresh.taubeg = 15;
   mouse.meta.thresh.taufin = 2;
   
   mouse.meta.thresh.tauday = 2;
   
   % Slope parameters
   mouse.meta.slope.maxbeg = 3;
   mouse.meta.slope.maxfin = 6;
   
   mouse.meta.slope.minbeg = 1;
   mouse.meta.slope.minfin = 5;
   
   mouse.meta.slope.taubeg = 15;
   mouse.meta.slope.taufin = 2;
   
   mouse.meta.slope.tauday = 2;

end

