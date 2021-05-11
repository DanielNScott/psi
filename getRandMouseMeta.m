function mouse = getRandMouseMeta()

   % Fixed parameters
   mouse.fixed.lapse.min = 0.1;
   mouse.fixed.lapse.max = 0.1;
   mouse.fixed.lapse.tau = 2;

   mouse.fixed.guess.min = 0.1;
   mouse.fixed.guess.max = 0.1;
   mouse.fixed.guess.tau = 2;
   
   % Parameters govorning by-day learning
   % Threshold
   ur = 0.05 + sort(rand(4,1)*(0.5-0.05),'descend');
   mouse.meta.thresh.maxbeg = ur(1);
   mouse.meta.thresh.maxfin = ur(2);
   
   mouse.meta.thresh.minbeg = ur(3)-0.02;
   mouse.meta.thresh.minfin = ur(4)-0.02;
   
   ur = rand(2,1);
   mouse.meta.thresh.taubeg = ur(2)*10 + ur(1)*10;
   mouse.meta.thresh.taufin = ur(2)*10;
   
   mouse.meta.thresh.tauday = 2 + rand()*4;
   
   % Slope parameters
   mouse.meta.slope.maxbeg = 3;
   mouse.meta.slope.maxfin = 6;
   
   mouse.meta.slope.minbeg = 1;
   mouse.meta.slope.minfin = 5;
   
   ur = rand(2,1);
   mouse.meta.slope.taubeg = ur(2)*5 + ur(1)*10;
   mouse.meta.slope.taufin = ur(2)*5;
   
   mouse.meta.slope.tauday = rand()*4;

end