function curve = getCurve(min, max, tau, nSteps, cond)

switch(cond)
   case 'down'
      curve = (max-min)*exp(-(0:nSteps-1)/tau) + min;
   case 'up'
      curve = (max-min)*(1-exp(-(0:nSteps-1)/tau)) + min;
end

end