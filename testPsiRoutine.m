close all;

% Plot
plots_on = 1;

% Number of days to run, trials per day, mice
n_mice = 1;
n_days = 6;
n_trials = 50;

% Regularization parameter for "Bayes" updates (0 = full bayes, 1 = no updates)
reg = 0.5;

% Forgetting parameter for helping deal with non-stationarity
fget = 0.00;

% Fraction catch trials
ct_mix = 0.0;

% In-mixture from method of constant stimuli
cs_mix = 0.0;

% Parameter grid bounds for the presentation
% algorithm. These are NOT mouse parameters.
bnds.threshs = [0.005, 0.5]; %[0.005, 0.5];
bnds.slopes  = [1,9];
bnds.lapses  = [0.1,0.1];
bnds.guess   = 0.1;

% Example data generating params
% true.thresh = 0.087;
% true.slope  = 3;
% true.lapse  = 0.05;
% true.guess  = 0;

% Old FIXED mouse parameters
% thresh.min = 0.08;
% thresh.max = 0.2;
% thresh.tau = 15;
%
% slope.min = 2;
% slope.max = 4;
% slope.tau = 15;

for mnum = 1:n_mice

   % New changing-by-day mouse parameters
   %mouse = getBasicMouseMeta();
   mouse = getRandMouseMeta();
   mouse = getMouse(n_days, n_trials, mouse);

   lapse = mouse.fixed.lapse;
   guess = mouse.fixed.guess;

   % Day 1 parameter priors for the presentation
   % algorithm. These are NOT mouse parameters.
   priors.mu_thresh = 0.25;
   priors.sd_thresh = 1.0;

   priors.mu_slope = 3.0;
   priors.sd_slope = 2.0;

   priors.a_lapse = 1;
   priors.b_lapse = 1;

   % Loop over sessions
   for day = 1:n_days

      % Use the values defined by our mouse
      thresh.curve = mouse.day(day).thresh.curve;
      slope.curve  = mouse.day(day).slope.curve;

      % When using fixed params defined above...
      %thresh.curve = getCurve(thresh, n_trials, 'down');
      %slope.curve  = getCurve(slope , n_trials, 'up');

      % Params which are always fixed (as of now)
      lapse.curve  = getCurve(lapse.min, lapse.max, lapse.tau, n_trials, 'down');
      guess.curve  = getCurve(guess.min, guess.max, guess.tau, n_trials, 'down');

      %
      const_slope  = bnds.slopes(1)  == bnds.slopes(2);
      const_thresh = bnds.threshs(1) == bnds.threshs(2);
      %const_slope  = mouse.day(day).slope.min  == mouse.day(day).slope.max;
      %const_thresh = mouse.day(day).thresh.min == mouse.day(day).thresh.max;
      if const_slope || const_thresh
         hist.posterior = zeros(40,1 ,n_trials);
      else
         hist.posterior = zeros(40,40,n_trials);
      end

      if day == 1
         psi = PsiRoutine(bnds.guess, bnds, priors);
         hist.posterior(:,:,1) = psi.posterior;
         t_beg = 2;
         clist = psi.threshs(1:5:end);
         ccounts = zeros(1,length(clist));
      else
         psi.intensities = [];
         psi.responses = [];
         psi.nr_observations = 0;
         t_beg = 1;
         ccounts = zeros(1,length(clist));
      end

      % Trial loop
      for t = t_beg:n_trials
         % Data generating parameters
         true.thresh = thresh.curve(t);
         true.slope  = slope.curve(t);
         true.lapse  = lapse.curve(t);
         true.guess  = guess.curve(t);

         % Get another stimulus
         rnum = rand();
         if rnum < ct_mix
            % Catch trial
            x = bnds.threshs(1);
            
         elseif rnum < ct_mix + cs_mix
            % Random stimulus from under-presented 
            % method of constant stimuli stims
            inds = find(ccounts == min(ccounts));
            inds = inds(randperm(length(inds)));
            x = clist(inds(1));
            
         else
            % Maximum entropy reduction suggestion
            x = psi.suggestIntensity(); 
            
         end
            
         

         % Get a response
         r = psi.queryObserver(x, true.thresh, true.slope, true.lapse, true.guess);

         % Update the psi posterior densities
         psi = psi.updatePosterior(x,r,reg, fget);

         % Save the posterior
         hist.posterior(:,:,t) = psi.posterior;
      end

      %%
      map_thresh = getMAP(hist.posterior, 2, psi.threshs);
      map_slope  = getMAP(hist.posterior, 1, psi.slopes);

      hist.recovery_thresh(mnum,:,day) = [mean(mouse.day(day).thresh.curve(10:end)), mean(map_thresh(10:end))];
      hist.recovery_slope(mnum,:,day)  = [mean(mouse.day(day).slope.curve(10:end) ), mean(map_slope(10:end)) ];

      if plots_on
         figure();
         set(gcf,'Position', [95, 42, 1666, 876])

         % Threshold estimates
         subplot(3,4,1)
         plot(psi.threshs, mean(hist.posterior(:,:,  1),2), '-o'); hold on;
         plot(psi.threshs, mean(hist.posterior(:,:,end),2), '-o')
         ylims = ylim();
         plot([true.thresh, true.thresh], ylims, '--k', 'LineWidth', 2);
         grid on
         title('Threshold Posterior Densities')
         xlabel('Threshold')
         ylabel('Density')
         
         subplot(3,4,2)
         imagesc(squeeze(mean(hist.posterior(:,:,:),2)))
         yticklabels(round(psi.threshs(yticks),4))
         xlabel('Trial')
         ylabel('Threshold')
         title('Posterior Threshold Density')
         grid on
         
         subplot(3,4,3)
         plot(1:n_trials, map_thresh, '-o'); hold on
         plot(1:n_trials, thresh.curve, '--')
         ylim([0,0.5])
         title('Threshold MAP Estimate')
         xlabel('Trial')
         ylabel('Estimate')
         grid on

         % Slope estimates
         subplot(3,4,5)
         plot(psi.slopes, mean(hist.posterior(:,:,  1),1), '-o'); hold on;
         plot(psi.slopes, mean(hist.posterior(:,:,end),1), '-o')
         ylims = ylim();
         plot([true.slope, true.slope], ylims, '--k', 'LineWidth', 2);
         grid on
         title('Slope Posterior Densities')
         xlabel('Slope')
         ylabel('Density')
         
         subplot(3,4,6)
         imagesc(squeeze(mean(hist.posterior(:,:,:),1)))
         yticks(5:5:40)
         yticklabels(round(psi.slopes(yticks),4))
         xlabel('Trial')
         ylabel('Slope')
         title('Posterior Slope Density')
         grid on

         subplot(3,4,7)
         plot(1:n_trials, map_slope, '-o'); hold on
         plot(1:n_trials, slope.curve, '--')
         ylim([0,6])
         title('Slope MAP Estimate')
         xlabel('Trial')
         ylabel('Estimate')
         grid on

         % Misc
         
         subplot(3,4,4)
         plot(t_beg:n_trials, psi.intensities, '-o'); hold on;
         plot(1:n_trials, thresh.curve, '--')
         ylim([0,0.5])
         title('Stimulus Presentations')
         xlabel('Trial')
         ylabel('Intensity')
         grid on
         
         sind = floor(find(hist.posterior(:,:,n_trials) == max(max(hist.posterior(:,:,n_trials))))/40)+1;
         tind = mod(find(hist.posterior(:,:,n_trials) == max(max(hist.posterior(:,:,n_trials)))),40);

         subplot(3,4,8)
         plot(psi.pints, squeeze(psi.pfunc(tind,sind,1,:)), '-o')
         hold on; 
         plot(psi.pints, psi.weibull(psi.pints, true.thresh, true.slope, true.lapse, true.guess), 'r--');
         ylim([0,1])
         grid on
         title('MAP Psychometric Function')
         legend({'Est', 'True'}, 'Location', 'SouthEast')
         
         subplot(3,4,9)
         imagesc(hist.posterior(:,:,1))
         title('Initial Posterior')
         xlabel('Thresh')
         ylabel('Slope')
         
         subplot(3,4,10)
         imagesc(hist.posterior(:,:,n_trials))
         title('Final Posterior')
         xlabel('Thresh')
         ylabel('Slope')
         suptitle(['\bf{Day ' num2str(day) ' Simulation Data}'])
         drawnow()

      end
   end
end

figure()
cmap = parula(n_days);
dlabels = cell(n_days,1);
for d = 1:n_days
   scatter(hist.recovery_thresh(:,1,d), hist.recovery_thresh(:,2,d), 'MarkerEdgeColor', cmap(d,:)); hold on;
   dlabels{d} = ['Day ', num2str(d)];
end
xlim([0,0.5])
ylim([0,0.5])
grid on
plot([0,0.5], [0,0.5], '--k')
title('Threshold Parameter Recovery')
xlabel('True Threshold')
ylabel('Estimate')
legend(dlabels, 'Location', 'NorthWest')
