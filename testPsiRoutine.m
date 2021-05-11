close all;

% Plot
plots_on = 0;

% Number of days to run, trials per day, mice
n_mice = 20;
n_days = 10;
n_trials = 50;

% Regularization parameter for "Bayes" updates
reg = 0.5;

% Parameter grid bounds for the presentation
% algorithm. These are NOT mouse parameters.
bnds.threshs = [0.005, 0.5];
bnds.slopes  = [1,6];
bnds.lapses  = [0.1,0.1];
bnds.guess   = 0.1;

% Example data generating params
% true.thresh = 0.087;
% true.slope  = 3;
% true.lapse  = 0.05;
% true.guess  = 0;

% We represent attentional state as a function which 
% starts low, rapidly attains a maximum, fluctuates,
% and slowly increases over trial time.

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
   mouse = getRandMouseMeta();
   mouse = getMouse(n_days, n_trials, mouse);

   lapse = mouse.fixed.lapse;
   guess = mouse.fixed.guess;

   % Day 1 parameter priors for the presentation
   % algorithm. These are NOT mouse parameters.
   priors.mu_thresh = 0.25;
   priors.sd_thresh = 1.0;

   priors.mu_slope = 3.0;
   priors.sd_slope = 1.0;

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
      if mouse.day(day).slope.min == mouse.day(day).slope.max
         hist.posterior = zeros(40,1 ,n_trials);
      else
         hist.posterior = zeros(40,40,n_trials);
      end

      if day == 1
         psi = PsiRoutine(bnds.guess, bnds, priors);
         hist.posterior(:,:,1) = psi.posterior;
         t_beg = 2;
      else
         psi.intensities = [];
         psi.responses = [];
         psi.nr_observations = 0;
         t_beg = 1;
      end

      % Trial loop
      for t = t_beg:n_trials
         % Data generating parameters            
         true.thresh = thresh.curve(t);
         true.slope  = slope.curve(t);
         true.lapse  = lapse.curve(t);
         true.guess  = guess.curve(t);

         % Get another stimulus
         x = psi.suggestIntensity();

         % Get a response
         r = psi.queryObserver(x, true.thresh, true.slope, true.lapse, true.guess);

         % Update the psi posterior densities
         psi = psi.updatePosterior(x,r,reg);

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
         set(gcf,'Position', [95         141        1058         833])

         subplot(3,2,1)
         plot(psi.threshs, mean(hist.posterior(:,:,  1),2), '-o'); hold on;
         plot(psi.threshs, mean(hist.posterior(:,:,end),2), '-o')
         ylims = ylim();
         plot([true.thresh, true.thresh], ylims, '--k', 'LineWidth', 2);
         grid on
         title('Threshold Posterior Densities')
         xlabel('Threshold')
         ylabel('Density')

         subplot(3,2,2)
         plot(t_beg:n_trials, psi.intensities, '-o'); hold on;
         plot(1:n_trials, thresh.curve, '--')
         ylim([0,0.5]) 
         title('Stimulus Presentations')
         xlabel('Trial')
         ylabel('Intensity')
         grid on

         subplot(3,2,3)
         imagesc(squeeze(mean(hist.posterior(:,:,:),2)))
         yticklabels(round(psi.threshs(yticks),4))
         xlabel('Trial')
         ylabel('Threshold')
         title('Posterior Threshold Density')
         grid on

         subplot(3,2,4)
         plot(1:n_trials, map_thresh, '-o'); hold on
         plot(1:n_trials, thresh.curve, '--')
         ylim([0,0.5]) 
         title('Threshold MAP Estimate')
         xlabel('Trial')
         ylabel('Estimate')
         grid on

         subplot(3,2,5)
         imagesc(squeeze(mean(hist.posterior(:,:,:),1)))
         yticks(5:5:40)
         yticklabels(round(psi.slopes(yticks),4))
         xlabel('Trial')
         ylabel('Slope')
         title('Posterior Slope Density')
         grid on

         subplot(3,2,6)
         plot(1:n_trials, map_slope, '-o'); hold on
         plot(1:n_trials, slope.curve, '--')
         ylim([0,6])
         title('Slope MAP Estimate')
         xlabel('Trial')
         ylabel('Estimate')
         grid on

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