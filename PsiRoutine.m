classdef PsiRoutine < handle
   
   properties(Hidden)
      xx;
      ss;
      tt;
      ll;
      possible_intensities;
      nx;
      nt;
      ns;
      nl;
   end
   
   properties(SetAccess = public)
      posterior
      intensities
      responses
      threshs
      slopes
      lapses
      guess
      %lapse
      pfunc
      nr_observations = 0;
   end
   
   methods
      function psi = PsiRoutine(guess_rate, bnds, priors)
         
         % Default bounds for parameters
         if nargin < 3
            bnds.threshs = [0.005, 0.5];
            bnds.slopes  = [-3, 3];
            bnds.lapses  = [0.005, 0.05];
         end
         
         % Grid resolution in parameter space
         n_bins_threshs = 40;
         n_bins_slopes  = 40;
         n_bins_lapses  = 10;
         
         % If the upper and lower bounds are the same, override the resolution
         if bnds.threshs(1) == bnds.threshs(2); n_bins_threshs = 1; end
         if bnds.slopes(1)  == bnds.slopes(2);  n_bins_slopes  = 1; end
         if bnds.lapses(1)  == bnds.lapses(2);  n_bins_lapses  = 1; end

         % Set up the parameter space limits
         psi.threshs = exp( linspace(log(bnds.threshs(1)), log(bnds.threshs(2)), n_bins_threshs) );         
         psi.slopes  = linspace(bnds.slopes(1), bnds.slopes(2), n_bins_slopes);
         psi.lapses  = linspace(bnds.lapses(1), bnds.lapses(2), n_bins_lapses);
         
         %psi.lapses = lapse_rate;
         psi.guess = guess_rate;
         %psi.lapse = lapse_rate;
         
         % Set up parameter space grid and prior density
         % Assume uniform priors on the threshold and slope, beta on lapse.
         [tt,ss,ll] = ndgrid(psi.threshs, psi.slopes, psi.lapses);
         
         p_thresh = normpdf(log(tt), log(priors.mu_thresh), priors.sd_thresh);
         p_thresh = p_thresh/sum(p_thresh(:));
         
         p_slope  = normpdf(ss, priors.mu_slope, priors.sd_slope);
         p_slope  = p_slope/sum(p_slope(:));
         
         p_lapse  = betapdf(ll, priors.a_lapse, priors.b_lapse);
         %p_lapse  = betapdf(ll/0.05, 1.5, 3) + 0.5;
         p_lapse  = p_lapse/sum(p_lapse(:));
         
         % Theta is the vector of slope, thresh, lapse
         p_theta = p_slope.*p_thresh.*p_lapse;
         p_theta = p_theta/sum(p_theta(:));
         
         %p_theta = ones(length(psi.threshs),length(psi.slopes),length(psi.lapses));
         
         psi.posterior = p_theta;
         psi.possible_intensities = linspace(0.01,0.5,50);
         
         % Number of intensities, thresholds, slopes, lapses
         psi.nx = length(psi.possible_intensities);
         psi.nt = length(psi.threshs);
         psi.ns = length(psi.slopes);
         psi.nl = length(psi.lapses);
         
         % Grid of 
         [psi.tt, psi.ss, psi.ll, psi.xx] = ...
            ndgrid(psi.threshs, psi.slopes, psi.lapses, psi.possible_intensities);
         
         % 
         psi.pfunc = psi.weibull(psi.xx, psi.tt, psi.ss, psi.ll);
      end
      
      
      % Function for suggesting a next stimulus
      function res = suggestIntensity(psi)
         % Get the expected entropy of ... after each possible
         % next stimulus presentation
         entropy = psi.expectedEntropy();
         
         % Find the point (stimulus) at which it is minimum
         [~,min_idx] = min(entropy);
         
         % Return the stimulus to present
         res = psi.possible_intensities(min_idx);
      end
      
      
      % Function for updating the posterior probability of the parameters
      function res = updatePosterior(psi, x, r, reg)
         
         % Save the intensity and response
         psi.intensities = [psi.intensities, x];
         psi.responses   = [psi.responses  , r];

         % Increment the observation counter
         psi.nr_observations = psi.nr_observations+1;
                  
         % Retrieve (return) the parameter grid
         [t,s,l] = ndgrid(psi.threshs, psi.slopes, psi.lapses);
         
         % Determine the likelihood of this response under each parameter set
         likelihood = psi.weibull(x,t,s,l);
         if(~r)
            likelihood = 1-likelihood;
         end
         
         % Update the posterior distribution
         % Note the regularization, so that this is not strictly Bayes.
         post = reg*psi.posterior + (1-reg)*(psi.posterior.*likelihood);
         psi.posterior = post./sum(post(:));
         
         % Return
         res = psi;
      end
  
      
      % Maximum-A-Posteriori parameter estimate
      function params = computeMapTheta(psi)
         
         dims = size(psi.posterior);
         
         ts_posterior = sum(psi.posterior,3);
         [~,max_idx] = max(ts_posterior(:));
         
         [it,is] = ind2sub(dims,max_idx);
         
         disp([it,is])
         params = [psi.threshs(it), psi.slopes(is)];
      end
      
      % Weibull function with given parameters
      function p = weibull(psi, x, t, s, l)
         p = psi.guess + (1-l-psi.guess).*(1-exp(-(x./t).^s));
      end   
      
      % Function for querying an observer with a known PF (for testing purposes)
      function res = queryObserver(psi, x, thresh, slope, lapse, guess)
         if(nargin < 6)
            guess = 0.5;
         end
         if(nargin < 5)
            lapse = 0.01;
         end
         if(nargin < 4)
            slope = 3.5;
         end
         if(nargin < 3)
            thresh = 0.087;
         end
         p   = psyWeib(x, thresh, slope, lapse, guess);
         res = binornd(1,p);
      end
   end   
   
   methods(Hidden)
      function  p = p_r_x(psi,prior)
         p =  sum(sum(sum(psi.pfunc.*repmat(prior,[1,1,1,psi.nx]),3),2),1);
      end
      
      function [p_theta_r1,p_theta_r0] = p_theta_x_r(psi,prior)
         p_r = psi.p_r_x(prior);
         p_theta_r1 = psi.pfunc.*repmat(prior,[1,1,1,psi.nx])./repmat(p_r,[psi.nt,psi.ns,psi.nl,1]);
         p_theta_r0 = (1-psi.pfunc).*repmat(prior,[1,1,1,psi.nx])./repmat((1-p_r),[psi.nt,psi.ns,psi.nl,1]);
      end
      
      % Get the entropy 
      function mu_H = expectedEntropy(psi)
         
         prior = psi.posterior;
         [p_theta_r1, p_theta_r0] = psi.p_theta_x_r(prior);
         p_r1  = squeeze(psi.p_r_x(prior));
         H_x_1 = squeeze(-sum(sum(sum(p_theta_r1.*log(p_theta_r1),3),2),1));
         H_x_0 = squeeze(-sum(sum(sum(p_theta_r0.*log(p_theta_r0),3),2),1));
         
         mu_H = (H_x_1.*p_r1)+(H_x_0.*(1-p_r1));
      end
   end
   
   
end

function p = psyWeib(x,t,s,lapse,guess)
p = guess+(1-lapse-guess).*(1-exp(-(x./t).^s));
end