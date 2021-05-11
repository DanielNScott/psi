function [alpha,beta]  = psyWeibFit(x,k,n)
	% [alpha,beta] mlPsyWeibFit(x,k,n)
	% Takes as parameters the dimension measured (x), the observed number
	% correct (k), and the total number of trials (n). It returns the
	% maximum-likelihood estimates [alpha,beta] for a fitting Weibull
	% function
	mu = mean(x);

	start_val = [mu;3.0];
	%%options = optimset('MaxFunEvals',1e5,'MaxIter',1e5,'TolFun',1e-4);
	params = fminsearch(@(val) weibLike(x,val,k,n),start_val);
	beta = params(2);
	alpha = params(1);
%%
end

function res2 = weibLike(x_vals,params,k,n)
	if(any(params<1e-4))
		res2 = 1e-100;
	else
		beta = max(0,params(2));
		alpha = max(0,params(1));
		weibull =inline('0.5+0.5*(1-exp(-(contrast./thresh).^slope))',...
			'contrast','thresh','slope');
        p_est = weibull(x_vals,alpha,beta);
        y = sum(-logbinopdf(k,n,p_est)); 
		res2 = y;
	end
end

function res3 = logbinopdf(k,n,p)
	logbinocoeff = gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1);
	res3 = logbinocoeff+k.*log(p)+(n-k).*log(1-p);
end