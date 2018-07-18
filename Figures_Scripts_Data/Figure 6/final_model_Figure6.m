%% Modeling ODEs with and without noise

tspan = [0:120]; % time is in hours
x0 = [0.12 0.1 0.1]; % initial values
R =.52; % 
stimend = 48;
c = 1; % c=1 for oscillating input, c=2 for continuous input
mu=0;
sigma=0.3;

numsamples = 30;
lastvalue = zeros(numsamples,3);
figure, hold on
diff=0;
for ind = 1:numsamples
 %   noise = lognrnd(mu,sigma);   % random lognormal noise
     noise=1;   % set noise =1 for no added noise
    [t,x] = ode45(@(t,x) model1(t,x,R,noise,stimend,c), tspan,x0);
    lastvalue(ind,1) = x(end,1);
    lastvalue(ind,2) = x(end,2);
    lastvalue(ind,3) = x(end,3);
    if lastvalue(ind,1) > 2
        diff=diff+1
    end
    figure(gcf)
    plot(t,x(:,1)) % plot PPARG(t)
    legend({'PPARG','CEBPA','FABP4'})
end
percentdiff=diff/numsamples
