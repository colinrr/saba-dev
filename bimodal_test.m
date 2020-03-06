

nruns = 100;
bc = zeros([nruns 1]);
pv = bc;

mu_std = 3.3; % Mean separation of means, in units of standard deviation
std_std = 1;  % Spread of standard deviations

devs = randn([nruns 1])*mu_std;
sig1 = ones([nruns 1]); %randn([nruns 1])+1;
sig2 = sig1; %randn([nruns 1])+1;

for ii=1:nruns
    mu1 = randn([1 1])+1;
    mu2 = randn([1 1])+1;
    
    b = sort([randn([128 1])*sig1(ii) ; randn([128 1])*sig2(ii) + devs(ii) ]);
    [~, pv(ii), ~,~]=HartigansDipSignifTest(b,500);
    bc(ii) = bimodalityCoeff(b);
    
end

[devs,I] = sort(devs);

figure
subplot(2,1,1)
plot(pv(I))
hold on
plot(bc(I))
plot([0 nruns],[5/9 5/9],'--k')

subplot(2,1,2)
plot(devs)
hold on
plot(sig1(I),'--')
plot(sig2(I),'--')
plot([0 nruns],[3.5 3.5],'--k')
plot([0 nruns],[-3.5 -3.5],'--k')
