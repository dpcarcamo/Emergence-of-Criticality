% Plot stored Stringer h/lambda fits versus population size.

load('stringerhldata.mat')

%% h/lambda trajectories
figure
hold on
for i = 1:45
    
    plot(hs(Nss(:,i)<1200,i), ls(Nss(:,i)<1200,i), '-o', 'MarkerSize', 10)
    
end

axis square
xlabel('h')
ylabel('\lambda')

[20,30,39]


%% h versus N

figure
hold on
for i = 1:45
    
    plot(Nss(Nss(:,i)<1200,i), hs(Nss(:,i)<1200,i), '-o', 'MarkerSize', 10)
    
end

axis square
xlabel('N')
ylabel('h')
xscale('log')



%% lambda versus N

figure
hold on
for i = 1:45
    
    plot(Nss(Nss(:,i)<1200,i), ls(Nss(:,i)<1200,i), '-o', 'MarkerSize', 10)
    
end

axis square
xlabel('N')
ylabel('\lambda')
xscale('log')
