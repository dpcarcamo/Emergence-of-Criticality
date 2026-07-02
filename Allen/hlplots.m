% Plot stored Allen h/lambda fits and compare against full-data theory.

load('Allenhldata.mat')

%% Clean placeholders and select plotted range

Nss(Nss==0)=nan;
hs(hs==0) = nan;

valid = Nss < 1200;

%% h/lambda trajectories with full-data theory point

figure
hold on
for i = 1:size(Nss,2)
    
    plot(hs(valid(:,i),i), ls(valid(:,i),i), '-', 'MarkerSize', 10)
    
end


mus = mean(squeeze(Data(1,:,:,:)),1, 'omitnan');
A = reshape(mus,[40,size(Nss,2)]);
A(:,isnan(A(1,:))) = [];
B = ~isnan(A);
Indices = arrayfun(@(x) find(B(:, x), 1, 'last'), 1:size(A, 2));
Values = arrayfun(@(x,y) A(x,y), Indices, 1:size(A, 2));
muf = Values;
chis = mean(squeeze(Data(2,:,:,:)),1, 'omitnan');
A = reshape(chis,[40,size(Nss,2)]);
A(:,isnan(A(1,:))) = [];
B = ~isnan(A);
Indices = arrayfun(@(x) find(B(:, x), 1, 'last'), 1:size(A, 2));
Values = arrayfun(@(x,y) A(x,y), Indices, 1:size(A, 2));
chif = mean(Values);
A = Nss;
A(:,isnan(A(1,:))) = [];
B = ~isnan(A);
Indices = arrayfun(@(x) find(B(:, x), 1, 'last'), 1:size(A, 2));
Values = arrayfun(@(x,y) A(x,y), Indices, 1:size(A, 2));
N = mean(Values);


u0 = sqrt(chif/N + muf.^2);

happrx = atanh( muf./u0)/(N.*u0);
lamapprx = atanh(u0)./u0;

plot(happrx,lamapprx,'r.', MarkerSize=20)


axis square
xlabel('h', FontSize=18)
ylabel('\lambda', FontSize=18)

%% h versus N

figure
hold on
for i = 1:size(Nss,2)
    
    plot(Nss(valid(:,i),i), hs(valid(:,i),i), '-', 'MarkerSize', 10)
    
end
plot(N,happrx,'r.', MarkerSize=20)
axis square
xlabel('N', FontSize=18)
ylabel('h', FontSize=18)
xscale('log')

%% lambda versus N

figure
hold on
for i = 1:size(Nss,2)
    
    plot(Nss(valid(:,i),i), ls(valid(:,i),i), '-', 'MarkerSize', 10)
    
end
plot(N,lamapprx,'r.', MarkerSize=20)
axis square
xlabel('N', FontSize=18)
ylabel('\lambda', FontSize=18)
xscale('log')

%% h versus chi

figure
hold on
for i = 1:size(Nss,2)
    
    plot(chis(1,valid(:,i),i), hs(valid(:,i),i), '-', 'MarkerSize', 10)
    
end
%plot(chif,happrx,'r.', MarkerSize=20)
axis square
xlabel('\chi', FontSize=18)
ylabel('h', FontSize=18)

%% lambda versus mu

figure
hold on
for i = 1:size(Nss,2)
    
    plot(mus(1,:,i), ls(:,i), '-', 'MarkerSize', 10)
    
end
plot(muf,lamapprx,'r.', MarkerSize=20)
axis square
xlabel('\mu', FontSize=18)
ylabel('\lambda', FontSize=18)

