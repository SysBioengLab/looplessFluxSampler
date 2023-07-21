% Plot figures.
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------

% Load data for figure 2, S1 and S2 (run example3.m previously)
load('samples\R0.mat',"R0")
load('samples\R2.mat',"R2")
load('samples\R10.mat',"R10")
load('samples\R12.mat',"R12")
load('samples\R0_HRB.mat',"R0_HRB")
load('samples\R2_HRB.mat',"R2_HRB")
load('samples\R10_HRB.mat',"R10_HRB")
load('samples\R12_HRB.mat',"R12_HRB")

% Load data for figure 3
load("samples\m1_ADSB.mat","m1_ADSB")          % model: e_coli_core
load("samples\m1_llACHRB.mat","m1_llACHRB")    % model: e_coli_core
load("samples\m2_ADSB.mat","m2_ADSB")          % model: iIT341
load("samples\m2_llACHRB.mat","m2_llACHRB")    % model: iIT341
load("samples\m3_ADSB.mat","m3_ADSB")          % model: iYO844
load("samples\m3_llACHRB.mat","m3_llACHRB")    % model: iYO844
load("samples\m4_ADSB.mat","m4_ADSB")          % model: iMM904
load("samples\m4_llACHRB.mat","m4_llACHRB")    % model: iMM904

%% Figure 2
figure()
tiledlayout(1,3)
% A: Mu
ax(1) = nexttile(1);
loglog(abs(R0.mu),abs(R0_HRB.mu), 'x', 'markersize', 6)
hold on
loglog(abs(R2.mu),abs(R2_HRB.mu), 'x', 'markersize', 6)
loglog(abs(R10.mu),abs(R10_HRB.mu), 'x', 'markersize', 6)
loglog(abs(R12.mu),abs(R12_HRB.mu), 'x', 'markersize', 6)
hold off
title('A');
ax(1).TitleHorizontalAlignment = 'Left';
xlabel('Flux mean (ADSB)')
ylabel('Flux mean (HR)')
legend('R0','R2','R10','R12','location','NorthWest','FontSize',8);
% B: Sigma
ax(2) = nexttile(2);
loglog(R0.sigma,R0_HRB.sigma, 'x', 'markersize', 6)
hold on
loglog(R2.sigma,R2_HRB.sigma, 'x', 'markersize', 6)
loglog(R10.sigma,R10_HRB.sigma, 'x', 'markersize', 6)
loglog(R12.sigma,R12_HRB.sigma, 'x', 'markersize', 6)
hold off
title('B');
ax(2).TitleHorizontalAlignment = 'Left';
xlabel('Flux standard deviation (ADSB)')
ylabel('Flux standard deviation (HR)')
legend('R0','R2','R10','R12','location','NorthWest','FontSize',8);
% C: PSRF
x1 = R0_HRB.R';
x2 = R0.R';
x3 = R2_HRB.R';
x4 = R2.R';
x5 = R10_HRB.R';
x6 = R10.R';
x7 = R12_HRB.R';
x8 = R12.R';

group = [ones(numel(x1),1);2*ones(numel(x2),1);3*ones(numel(x3),1);4*ones(numel(x4),1);...
            5*ones(numel(x5),1);6*ones(numel(x6),1);7*ones(numel(x7),1);8*ones(numel(x8),1)];

positions = [.85,1.15,1.85,2.15,2.85,3.15,3.85,4.15];

ax(3) = nexttile(3);
boxplot([x1;x2;x3;x4;x5;x6;x7;x8],group,'positions', positions,'symbol','');
bar(positions([1,3,5,7]),[mean(x1),mean(x3),mean(x5),mean(x7)],.15,'FaceColor',[0,0.45,0.74]);
hold on
bar(positions([2,4,6,8]),[mean(x2),mean(x4),mean(x6),mean(x8)],.15,'FaceColor',[1,0.41,0.16]);

p1 = prctile(x1,[2.5,97.5]);
p2 = prctile(x2,[2.5,97.5]);
p3 = prctile(x3,[2.5,97.5]);
p4 = prctile(x4,[2.5,97.5]);
p5 = prctile(x5,[2.5,97.5]);
p6 = prctile(x6,[2.5,97.5]);
p7 = prctile(x7,[2.5,97.5]);
p8 = prctile(x8,[2.5,97.5]);

errorbar(positions,[mean(x1),mean(x2),mean(x3),mean(x4),mean(x5),mean(x6),mean(x7),mean(x8)],...
    [p1(1),p2(1),p3(1),p4(1),p5(1),p6(1),p7(1),p8(1)]-[mean(x1),mean(x2),mean(x3),mean(x4),mean(x5),mean(x6),mean(x7),mean(x8)],...
    [p1(2),p2(2),p3(2),p4(2),p5(2),p6(2),p7(2),p8(2)]-[mean(x1),mean(x2),mean(x3),mean(x4),mean(x5),mean(x6),mean(x7),mean(x8)],'.k')

set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) mean(positions(5:6)) mean(positions(7:8))])
set(gca,'xticklabel',{'R0','R2','R10','R12'})

legend('HR','ADSB','location','NorthWest','FontSize',8);
axis([.5,4.5,0,4])
title('C');
ax(3).TitleHorizontalAlignment = 'Left';
xlabel('Model')
ylabel('Potential Scale Redution Factor (psrf)')

%% Figure 3

figure()
tiledlayout(1,2)
% Time per Neff (left)
x1 = m1_llACHRB.samplingTime./m1_llACHRB.Neff';
x2 = m1_ADSB.samplingTime./m1_ADSB.Neff';
x3 = m2_llACHRB.samplingTime./m2_llACHRB.Neff';
x4 = m2_ADSB.samplingTime./m2_ADSB.Neff';
x5 = m3_llACHRB.samplingTime./m3_llACHRB.Neff';
x6 = m3_ADSB.samplingTime./m3_ADSB.Neff';
x7 = m4_llACHRB.samplingTime./m4_llACHRB.Neff';
x8 = m4_ADSB.samplingTime./m4_ADSB.Neff';

group = [ones(numel(x1),1);2*ones(numel(x2),1);3*ones(numel(x3),1);4*ones(numel(x4),1);...
            5*ones(numel(x5),1);6*ones(numel(x6),1);7*ones(numel(x7),1);8*ones(numel(x8),1)];

positions = [.85,1.15,1.85,2.15,2.85,3.15,3.85,4.15];

ax(1) = nexttile(1);
boxplot(([x1;x2;x3;x4;x5;x6;x7;x8]),group,'positions', positions,'symbol','');

set(gca,'YScale','log')
set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) mean(positions(5:6)) mean(positions(7:8))])
set(gca,'xticklabel',{'E. coli core','iIT341','iYO844','iMM904'})

color = [1,0.41,0.16;...
         0,0.45,0.74;...
         1,0.41,0.16;...
         0,0.45,0.74;...
         1,0.41,0.16;...
         0,0.45,0.74;...
         1,0.41,0.16;...
         0,0.45,0.74];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.9);
end
c = get(gca,'Children');

hleg1 = legend(c(1:2),'ll-ACHRB','ADSB');
legend('Location','northwest')
axis([.5,4.5,10^-6,10^1])
title('A');
ax(1).TitleHorizontalAlignment = 'Left';
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
xlabel('Model')
ylabel('Time per effective sample')

% Potential scale reduction factor (right)
x1 = m1_llACHRB.R';
x2 = m1_ADSB.R';
x3 = m2_llACHRB.R';
x4 = m2_ADSB.R';
x5 = m3_llACHRB.R';
x6 = m3_ADSB.R';
x7 = m4_llACHRB.R';
x8 = m4_ADSB.R';

group = [ones(numel(x1),1);2*ones(numel(x2),1);3*ones(numel(x3),1);4*ones(numel(x4),1);...
            5*ones(numel(x5),1);6*ones(numel(x6),1);7*ones(numel(x7),1);8*ones(numel(x8),1)];

positions = [.85,1.15,1.85,2.15,2.85,3.15,3.85,4.15];

ax(2) = nexttile(2);
boxplot([x1;x2;x3;x4;x5;x6;x7;x8],group,'positions', positions,'symbol','');
bar(positions([1,3,5,7]),[mean(x1),mean(x3),mean(x5),mean(x7)],.15,'FaceColor',[0,0.45,0.74]);
hold on
bar(positions([2,4,6,8]),[mean(x2),mean(x4),mean(x6),mean(x8)],.15,'FaceColor',[1,0.41,0.16]);

p1 = prctile(x1,[2.5,97.5]);
p2 = prctile(x2,[2.5,97.5]);
p3 = prctile(x3,[2.5,97.5]);
p4 = prctile(x4,[2.5,97.5]);
p5 = prctile(x5,[2.5,97.5]);
p6 = prctile(x6,[2.5,97.5]);
p7 = prctile(x7,[2.5,97.5]);
p8 = prctile(x8,[2.5,97.5]);

errorbar(positions,[mean(x1),mean(x2),mean(x3),mean(x4),mean(x5),mean(x6),mean(x7),mean(x8)],...
    [p1(1),p2(1),p3(1),p4(1),p5(1),p6(1),p7(1),p8(1)]-[mean(x1),mean(x2),mean(x3),mean(x4),mean(x5),mean(x6),mean(x7),mean(x8)],...
    [p1(2),p2(2),p3(2),p4(2),p5(2),p6(2),p7(2),p8(2)]-[mean(x1),mean(x2),mean(x3),mean(x4),mean(x5),mean(x6),mean(x7),mean(x8)],'.k')

set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) mean(positions(5:6)) mean(positions(7:8))])
set(gca,'xticklabel',{'E. coli core','iIT341','iYO844','iMM904'})

legend('ll-ACHRB','ADSB','location','NorthWest','FontSize',8);
axis([.5,4.5,0,4])
title('B');
ax(2).TitleHorizontalAlignment = 'Left';
xlabel('Model')
ylabel('Potential Scale Redution Factor (psrf)')

%% Figure S1
% Calculate standarized differences for mu and sigma
x1 = R0_HRB.mu - R0.mu;
x2 = R2_HRB.mu - R2.mu;
x3 = R10_HRB.mu - R10.mu;
x4 = R12_HRB.mu - R12.mu;
x5 = R0_HRB.sigma - R0.sigma;
x6 = R2_HRB.sigma - R2.sigma;
x7 = R10_HRB.sigma - R10.sigma;
x8 = R12_HRB.sigma - R12.sigma;

mu1 = (x1 - mean(x1))/ std(x1);
mu2 = (x2 - mean(x2))/ std(x2);
mu3 = (x3 - mean(x3))/ std(x3);
mu4 = (x4 - mean(x4))/ std(x4);
sigma1 = (x5 - mean(x5))/ std(x5);
sigma2 = (x6 - mean(x6))/ std(x6);
sigma3 = (x7 - mean(x7))/ std(x7);
sigma4 = (x8 - mean(x8))/ std(x8);

p1 = prctile(mu1,[2.5,97.5]);
p2 = prctile(mu2,[2.5,97.5]);
p3 = prctile(mu3,[2.5,97.5]);
p4 = prctile(mu4,[2.5,97.5]);
p5 = prctile(sigma1,[2.5,97.5]);
p6 = prctile(sigma2,[2.5,97.5]);
p7 = prctile(sigma3,[2.5,97.5]);
p8 = prctile(sigma4,[2.5,97.5]);

% Plot
MU = [mu1' mu2' mu3' mu4'];
SIGMA = [sigma1' sigma2' sigma3' sigma4'];
grp = [zeros(1,80),ones(1,86),2*ones(1,86),3*ones(1,86)];

figure()
tiledlayout(1,2)
ax(1) = nexttile(1);
boxplot(MU,grp,'Labels',{'R0','R2','R10','R12'},'colors','k','symbol','','whisker',0)
title('A');
xlabel('Models')
ylabel('Mean flux differences (standarized)')
ax(1).TitleHorizontalAlignment = 'Left';

hold on

errorbar([1,2,3,4],[mean(mu1),mean(mu2),mean(mu3),mean(mu4)],...
    [p1(1),p2(1),p3(1),p4(1)]-[mean(mu1),mean(mu2),mean(mu3),mean(mu4)],...
    [p1(2),p2(2),p3(2),p4(2)]-[mean(mu1),mean(mu2),mean(mu3),mean(mu4)],'.k')

hold off

ax(2) = nexttile(2);
boxplot(SIGMA,grp,'Labels',{'R0','R2','R10','R12'},'colors', 'k','symbol','','whisker',0)
xlabel('Models')
ylabel('Standard deviation of flux differences (standarized)')
title('B')
ax(2).TitleHorizontalAlignment = 'Left';
hold on
errorbar([1,2,3,4],[mean(sigma1),mean(sigma2),mean(sigma3),mean(sigma4)],...
    [p5(1),p6(1),p7(1),p8(1)]-[mean(sigma1),mean(sigma2),mean(sigma3),mean(sigma4)],...
    [p5(2),p6(2),p7(2),p8(2)]-[mean(sigma1),mean(sigma2),mean(sigma3),mean(sigma4)],'.k')
hold off

%% Figure S2

% Time per Neff
x1 = R0_HRB.samplingTime./R0_HRB.Neff';
x2 = R0.samplingTime./R0.Neff';
x3 = R2_HRB.samplingTime./R2_HRB.Neff';
x4 = R2.samplingTime./R2.Neff';
x5 = R10_HRB.samplingTime./R10_HRB.Neff';
x6 = R10.samplingTime./R10.Neff';
x7 = R12_HRB.samplingTime./R12_HRB.Neff';
x8 = R12.samplingTime./R12.Neff';

group = [ones(numel(x1),1);2*ones(numel(x2),1);3*ones(numel(x3),1);4*ones(numel(x4),1);...
            5*ones(numel(x5),1);6*ones(numel(x6),1);7*ones(numel(x7),1);8*ones(numel(x8),1)];

positions = [.85,1.15,1.85,2.15,2.85,3.15,3.85,4.15];

figure()
boxplot(([x1;x2;x3;x4;x5;x6;x7;x8]),group,'positions', positions,'symbol','');

set(gca,'YScale','log')
set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) mean(positions(5:6)) mean(positions(7:8))])
set(gca,'xticklabel',{'R0','R2','R10','R12'})

color = [1,0.41,0.16;...
         0,0.45,0.74;...
         1,0.41,0.16;...
         0,0.45,0.74;...
         1,0.41,0.16;...
         0,0.45,0.74;...
         1,0.41,0.16;...
         0,0.45,0.74];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.9);
end
c = get(gca,'Children');

hleg1 = legend(c(1:2),'HR','ADSB');
legend('Location','northwest')
axis([.5,4.5,10^-6,10^1])
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
xlabel('Model')
ylabel('Time per effective sample')

