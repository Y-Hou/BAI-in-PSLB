% To draw the figures, please run naive_DBAI and naive_DBAIbeta.m first.
% Then run this code.
close all;

% load('processing.mat',"table0_8","table1_0","table1_2"); % processing.mat contains all results for misspecified Lmax experiments 
load('trial_2.mat',"table1_0");
load('trial_1.mat',"table0_8");
load('trial_3.mat',"table1_2");
n_epsilon = 12;
number_of_exp = 20;

table_PSBAI = table1_0;

output = 0; % set it to be 1 if needs to save the figures.

% color
color_PSBAI = [0 0.4470 0.7410];
color_naive = [0.8500 0.3250 0.0980];
color_DBAI = [0.9290 0.6940 0.1250];
color_DBAIbeta = [0.4940 0.1840 0.5560];
color_alpha = [0.4660 0.6740 0.1890];
color_xi =[0.6350 0.0780 0.1840];

%% Number of change points
%DBAI
average = mean(table_DBAI,3);
ncp_DBAI = average(1,:);
stand_dv = std (table_DBAI,0,3);
ncp_DBAI_std = stand_dv(1,:);
%DBAI_beta
average = mean(table_DBAIbeta,3);
ncp_DBAIbeta = average(1,:);
stand_dv = std (table_DBAIbeta,0,3);
ncp_DBAIbeta_std = stand_dv(1,:);
%PSBAI
ncp_PSBAI = zeros(n_epsilon,1);
ncp_PSBAI_std = zeros(n_epsilon,1);
failure_count = 0;

% the sample complexity for the failed experiments are excluded.
for i = 1:n_epsilon
    bestarm = table_PSBAI(2,i,:);
    ncp_i = table_PSBAI(4,i,:);
    failure_count = failure_count + (number_of_exp-sum(bestarm==1)); % for large epsilon, arm 3, even arm 2, maybe an epsilon best arm
    valid_ncp_i = ncp_i(bestarm==1);
    ncp_PSBAI(i) = mean(valid_ncp_i);
    ncp_PSBAI_std(i) = std(valid_ncp_i);
end

%plot figures
f1 = figure('Position',[0 0 800 400]);
set(gcf,'PaperPositionMode','auto');
set(gcf,'Visible', 'on');
n_epsilon = 12;
epsilon_range = 0.03*ones(n_epsilon,1);
for i=1:n_epsilon
    epsilon_range(i) = epsilon_range(i)*1.35^(i);
end
epsilon_range = 1./(0.0703+epsilon_range).^2;


errorbar(epsilon_range,ncp_DBAI,ncp_DBAI_std,'-o','Color', color_DBAI);
hold on;
errorbar(epsilon_range,ncp_DBAIbeta,ncp_DBAIbeta_std,'-^','Color',color_DBAIbeta);
hold on;
errorbar(epsilon_range,ncp_PSBAI,ncp_PSBAI_std,'b-*','Color',color_PSBAI);
hold on;
ylim([0, 900])
xlim([0, 82])

rectangle('Position',[0.1 0.1 7 40],'LineStyle', '-','LineWidth', 1)
x = [7, 58];
y = [0, 45];
line(x, y, 'LineStyle', ':','LineWidth', 1,'Color', 'k');
x = [7, 58];
y = [40, 275];
line(x, y, 'LineStyle', ':','LineWidth', 1,'Color', 'k');

ylabel('Number of Context Samples');
xlabel('Squared relaxed mean gap $1/(\Delta_{\min}+\varepsilon)^2$','Interpreter','latex');
leg = legend('D$\varepsilon$BAI','D$\varepsilon$BAI$_\beta$','PS$\varepsilon$BAI$^+$','Interpreter','latex','Orientation','horizontal');
set(gca, 'Fontname', 'Times New Roman','FontSize',18);
str=['ContextSample','.eps'];

% the small figure
axes('Position',[.68 .17 .21 .2])
box on
sub = 8:12;
errorbar(epsilon_range(sub),ncp_DBAI(sub),ncp_DBAI_std(sub),'-o','Color', color_DBAI);
hold on;
errorbar(epsilon_range(sub),ncp_DBAIbeta(sub),ncp_DBAIbeta_std(sub),'-^','Color', color_DBAIbeta);
hold on;
errorbar(epsilon_range(sub),ncp_PSBAI(sub),ncp_PSBAI_std(sub),'-*','Color', color_PSBAI);
hold on;
xlim([0.5, 6.5]);
ylim([0.5, 30]);

% naive
average = mean(table_PSBAI,3);
ncp_naive = average(8,:);
stand_dv = std (table_PSBAI,0,3);
ncp_naive_std = stand_dv(8,:);
axes('Position',[.17 .6 .3 .3])
box on
errorbar(epsilon_range,ncp_naive,ncp_naive_std,'->','Color', color_naive);
ylim([10^3, 3*10^5])
xlim([0, 82])
set(gca,'Yscale','log');
legend('N$\varepsilon$BAI','Fontname', 'Times New Roman','FontSize',18,'Interpreter','latex','Location','southeast');

if output == 1
    exportgraphics(f1,str);
end
%% Sample complexity for misspecified Lmax
% PSBAI
% 1
complexity_PSBAI = zeros(n_epsilon,1);
complexity_PSBAI_std = zeros(n_epsilon,1);
% 0.8
table_PSBAI0_8 = table0_8;
complexity_PSBAI0_8 = zeros(n_epsilon,1);
complexity_PSBAI_std0_8 = zeros(n_epsilon,1);
% 1.2
table_PSBAI1_2 = table1_2;
complexity_PSBAI1_2 = zeros(n_epsilon,1);
complexity_PSBAI_std1_2 = zeros(n_epsilon,1);

% naive
complexity_naive = zeros(n_epsilon,1);
complexity_naive_std = zeros(n_epsilon,1);

% the sample complexity for the failed experiments are excluded.
for i = 1:n_epsilon
    bestarm = table_PSBAI(2,i,:);
    complexity_PSBAI_i = table_PSBAI(1,i,:);
    valid_complexity_PSBAI_i = complexity_PSBAI_i(bestarm==1);
    complexity_PSBAI(i) = mean(valid_complexity_PSBAI_i);
    complexity_PSBAI_std(i) = std(valid_complexity_PSBAI_i);

    bestarm = table_PSBAI0_8(2,i,:);
    complexity_PSBAI_i = table_PSBAI0_8(1,i,:);
    valid_complexity_PSBAI_i = complexity_PSBAI_i(bestarm==1);
    complexity_PSBAI0_8(i) = mean(valid_complexity_PSBAI_i);
    complexity_PSBAI_std0_8(i) = std(valid_complexity_PSBAI_i);

    bestarm = table_PSBAI1_2(2,i,:);
    complexity_PSBAI_i = table_PSBAI1_2(1,i,:);
    valid_complexity_PSBAI_i = complexity_PSBAI_i(bestarm==1);
    complexity_PSBAI1_2(i) = mean(valid_complexity_PSBAI_i);
    complexity_PSBAI_std1_2(i) = std(valid_complexity_PSBAI_i);

    bestarm_naive = table_PSBAI(6,i,:);
    complexity_naive_i = table_PSBAI(5,i,:);
    valid_complexity_naive_i = complexity_naive_i(bestarm==1);
    complexity_naive(i) = mean(valid_complexity_naive_i);
    complexity_naive_std(i) = std(valid_complexity_naive_i);
end

clear bar;
tc=[complexity_PSBAI0_8,complexity_PSBAI,complexity_PSBAI1_2,complexity_naive];
f2 = figure('Position',[0 0 800 400]);
set(gcf,'PaperPositionMode','auto');
set(gcf,'Visible', 'on');
n_epsilon = length(complexity_naive_std);

% Manually set colors for each bar series
colors = [color_alpha; 
          color_PSBAI; 
          color_xi;
          color_naive];

b= bar(1:n_epsilon,tc);
for k = 1:4
    b(k).FaceColor = colors(k, :);
end

hold on;
err=[complexity_PSBAI_std0_8,complexity_PSBAI_std,complexity_PSBAI_std1_2,complexity_naive_std];
[ngroups,nbars] = size(tc);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',tc,err,'k','linestyle','none');
set(gca,'Yscale','log','Fontname', 'Times New Roman','FontSize',18);
xlabel('Instance $k$','Interpreter','latex');
ylabel('Sample complexity', 'Interpreter','latex');
ylim([10^5,1.1*10^10]);
legend('PS$\varepsilon$BAI$^+$ w/ $\nu=0.8$','PS$\varepsilon$BAI$^+$','PS$\varepsilon$BAI$^+$ w/ $\nu=1.2$','N$\varepsilon$BAI','Interpreter','latex','Location','northeast');
hold off
str=['SampleComplexityComparison_misspe_L','.eps'];
grid;
if output == 1
    exportgraphics(f2,str);
end

%% Sample complexity
clear bar;
tc=[complexity_PSBAI,complexity_naive];
f3 = figure('Position',[0 0 800 400]);
set(gcf,'PaperPositionMode','auto');
set(gcf,'Visible', 'on');
n_epsilon = length(complexity_naive_std);

% Manually set colors for each bar series
colors = [color_PSBAI;  
          color_naive];

b= bar(1:n_epsilon,tc);
for k = 1:2
    b(k).FaceColor = colors(k, :);
end

hold on;
err=[complexity_PSBAI_std,complexity_naive_std];
[ngroups,nbars] = size(tc);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',tc,err,'k','linestyle','none');
set(gca,'Yscale','log','Fontname', 'Times New Roman','FontSize',18);
xlabel('Instance $k$','Interpreter','latex');
ylabel('Sample complexity', 'Interpreter','latex');
ylim([10^5,1.1*10^10]);
legend('PS$\varepsilon$BAI$^+$','N$\varepsilon$BAI','Interpreter','latex','Location','northeast');

hold off
str=['SampleComplexityCompar','.eps'];
grid;
if output == 1
    exportgraphics(f3,str);
end
