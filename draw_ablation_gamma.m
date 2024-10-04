%% Run this code to get the plot for the experiment of abaltion on gamma.
clear;
close all;

output = 1; % set it to be 1 if needs to save the figures.

%% obtain the results for gamma=6 from the previous results
load("trial_2.mat","table1_0");
table_PSBAI = table1_0(:,5:12,:);

%% Aggregate all the results
gamma_set = [2,3,6,12];
n_gamma = length(gamma_set);
table_collector = {};
table_collector{3} = table_PSBAI;

for i_gamma=1:n_gamma
    if i_gamma == 3
        continue;
    end
    filename = strcat('trial_gamma_',num2str(i_gamma),'.mat');
    load(filename,"table");
    disp(filename)
    table_collector{i_gamma} = table;
end

%%
n_epsilon = 8;
number_of_exp = 20;

% color
color_PSBAI = [0 0.4470 0.7410];
color_naive = [0.8500 0.3250 0.0980];
color_DBAI = [0.9290 0.6940 0.1250];
color_DBAIplus = [0.4940 0.1840 0.5560];
color_alpha = [0.4660 0.6740 0.1890];
color_xi =[0.6350 0.0780 0.1840];

%% Sample complexity for different gamma

complexity_collector = {};
complexity_collector_std = {};

complexity_collector_naive = {};
complexity_collector_naive_std = {};

for index = 1:n_gamma
    complexity_PSBAI = zeros(n_epsilon,1);
    complexity_PSBAI_std = zeros(n_epsilon,1);

    complexity_naive = zeros(n_epsilon,1);
    complexity_naive_std = zeros(n_epsilon,1);

    table_PSBAI = table_collector{index};
% the sample complexity for the failed experiments are excluded.
    for i = 1:n_epsilon
        bestarm = table_PSBAI(2,i,:);
        complexity_PSBAI_i = table_PSBAI(1,i,:);
        valid_complexity_PSBAI_i = complexity_PSBAI_i(bestarm==1);
        complexity_PSBAI(i) = mean(valid_complexity_PSBAI_i);
        complexity_PSBAI_std(i) = std(valid_complexity_PSBAI_i);
        
        bestarm = table_PSBAI(6,i,:);
        complexity_naive_i = table_PSBAI(5,i,:);
        valid_complexity_naive_i = complexity_naive_i(bestarm==1);
        complexity_naive(i) = mean(valid_complexity_naive_i);
        complexity_naive_std(i) = std(valid_complexity_naive_i);
    end
    complexity_collector{index} = complexity_PSBAI;
    complexity_collector_std{index} = complexity_PSBAI_std;
    complexity_collector_naive{index} = complexity_naive;
    complexity_collector_naive_std{index} = complexity_naive_std;
end

clear bar;
tc = [];
for index = 1:n_gamma
    tc=[tc,complexity_collector{index}];
end

tc=[tc,complexity_collector_naive{3}];

f2 = figure('Position',[0 0 800 400]);
set(gcf,'PaperPositionMode','auto');
set(gcf,'Visible', 'on');

% Manually set colors for each bar series
color_PSBAI_gamma = [   [0.9290 0.6940 0.1250];
                        [0.4940 0.1840 0.5560];
                        color_PSBAI;
                        [0.6350 0.0780 0.1840]
                    ];

colors = [color_PSBAI_gamma;
            color_naive
          ];

b= bar(1+4:n_epsilon+4,tc);
for k = 1: n_gamma+1
    b(k).FaceColor = colors(k, :);
end
hold on;
err = [];
for index = 1:n_gamma
    err=[err,complexity_collector_std{index}];
end
err=[err,complexity_collector_naive_std{3}];

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
lgd = legend('PS$\varepsilon$BAI$^+$ w/$w=3\tilde{w}$' ...
    ,'PS$\varepsilon$BAI$^+$ w/$w=2\tilde{w}$' ...
    ,'PS$\varepsilon$BAI$^+$ w/$w=\tilde{w}$' ...
    , 'PS$\varepsilon$BAI$^+$ w/$w=\tilde{w}/2$'...
    ,'N$\varepsilon$BAI' ...
    ,'Interpreter','latex','Location','northeast');
lgd.NumColumns = 2;
hold off
str=['SampleComplexityComparison_ablation_w','.eps'];
% str=['SampleComplexityComparison_ablation_w','.pdf'];
grid;
if output == 1
    exportgraphics(f2,str);
end