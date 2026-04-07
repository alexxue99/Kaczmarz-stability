function [] = graph_test_suite(n, cond_num, err,y,save,varargin)
set(groot, 'defaultLineLineWidth', 4);
set(groot, 'defaultAxesLineWidth', 1);    

id = 'MATLAB:handle_graphics:exceptions:SceneNode';
warning('off', id);
names = {"poly", "exp", "highrank", "harmonic"};

% old names - use the following line if trying to graph the old entire test
% suite consisting of 9 matrices
% names = {"poly", "exp", "cluster1", "cluster2", "lowrank1", "lowrank2", "highrank", "linear", "harmonic"};

labels = {"base", "base refinement", "accelerated", "accelerated refinement"};
shapes = {':','--',':', '--'};
colors = {'r','b','g','black'};
idx = round(linspace(1, numel(y), 200));
for i = 1:size(err, 1)
    fig=figure(i);
    fig.Position(3:4) = [400, 300]; 
    set(gcf, 'Name', names{i}, 'NumberTitle', 'off');
    
    for j = 1:numel(varargin)
        e = varargin{j};
        loglog(y(idx), e(i, idx), shapes{j}, 'Color', colors{j},'DisplayName', names{i} + " " + labels{j}); hold on;
    end 
    yline(err(i),"--", 'DisplayName', sprintf('A\\b'), 'Color', [0 100/256 50/256], 'LineWidth', 4);
    %legend show;
    %disp(names{i} + " " + err(i))

    if cond_num == 1e5
        ylim([1e-5*7, 2]);
        yticks([1e-4, 1e-3, 1e-2, 1e-1, 1]);
        yticklabels({'10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}', '10^{0}'});
    end

    if cond_num == 1e4
        ylim([1e-6*7, 2]);
        yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1]);
        yticklabels({'10^{-5}', '10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}', '10^{0}'});
    end

    set(gca, 'FontWeight', 'bold', 'FontSize', 20, 'TickLabelInterpreter', 'tex');
    xlabel('Iteration', 'FontSize', 30);
    ylabel('Error', 'FontSize', 30);

    xlim([1e6, 1e9]);
    xticks([1e6, 1e7, 1e8, 1e9]);
    cond_num_str = strrep(sprintf('%.0e', cond_num), '+', '');
    if save
        saveas(fig, "figs/"+names{i}+"_"+ n + "_" + cond_num_str + ".png");
        saveas(fig, "figs/"+names{i}+"_"+ n + "_" + cond_num_str + ".fig");
    end
end
end