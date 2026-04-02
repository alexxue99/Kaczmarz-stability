function [] = graph_test_suite(n, cond_num, err,y,save,varargin)
id = 'MATLAB:handle_graphics:exceptions:SceneNode';
warning('off', id);
names = {"poly", "exp", "lowrank1", "highrank", "harmonic"};

% old names - use the following line if trying to graph the old entire test
% suite consisting of 9 matrices
% names = {"poly", "exp", "cluster1", "cluster2", "lowrank1", "lowrank2", "highrank", "linear", "harmonic"};

labels = {"base", "base refinement", "accelerated", "accelerated refinement"};
shapes = {'--',':','-', '-'};
colors = {'r','b','g','black'};
idx = round(linspace(1, numel(y), 200));
for i = 1:size(err, 1)
    fig=figure(i);
    set(gcf, 'Name', names{i}, 'NumberTitle', 'off');
    
    maximum = 0;
    minimum = 1e6;
    for j = 1:numel(varargin)
        e = varargin{j};
        loglog(y(idx), e(i, idx), shapes{j}, 'Color', colors{j},'DisplayName', names{i} + " " + labels{j}); hold on;
        maximum = max([maximum e(i, :)]);
        minimum = min([minimum e(i, :)]);
        %disp(maximum)
        %disp(minimum)
    end 
    yline(err(i),"--", 'DisplayName', sprintf('A\\b'));
    xlabel('Iteration'); ylabel('Error');
    legend show;
    %disp(names{i} + " " + err(i))
    maximum = max([maximum err(i, :)]);
    minimum = min([minimum err(i, :)]);
    ylim([minimum/2, maximum*2]);
    cond_num_str = strrep(sprintf('%.0e', cond_num), '+', '');
    if save
        saveas(fig, "figs/"+names{i}+"_"+ n + "_" + cond_num_str + ".png");
        saveas(fig, "figs/"+names{i}+"_"+ n + "_" + cond_num_str + ".fig");
    end
end
end