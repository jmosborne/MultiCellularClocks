clear; close all;

addpath ../../ChasteNumericalMethods/MatlabScripts/

base_dir = '../../../testoutput/MultiCellularClocks/ClocksCompetition/'

num_sims = 10;

cell_numbers = [16, 0;
    0, 16;
    8, 8];

for case_index = 1:size(cell_numbers,1)

    figure
    hold on

    for i=1:num_sims

        data{i} = load([base_dir,'Case', num2str(case_index-1),'_Cycling', num2str(cell_numbers(case_index,1)),'_NonCycling',num2str(cell_numbers(case_index,2)), '/', num2str(i-1),'/results_from_time_0/results.non_cycling_cell_counts'])

        
        plot(data{i}(:,1),data{i}(:,2),'LineWidth',0.2,'LineStyle',':','Color',[1 0.7 0.7]);
        plot(data{i}(:,1),data{i}(:,3),'LineWidth',0.2,'LineStyle',':','Color',[0.7 0.7 1]);

        if i==1
            mean_data = data{i};
        else
            mean_data = mean_data + data{i};
        end
    end

    mean_data = mean_data/num_sims;



    % Plot the loaded data

    total_mean = mean_data(:,2) + mean_data(:,3);
    P0 = plot(mean_data(:,1),total_mean,'k','LineWidth',6.0);
    P1 = plot(mean_data(:,1),mean_data(:,2),'r','LineWidth',6.0);
    P2 = plot(mean_data(:,1),mean_data(:,3),'b','LineWidth',6.0);
    xlabel('Time');
    ylabel('Cell Counts');
    ylim([0,600]);
    title('Mean Cell Counts Over Time');
    legend([P0,P1,P2],{'Total','Non Cycing','Cycling'},'Location','northwest')
    SaveAsPngEpsAndFig(-1,['ClocksMonolayer_Cycling',num2str(cell_numbers(case_index,1)),'_NonCycling',num2str(cell_numbers(case_index,2))], 12, 7/5, 12)

end