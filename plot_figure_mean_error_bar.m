function plot_figure_mean_error_bar(data_output,varargin)

%PLOT_FIGURE_MEAN_ERROR_BAR computes the mean of the data and standard error to plot data
%% Parse parameters
p = inputParser;
addRequired(p,'data_output',@isstruct);
addParamValue(p,'pathOut','./',@isstr);
addParamValue(p,'file_append',date,@isstr);

parse(p,data_output,varargin{:});

pathOut=p.Results.pathOut;
file_append=p.Results.file_append;

%% Compute all set point values

for iStrain=1:length(data_output)
    
    mean_data(iStrain)=mean(data_output(iStrain).values);
    standard_error(iStrain)=std(data_output(iStrain).values)./sqrt(length(data_output(iStrain).values));
    
end

%% Plot data

[meanDataSorted,idx]=sort(mean_data);
medianDataSort=median(meanDataSorted);

namesStrains={data_output.strain};
hfig=figure('Position',[ 328   198   865   397]);
hold all;

errorbar(meanDataSorted,standard_error(idx),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',8)
ylim([-9 -3])

xticklabel_rotate(1:length(data_output),45,namesStrains(idx))
set(gca,'box','off')


Set_fig_RE(hfig,9,9,20);

filename=[pathOut 'Mean_error_bar_' file_append];
export_fig(filename,'-pdf',  '-transparent', '-nocrop')
