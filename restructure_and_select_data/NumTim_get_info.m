function [subjNames,Hemispheres,Run_1, Run_2] = NumTim_get_info(NumTim_data, DT_runs)
%% gets subject & hemisphere names and stores DT-run information in struct for this run (incl maps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% NumTim_data: made with the NumTim_load_data script. Contains all values
%       that you could need for all following analyses
% DT_runs: the two data runs you want to compare (options: TimingAll, TimingEven,
%       TimingOdd, NumerosityAll, NumerosityEven, NumerosityOdd)
%
% Output
% subjNames: cells with all subject names in NumTim_data
% Hemispheres: cells with all Hemispheres in NumTim_data
% Run_1: struct with info about DT_run_1 (e.g., DTrun = TimingAll;
%           Type = Timing; Run = All; Map = {collection of timing maps})
% Run_2: struct with info about DT_run_2 (e.g., DTrun = TimingAll;
%           Type = Timing; Run = All; Map = {collection of timing maps})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% get names general variables
subjNames = fieldnames(NumTim_data);
mapNames = fieldnames(NumTim_data.S6); % S6 has all maps
Hemispheres = fieldnames(NumTim_data.S6.(mapNames{1}));

%% select info for specific comparison
all_index = find(contains(mapNames,"All"));
grayMap = 1;
if length(all_index) > 1 
    all_index = all_index(1);
    grayMap = 0;
end

TimingMapNames = mapNames(1:all_index);
NumerosityMapNames = mapNames(all_index+1:end);
if grayMap == 1
    NumerosityMapNames = [NumerosityMapNames; mapNames(all_index)];
end


for run = 1:length(DT_runs)
    eval(strcat('Run_',char(num2str(run)),'.DTrun = DT_runs{',char(num2str(run)),'};'));

    T_or_N = split(eval(strcat('Run_',char(num2str(run)),'.DTrun')),{'All','Even','Odd','Halves','RM'});
    A_O_or_E = split(eval(strcat('Run_',char(num2str(run)),'.DTrun')),{'Timing','Numerosity'});
    eval(strcat('Run_',char(num2str(run)),'.Type = T_or_N{1};'));
    eval(strcat('Run_',char(num2str(run)),'.Run = A_O_or_E{2};'));
    
    if eval(strcat('Run_',char(num2str(run)),'.Type == "Timing"'))
        eval(strcat('Run_',char(num2str(run)),'.Names = TimingMapNames;'));
    elseif eval(strcat('Run_',char(num2str(run)),'.Type == "Numerosity"'))
        eval(strcat('Run_',char(num2str(run)),'.Names = NumerosityMapNames;'));
    end
end