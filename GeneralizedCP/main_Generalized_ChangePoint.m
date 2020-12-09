function [G, MDL, states, eff, eff_fit, breaks, output, records, excluded] = main_Generalized_ChangePoint(eff)
%% The main function of the generalized change point algorithm
% Input: 
%       single 1D trace (eff) or multiple traces selected in dialog
% Output:
%       G: structure, recording all the optimum clusterings under different
%   number of states
%       MDL: the minimum discription length, used to determine the optimum
%   number of states
%       states: the fitting based on the optimum number of states
%       eff: group the traces of different data together, if use input eff
%   trace, the output will be identical to the input
%       eff_fit: the fitting of all the feasible number of states, up to 30
%       breaks: recording the separations among different traces
%       output: recording several important parameters for potential usage
%       records: structure, recording the analysis of each loaded trace,
%   and also recording the location of each trace in the output eff
%       excluded: structure, recording all the traces not being used

%% step 1: loading the traces and change-points detection
[filelist,datadir] = uigetfile({'*.mat';'*.dat'},'multiselect','on');
groups = [];
breaks = [];
records = struct([]);
excluded = struct([]);
% make a judgement if use input trace or use loaded traces
if iscell(filelist)
    numfiles = numel(filelist);
    eff = [];
elseif isstr(filelist)
    numfiles = 1;
    eff = [];
else % use the input eff and detect the change point
    numfiles = 0;
    sd = w1_noise(diff(eff))/sqrt(2);% estimate the noise level
    points = main_change_point_wavelet(eff);% change points detection
    T = numel(eff);
    groups = [1, points+1; points, T];% record the starts and ends of each segments
    breaks(end+1) = T;
end
% load each trace and detect change points
for n = 1 : numfiles
    if numfiles == 1
        fullname = strcat(datadir, filelist);
    else
        fullname = strcat(datadir, filelist{n});
    end
    data = importdata(fullname);
    temp = data.obsf;% for FRET trace only
    %temp = data.denf;
    T1 = numel(eff)+1;
    sd(n) = w1_noise(diff(temp))/sqrt(2);% estimate the noise level
    points = main_change_point_wavelet(temp);% change points detection
    if numel(points) < 300
        eff = [eff, temp];% group traces together
        T2 = numel(eff);
        breaks(end+1) = T2;
        groups = cat(2,groups,[T1, points+T1; points+T1-1, T2]);
        % record all the information for each trace, specially for our FRET
        % data
        records(end+1).filename = fullname;
        records(end).obsf = data.obsf;
        records(end).denf = data.denf;
        records(end).acc = data.accf_bf;
        records(end).don = data.donf_bf;
        records(end).loc = [T1, T2];
    else% too many change-points, exclude this trace
        excluded(end+1).filename = fullname;
        excluded(end).obsf = data.obsf;
        excluded(end).denf = data.denf;
        excluded(end).acc = data.accf_bf;
        excluded(end).don = data.donf_bf;
    end
end
sd = max(sd);% use the maximum noise level among these traces as the global noise level

%% step 2 and 3: clustering the segments and calculate MDL
[G, Ij, Tj] = clustering_GCP(eff, groups);
G = G(end:-1:1);% flip the G
n_mdl = min(30, numel(G));% calculate up to 30 states
MDL = zeros(1,n_mdl);
eff_fit = zeros(n_mdl, numel(eff));
for i = 1:n_mdl;
    [MDL(i), eff_fit(i,:)] = MDL_piecewise(Ij, Tj, G(i), eff, groups, sd, breaks);
end
[~, q] = min(MDL);%now the BIC is actually MDL

output.Ij = Ij;
output.Tj = Tj;
output.sd = sd;
output.groups = groups;

figure
plot(MDL)

disp([num2str(q) ' states were used']);
states = eff_fit(q,:);
if ~isempty(records)
    for n = 1:numel(records)
        records(n).states = states(records(n).loc(1):records(n).loc(2));
    end
end
figure
plot(eff)
hold on
plot(states,'r')
hold off