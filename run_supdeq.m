% SUpDEq time-aligned SH interpolation (onset-based alignment, OBTA) with
% magnitude correction off and on, for every file written by `python mca.py export`.
% The wall-clock time of each call is written to results/supdeq_times.csv.
root = fullfile(pwd, 'external', 'SUpDEq');
addpath(root, genpath(fullfile(root, 'thirdParty', 'AKtools')), ...
        genpath(fullfile(root, 'thirdParty', 'SOFiA R13_MIT-License')), ...
        fullfile(pwd, 'external', 'AuditoryToolbox'));  % Slaney's ERB filters, needed for MCA

files = dir(fullfile('results', 'mca', 'in', '*.mat'));
outFolder = fullfile('results', 'mca', 'out');
if ~exist(outFolder, 'dir'), mkdir(outFolder); end

seconds = zeros(numel(files), 2);
for k = 1:numel(files)
    in = load(fullfile(files(k).folder, files(k).name));
    hrtfs = struct('HRTF_L', in.HRTF_L, 'HRTF_R', in.HRTF_R, 'f', in.f, ...
                   'samplingGrid', in.samplingGrid, 'Nmax', double(in.Nmax), 'FFToversize', 1);
    if k == 1
        supdeq_interpHRTF(hrtfs, in.ipGrid, 'OBTA', 'SH', inf);  % warm-up call, not timed
    end
    tic; off = supdeq_interpHRTF(hrtfs, in.ipGrid, 'OBTA', 'SH', nan); seconds(k, 1) = toc;
    tic; on = supdeq_interpHRTF(hrtfs, in.ipGrid, 'OBTA', 'SH', inf); seconds(k, 2) = toc;
    HRTF_L_off = off.HRTF_L; HRTF_R_off = off.HRTF_R;
    HRTF_L_on = on.HRTF_L; HRTF_R_on = on.HRTF_R;
    save(fullfile(outFolder, files(k).name), 'HRTF_L_off', 'HRTF_R_off', 'HRTF_L_on', 'HRTF_R_on');
    fprintf('%d/%d %s\n', k, numel(files), files(k).name);
end
writetable(table({files.name}', seconds(:, 1), seconds(:, 2), ...
                 'VariableNames', {'file', 'obta_s', 'mca_s'}), fullfile('results', 'supdeq_times.csv'));
