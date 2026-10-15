% SUpDEq time-aligned SH interpolation (onset-based alignment, OBTA) with
% magnitude correction off and on, for every file written by `python mca.py export`.
root = fullfile(pwd, 'external', 'SUpDEq');
addpath(root, genpath(fullfile(root, 'thirdParty', 'AKtools')), ...
        genpath(fullfile(root, 'thirdParty', 'SOFiA R13_MIT-License')), ...
        fullfile(pwd, 'external', 'AuditoryToolbox'));  % Slaney's ERB filters, needed for MCA

files = dir(fullfile('results', 'mca', 'in', '*.mat'));
outFolder = fullfile('results', 'mca', 'out');
if ~exist(outFolder, 'dir'), mkdir(outFolder); end

for k = 1:numel(files)
    in = load(fullfile(files(k).folder, files(k).name));
    hrtfs = struct('HRTF_L', in.HRTF_L, 'HRTF_R', in.HRTF_R, 'f', in.f, ...
                   'samplingGrid', in.samplingGrid, 'Nmax', double(in.Nmax), 'FFToversize', 1);
    off = supdeq_interpHRTF(hrtfs, in.ipGrid, 'OBTA', 'SH', nan);
    on = supdeq_interpHRTF(hrtfs, in.ipGrid, 'OBTA', 'SH', inf);
    HRTF_L_off = off.HRTF_L; HRTF_R_off = off.HRTF_R;
    HRTF_L_on = on.HRTF_L; HRTF_R_on = on.HRTF_R;
    save(fullfile(outFolder, files(k).name), 'HRTF_L_off', 'HRTF_R_off', 'HRTF_L_on', 'HRTF_R_on');
    fprintf('%d/%d %s\n', k, numel(files), files(k).name);
end
