function reconDir = saveRecon(sampDir, volume, varNames, options)
%% saveRecon:   Save reconstruction with variables, metadata, and text 
%%              files of the code.
% Arguments:
%   sampDir     - Path to sample directory
%   volume      - Reconstructed volume for saving
%   varNames    - Cell array of variable names for saving
%   options     - Struct of options for saving (see below)
% Results:
%   reconDir    - Path to save directory
%
%   Data is save to the following structure:
%       SampleDir/Reconstructions/[ID]_Recon/vars/
%       SampleDir/Reconstructions/[ID]_Recon/metadata/
%       SampleDir/Reconstructions/[ID]_Recon/code/
%
% Options: 
%     'Format'           - 'mat' or 'h5'
%     'IDMode'           - 'time' or 'number'
%                          (increments based on existing *_Recon folders)
%     'VolumeName'       - name to use for the volume dataset/field
%     'Metadata'         - struct of extra metadata fields to record (e.g.
%                          struct('Notes','full NA test','NikonObjective','20x_0.8NA'))
%     'CodeFiles'        - cell array of function names or file paths to save
%                          to code/ INSTEAD OF auto-tracing dependencies, e.g.
%                          {'MSBP_fwd.m','BPM_upd.m'}. The calling file itself
%                          is always included in addition to these. Bare
%                          names (with or without .m) are resolved via
%                          `which`. Use this when auto-detection pulls in
%                          more files than you want archived.
%     'ExcludeCodeFiles' - cell array of function/file names to drop from the
%                          AUTO-TRACED dependency list (ignored if 'CodeFiles'
%                          is given). Matched by filename, extension optional,
%                          e.g. {'someUtilityFunction'}.
%
% Example call:
%       saveRecon( 'D:\Data\Sample01', ...
%                  volume, ...
%                  {'r','p'}, ...
%                  'Format','h5', ...
%                  'IDMode','number', ...
%                  'Metadata', struct('Notes','Example notes...'), ...
%                  'CodeFiles', {'msFwd.mat, msBwd.mat});
%
% Notes:
%     - Complex-valued volumes saved as h5 are split into separate
%       '<name>_real' and '<name>_imag' datasets (HDF5 has no native
%       complex type). An attribute IsComplex=1 is added so this is
%       documented in the file itself.
%     - Code archiving uses matlab.codetools.requiredFilesAndProducts on
%       the calling file, so any subfunctions that live in their own .m
%       files (e.g. BPM_upd.m called by MSBP_fwd.m) are captured
%       automatically.
    arguments
        sampDir (1,:) char
        volume
        varNames cell = {}
        options.Format (1,:) char {mustBeMember(options.Format,{'mat','h5'})} = 'mat'
        options.IDMode (1,:) char {mustBeMember(options.IDMode,{'time','number'})} = 'time'
        options.VolumeName (1,:) char = 'reconObj'
        options.Metadata struct = struct()
        options.CodeFiles cell = {}
        options.ExcludeCodeFiles cell = {}
    end

    % Set up folder structure
    reconsDir = fullfile(sampDir, 'Reconstructions');
    if ~exist(reconsDir, 'dir')
        mkdir(reconsDir);
    end

    ID = generateID(reconsDir, options.IDMode);
    reconDir = fullfile(reconsDir, [ID '_Recon']);
    varsDir     = fullfile(reconDir, 'vars');
    metadataDir = fullfile(reconDir, 'metadata');
    codeDir     = fullfile(reconDir, 'code');

    mkdir(reconDir);
    mkdir(varsDir);
    mkdir(metadataDir);
    mkdir(codeDir);

    % Pull requested variables from the caller's workspace
    %   Large GPU array recon is kept separate from the small named variables
    extraData = struct();
    missingVars = {};
    for i = 1:numel(varNames)
        name = varNames{i};
        try
            extraData.(name) = evalin('caller', name);
        catch
            missingVars{end+1} = name; %#ok<AGROW>
        end
    end
    if ~isempty(missingVars)
        warning('saveReconstruction:missingVars', ...
            'Could not find these variables in the caller workspace: %s', ...
            strjoin(missingVars, ', '));
    end

    % Gather gpuArray named variables back to host memory
    gatheredFields = {};
    fnames = fieldnames(extraData);
    for i = 1:numel(fnames)
        f = fnames{i};
        if isa(extraData.(f), 'gpuArray')
            extraData.(f) = gather(extraData.(f));
            gatheredFields{end+1} = f; %#ok<AGROW>
        end
    end
    if isa(volume, 'gpuArray')
        gatheredFields{end+1} = options.VolumeName;
    end
    if ~isempty(gatheredFields)
        fprintf('Gathered %d gpuArray variable(s) to host memory before saving: %s\n', ...
            numel(gatheredFields), strjoin(gatheredFields, ', '));
    end

    % Save variables
    switch options.Format
        case 'mat'
            saveVarsMat(varsDir, volume, options.VolumeName, extraData);
        case 'h5'
            saveVarsH5(varsDir, volume, options.VolumeName, extraData);
    end

    % Identify the calling file and archive code
    st = dbstack('-completenames');
    if numel(st) >= 2
        callerFile = st(2).file;
    else
        callerFile = '';
        warning('saveReconstruction:noCaller', ...
            'saveReconstruction was called from the command line; skipping code archive.');
    end
    savedCodeFiles = saveCode(codeDir, callerFile, options.CodeFiles, options.ExcludeCodeFiles);

    % Save metadata
    meta = options.Metadata;
    meta.ID = ID;
    meta.Timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    meta.MATLABRelease = version;
    meta.Hostname = getHostname();
    meta.CallerFile = callerFile;
    meta.SavedVariableNames = [{options.VolumeName}, fieldnames(extraData)'];
    meta.GatheredFromGPU = gatheredFields;
    meta.SavedCodeFiles = savedCodeFiles;
    meta.Format = options.Format;

    saveMetadata(metadataDir, meta);

    fprintf('Reconstruction saved to:\n  %s\n', reconDir);
end

%% Helper Functions
function ID = generateID(reconsDir, mode)
%% generateID: Define ID for folder naming
% Arguments:
%   reconsDir   -   Path to reconstruction directory
%   mode        -   Type of naming structure
% Returns:
%   ID          -   String holding ID value
    switch mode
        case 'time'
            ID = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        case 'number'
            d = dir(reconsDir);
            d = d([d.isdir]);
            names = {d.name};
            isRecon = ~cellfun(@isempty, regexp(names, '_Recon$', 'once'));
            n = sum(isRecon);
            ID = sprintf('%03d', n + 1);
    end
end

function hostname = getHostname()
%% getHostname: Get hostname for metadata
% Returns:
%   hostname    -   Hostname char array
    hostname = getenv('COMPUTERNAME');
    if isempty(hostname)
        hostname = getenv('HOSTNAME');
    end
    if isempty(hostname)
        hostname = 'unknown';
    end
end

function saveVarsMat(varsDir, volume, volumeName, extraData)
%% saveVarsMat: Gather and save volume as mat
% Arguments:
%   varsDir     -   Path to variable directory
%   volume      -   Reconstruction volume
%   volumeName  -   Name of volume
%   extraData   -   Additional data to save along volume
    matFile = fullfile(varsDir, 'vars.mat');
    if isa(volume, 'gpuArray')
        volume = gather(volume);
    end
    s = extraData;
    s.(volumeName) = volume;
    save(matFile, '-struct', 's', '-v7.3');
end

function saveVarsH5(varsDir, volume, volumeName, extraData)
%% saveVarsH5: Gather and save volume as h5
% Arguments:
%   varsDir     -   Path to variable directory
%   volume      -   Reconstruction volume
%   volumeName  -   Name of volume
%   extraData   -   Additional data to save along volume
    h5File = fullfile(varsDir, 'vars.h5');
    if exist(h5File, 'file')
        delete(h5File);
    end

    % Volume gets a dedicated, memory-conscious path.
    writeVolumeH5(h5File, volumeName, volume);

    % Small named variables use the generic path.
    names = fieldnames(extraData);
    for i = 1:numel(names)
        writeVarH5(h5File, names{i}, extraData.(names{i}));
    end
end

function writeVolumeH5(h5File, name, volume)
%% writeVolumeH5: Write a volume to HDF5.
% Arguments:
%   h5file      -   Path to h5 file directory
%   name        -   Name of h5 file
%   volume      -   Volume to save
    dsetBase = ['/' name];

    if ~isreal(volume)
        realPart = gather(real(volume));
        h5create(h5File, [dsetBase '_real'], size(realPart), 'Datatype', class(realPart));
        h5write(h5File, [dsetBase '_real'], realPart);
        clear realPart;

        imagPart = gather(imag(volume));
        h5create(h5File, [dsetBase '_imag'], size(imagPart), 'Datatype', class(imagPart));
        h5write(h5File, [dsetBase '_imag'], imagPart);
        clear imagPart;

        h5writeatt(h5File, [dsetBase '_real'], 'IsComplex', uint8(1));
    else
        hostVol = gather(volume);
        h5create(h5File, dsetBase, size(hostVol), 'Datatype', class(hostVol));
        h5write(h5File, dsetBase, hostVol);
        clear hostVol;
    end
end

function writeVarH5(h5File, name, val)
%% writeVarH5: Write non-volume vars to h5
% Arguments:
%   h5file      -   Path to h5 file directory
%   name        -   Name of h5 file
%   val         -   Variable to save
%
%       Numeric values are stored normally with the exception of complex 
%       data; this is split across real/imag. 
%       Non-numeric (struct/cell/char/string) are JSON-encoded as key-value 
%       pairs.
    if isnumeric(val) && ~isempty(val)
        dsetBase = ['/vars/' name];

        if ~isreal(val)
            h5create(h5File, [dsetBase '_real'], size(val), 'Datatype', class(real(val)));
            h5write(h5File, [dsetBase '_real'], real(val));
            h5create(h5File, [dsetBase '_imag'], size(val), 'Datatype', class(imag(val)));
            h5write(h5File, [dsetBase '_imag'], imag(val));
            h5writeatt(h5File, [dsetBase '_real'], 'IsComplex', uint8(1));
        else
            h5create(h5File, dsetBase, size(val), 'Datatype', class(val));
            h5write(h5File, dsetBase, val);
        end
    else
        % Store as a JSON string attribute
        try
            jsonStr = jsonencode(val);
        catch
            jsonStr = sprintf('<could not encode variable "%s" of class %s>', name, class(val));
        end
        h5writeatt(h5File, '/', ['vars_' name '_json'], jsonStr);
    end
end

function savedCodeFiles = saveCode(codeDir, callerFile, codeFilesOverride, excludeCodeFiles)
%% saveCode: Save code as .txt files for version control
% Arguments:
%   codeDir             -   Path to code storage directory
%   callerFile          -   File where saveRecon is called
%   codeFilesOverride   -   List of files to explicitly save (Does not save
%                           all dependencies
%   excludeCodeFiles    -   List of files to exclude from saving (Saves all
%                           but these dependencies)
%
%       Numeric values are stored normally with the exception of complex 
%       data; this is split across real/imag. 
%       Non-numeric (struct/cell/char/string) are JSON-encoded as key-value 
%       pairs.
    savedCodeFiles = {};

    if ~isempty(codeFilesOverride)
        % Manual mode: Save caller file and the files listed
        deps = {};
        if ~isempty(callerFile) && exist(callerFile, 'file')
            deps{end+1} = callerFile;
        end
        for i = 1:numel(codeFilesOverride)
            resolved = resolveCodeFile(codeFilesOverride{i});
            if isempty(resolved)
                warning('saveReconstruction:fileNotFound', ...
                    'Could not locate code file "%s" on the MATLAB path; skipping.', ...
                    codeFilesOverride{i});
            else
                deps{end+1} = resolved; %#ok<AGROW>
            end
        end
        deps = unique(deps, 'stable');
    else
        % Automatic mode: Trace and save all dependencies of the caller file
        if isempty(callerFile) || ~exist(callerFile, 'file')
            return;
        end
        try
            deps = matlab.codetools.requiredFilesAndProducts(callerFile);
        catch ME
            warning('saveReconstruction:depTrace', ...
                'Could not trace dependencies (%s). Saving only the main file.', ME.message);
            deps = {callerFile};
        end

        % Filter out anything in excludeCodeFiles
        if ~isempty(excludeCodeFiles)
            excludeNames = cellfun(@(x) localBaseName(x), excludeCodeFiles, 'UniformOutput', false);
            keep = true(size(deps));
            for i = 1:numel(deps)
                if any(strcmpi(localBaseName(deps{i}), excludeNames))
                    keep(i) = false;
                end
            end
            deps = deps(keep);
        end
    end

    for i = 1:numel(deps)
        srcFile = deps{i};
        [~, fname, fext] = fileparts(srcFile);
        destFile = fullfile(codeDir, [fname fext '.txt']);
        try
            copyfile(srcFile, destFile);
            savedCodeFiles{end+1} = [fname fext]; %#ok<AGROW>
        catch ME
            warning('saveReconstruction:copyFail', ...
                'Could not copy %s: %s', srcFile, ME.message);
        end
    end

    function name = localBaseName(f)
        % Filename without path or extension
        [~, name] = fileparts(f);
    end
end

function resolved = resolveCodeFile(f)
%% resolveCodeFile: Accepts a path, function name, or .m file and returns a 
%%                  full path or '' if not found.
% Arguments:
%   f   - File name as char array
    name = f;
    if endsWith(name, '.m')
        name = name(1:end-2);
    end
    resolved = which(name);
    if ~isempty(resolved)
        return;
    end

    if exist(f, 'file') == 2
        fi = dir(f);
        if ~isempty(fi)
            resolved = fullfile(fi(1).folder, fi(1).name);
            return;
        end
    end

    resolved = '';
end

function saveMetadata(metadataDir, meta)
%% saveMetaData: Save a .mat and .txt metadata file
% Arguments:
%   metadataDir -   Path to metadata directory
%   meta        -   Metadata as struct
    save(fullfile(metadataDir, 'metadata.mat'), '-struct', 'meta');

    txtFile = fullfile(metadataDir, 'metadata.txt');
    fid = fopen(txtFile, 'w');
    fields = fieldnames(meta);
    for i = 1:numel(fields)
        f = fields{i};
        v = meta.(f);
        if iscell(v)
            v = strjoin(v, ', ');
        elseif isnumeric(v)
            v = mat2str(v);
        elseif isstruct(v)
            try
                v = jsonencode(v);
            catch
                v = '<struct>';
            end
        end
        fprintf(fid, '%s: %s\n', f, v);
    end
    fclose(fid);
end