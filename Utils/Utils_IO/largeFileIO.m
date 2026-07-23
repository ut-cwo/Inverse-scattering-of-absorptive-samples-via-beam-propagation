classdef largeFileIO < handle
%% largeFileIO: HDF5 chunked file I/O
%       Synchronous acquisition and async writing and reading
%
% Async Write Usage:
%   fileIO = largeFileIO(dataPath, fileName, dims, saveParams)
%   fileIO.bindWriter(pool)         % Bind async writer to queue
%   fileIO.sendChunk(chunk,lamIdx)  % Send chunk for write
%   fileIO.finalizeWrite()          % Write Metadata
%
%   (Static) largeFileIO.writerChunkWorker(fullPath, dims, data)  % Write Loop
%
% Serial Write Usage:
%   fileIO.writeAllChunks(imageMatrix)      % Optional full block writing
%
% Read Usage:
%   [imData, params] = fileIO.read()
%
%   (Static) [imData, params] = largeFileIO.readAllChunks(dataPath, fileName)
properties (SetAccess = private, GetAccess = public)
    FullPath
    FileName
    SaveParams
    Dims
    ChunkDims
    Queue           % DataQueue - main thread sends chunks to workers
    ChunksWritten   % counter incremented on each ack
end

%% Methods
methods
    function fileIO = largeFileIO(dataPath, fileName, dims, saveParams)
        %% largeFileIO: Constructor for largeFileIO class
        % Arguments:
        %   dataPath   - Output directory
        %   fileName   - Base file name (no extension)
        %   dims       - Full dataset dimensions [H x W x N x nLambdas]
        %   saveParams - Metadata struct
        arguments
            dataPath    (1,:) char
            fileName    (1,:) char
            dims
            saveParams
        end

        if ~exist(dataPath, 'dir')
            mkdir(dataPath);
        end
        
        if nargin > 3
            fileIO.SaveParams = saveParams;
            save(fullfile(dataPath, 'params.mat'), 'saveParams');
            fprintf('Saved acquisiton parameters.\n')
        else
            fprintf('Did not save acquisiton parameters.\n')
        end

        fileIO.FileName = [fileName '.h5'];
        fileIO.FullPath = fullfile(dataPath, fileIO.FileName);

        assert(length(dims) == 4, 'Incremental write requires 4D data [H x W x N x n Lambdas].');
        chunkDims = [dims(1), dims(2), dims(3), 1];
        fileIO.Dims = dims;
        fileIO.ChunkDims = chunkDims;

        fileIO.Queue        = parallel.pool.DataQueue;
        fileIO.ChunksWritten = 0;

        h5create(fileIO.FullPath, '/ChunkedImage', dims, ...
            'Datatype',  'uint8', ...
            'ChunkSize', chunkDims, ...
            'Deflate',   0);

        fprintf('HDF5 file created: %s\n', fileIO.FullPath);
        fprintf('Dimensions: [%s], Chunk: [%s]\n\n', ...
            num2str(dims), num2str(chunkDims));
    end

    function writeChunk(fileIO, chunk, lamIdx)
        %% writerChunkWorker: Write one chunk, called by afterEach
        % Arguments:
        %   chunk - Directory containing the file
        %   lamIdx - Base file name (no extension)
        arguments
            fileIO
            chunk       
            lamIdx      
        end
        start  = [1 1 1 lamIdx];
        count  = [fileIO.Dims(1), fileIO.Dims(2), fileIO.Dims(3), 1];
    
        tWrite = tic;

        h5write(fileIO.FullPath, '/ChunkedImage', chunk, start, count);

        writeTime = toc(tWrite);
        writeTime_m = floor(writeTime / 60);
        writeTime_s = mod(writeTime, 60);
        fprintf('Chunk %d write time: %dm %.2fs\n', lamIdx, writeTime_m, writeTime_s);
    end

    function finalizeWrite(fileIO)
        %% finalizeWrite: Write metadata after all chunks confirmed written
        h5writeatt(fileIO.FullPath, '/', 'dtype',   'uint8');
        h5writeatt(fileIO.FullPath, '/', 'created', char(datetime));
        fprintf('HDF5 finalized: %s\n', fileIO.FullPath);
    end

    function writeAllChunks(fileIO, imageMatrix)
        %% writeAllChunks: Write full image matrix in one call (blocking)
        % Arguments:
        %   imageMatrix - [H x W x N x nLambdas] uint8
        arguments
            fileIO
            imageMatrix     (:,:,:,:) uint8
        end

        imageSize_gb = numel(imageMatrix) / 1e9;
        fprintf('Starting image write (%.2f GB)...\n', imageSize_gb)
        tic;

        h5write(fileIO.FullPath, '/ChunkedImage', imageMatrix);
        fileIO.finalizeWrite();

        writeTime   = toc;
        writeTime_m = floor(writeTime / 60);
        writeTime_s = mod(writeTime, 60);
        fprintf('Finished in %dm %.2fs\n', writeTime_m, writeTime_s);
    end

    function [imData, params] = read(fileIO, verbose)
        %% read: Read the file this object points to
        % Arguments:
        %   verbose - If true, display h5 metadata
        arguments
            fileIO
            verbose     logical = false
        end

        [imData, params] = largeFileIO.readAllChunks(...
            fileparts(fileIO.FullPath), fileIO.FileName, verbose);
    end
end

%% Static Methods
methods (Static)
    function writeChunkWorker(fullPath, dims, chunk, lamIdx)
        %% writeChunkWorker: Static worker for parfeval
        % Arguments:
        %   fullPath - Directory containing the file
        %   dims - Base file name (no extension)
        %   chunk  - (optional) If true, display h5 metadata
        %   lamIdx - 
        arguments
            fullPath    (1,:) char
            dims        
            chunk
            lamIdx
        end

        assert(isa(chunk, 'uint8'), 'largeFileIO: Chunk must be uint8.');

        start = [1 1 1 lamIdx];
        count = [dims(1), dims(2), dims(3), 1];
        tWrite = tic;

        h5write(fullPath, '/ChunkedImage', chunk, start, count);

        fprintf('Lambda %d written in %.2fs\n', lamIdx, toc(tWrite));
    end

    function [imData, params] = readAllChunks(dataPath, fileName, verbose)
        %% readAllChunks: Read a chunked HDF5 file from disk
        % Arguments:
        %   dataPath - Directory containing the file
        %   fileName - Base file name (no extension)
        %   verbose  - (optional) If true, display h5 metadata
        % Results:
        %   imData
        %   params
        arguments
            dataPath    (1,:) char
            fileName    (1,:) char
            verbose     logical = false
        end

        fprintf('Loading acquisition parameters...\n')
        params = load(fullfile(dataPath, 'params.mat'));

        [~,~,ext] = fileparts(fileName);
        ext = lower(ext);
        if ~ext; fullPath = fullfile(dataPath, [fileName '.h5']);
        else; fullPath = fullfile(dataPath, fileName); end

        if verbose
            fprintf('HDF5 Metadata:\n')
            h5disp(fullPath);
        end

        fprintf('Reading chunked file...\n')
        tic;

        datasetName = '/ChunkedImage';
        output      = h5read(fullPath, datasetName);

        assert(isa(output, 'uint8'), 'largeFileIO: Dataset %s is not uint8.', datasetName);
        imData = output;

        readTime   = toc;
        readTime_m = floor(readTime / 60);
        readTime_s = mod(readTime, 60);
        fprintf('Read complete in %dm %.2fs\n', readTime_m, readTime_s);
    end
end
end