% interrogation from netcdf sample files

% chains
chains = [5];
% burn-in
burn_in = [400000, 400000, 400000, 400000,  400000, 400000, 400000, 400000];
% thin
thin = 100;
% netcdf dir
dir = 'Results';
% nsamples
nsamples = 1200000;

% load interrogation points
% query = load('query.dat');
query = [2.8 9];

nthin = 0;
vs = zeros(length(chains)*(nsamples-min(burn_in))/thin,size(query,1));
stepType = zeros(nsamples,1);
accept = zeros(nsamples,1);
index = zeros(nsamples,1);
ncells = zeros(nsamples,1);
misfit = zeros(nsamples,1);
umisfit = zeros(nsamples,1);
likelihood = zeros(nsamples,1);
coord = zeros(2,nsamples);
values = zeros(3,nsamples);
noise0 = zeros(nsamples,1);
noise1 = zeros(nsamples,1);
% loop over chains
for i = chains
    % read initial points and parameters from netcdf file
    initialSample_ = dlmread(['../' num2str(i) '/' dir '/ini/InitialSample_' num2str(i) '.txt']);
    ncells = initialSample_(1,1);
    points = initialSample_(2:ncells+1,1:2);
    parameters = initialSample_(2:ncells+1,3:5);
    points = points'; parameters = parameters';
    ncfile = ['../' num2str(i) '/' dir '/samples_' num2str(i) '.out'];
    fileID = fopen(ncfile,'r');
    disp(ncfile);
    
    % update model sample by sample
    for j = 1 : nsamples
        % read one sample
        stepType(j)=fread(fileID,1,'int');
        accept(j) = fread(fileID,1,'char');
        index(j) = fread(fileID,1,'int');
        ncells(j) = fread(fileID,1,'int64');
        misfit(j) = fread(fileID,1,'double');
        umisfit(j) = fread(fileID,1,'double');
        likelihood(j) = fread(fileID,1,'double');
        coord(:,j) = fread(fileID,2,'double');
        values(:,j) = fread(fileID,3,'double');
        noise0(j) = fread(fileID,1,'double');
        noise1(j) = fread(fileID,1,'double');
        
        if(accept(j))
            switch stepType(j)
                case 1
                    points(:,ncells(j)) = coord(:,j);
                    parameters(:,ncells(j)) = values(:,j);
                case 2
                    points_copy = points;
                    parameters_copy = parameters;
                    points(:,index(j):ncells(j)) = points_copy(:,index(j)+1:ncells(j)+1);
                    parameters(:,index(j):ncells(j)) = parameters_copy(:,index(j)+1:ncells(j)+1);
                case 3
                    points(:,index(j)) = coord(:,j);
                case 4
                    parameters(:,index(j)) = values(:,j);
            end
        end
        % if after burn-in, sampling
        if( j > burn_in(i) && mod(j-burn_in(i),thin) == 0)
            nthin = nthin + 1;
            % extract parameters for query point
            % create a kd tree
            ns = createns(points(:,1:ncells(j))','nsmethod','kdtree');
            % query
            idx = knnsearch(ns,query,'k',1);
            vs(nthin,:) = parameters(1,idx);
%             query_pm{nthin} = parameters(:,idx);
        end
            
    end
    fclose(fileID);
end

%% plot
figure; plot(vs,'.');
%ylim([0.415 0.44]);
%for i = 1 : length(vs)
    figure;
    histogram(vs,20);
    axis([0.415 0.44 0 400]);
%    print(['hist' num2str(i)],'-dpng');
%end