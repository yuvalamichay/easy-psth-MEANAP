% %% VERSION 4
% % Adds bslshuff to outputs and supports a 1ms time step
% function [kern,hcoeffs,hcoeffs2D,bslshuff] = hcoeff(peristim,rawdata,mfr,nstims,rsp,fast,time_step_ms)
% 
% % *********************************************************************************************
% % *********************************************************************************************
% % *******   This code was authored by Michael Hill (buffalohill@me.com), it may not     *******
% % *******   be used as a whole or in parts without the written permission of Michael    *******
% % *******   Hill and any use of this code (in whole or in parts) must reference the     *******
% % *******   paper in which the h-coefficient method was first introduced (contact me    *******
% % *******   to inquire about appropriate referencing formats for this code).             ******
% % *********************************************************************************************
% % *********************************************************************************************
% 
% dbstop if error
% 
% % --- Set default for time_step_ms if not provided (for standalone use) ---
% if nargin < 7
%     time_step_ms = 1; % Default to 1 ms resolution
% end
% 
% % --- HARDCODED TIME STEP ---
% % This code is optimized for a 1ms time step
% time_step_ms = 1;  % Force 1ms time step
% time_step_s = time_step_ms / 1000;
% 
% %% Generate the shuffled responses
% rspdur = rsp(1,2)-rsp(1,1);                                                 % Response duration in ms
% bslshuff = randshuff(nstims,rspdur,rawdata,mfr,fast);                       % Generate your shuffled responses
% 
% %% Generate a smooth PSTH from your response
% corrkern = ones(nstims,1)*NaN;
% for j = 1:nstims                                                            % Prep for the adaptive kernel analysis
%     if j == 1
%         forkern = peristim{1,1};
%     else
%         forkern = cat(2, forkern, peristim{j,1});
%     end
%     corrkern(j,1) = length(peristim{j,1});                                  % Preparation of the factor by which you need to multipy the ssvkernel output, as it is normalized to the mean firing rate
% end
% 
% % Create time vector with 1ms resolution
% psth_time_vector_s = (rsp(1)/1000) : (time_step_s) : (rsp(2)/1000);
% [kern(2,:),kern(1,:)] = ssvkernel(forkern.*10^(-3), psth_time_vector_s);     % Smooth your PSTH
% 
% kern(2,:) = kern(2,:) * (squeeze(nanmean(corrkern)));                       % Reverse the normaliztion applied by ssvkernel and transform to Hz.
% if isempty(kern(2,~isnan(kern(2,:))))                                       % Check for the presence of any firing rate whatsoever
%     kern(2,:) = zeros(1,length(kern));
% end
% 
% %% Find your local extrema
% % Find index of the response window in the PSTH
% fromrsp = 1;
% torsp = length(psth_time_vector_s);
% indmax = find(kern(2,fromrsp:torsp) == max(kern(2,fromrsp:torsp)),1,'first') + (fromrsp - 1);
% indmin = find(kern(2,fromrsp:torsp) == min(kern(2,fromrsp:torsp)),1,'first') + (fromrsp - 1);
% if isempty(indmax) || isempty(indmin)
%     indmax = NaN;
%     indmin = NaN;
% end
% 
% %% Calculate the h-coefficient
% [hcoeff_max,hcoeff_min,hcoeff_max2D,hcoeff_min2D] = hcoeff_gen(bslshuff,kern(2,:),indmin,indmax,mfr);
% hcoeffs = [hcoeff_max,hcoeff_min];
% hcoeffs2D = [hcoeff_max2D;hcoeff_min2D];
% 
% end
% 
% % --- SUB-FUNCTION: randshuff ---
% function [bslshuff] = randshuff(nstims,rspdur,rawdata,mfr,fast,nshuff)
% dbstop if error
% if nargin < 5, fast = 1; end
% if nargin < 6, nshuff = 100; end
% 
% % Optimized for 1ms time step
% time_step_ms = 1;
% time_step_s = time_step_ms / 1000;
% 
% kernshuff = randperm(nshuff); shuffsy = zeros(nstims,1) * NaN;
% nsegs = floor(((floor(rawdata(end))-2000)-(ceil(rawdata(1,1))))/(rspdur*10));
% if fast
%     % Use ceil to ensure we process all shuffles (e.g., 100/10 = 10 loops)
%     for i = 1:ceil(nshuff/10)
%         clc; disp([num2str(((i-1)*10)+1),' to ',num2str(i*10),' out of ',num2str(nshuff),' shuffles being processed...'])
%         permind1 = randperm(nsegs*10); permind2 = randperm(nsegs);
%         for j = 1:nstims
%             shuffsy(j,1) = length(rawdata((rawdata > ((permind1(1,j)-1)*rspdur) & (rawdata <= (((permind1(1,j)-1)*rspdur)+rspdur)))));
%             if j == 1, forkern = (rawdata((rawdata > ((permind2(1,j)-1)*(rspdur*10))) & (rawdata <= ((permind2(1,j)-1)*(rspdur*10))+(rspdur*10)+2000)) - ((permind2(1,j)-1)*(rspdur*10)));
%             else, forkern = cat(2,forkern,(rawdata((rawdata > ((permind2(1,j)-1)*(rspdur*10))) & (rawdata <= ((permind2(1,j)-1)*(rspdur*10))+(rspdur*10)+2000)) - ((permind2(1,j)-1)*(rspdur*10)))); end
%         end
% 
%         % Use 1ms time step for baseline kernel
%         prekern_time_vec_s = (time_step_s:time_step_s:(((rspdur*10)+2000)/1000));
%         [prekern{i}(2,:),prekern{i}(1,:)] = ssvkernel(forkern.*10^(-3), prekern_time_vec_s);
% 
%         prekern{i}(2,:) = prekern{i}(2,:) .* 10 .* (squeeze(nanmean(shuffsy))*(1000/rspdur));
%         for j = 1:10
%             k = ((i-1)*10)+j;
%             if k > nshuff, continue; end % Don't create more shuffles than requested
% 
%             bslshuff.kern{kernshuff(k)}(1,:) = (time_step_s:time_step_s:rspdur/1000);
% 
%             % Calculate indices for slicing the pre-generated kernel
%             start_idx = round((j-1)*(rspdur/10) / time_step_s) + 1;
%             end_idx = round(j*(rspdur/10) / time_step_s);
% 
%             % Ensure indices are within bounds
%             start_idx = max(1, start_idx);
%             end_idx = min(length(prekern{i}(2,:)), end_idx);
% 
%             % Ensure the resulting vector has the correct length
%             temp_kern = prekern{i}(2,start_idx:end_idx);
%             target_len = length(bslshuff.kern{kernshuff(k)}(1,:));
%             if length(temp_kern) > target_len
%                 temp_kern = temp_kern(1:target_len);
%             elseif length(temp_kern) < target_len
%                 temp_kern(end+1:target_len) = 0; % Pad with zeros if needed
%             end
%             bslshuff.kern{kernshuff(k)}(2,:) = temp_kern;
%         end
%     end
% end
% clc
% bslshuffmaxes = zeros(10001,nshuff); bslshuffminses = bslshuffmaxes;
% for a = 1:nshuff
%     kern = (bslshuff.kern{a}(2,:)./mfr)-1;
%     indmax = find(kern == max(kern),1,'first');
%     if ~isempty(indmax) && kern(1,indmax) > 0
%         firstdec = norm(floor(norm(kern(1,indmax))*10)/10);
%         u = indmax(1); v = u;
%         for b = 1:round(firstdec*10)+1
%             while u > 1 && kern(1,u) > roundn(firstdec-((b-1)*0.1),-1), u = u-1; end
%             while v < length(kern) && kern(1,v) > roundn(firstdec-((b-1)*0.1),-1), v = v+1; end
%             % Modified: Use 0.001 for 1ms time step (instead of 0.01 for 10ms)
%             bslshuffmaxes(end-(round(firstdec*10)+1)+b,a) = trapz((kern(1,u:v))-roundn(firstdec-((b-1)*0.1),-1)).*0.001;
%         end
%     end
%     indmin = find(kern == min(kern),1,'first');
%     if ~isempty(indmax) && kern(1,indmin) < 0
%         firstdec = norm(floor(norm(kern(1,indmin))*10)/10);
%         u = indmin(1); v = u;
%         for b = 1:round(firstdec*10)+1
%             while u > 1 && kern(1,u) < -roundn(firstdec-((b-1)*0.1),-1), u = u-1; end
%             while v < length(kern) && kern(1,v) < -roundn(firstdec-((b-1)*0.1),-1), v = v+1; end
%             % Modified: Use 0.001 for 1ms time step (instead of 0.01 for 10ms)
%             bslshuffminses(end-(round(firstdec*10)+1)+b,a) = norm(trapz(kern(1,u:v) + roundn(firstdec-((b-1)*0.1),-1)).*0.001);
%         end
%     end
% end
% bslshuffmaxes = bslshuffmaxes(2:end,:) - bslshuffmaxes(1:end-1,:);
% bslshuffminses = bslshuffminses(2:end,:) - bslshuffminses(1:end-1,:);
% bslshuff.maxes = max(bslshuffmaxes,[],2);
% bslshuff.minses = max(bslshuffminses,[],2);
% end
% 
% function [hcoefficient_max,hcoefficient_min,hcoefficient_max2D,hcoefficient_min2D] = hcoeff_gen(bslshuff,kern,indmin,indmax,normit)
% maxes = zeros(10001,1); minses = maxes;
% kern = (kern/normit)-1;
% 
% % Optimized for 1ms time step
% time_step_ms = 1;
% time_step_s = time_step_ms / 1000;
% 
% if ~isnan(indmax) && kern(1,indmax) > 0
%     firstdec = norm(floor(norm(kern(1,indmax))*10)/10);
%     u = indmax(1); v = u;
%     for b = 1:round(firstdec*10)+1
%         while u > 1 && kern(1,u) > roundn(firstdec-((b-1)*0.1),-1), u = u-1; end
%         while v < length(kern) && kern(1,v) > roundn(firstdec-((b-1)*0.1),-1), v = v+1; end
%         % Modified: Use 0.01 for 1ms time step (instead of 0.1 for 10ms)
%         % Note: The factor here should be 10x larger than in randshuff
%         maxes(end-(round(firstdec*10)+1)+b,1) = trapz((kern(1,u:v))-roundn(firstdec-((b-1)*0.1),-1)).*0.01;
%     end
% end
% if ~isnan(indmin) && kern(1,indmin) < 0
%     firstdec = norm(floor(norm(kern(1,indmin))*10)/10);
%     u = indmin(1); v = u;
%     for b = 1:round(firstdec*10)+1
%         while u > 1 && kern(1,u) < -roundn(firstdec-((b-1)*0.1),-1), u = u-1; end
%         while v < length(kern) && kern(1,v) < -roundn(firstdec-((b-1)*0.1),-1), v = v+1; end
%         % Modified: Use 0.01 for 1ms time step (instead of 0.1 for 10ms)
%         minses(end-(round(firstdec*10)+1)+b,1) = norm(trapz(kern(1,u:v)+roundn(firstdec-((b-1)*0.1),-1)).*0.01);
%     end
% end
% maxes = maxes(2:end,1) - maxes(1:end-1,:);
% minses = minses(2:end,1) - minses(1:end-1,:);
% hcoefficient_max = sum(maxes > bslshuff.maxes)/sum(bslshuff.maxes > 0);
% hcoefficient_min = sum(minses > bslshuff.minses)/sum(bslshuff.minses > 0);
% prov = (nnz(maxes)-nnz(bslshuff.maxes));
% if prov > 0, hcoefficient_max2D(1,1) = prov; else, hcoefficient_max2D(1,1) = 0; end
% hcoefficient_max2D(1,2) = sum(maxes > bslshuff.maxes) - hcoefficient_max2D(1,1);
% hcoefficient_max2D(1,3) = sum(bslshuff.maxes > 0);
% prov = (nnz(minses)-nnz(bslshuff.minses));
% if prov > 0, hcoefficient_min2D(1,1) = prov; else, hcoefficient_min2D(1,1) = 0; end
% hcoefficient_min2D(1,2) = sum(minses > bslshuff.minses) - hcoefficient_min2D(1,1);
% hcoefficient_min2D(1,3) = sum(bslshuff.minses > 0);
% end

%% VERSION 6
% Adds bslshuff to outputs and supports a 1ms time step
function [kern,hcoeffs,hcoeffs2D,bslshuff] = hcoeff(peristim,rawdata,mfr,nstims,rsp,fast,time_step_ms)

% *********************************************************************************************
% *********************************************************************************************
% *******   This code was authored by Michael Hill (buffalohill@me.com), it may not     *******
% *******   be used as a whole or in parts without the written permission of Michael    *******
% *******   Hill and any use of this code (in whole or in parts) must reference the     *******
% *******   paper in which the h-coefficient method was first introduced (contact me    *******
% *******   to inquire about appropriate referencing formats for this code).             ******
% *********************************************************************************************
% *********************************************************************************************

dbstop if error

% --- Set default for time_step_ms if not provided (for standalone use) ---
if nargin < 7
    time_step_ms = 1; % Default to 1 ms resolution
end

% --- HARDCODED TIME STEP ---
% This code is optimized for a 1ms time step
time_step_ms = 1;  % Force 1ms time step
time_step_s = time_step_ms / 1000;

%% Generate the shuffled responses
rspdur = rsp(1,2)-rsp(1,1);                                                 % Response duration in ms
bslshuff = randshuff(nstims,rspdur,rawdata,mfr,fast);                       % Generate your shuffled responses

%% Generate a smooth PSTH from your response
corrkern = ones(nstims,1)*NaN;
for j = 1:nstims                                                            % Prep for the adaptive kernel analysis
    if j == 1
        forkern = peristim{1,1};
    else
        forkern = cat(2, forkern, peristim{j,1});
    end
    corrkern(j,1) = length(peristim{j,1});                                  % Preparation of the factor by which you need to multipy the ssvkernel output, as it is normalized to the mean firing rate
end

% Create time vector with 1ms resolution, starting at 0
psth_time_vector_s = (0) : time_step_s : ((rsp(2)-rsp(1))/1000);
[kern(2,:),kern(1,:)] = ssvkernel(forkern.*10^(-3), psth_time_vector_s);     % Smooth your PSTH

kern(2,:) = kern(2,:) * (squeeze(nanmean(corrkern)));                       % Reverse the normaliztion applied by ssvkernel and transform to Hz.
if isempty(kern(2,~isnan(kern(2,:))))                                       % Check for the presence of any firing rate whatsoever
    kern(2,:) = zeros(1,length(kern));
end

%% Find your local extrema
% Find index of the response window in the PSTH
fromrsp = 1;
torsp = length(psth_time_vector_s);
indmax = find(kern(2,fromrsp:torsp) == max(kern(2,fromrsp:torsp)),1,'first') + (fromrsp - 1);
indmin = find(kern(2,fromrsp:torsp) == min(kern(2,fromrsp:torsp)),1,'first') + (fromrsp - 1);
if isempty(indmax) || isempty(indmin)
    indmax = NaN;
    indmin = NaN;
end

%% Calculate the h-coefficient
[hcoeff_max,hcoeff_min,hcoeff_max2D,hcoeff_min2D] = hcoeff_gen(bslshuff,kern(2,:),indmin,indmax,mfr);
hcoeffs = [hcoeff_max,hcoeff_min];
hcoeffs2D = [hcoeff_max2D;hcoeff_min2D];

end

% --- SUB-FUNCTION: randshuff ---
function [bslshuff] = randshuff(nstims,rspdur,rawdata,mfr,fast,nshuff)
dbstop if error
if nargin < 5, fast = 1; end
if nargin < 6, nshuff = 500; end

% Optimized for 1ms time step
time_step_ms = 1;
time_step_s = time_step_ms / 1000;

kernshuff = randperm(nshuff); shuffsy = zeros(nstims,1) * NaN;
nsegs = floor(((floor(rawdata(end))-2000)-(ceil(rawdata(1,1))))/(rspdur*10));
if fast
    % Use ceil to ensure we process all shuffles (e.g., 100/10 = 10 loops)
    for i = 1:ceil(nshuff/10)
        clc; disp([num2str(((i-1)*10)+1),' to ',num2str(i*10),' out of ',num2str(nshuff),' shuffles being processed...'])
        permind1 = randperm(nsegs*10); permind2 = randperm(nsegs);
        for j = 1:nstims
            shuffsy(j,1) = length(rawdata((rawdata > ((permind1(1,j)-1)*rspdur) & (rawdata <= (((permind1(1,j)-1)*rspdur)+rspdur)))));
            if j == 1, forkern = (rawdata((rawdata > ((permind2(1,j)-1)*(rspdur*10))) & (rawdata <= ((permind2(1,j)-1)*(rspdur*10))+(rspdur*10)+2000)) - ((permind2(1,j)-1)*(rspdur*10)));
            else, forkern = cat(2,forkern,(rawdata((rawdata > ((permind2(1,j)-1)*(rspdur*10))) & (rawdata <= ((permind2(1,j)-1)*(rspdur*10))+(rspdur*10)+2000)) - ((permind2(1,j)-1)*(rspdur*10)))); end
        end
        
        % Use 1ms time step for baseline kernel, starting at 0
        % FIXED: Removed the subtraction of time_step_s from the end time
        prekern_time_vec_s = (0:time_step_s:(((rspdur*10)+2000)/1000));
        [prekern{i}(2,:),prekern{i}(1,:)] = ssvkernel(forkern.*10^(-3), prekern_time_vec_s);
        
        prekern{i}(2,:) = prekern{i}(2,:) .* 10 .* (squeeze(nanmean(shuffsy))*(1000/rspdur));
        for j = 1:10
            k = ((i-1)*10)+j;
            if k > nshuff, continue; end % Don't create more shuffles than requested
            
            % FIXED: Removed the subtraction of time_step_s from the end time
            bslshuff.kern{kernshuff(k)}(1,:) = (0:time_step_s:rspdur/1000);
            
            % Calculate indices for slicing the pre-generated kernel
            start_idx = round((j-1)*(rspdur/10) / time_step_s) + 1;
            end_idx = round(j*(rspdur/10) / time_step_s);
            
            % Ensure indices are within bounds
            start_idx = max(1, start_idx);
            end_idx = min(length(prekern{i}(2,:)), end_idx);
            
            % Ensure the resulting vector has the correct length
            temp_kern = prekern{i}(2,start_idx:end_idx);
            target_len = length(bslshuff.kern{kernshuff(k)}(1,:));
            if length(temp_kern) > target_len
                temp_kern = temp_kern(1:target_len);
            elseif length(temp_kern) < target_len
                temp_kern(end+1:target_len) = 0; % Pad with zeros if needed
            end
            bslshuff.kern{kernshuff(k)}(2,:) = temp_kern;
        end
    end
end
clc
bslshuffmaxes = zeros(10001,nshuff); bslshuffminses = bslshuffmaxes;
for a = 1:nshuff
    kern = (bslshuff.kern{a}(2,:)./mfr)-1;
    indmax = find(kern == max(kern),1,'first');
    if ~isempty(indmax) && kern(1,indmax) > 0
        firstdec = norm(floor(norm(kern(1,indmax))*10)/10);
        u = indmax(1); v = u;
        for b = 1:round(firstdec*10)+1
            while u > 1 && kern(1,u) > roundn(firstdec-((b-1)*0.1),-1), u = u-1; end
            while v < length(kern) && kern(1,v) > roundn(firstdec-((b-1)*0.1),-1), v = v+1; end
            % Modified: Use 0.001 for 1ms time step (instead of 0.01 for 10ms)
            bslshuffmaxes(end-(round(firstdec*10)+1)+b,a) = trapz((kern(1,u:v))-roundn(firstdec-((b-1)*0.1),-1)).*0.001;
        end
    end
    indmin = find(kern == min(kern),1,'first');
    if ~isempty(indmax) && kern(1,indmin) < 0
        firstdec = norm(floor(norm(kern(1,indmin))*10)/10);
        u = indmin(1); v = u;
        for b = 1:round(firstdec*10)+1
            while u > 1 && kern(1,u) < -roundn(firstdec-((b-1)*0.1),-1), u = u-1; end
            while v < length(kern) && kern(1,v) < -roundn(firstdec-((b-1)*0.1),-1), v = v+1; end
            % Modified: Use 0.001 for 1ms time step (instead of 0.01 for 10ms)
            bslshuffminses(end-(round(firstdec*10)+1)+b,a) = norm(trapz(kern(1,u:v) + roundn(firstdec-((b-1)*0.1),-1)).*0.001);
        end
    end
end
bslshuffmaxes = bslshuffmaxes(2:end,:) - bslshuffmaxes(1:end-1,:);
bslshuffminses = bslshuffminses(2:end,:) - bslshuffminses(1:end-1,:);
bslshuff.maxes = max(bslshuffmaxes,[],2);
bslshuff.minses = max(bslshuffminses,[],2);
end

function [hcoefficient_max,hcoefficient_min,hcoefficient_max2D,hcoefficient_min2D] = hcoeff_gen(bslshuff,kern,indmin,indmax,normit)
maxes = zeros(10001,1); minses = maxes;
kern = (kern/normit)-1;

% Optimized for 1ms time step
time_step_ms = 1;
time_step_s = time_step_ms / 1000;

if ~isnan(indmax) && kern(1,indmax) > 0
    firstdec = norm(floor(norm(kern(1,indmax))*10)/10);
    u = indmax(1); v = u;
    for b = 1:round(firstdec*10)+1
        while u > 1 && kern(1,u) > roundn(firstdec-((b-1)*0.1),-1), u = u-1; end
        while v < length(kern) && kern(1,v) > roundn(firstdec-((b-1)*0.1),-1), v = v+1; end
        % Modified: Use 0.01 for 1ms time step (instead of 0.1 for 10ms)
        % Note: The factor here should be 10x larger than in randshuff
        maxes(end-(round(firstdec*10)+1)+b,1) = trapz((kern(1,u:v))-roundn(firstdec-((b-1)*0.1),-1)).*0.01;
    end
end
if ~isnan(indmin) && kern(1,indmin) < 0
    firstdec = norm(floor(norm(kern(1,indmin))*10)/10);
    u = indmin(1); v = u;
    for b = 1:round(firstdec*10)+1
        while u > 1 && kern(1,u) < -roundn(firstdec-((b-1)*0.1),-1), u = u-1; end
        while v < length(kern) && kern(1,v) < -roundn(firstdec-((b-1)*0.1),-1), v = v+1; end
        % Modified: Use 0.01 for 1ms time step (instead of 0.1 for 10ms)
        minses(end-(round(firstdec*10)+1)+b,1) = norm(trapz(kern(1,u:v)+roundn(firstdec-((b-1)*0.1),-1)).*0.01);
    end
end
maxes = maxes(2:end,1) - maxes(1:end-1,:);
minses = minses(2:end,1) - minses(1:end-1,:);
hcoefficient_max = sum(maxes > bslshuff.maxes)/sum(bslshuff.maxes > 0);
hcoefficient_min = sum(minses > bslshuff.minses)/sum(bslshuff.minses > 0);
prov = (nnz(maxes)-nnz(bslshuff.maxes));
if prov > 0, hcoefficient_max2D(1,1) = prov; else, hcoefficient_max2D(1,1) = 0; end
hcoefficient_max2D(1,2) = sum(maxes > bslshuff.maxes) - hcoefficient_max2D(1,1);
hcoefficient_max2D(1,3) = sum(bslshuff.maxes > 0);
prov = (nnz(minses)-nnz(bslshuff.minses));
if prov > 0, hcoefficient_min2D(1,1) = prov; else, hcoefficient_min2D(1,1) = 0; end
hcoefficient_min2D(1,2) = sum(minses > bslshuff.minses) - hcoefficient_min2D(1,1);
hcoefficient_min2D(1,3) = sum(bslshuff.minses > 0);
end