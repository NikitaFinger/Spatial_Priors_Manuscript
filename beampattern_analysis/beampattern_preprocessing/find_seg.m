
function seg_idx = find_seg(pos,sm_len)
% Find continous segment in the trajectories with small gaps filled
notnan = ~isnan(pos(:,1));
idx_nan = find(diff(notnan)~=0)+1;

g = normpdf(-sm_len:sm_len,0,sm_len);
w = conv(double(notnan),g,'same');
idx = find(diff(w~=0)~=0)+1;
idx_up = find(diff(w~=0)>0)+1;
idx_dn = find(diff(w~=0)<0)+1;

if ~isempty(idx)
    [~,iconv] = min(abs(repmat(idx',length(idx_nan),1)-repmat(idx_nan,1,length(idx))),[],1);
    if isempty(idx_up) && ~isempty(idx_dn)
        seg_idx = [1;idx_nan(iconv)];
    elseif ~isempty(idx_up) && isempty(idx_dn)
        seg_idx = [idx_nan(iconv);size(pos,1)];
    else
        if idx_up(1)~=idx(1)  % if the first up edge not at the beginning
            seg_idx = [1;idx_nan(iconv)];
        else
            seg_idx = idx_nan(iconv);
        end
    end
else
    seg_idx = [1,size(pos,1)];
end
if mod(length(seg_idx),2)~=0  % if the last down edge not at the end
    seg_idx = [seg_idx;length(notnan)];
end
seg_idx = reshape(seg_idx,2,[])';