function v_sm = sm_track(v,sm_len,seg_idx)
v_sm = nan(size(v));
for iS=1:size(seg_idx,1)
    v_sm(seg_idx(iS,1):seg_idx(iS,2),1) = smooth(v(seg_idx(iS,1):seg_idx(iS,2),1),sm_len);
    v_sm(seg_idx(iS,1):seg_idx(iS,2),2) = smooth(v(seg_idx(iS,1):seg_idx(iS,2),2),sm_len);
    v_sm(seg_idx(iS,1):seg_idx(iS,2),3) = smooth(v(seg_idx(iS,1):seg_idx(iS,2),3),sm_len);
end