function mic_temp = mic_nan_convert(mic,mic_idx,n_channels)
    mic_temp = NaN(n_channels,size(mic,2));
    mic_temp(mic_idx,:) = mic;
end