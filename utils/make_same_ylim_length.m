function make_same_ylim_length(axes_hs)
ylims = [];
for i = 1:length(axes_hs)
    ylims(i,:) = get(axes_hs(i), 'ylim');
end

lens = ylims(:,2) - ylims(:,1);
max_len = max(lens);

for i = 1:length(axes_hs)
    pad_sz = (max_len - lens(i)) / 2;
    new_ylim = [ylims(i,1) - pad_sz, ylims(i,2) + pad_sz];
    ylim(axes_hs(i), new_ylim);
end

end