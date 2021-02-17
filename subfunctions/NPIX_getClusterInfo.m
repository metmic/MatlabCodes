function info = NPIX_getClusterInfo(dire)

cluster_info = tdfread ([dire, '\cluster_info.tsv']); % [DESKTOP]
%cluster_info = tdfread ([dire, '/cluster_info.tsv']); % [LAPTOP]

sID = false(size(cluster_info.id));
for i = 1:length(cluster_info.sh)
%for i = 1:length(cluster_info.shank)
    if strcmp(cluster_info.group(i,:), 'good ') % high contrast data 20191011 and all more recent versions
    %if strcmp(cluster_info.group(i,:), 'good    ') % low contrast data 20191011 and other earlier softwear sorting?
    %if strcmp(cluster_info.group(i,:), 'good')
        sID(i) = true;
    end
end


info.id = cluster_info.id(sID);
info.channel = cluster_info.ch(sID);
%info.channel = cluster_info.channel(sID);
info.depth = cluster_info.depth(sID);
info.firing_rate = cluster_info.fr(sID);
%info.firing_rate = cluster_info.firing_rate(sID);
info.n_spikes = cluster_info.n_spikes(sID);

if isfield(cluster_info,'sorted')
    info.sorted = cluster_info.sorted(sID);
end

end