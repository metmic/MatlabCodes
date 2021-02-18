function info = NPIX_getClusterInfo(dire)
% set correct delimiter
if ispc
    pLoc = '\';
elseif ismac
    pLoc = '/';
end

cluster_info = tdfread ([dire, pLoc, 'cluster_info.tsv']);

sID = false(size(cluster_info.id));
for i = 1:length(cluster_info.sh)
    if strcmp(cluster_info.group(i,:), 'good ')
        sID(i) = true;
    end
end

info.id = cluster_info.id(sID);
info.channel = cluster_info.ch(sID);
info.depth = cluster_info.depth(sID);
info.firing_rate = cluster_info.fr(sID);
info.n_spikes = cluster_info.n_spikes(sID);

if isfield(cluster_info,'sorted')
    info.sorted = cluster_info.sorted(sID);
end

end