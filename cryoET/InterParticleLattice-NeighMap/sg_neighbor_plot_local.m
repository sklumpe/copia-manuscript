%% sg_neighbor_plot_local
% A function for reading in motivelist and generating a neighbor plot for
% a motivelist. When calculating the neighbor plot, only inter-subunit
% distances within each object is considered.
%
% WW 08-2018

%% Inputs


% Input motl
motl_nameWk = 'lists/pos_1-patch.star';
distMapName='lists/pos_1-distPairPatch.em';
classes=[1 1;1 5; 5 5];% 5 3; 1 3];
tomoID=1;

% Output root name
output_root = 'maps/neigh/';

% Plot parameters
%boxsize = 210;% %210; %200
boxsize = 200;% %210; %200
scaling = 1; %0.2

d=dir(motl_nameWk);
for i=1:length(d)
    genNeighMap([d(i).folder filesep d(i).name],boxsize,scaling,[output_root strrep(d(i).name,'.star','.mrc') ]);
end

class2Sym=[1 3 5];
distmat=tom_emread(distMapName);
distmat=distmat.Value;

for ic=1:size(classes,1)
        vT=[distmat(:,6)==tomoID]';
        clV=sort([distmat(:,4) distmat(:,5)],2);
        vc1=class2Sym(clV(:,1))==classes(ic,1);
        vc2=class2Sym(clV(:,2))==classes(ic,2);
        %vc1=class2Sym(distmat(:,4))==classes(ic,1);
        %vc2=class2Sym(distmat(:,5))==classes(ic,2);
        idx=find(vT.*vc1.*vc2);
        outputName=[output_root filesep 'neighConn-' num2str(classes(ic,1)) '-'  num2str(classes(ic,2)) '.mrc'];
        genConnectNeighMap(distmat(idx,:),boxsize,scaling,outputName);
end

function genConnectNeighMap(distmat,boxsize,scaling,outputName)
    
    cen = floor(boxsize/2)+1;
    pos1=distmat(:,8:10);
    pos2=distmat(:,11:13);
    posInMap=round((pos2-pos1)+cen);
    
    nplot=zeros(boxsize,boxsize,boxsize);
    for i=1:size(posInMap,1)
        nplot(posInMap(i,1),posInMap(i,2),posInMap(i,3))=nplot(posInMap(i,1),posInMap(i,2),posInMap(i,3))+1;
    end
    tom_mrcwrite(nplot,'name',outputName);

    
end

function genNeighMap(motl_name,boxsize,scaling,outputName)

% Read motl
motl = tom_starread(motl_name);
%st=tom_starread(star);

% Center of plot
cen = floor(boxsize/2)+1;

% Distance cutoff
d_cut = floor((boxsize-1)/2);


% Generate plots

% Initialize plot
nplot=zeros(boxsize,boxsize,boxsize);
% Parse object parameters
n_pos = length(motl);
pos = zeros(3,n_pos);
pos(1,:) = [motl(:).rlnCoordinateX];
pos(2,:) = [motl(:).rlnCoordinateY];
pos(3,:) = [motl(:).rlnCoordinateZ];
pos = pos.*scaling;


% Determine neighbors for each position
for j = 1:n_pos

    % Calculate Euclidean distances
    dist = sg_pairwise_dist(pos(:,j),pos);

    % Threshold by distance cutoff
    d_idx = dist <= d_cut;

    % Parse positions and shift to center
    temp_pos = pos(:,d_idx) - repmat(pos(:,j),[1,sum(d_idx)]);

    % Generate rotation matrix
    rmat=eye(3);

    % Rotate positions
    rpos = round(rmat*temp_pos)+cen;

    % Add positions to map
    for k = 1:size(rpos,2)
        nplot(rpos(1,k),rpos(2,k),rpos(3,k)) = nplot(rpos(1,k),rpos(2,k),rpos(3,k)) + 1;
    end

end

% Normalize central peak
% nplot(cen,cen,cen) = 0;
% nplot(cen,cen,cen) = max(nplot(:));
nplot = nplot./nplot(cen,cen,cen);
nplot(cen,cen,cen) = 0;
nplot(cen,cen,cen) = max(nplot(:));

% Write output
sg_volume_write(outputName,nplot);


end
