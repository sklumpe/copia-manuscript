function findUnitCell()


%neighMap to pos
neighCentName='lists/pos_1patch-Neigh.star';
pixS=5.9; 
neigCentscaleFact=1;
neighCentMapThr=0.02;
neighCentMap='maps/neigh/pos_1-patch.mrc';
centerPos='lists/pos_1-patch.star';
colors{1}=[1 0.922 0.701]; %C1==>C5 %flax
colors{2}=[1 0.922 0.701]; %C1==>C5 %flax
colors{3}=[1 0.188 0.1843]; %c5==>c5 %red
colors{4}=[0.466 0.655 0.850]; %c1==>c1 %blue
colors{5}=[0.466 0.655 0.850]; %c1==>c1 %blue


%unit cell
ucell.refVectId=3;
ucell.nr4Unit=12;
ucell.final4Unit=[1 3 7]; %min angle
ucell.outputName='uCell.mat';
ucell.vis='vis/ucell/pos_1-patch-lattice.bild';
ucell.vis2='vis/ucell/ucellOnNeighMap.bild';

%unit cell Extended
ucellext.refVectId=5;
ucellext.nr4Unit=12;
ucellext.final4Unit=[5 10]; %min angle
ucellext.outputName='uCell/uCellExt.mat';
ucellext.vis='vis/ucell/pos_1-patch-latticeExt.bild';
ucellext.vis2='vis/ucell/ucellOnNeighMapExt.bild';

neighMap2Pos(neighCentMap,centerPos,neighCentMapThr,neigCentscaleFact,neighCentName);
%minimal unitcell
calcUnitCell(ucell,neighCentName,pixS);
visUnitCell(ucell.outputName,centerPos,-1,colors(1:3),2,ucell.vis);
visUnitCell(ucell.outputName,[101 101 101],-1,colors(1:3),1,ucell.vis2);

%extended unitcell
calcUnitCell(ucellext,neighCentName,pixS);
visUnitCell(ucellext.outputName,centerPos,-1,colors(4:5),2,ucellext.vis);
visUnitCell(ucellext.outputName,[101 101 101],-1,colors(4:5),1,ucellext.vis2);



function visUnitCell(uVectName,baseList,idBase,colors,thick,outputName)

%uVectName='uCell.mat';

%baseList='patch.star';
%idBase=[36 107 306];
%idBase=1:51;

%gridOut='visGrid/latticePatch.bild';
%gridOut='visGrid/latticeUnit.bild';

load(uVectName);

if (ischar(baseList))
    st=tom_starread(baseList);
else
    st=baseList;
end

if (idBase==-1)
    idBase=1:length(st);
end

fid=fopen(outputName,'wt');
for i=1:length(idBase)
    ind=idBase(i);
    if (isstruct(st))
        pos1=[st(ind).rlnCoordinateX st(ind).rlnCoordinateY st(ind).rlnCoordinateZ];
    else
        pos1=st;
    end

    pos1_1=pos1+ uVect.vectOrg(1,:);
    pos1_2=pos1+ uVect.vectOrg(2,:);
    
    %fprintf(fid,'.color 1 0 0\n');
    fprintf(fid,'.color %s\n',num2str(colors{1}));
    fprintf(fid,'.arrow %f %f %f %f %f %f %f\n',pos1(1),pos1(2),pos1(3),pos1_1(1),pos1_1(2),pos1_1(3),thick);
    fprintf(fid,'.color %s\n',num2str(colors{2}));
    fprintf(fid,'.arrow %f %f %f %f %f %f %f\n',pos1(1),pos1(2),pos1(3),pos1_2(1),pos1_2(2),pos1_2(3),thick);
    
    if (size(uVect.vectOrg,1)>2)
        pos1_3=pos1+ uVect.vectOrg(3,:);
        fprintf(fid,'.color %s\n',num2str(colors{3}));
        fprintf(fid,'.arrow %f %f %f %f %f %f %f\n',pos1(1),pos1(2),pos1(3),pos1_3(1),pos1_3(2),pos1_3(3),thick);
    end
end
fclose(fid);



function calcUnitCell(ucell,neighCentName,pixS)

st=tom_starread(neighCentName);

for i=1:length(st)
    pos=[st(i).rlnCoordinateX st(i).rlnCoordinateY st(i).rlnCoordinateZ];
    vn(i,:)=pos;
    dist(i)=norm(vn(i,:)).*pixS;
    all(i,:)=[vn(i,:) dist(i) sign(vn(i,1)).*sign(vn(i,2)).*sign(vn(i,3))];
end
id=find(all(:,4)>0);
all=all(id,:);

allF=all(1:ucell.nr4Unit,:);
stF=st(1:ucell.nr4Unit);
ref=allF(ucell.refVectId,1:3);

for i=1:size(allF)
    v1=allF(i,1:3);
    ang(i)=atan2d(norm(cross(v1,ref)),dot(v1,ref)); 
    allF(i,6)= ang(i);
end
disp(' ');

unitFin=allF(ucell.final4Unit,:);

if size(unitFin,1)==3
    cmp={1 2; 1 3; 2 3 };

    for i=1:size(unitFin,1)
        v1=unitFin(cmp{i,1},1:3);
        v2=unitFin(cmp{i,2},1:3);
        angU(i)=atan2d(norm(cross(v1,v2)),dot(v1,v2));
    end

else
    v1=unitFin(1,1:3);
    v2=unitFin(2,1:3);
    angU(1)=atan2d(norm(cross(v1,v2)),dot(v1,v2));
    cmp={1 2};
end

uVect.vectOrg=unitFin(:,1:3);
uVect.len=unitFin(:,4);
uVect.ang=angU;
uVect.angIdx=cmp;
save(ucell.outputName,'uVect');




function neighMap2Pos(nm_name,dummyName,thr,scaleFact,outputName)

% nm_name='maps/neigh_patch.mrc';
% thr=0.02;
% scaleFact=1;
% outputName='lists/neig_patch_5.92.star';
% dummyName='lists/pos_1-patch.star';

v=tom_mrcread(nm_name); v=v.Value;
mid=floor(size(v)./2)+1;
dummy=tom_starread(dummyName);

vb=v>thr;
mt=v;%.*vb;
for i=1:13
    [pos(i,:),val(i),mt]=tom_peak(mt,20);
end

for i=1:13
    tmp=dummy(1);
    tmp.rlnCoordinateX=(pos(i,1)-mid(1)).*scaleFact;
    tmp.rlnCoordinateY=(pos(i,2)-mid(2)).*scaleFact;
    tmp.rlnCoordinateZ=(pos(i,3)-mid(3)).*scaleFact;
    stMap(i)=tmp;
end
tom_starwrite(outputName,stMap);


% figure; histogram(allAng);
%                
% disp(' ');