function cleanPicklists()

%st.io.outputFolder='runs/testHT3/';
%st.io.volumeList='vols/f2R58GP1Cut.star';

st.io.outputFolder='runs/f2All/';
st.io.volumeList='vols/f2R58.star'; %vols/f2R58tomo3Cut.star';
st.io.cmbMotlName='runs/f2All/cmbLists/motlAll.em';
st.io.tag={'C1';'C3';'C5'};
st.io.szInfo=[272 272 272];
st.io.tomoNumToVolName={'Position_1.mrc','Position_36.mrc','Position_37.mrc'};
st.clean.do=1;

%% Cross correlatin cleaning
st.clean.minCCGlobal=0.10;  %0.13; 
st.clean.minCCLocalDiff=3.0;  %2.5; %in std
st.clean.minCCLocalRadius=33;

%% Distance cleaning
st.clean.minDistNeighbour=20.1;  % in pixel maybe 20
st.clean.maxDistCenterABS=77; %76; %in pixel mean~75 %check for elliptical
st.clean.maxDistCenterLocalDiff=2.1; %in std

%% Anguglar Cleaning
st.clean.maxAngDistToNormalABS=40; %28;  %25 suggested by hist%in deg

st.reg.do=1;
st.reg.neighRadius=35;

removeOutlier(st.io,st.clean);
regNeigh(st.io,st.reg);


function regNeigh(io,reg)

motlAll=tom_emread(strrep(io.cmbMotlName,'.em','_clean.em')); motlAll=motlAll.Value;

starAll=tom_starread(io.volumeList);
allFeat=motlAll(6,:);
featU=unique(allFeat);
allMotlReg=[];
for i=1:length(featU)
    idx=find(allFeat==featU(i));
    motlOneFeat=motlAll(:,idx);
    motlRegOneFeat=correctByNeighClass(motlOneFeat,reg);
    allMotlReg=cat(2,allMotlReg,motlRegOneFeat);
end
tom_emwrite(strrep(io.cmbMotlName,'.em','_reg.em'),allMotlReg);
%motlToStar(motlName,io)
motlToStar(strrep(io.cmbMotlName,'.em','_reg.em'),io);



function motlClean=correctByNeighClass(motl,reg)

motlClean=motl;
for i=1:size(motl,2)
     [~,motlNeigh,dist]=findNeighbours(motl,i,reg.neighRadius);
    nrNeigh=size(motlNeigh,2);
    meanDist(i)=mean(dist);
    stdDist(i)=std(dist);
    nrC1=length(find(motlNeigh(20,:)==1));
    nrC3=length(find(motlNeigh(20,:)==2));
    nrC5=length(find(motlNeigh(20,:)==3));
   
    scoreC1=abs(nrC1-3)+abs(nrC3-2)+abs(nrC5-1);
    scoreC3=abs(nrC1-6)+nrC3+nrC5;
    scoreC5=abs(nrC1-5)+nrC3+nrC5;
    [a,b]=min([scoreC1 scoreC3 scoreC5]);
     if (a==0 && b~=motlClean(20,i) )
         if (b==3)
           if (meanDist(i)<26)
                motlClean(20,i)=b;
           end
        else
             motlClean(20,i)=b;
        end
    end
%     if (nrC5==0 && nrC3==0 && meanDist(i)>26 && nrNeigh>2)
%         motlClean(20,i)=2;
%     end
%     if (nrC5==0 && nrC3==0 && meanDist(i)<25 && nrNeigh>3)
%         motlClean(20,i)=3;
%     end

    %disp([num2str(i) ' nrNeigh ' num2str(nrNeigh) ' mean Dist: '  num2str(meanDist(i)) ' std Dist: ' num2str(stdDist(i))]);
  %  tom_emwrite('tmpDebug/neigh.em',motlNeigh);
    if (motlClean(20,i)==2);
       % disp(' ');
    end
        %disp(' ');
end



function removeOutlier(io,clean)


motlAll=tom_emread(io.cmbMotlName); motlAll=motlAll.Value;
starAll=tom_starread(io.volumeList);
allFeat=motlAll(6,:);
featU=unique(allFeat);
allMotlClean=[];

for i=1:length(featU) %110:110%length(featU)
    idx=find(allFeat==featU(i));
    motlOneFeat=motlAll(:,idx);
    cent=[starAll(i).rlnCoordinateX starAll(i).rlnCoordinateY starAll(i).rlnCoordinateZ];
    tom_emwrite(strrep(io.cmbMotlName,'.em','_B4clean.em'),motlOneFeat);
    motlOneFeatClean=cleanByCC(motlOneFeat,clean);
    %tom_emwrite(strrep(io.cmbMotlName,'.em','_B4Distclean.em'),motlOneFeatClean);
    motlOneFeatClean=cleanByDist(motlOneFeatClean,clean,cent);
    
    motlOneFeatClean=cleanByAngle(motlOneFeatClean,clean,cent);

    allMotlClean=cat(2,allMotlClean,motlOneFeatClean);
end
tom_emwrite(strrep(io.cmbMotlName,'.em','_clean.em'),allMotlClean);
motlToStar(strrep(io.cmbMotlName,'.em','_clean.em'),io);


function motlClean=cleanByAngle(motl,clean,cent)

if (size(motl,2)==0)
    motlClean=motl;
    return;
end

motlNormal=normalvec(motl,cent);

for i=1:size(motl,2)
    angDist(i)=tom_angular_distance(motl(17:19,i)',motlNormal(17:19,i)');
end
%figure; histogram(angDist,60);
use=angDist<clean.maxAngDistToNormalABS;
motlClean=motl(:,use);


function motlToStar(motlName,io)

dummy=tom_starread(io.volumeList);
motl=tom_emread(motlName); motl=motl.Value;

allCl=motl(20,:);
allClU=unique(allCl);                       

for icl=1:length(allClU)
    idx=find(allCl==allClU(icl));
    zz=1;
    for i=idx
        [~,angZYZ]=tom_eulerconvert_xmipp(motl(17,i),motl(18,i),motl(19,i),'tom2xmipp');
        newStar(zz)=dummy(1);
        newStar(zz).rlnAngleRot=angZYZ(1);
        newStar(zz).rlnAngleTilt=angZYZ(2);
        newStar(zz).rlnAnglePsi=angZYZ(3);
        newStar(zz).rlnCoordinateX=motl(8,i);
        newStar(zz).rlnCoordinateY=motl(9,i);
        newStar(zz).rlnCoordinateZ=motl(10,i);
        if (isfield(newStar,'rlnRandomSubset'))
            newStar(zz).rlnRandomSubset=round(rand(1)+1);
        end
        newStar(zz).rlnMicrographName=strrep(newStar(zz).rlnMicrographName,'.tomostar','.mrc');
         newStar(zz).rlnMicrographName=io.tomoNumToVolName{motl(5,i)};
        zz=zz+1;
    end
    tom_starwrite(strrep(motlName,'.em',['_' io.tag{icl} '.star']),newStar);
    tom_emwrite(strrep(motlName,'.em',['_' io.tag{icl} '.em']),motl(:,idx));
    sg_motl_av3_to_stopgap(strrep(motlName,'.em',['_' io.tag{icl} '.em']),strrep(motlName,'.em',['_' io.tag{icl} 'SG.star']));


    clear('newStar');
end






function motlClean=cleanByDist(motl,clean,cent)

allPos=motl(8:10,:);

diff=allPos-repmat(cent',1,size(allPos,2));
dist=vecnorm(diff);
idx=find(dist<clean.maxDistCenterABS);
motlClean=motl(:,idx);
%figure; histogram(dist,50);
%figure; plot3(allPos(1,:),allPos(2,:),allPos(3,:),'ro');
%hold on; plot3(cent(1),cent(2),cent(3),'g+');
 motlClean=cleanInterDist(motlClean,clean.minDistNeighbour,0);




function motlClean=cleanInterDist(motl,cutOff,vis)

allPos=motl(8:10,:)';
Y = pdist(allPos);
if (isempty(find(Y<cutOff)))
    motlClean=motl;    
    return;
end

Z = linkage(Y);
if (vis==1)
    figure; 
end
[st.groups,st.cmap,st.groupIdx,st.ColorThreshold]=tom_dendrogram(Z,cutOff,size(allPos,1),vis);


if (vis==1)
    figure; hold on;
    for i=1:length(st.groups)
        id=st.groups(i).members;
        if (id==-1)
            continue;
        end
        plot3(allPos(id,1),allPos(id,2),allPos(id,3),'ro','Color',st.groups(i).color); axis image;
    end
    hold off;
end

idToRemove=[];
for i=1:length(st.groups)
        id=st.groups(i).members;
        if (length(id)==1 || st.groups(i).id==0)
            continue;
        end
        ccOneCluser=motl(1,id);
        [~,idcc]=sort(ccOneCluser);
        idToRemove=cat(1,idToRemove,id(idcc(1:end-1)));
end
idxkeep=ismember(1:size(motl,2),idToRemove)==0;
motlClean=motl(:,idxkeep);






function motlClean=cleanByCC(motl,clean)

%motlClean=motl;
allCC=motl(1,:);
v=allCC>clean.minCCGlobal;
motlClean=motl(:,v);


for i=1:size(motlClean,2)
    [~,motlNeigh]=findNeighbours(motlClean,i,clean.minCCLocalRadius);
    ccCenter=motlClean(1,i);
    ccNeigh=motlNeigh(1,:);
    meaCC=mean(ccNeigh);
    stdCC=std(ccNeigh);
    thrLocal=meaCC-(stdCC.*clean.minCCLocalDiff);
    if (ccCenter>thrLocal)
        v(i)=1;
    else
        v(i)=0;
    end
end

motlClean=motlClean(:,v);




function  [idxNeigh,motlNeigh,dist]=findNeighbours(motl,indCent,rad)

posCent=motl(8:10,indCent);
allPos=motl(8:10,:);

diff=allPos-repmat(posCent,1,size(allPos,2));
dist=vecnorm(diff);

idxNeigh=find((dist<rad).*(dist>0));

motlNeigh=motl(:,idxNeigh);
dist=dist(idxNeigh);



