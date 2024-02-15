function filterByDist()
 

starName='../genList/motlDRSh.star';
outputFolder='connect';
pixS=2.24;
cutOff=16.*pixS;
% bins for hist
calc_distMat=1;
avgCall='mpirun -np 25 relion_reconstruct_mpi --i particles.star --o vol.mrc --maxres 35 --3d_rot --ctf --sym c1';
scaleF4Sub=0.37966;

%% Code
splitStarFiles([outputFolder '/classAvg/*/particles.star'],scaleF4Sub);

motl=tom_starread(starName);
warning off; mkdir(outputFolder);warning on;

if (calc_distMat==1)
    tomoNrU=unique({motl(:).rlnMicrographName});
    distPairMat=[];
    figure; hold on;
    for i=1:length(tomoNrU)
        idx=find(ismember({motl(:).rlnMicrographName},tomoNrU{i}));
        disp(['found:' num2str(length(idx)) ]);
        distPairMatOneT=calcAlldist(motl(idx));
        distPairMatOneT(:,2)=idx(distPairMatOneT(:,2));
        distPairMatOneT(:,3)=idx(distPairMatOneT(:,3));
        distPairMat=cat(1,distPairMat,distPairMatOneT);
    end
    tom_emwrite('distPair.em',distPairMat);
else
   distPairMat=tom_emread('distPair.em');
   distPairMat=distPairMat.Value;
end

genStatistics(distPairMat,motl,outputFolder,pixS,cutOff);
filterStarFiles(distPairMat,motl,outputFolder,pixS,cutOff,avgCall);



function distCutCl=cleanDistList(distCutCl)

% g=tom_emread('../distPair.em');distPairMat=g.Value;
% pixS=2.24;
% dist=distPairMat(:,1).*pixS;
% cutOff=16.*pixS;
% idxC=find(dist<cutOff);
% distCutCl=distPairMat(idxC,:);

distCutCl=testAndFix(distCutCl);
testAndFix(distCutCl);


function distCutCl=testAndFix(distCutCl)
idBad=-1;
allInd1=distCutCl(:,2);
allInd2=distCutCl(:,3);
verr=zeros(size(distCutCl,1),1);
for i=1:size(distCutCl,1)
    ind2=distCutCl(i,3);
    idx=find(allInd1==ind2);
    dist1=distCutCl(i,1);
    if (length(idx)~=1)
        verr(i)=1;
    end
    ind1P=distCutCl(idx,2);
    if (ind2~=ind1P)
        verr(i)=1;
    end
    idx=find(allInd2==ind1P);
    if (length(idx)~=1)
        verr(i)=1;
        idBad=idx(find(distCutCl(idx,1)~=dist1));
    end
end
disp(['found ' num2str(length(find(verr==1))./2 ) ' errors ']);

if (idBad>1)
    disp('fixing ');
    idxGood=find(([1:size(distCutCl,1)]==idBad)==0);
    distCutCl=distCutCl(idxGood,:);
end

function splitStarFiles(wkStarFiles,scaleF)

d=dir(wkStarFiles);
for i=1:length(d)
    st=tom_starread([d(i).folder filesep d(i).name] );
    tomoAll={st(:).rlnMicrographName};
    tomoU=unique(tomoAll);
    for ii=1:length(tomoU)
        idx=find(ismember(tomoAll,tomoU{ii}));
        sttmp=st(idx);
        sttmp(1).Header=st(1).Header;
        for iii=1:length(sttmp)
            sttmp(iii).rlnCoordinateX=sttmp(iii).rlnCoordinateX.*scaleF;
            sttmp(iii).rlnCoordinateY=sttmp(iii).rlnCoordinateY.*scaleF;
            sttmp(iii).rlnCoordinateZ=sttmp(iii).rlnCoordinateZ.*scaleF;
            sttmp(iii).rlnDetectorPixelSize=sttmp(iii).rlnDetectorPixelSize.*1/scaleF;
        end
        baseout=strrep(tomoU{ii},'.mrc','');
        outputFold=[d(i).folder filesep 'splitbyTomo' ];
        warning off; mkdir(outputFold); warning on;
        tom_starwrite([outputFold filesep baseout '_' num2str(sttmp(iii).rlnDetectorPixelSize) '.star'],sttmp);
    end
end


function filterStarFiles(distPairMat,motl,outputFolder,pixS,cutOff,avgCall)

warning off; mkdir([outputFolder filesep 'classAvg']); warning on;
class2Sym=[1 3 5];

idxC=find(distPairMat(:,1).*pixS<cutOff);
distCutCl=distPairMat(idxC,:);
distCutCl=testAndFix(distCutCl);
for i=1:size(distCutCl,1)
    clVect(i,:)=distCutCl(i,4:5);
end
[Mu] = unique(clVect, 'rows', 'stable');

for i=1:size(Mu,1)
    idSM=find( (distCutCl(:,4)==Mu(i,1)) .*(distCutCl(:,5)==Mu(i,2)) );
    %if (Mu(i,1)==Mu(i,2))
        ind=distCutCl(idSM,2);%  unique(cat(1,distPairMat(idSM,2),distPairMat(idSM,3)));
    %else
        ind=distCutCl(idSM,2);
    %end
    fSt=motl(ind);
    fSt(1).Header=motl.Header;
    clTag=[num2str(class2Sym(Mu(i,1))) '-' num2str(class2Sym(Mu(i,2)))];
    clFold=[outputFolder filesep 'classAvg' filesep clTag];
    warning off; mkdir([outputFolder filesep 'classAvg' filesep clTag]); warning on;
    tom_starwrite([clFold filesep 'particles.star'],fSt);
    call=['cd ' clFold ';' avgCall];
    unix(call);
end


disp(' ');


function genStatistics(distPairMat,motl,outputFolder,pixS,cutOff)

class2Sym=[1 3 5];

warning off; mkdir([outputFolder filesep 'distPlot']);
nrBin=800; 
dist=distPairMat(:,1).*pixS;
idx=find(dist<50.*pixS); %filter 4Disp
h=figure; histogram(dist(idx),nrBin); 
saveas(h,[outputFolder filesep 'distPlot/allDist.png']);
close(h);
csvwrite([outputFolder filesep 'distPlot/allDist.csv'],dist(idx));
idxC=find(dist<cutOff);


%% Count classes 
warning off; mkdir([outputFolder filesep 'classdist']); warning on;
distPairMatCut=cleanDistList(distPairMat(idxC,:));
tom_emwrite('distPairCut.em',distPairMatCut);
distCutCl=distPairMatCut(:,4:5);


fid=fopen([outputFolder filesep 'distPlot/contToFreeRatio.txt'],'wt');
fprintf(fid,'%d of %d == %f%\n',size(distCutCl,1),size(distPairMat,1),(length(idxC)/size(distPairMat,1)).*100);
fprintf('Connected Pores %d of Tot: Pores %d == %f  \n',size(distCutCl,1),size(distPairMat,1),(length(idxC)/size(distPairMat,1)).*100);
fclose(fid);

for i=1:size(distCutCl,1)
    distCutCl(i,:)=sort(distCutCl(i,:));
end

[Mu,ia,ic] = unique(distCutCl, 'rows', 'stable');          
h = accumarray(ic, 1);                             
res=[class2Sym(Mu) h (h/size(distCutCl,1))*100];
[~,id]=sortrows(res(:,1:2));
res=res(id,:);
[resRandABS,resRandRel]=genRandomCase(motl,size(distCutCl,1)-1,outputFolder);

zz=0;
for i=1:size(res,1)
    zz=zz+1;
    res4Plot(i,1)=res(i,4);
    res4Plot(i,2)=resRandRel(i,3);
    label{i}=strrep(num2str(res(i,1:2)),' ','-');
end
h=figure; bar(res4Plot);set(gca, 'XTickLabels', label);
ngroups = size(res4Plot, 1);
nbars = size(res4Plot, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
hold on;
for i = 2:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    hr=errorbar(x,res4Plot(:,i), resRandRel(:,5), '.');
    set(hr,'LineWidth',2.2);set(hr,'CapSize',22)
end
hold off;
legend({'measured';'random';'3std error'});
set(gcf,'Position',[1000         552        2199        1495]);
saveas(h,[outputFolder filesep 'classdist/cldist.png']);

res4Mat=res;
res4Mat(:,5)=resRandABS(:,3);
res4Mat(:,6)=resRandRel(:,3);
res4Mat(:,7)=resRandRel(:,5);

C = [{'c1' 'c2' 'nr abs-Mea' 'rel mea %' 'nr absRand' 'rel Rand %' 'rel error Rand 3std%'}; num2cell(res4Mat)];
writecell(C,[outputFolder filesep 'classdist/cldist.csv']);


% function filterStar()
% idTmp=find((distPairMat(:,4)==1) .*(distPairMat(:,5)==1));
% figure; hist(distPairMat(idTmp,1),nrBin);
% 
% idSM=find((distPairMat(:,1)>1).* (distPairMat(:,1)<cutOff).* (distPairMat(:,4)==1) .*(distPairMat(:,5)==1));
% ind=distPairMat(idSM,2);
% 
% %ind=unique(cat(1,distPairMat(idSM,2),distPairMat(idSM,3)));
% disp(['found:' num2str(length(ind)) ]);
% fSt=motl(ind);resStatRel(i,3:end)
% fSt(1).Header=motl.Header;
% tom_starwrite('filtered2.star',fSt);
% 
% 
% figure; hist(distPairMat(:,1),nrBin);
% 
% disp(' ');

function [resStatABS,resStatRel]=genRandomCase(motl,subSize,outputFolder)

%st=tom_starread('../../genList/motlDRSh.star');
%subSize=size(distPairMat)-1;
nrSim=1000;
class2Sym=[1 3 5];
clVectAll=[motl(:).rlnClassNumber];

for is=1:nrSim
    r=randperm(length(clVectAll));
    classVect=clVectAll(r(1:subSize));
    classVect=cat(2,classVect,classVect); r=randperm(length(classVect)); classVectR=classVect(r);
    %classVectR=classVect;
    classVectRPair=reshape(classVectR,[round(length(classVectR)./2) 2]);
    for i=1:size(classVectRPair,1)
        classVectRPair(i,:)=sort(classVectRPair(i,:));
    end
    [Mu,ia,ic] = unique(classVectRPair, 'rows', 'stable');
    h = accumarray(ic, 1);
    resTmp=[class2Sym(Mu) h (h/length(classVectRPair))*100];
    [~,id]=sortrows(resTmp(:,1:2));
    res(:,:,is)=resTmp(id,:);
end

%round(res)
absNr=reshape(res(:,3,:),[size(res,1) size(res,3)]);
for i=1:size(absNr,1)
    OneCl=absNr(i,:);
    resStatABS(i,:)=[res(i,1:2,1) mean(OneCl) std(OneCl) 3*std(OneCl)  min(OneCl) max(OneCl) ];
    resStatRel(i,:)=[res(i,1:2,1) mean(OneCl)./subSize std(OneCl)./subSize 3*std(OneCl)./subSize  min(OneCl)./subSize max(OneCl)./subSize ];
    resStatRel(i,3:end)=resStatRel(i,3:end)*100;
end

disp(' ');


function distPairMat=calcAlldist(motl)

distPairMat=zeros(length(motl),5);
allCoords=[motl(:).rlnCoordinateX; motl(:).rlnCoordinateY; motl(:).rlnCoordinateZ]';

for i=1:length(motl)
    d=sg_pairwise_dist(allCoords(i,:)',allCoords');
    d(i)=1000000;
    [dist,ind]=min(d);
    distPairMat(i,:)=[dist i ind motl(i).rlnClassNumber motl(ind).rlnClassNumber];
    %clear('d');
end


function D = sg_pairwise_dist(X, Y)
% A function to determine the pairwise distance between two n-dimensional
% arrays. I found it after googling "matlab pairwise distance". It should
% be substantially faster than the matlab function pdist.
%
% A fix had to be implemented as the bsxfun changed in the matlab 8.x
% and produces rounding errors close to zero. This results in some numbers 
% becoming negative, making the distance array a complex array. This fix
% rounds all negative values to zero. 
%
% WW 02-2016


% Squares of the distances
sq_dist = (bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y));

% Find negative values
neg = sq_dist < 0;
% Set negative values to zero
sq_dist(neg) = 0;

D = sq_dist.^0.5;

