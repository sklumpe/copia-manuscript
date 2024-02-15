function mergeCenterToDistpairs()

dispPairList='distPairCut.em';
dispPariOrgList='../genList/motlDR.star';
pixsPair=2.24;
mic2class={'Position_1.mrc';'Position_36.mrc';'Position_37.mrc'};

orgPatchList='../../../../motlAll_clean.em';
centList2Patch='../../../../../../../../relionM/MRun4PYTOM/partsMP295B272ABS.star';
pixsOrgPatch=2.95000;

targetPixelSize=5.9;
outputList='distPairCutCent.em';



%% Code
pairMat=tom_emread(dispPairList); pairMat=pairMat.Value;
orgList=tom_starread(dispPariOrgList);
for i=1:length(orgList)
    scf=pixsPair./targetPixelSize;
    orgList(i).rlnCoordinateX=orgList(i).rlnCoordinateX.*scf;
    orgList(i).rlnCoordinateY=orgList(i).rlnCoordinateY.*scf;
    orgList(i).rlnCoordinateZ=orgList(i).rlnCoordinateZ.*scf;
end
centListP=tom_starread(centList2Patch);
for i=1:length(centListP)
    scf=pixsOrgPatch./targetPixelSize;
    centListP(i).rlnCoordinateX=centListP(i).rlnCoordinateX.*scf;
    centListP(i).rlnCoordinateY=centListP(i).rlnCoordinateY.*scf;
    centListP(i).rlnCoordinateZ=centListP(i).rlnCoordinateZ.*scf;
end

orgPMotl=tom_emread(orgPatchList); orgPMotl=orgPMotl.Value;
orgPMotl(8:10,:)=orgPMotl(8:10,:).*pixsOrgPatch/targetPixelSize;


cOrgPL=orgPMotl(8:10,:); 
tnOrgPl=orgPMotl(5,:);
onOrgPl=orgPMotl(6,:);

for i=1:size(pairMat,1)
    cPair=[orgList(pairMat(i,2)).rlnCoordinateX orgList(pairMat(i,2)).rlnCoordinateY orgList(pairMat(i,2)).rlnCoordinateZ];
    id=find(ismember(mic2class,orgList(pairMat(i,2)).rlnMicrographName));
    idxOrgPl=find(tnOrgPl==id);
    d=sg_pairwise_dist(cPair',cOrgPL(:,idxOrgPl));
    [vald,idxOrgMatch]=min(d);
    pairMat(i,6)=id;
    pairMat(i,7)=onOrgPl(idxOrgPl(idxOrgMatch));
    pairMat(i,8:10)=[centListP(pairMat(i,7)).rlnCoordinateX centListP(pairMat(i,7)).rlnCoordinateY centListP(pairMat(i,7)).rlnCoordinateZ];
    pairMat(i,11:13)=cPair;
    pairMat(i,14:16)=[orgList(pairMat(i,2)).rlnAngleRot orgList(pairMat(i,2)).rlnAngleTilt orgList(pairMat(i,2)).rlnAnglePsi];
    [~,euler_tom]=tom_eulerconvert_xmipp(pairMat(i,14),pairMat(i,15),pairMat(i,16));
    pairMat(i,17:19)=euler_tom;
end
tom_emwrite(outputList,pairMat);

