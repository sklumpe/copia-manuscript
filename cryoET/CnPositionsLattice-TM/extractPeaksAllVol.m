function extractPeaksAllVol()

st.io.outputFolder='runs/f2All/';
st.io.volumeList='vols/f2R58.star';
st.io.templates{1}='ref/bin1/c1refined128.mrc';
st.io.templates{2}='ref/bin1/c3refined128.mrc';
st.io.templates{3}='ref/bin1/c5refined128.mrc';
st.io.tag={'C1';'C3';'C5'};
st.io.szInfo=[272 272 272];

st.extract.do=1;
st.extract.numPeaksCn=[60 20 12 92]; %[60 20 12 92]; %[60 20 12 2];
st.extract.pRemoveRad=25; %20
st.extract.mask='masks/bin1/mask4Extract2.mrc';  %'masks/bin1/maskExtracFull.mrc'; %'masks/bin1/mask4Extract2.mrc';
st.extract.angList='angles/a_7.em';

st.cmb.do=1;
st.cmb.motlName='cmbmotl.em';

if (st.extract.do)
    star=tom_starread(st.io.volumeList);
    dstar=tom_extractData(star);
    waitbar=tom_progress(length(star),['Extracting peaks ' num2str(length(star)) ' vols']);
    parfor i=1:length(star)
        extractPeaksOneSubVol(st,star(i),dstar.label.tomoID(i));
        waitbar.update();
    end
    waitbar.close();
else
    disp('skipping peak extraction');
end


combineLists(st.cmb,st.io)



function combineLists(cmb,io)

if (cmb.do==0)
    disp('skipping motl combination');
    return
end

star=tom_starread(io.volumeList);
dstar=tom_extractData(star);

waitbar=tom_progress(length(star),['combining lists ' num2str(length(star)) ' ']);
for i=1:length(star)
    [~,volName]=fileparts(star(i).rlnImageName);
    Orig=[star(i).rlnCoordinateX star(i).rlnCoordinateY star(i).rlnCoordinateZ];
    listname=[io.outputFolder filesep volName filesep cmb.motlName];
    list{i}=applyOrigToList(listname,round([Orig - (io.szInfo./2)]),dstar.label.tomoID(i),i);
    waitbar.update();
end
waitbar.close();

listCmb=[];
for i=1:length(list)
    listCmb=cat(2,listCmb,list{i});
end

warning off; mkdir([io.outputFolder filesep 'cmbLists']); warning on;
tom_emwrite([io.outputFolder filesep 'cmbLists' filesep 'motlAll.em'],listCmb);
disp('use dstar for mappin tomoId to tomoname')

disp(' ');


function list=applyOrigToList(listname,Orig,tomoNum,featNum)

 if (tom_isemfile(listname))
    list=tom_emread(listname);
    list=list.Value;
    for i=1:size(list,2)
        list(5,i)=tomoNum;
        list(6,i)=featNum;
        list(8,i)=list(8,i)+Orig(1);
        list(9,i)=list(9,i)+Orig(2);
        list(10,i)=list(10,i)+Orig(3);
    end
 end




function extractPeaksOneSubVol(st,star,tomoNum)


[~,volName]=fileparts(star.rlnImageName);
outFold=[st.io.outputFolder filesep volName];
for i=1:length(st.io.templates)
     [~,b,c]=fileparts(st.io.templates{i});
    st.extract.scName=[outFold filesep st.io.tag{i} filesep 'scores_' b '.em'];
    st.extract.angName=[outFold filesep st.io.tag{i} filesep 'angles_' b '.em'];
    st.extract.tomoNum=tomoNum;
    st.extract.NrPeaks=st.extract.numPeaksCn(i); 
    st.extract.outputName=[outFold filesep st.io.tag{i} 'motl' '.em'];
    extractPeaks(st.extract);
end

mergeCCF(st,outFold);

st.extract.scName=[outFold filesep 'cmb' filesep 'scores_cmb.em'];
st.extract.angName=[outFold filesep 'cmb' filesep 'angles_cmb.em'];
st.extract.clName=[outFold filesep 'cmb' filesep 'classes_cmb.em'];
st.extract.NrPeaks=st.extract.numPeaksCn(4); %20; 12
st.extract.outputName=[outFold filesep 'cmb' 'motl' '.em'];
extractPeaks(st.extract);


function mergeCCF(st,outFold)

scCmb=ones(st.io.szInfo).*-1;
angCmb=zeros(st.io.szInfo);
classCmb=zeros(st.io.szInfo);
for i=1:length(st.io.templates)
     [~,b,c]=fileparts(st.io.templates{i});
    scName=[outFold filesep st.io.tag{i} filesep 'scores_' b '.em'];
    angName=[outFold filesep st.io.tag{i} filesep 'angles_' b '.em'];
    scAct=tom_emread(scName);scAct=scAct.Value;
    angAct=tom_emread(angName);angAct=angAct.Value;
    ind=find(scAct>scCmb);
    scCmb(ind)=scAct(ind);
    angCmb(ind)=angAct(ind);
    classCmb(ind)=i;
end
warning off; mkdir([outFold filesep 'cmb']);  warning on;
tom_emwrite([outFold filesep filesep 'cmb' filesep 'scores_cmb.em'],single(scCmb));
tom_emwrite([outFold filesep 'cmb' filesep 'angles_cmb.em'],single(angCmb));
tom_emwrite([outFold filesep 'cmb' filesep 'classes_cmb.em'],single(classCmb));


function extractPeaks(st)

% if (nargin==0)
%     st.scName='outputFold6Real/scores_c5refined64.em';
%     st.angName='outputFold6Real/angles_c5refined64.em';
%     st.maskE='masks/bin2/mask4Extract2.mrc';
%     st.angListName='/fs/pool/pool-bmapps/hpcl8/app/soft/PYTOM/19-04-23/pytom/pytom/angles/angleLists/angles_07_45123.em';
%     st.tomoNum=1;
%     st.NrPeaks=12; %20; 12
%     st.PeakRad=25;
%     st.outputName='outputFold6Real/motl';
% end

vsc=tom_emread(st.scName); vsc=vsc.Value;
vang=tom_emread(st.angName); vang=vang.Value;
angL=tom_emreadc(st.angList); angL=angL.Value;
mE=tom_mrcread(st.mask); mE=mE.Value;
if (isfield(st,'clName'))
    vcl=tom_emread(st.clName); vcl=vcl.Value;
end

list=zeros(20,st.NrPeaks);


tmpSc=vsc.*mE;
for i=1:st.NrPeaks
    [pos,val,tmpSc]=tom_peak(tmpSc,st.pRemoveRad);
    angI=vang(pos(1),pos(2),pos(3))+1;
    ang=angL(:,angI).*(180/pi);
    list(1,i)=val;
    list(4,i)=i;
    list(5,i)=st.tomoNum;
    list(8,i)=pos(1);
    list(9,i)=pos(2);
    list(10,i)=pos(3);
    list(17,i)=ang(1);
    list(18,i)=ang(2);
    list(19,i)=ang(3);
    if (isfield(st,'clName'))
        cl=vcl(pos(1),pos(2),pos(3));
        list(20,i)=cl;
    end
end
tom_emwrite(st.outputName,list);





