function prepareTemplateMatching()

st.io.outputFolder='runs/f2All/';
st.io.volumeList='vols/f2R58.star';
st.io.templates{1}='ref/bin1/c1refined128.mrc';
st.io.templates{2}='ref/bin1/c3refined128.mrc';
st.io.templates{3}='ref/bin1/c5refined128.mrc';
st.io.masks{1}='masks/bin1/mBC1.mrc';
st.io.masks{2}='masks/bin1/mBC3.mrc';
st.io.masks{3}='masks/bin1/mBC5.mrc';
st.io.tag={'C1';'C3';'C5'};
st.io.szInfo=[272 272 272];

%% Matching
st.match.jobtemplate='jobTemplates/jobTM.xml';
st.match.angList='angles/a_7.em';
st.match.wedgeInfo.ang1=32;%32;
st.match.wedgeInfo.ang2=32;%32;
st.match.pytom='/fs/pool/pool-bmapps/hpcl8/app/soft/PYTOM/21-09-2023/conda3/envs/pytom_env//bin/localization.py';

%% Computing
comp.partion='hpcl9';
comp.nodes=3;
comp.GpuPerNode=4;
comp.qsubTempl='jobTemplates/qsubTempl.sh';
%comp.comandList=[p.tmRootFolder 'list.txt'];


%% Code


%prepare jobs
star=tom_starread(st.io.volumeList);
allJobs=[];
for i=1:length(star)
     jobsOneVol=prepareJobs(st,star(i));
     allJobs=cat(2,allJobs,jobsOneVol);
end

%write submission scripts
smFold=[st.io.outputFolder filesep 'submission/'];
smAllFile=[smFold '/submitAll.sh'];
warning off;  mkdir(smFold); warning on;
fid=fopen(smAllFile,'wt');
packages=tom_calc_packages(comp.nodes,length(allJobs));
for i=1:size(packages,1)
    subFile=[smFold '/pyTM_'  num2str(i) '.sh'];
    writeJoblistOneNode(subFile,allJobs(packages(i,1):packages(i,2)),comp.GpuPerNode,comp.qsubTempl);
    fprintf(fid,'%s\n',['sbatch ' subFile]);
end
fclose(fid);
unix(['chmod +x ' smAllFile]);
disp(['execute: ' smAllFile ]);


function writeJoblistOneNode(subFile,jobs,nrGpus,qsubTempl)

baseC=['cat ' qsubTempl ' | awk ''{'];
logC=['gsub("xxx","'  strrep(subFile,'.sh','') '");'];
endC=['print $0}'' > ' subFile];
call=[baseC logC endC];
unix(call);

pidAll=[];
fid=fopen(subFile,'a');
for i=1:length(jobs)
    gNr=mod(i,nrGpus);
    pidSingle=['pid' num2str(gNr) '=$!'];
    pidAll=[pidAll ' $pid' num2str(gNr) ];
    stmp=strsplit(jobs{i},'-j ');
    logS=[' &> ' strrep(stmp{2},'.xml','pytomTM.log')];
    fprintf(fid,'%s\n',[jobs{i} ' -g ' num2str(gNr) ' '  logS ' & ' pidSingle]);
     if (gNr==0 || i==length(jobs))
        fprintf(fid,'wait%s\n',pidAll);
        pidAll=[];
     end
end
fclose(fid);


function job=prepareJobs(st,star)

[~,volName]=fileparts(star.rlnImageName);
outFold=[st.io.outputFolder filesep volName];
warning off; mkdir(outFold); warning on;

for i=1:length(st.io.templates)
    tm.template=st.io.templates{i};
    tm.mask=st.io.masks{i};
    tm.angles=st.match.angList;
    tm.jobTemplate=st.match.jobtemplate;
    wedgeInfo=st.match.wedgeInfo;
    szInfo=st.io.szInfo; 
    jobOut=['job' st.io.tag{i} '.xml'];
    warning off; mkdir([outFold filesep st.io.tag{i}]); warning on;
    jobFile=tomoman_pytomTMprepare([outFold filesep st.io.tag{i}],tm,star.rlnImageName,wedgeInfo,szInfo,jobOut);
    job{i}=[st.match.pytom ' -j ' jobFile ];
end





