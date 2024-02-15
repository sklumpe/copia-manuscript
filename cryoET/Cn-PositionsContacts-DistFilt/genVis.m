function genVis()

distMatName='distPairCutCent.em';
outputFold='vis/lists/';
radiusCenter=20;
tomoId=[1 2 3] ;
classes=[1 1;1 5;5 1;5 5; 3 5; 5 3];
colCon{1}=[1 0 0]; colCon{2}=[0 0.9 0.1];  colCon{3}=[0 1 0]; 
colCon{4}=[0 0 1]; colCon{5}=[0 1 1];  colCon{6}=[0 1 0.9]; 
class2Sym=[1 3 5];
tomoNum2Pos={'pos_1';'pos_36';'pos_37'};

%% 

distmat=tom_emread(distMatName);
distmat=distmat.Value;


for i=1:length(tomoId)
    mkdir([outputFold tomoNum2Pos{tomoId(i)}]);
    for ic=1:size(classes,1)
        vT=[distmat(:,6)==tomoId(i)]';
        vc1=class2Sym(distmat(:,4))==classes(ic,1);
        vc2=class2Sym(distmat(:,5))==classes(ic,2);
        idx=find(vT.*vc1.*vc2);
        nameF=[tomoNum2Pos{tomoId(i)} '==>pToCent' num2str(classes(ic,1)) '-' num2str(classes(ic,2)) '.bild'];
        fid=fopen([outputFold tomoNum2Pos{tomoId(i)} filesep nameF],'wt');
        fprintf(fid,'.color %s\n',num2str(colCon{ic}));
        motl=zeros(20,length(idx));
        for ii=1:length(idx)
            motl(4,ii)=ii;
            motl(5,ii)=tomoId(i);
            motl(6,ii)=distmat(idx(ii),7);
            motl(8:10,ii)=distmat(idx(ii),11:13);
            motl(17:19,ii)=distmat(idx(ii),17:19);
            motl(20,ii)=distmat(idx(ii),7);
            pos2=distmat(idx(ii),11:13);
            pos1=distmat(idx(ii),8:10);
            fprintf(fid,'.arrow %f %f %f %f %f %f 2\n',pos1(1),pos1(2),pos1(3),pos2(1),pos2(2),pos2(3));
        end
        fclose(fid);
        nameF=[tomoNum2Pos{tomoId(i)} '==>pos' num2str(classes(ic,1)) '-' num2str(classes(ic,2)) '.em'];
        tom_emwrite([outputFold tomoNum2Pos{tomoId(i)} filesep nameF],motl);
        %gen bild connections
        nameF=[tomoNum2Pos{tomoId(i)} '==>pair' num2str(classes(ic,1)) '-' num2str(classes(ic,2)) '.bild'];
        fid=fopen([outputFold tomoNum2Pos{tomoId(i)} filesep nameF],'wt');
        fprintf(fid,'.color %s\n',num2str(colCon{ic}));
        for ii=1:length(idx)
             pos1=distmat(idx(ii),11:13);
             vid1=distmat(:,2)==distmat(idx(ii),3);
             vid2=distmat(:,3)==distmat(idx(ii),2);
             idCon=find(vid1.*vid2);
             pos2=distmat(idCon,11:13);
            fprintf(fid,'.arrow %f %f %f %f %f %f 2\n',pos1(1),pos1(2),pos1(3),pos2(1),pos2(2),pos2(3));
        end
        fclose(fid);
    end
end
