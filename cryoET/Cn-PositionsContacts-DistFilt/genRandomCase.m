function resStat=genRandomCase()

st=tom_starread('../../genList/motlDRSh.star');
subSize=5622;
nrSim=1000;
class2Sym=[1 3 5];
clVectAll=[st(:).rlnClassNumber];

for is=1:nrSim
    r=randperm(length(clVectAll));
    classVect=clVectAll(r(1:subSize));

    %r=randperm(length(classVect));
    classVectR=classVect;
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
    resStat(i,:)=[res(i,1:2,1) mean(OneCl) std(OneCl) 3*std(OneCl)  min(OneCl) max(OneCl) ];
end


