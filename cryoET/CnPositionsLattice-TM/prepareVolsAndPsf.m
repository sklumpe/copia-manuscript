function prepareVolsAndPsf()


volStar='../relionM/MRun4PYTOM/partsMP295B272ABS.star';
psfStar='../relionM/MRun4PYTOM/partsMP295B128ABS.star';
outputBase='/fs/pool/pool-drosophila/fbeck/batch1/tomoman/warpV1/TM-Neigh/vols/f2R58/vf_';
filtCer=2;

mask4Psf=tom_spheremask(ones(128,128,128),57,4); %org 58
stVol=tom_starread(volStar);
stPsf=tom_starread(psfStar);
stNew=stVol;

warning off; mkdir(fileparts(outputBase)); warning on;
waitbar=tom_progress(length(stVol),'Test Parfor Run'); 
for i=1:length(stVol)
    cf=tom_mrcread(stVol(i).rlnImageName);
    cf=cf.Value;
    cf=tom_filter(cf,filtCer);
    %cf=tom_apply_bfactor(cf,2.95,+350,1,[7.5 Inf]);
    tom_mrcwrite(cf,'name',[outputBase num2str(sprintf( '%04d', i)) '.mrc']);
    stNew(i).rlnImageName=[outputBase num2str(sprintf( '%04d', i)) '.mrc'];
    
    psf=tom_mrcread(stPsf(i).rlnCtfImage);
    psf=psf.Value;
    psf=psf.*mask4Psf;
    tom_mrcwrite(psf,'name',[outputBase num2str(sprintf( '%04d', i)) '.mrc.psf']);
    stNew(i).rlnCtfImage=[outputBase num2str(sprintf( '%04d', i)) '.mrc.psf']; 
    waitbar.update();
end
waitbar.close; 
tom_starwrite([fileparts(outputBase) '.star'],stNew);   

