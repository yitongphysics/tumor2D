function sp_dis(sim_name,aN)
%% voronoi tesselation

pathname = append(pwd,'/');
sim_name = string(sim_name);


FEND = 1000;

fstr = append(pathname,sim_name,'.pos');
finfo = dir(fstr);
fid = fopen(fstr);

sp = zeros(FEND,aN);
%loop over frames
for ii = 1:1:FEND
    % read in data
    tumorConfigData = readTumor2DInterface(fid);
    if tumorConfigData.NCELLS==0
        break;
    end
    a  = tumorConfigData.a; a = a(end-aN+1:end);
    p = tumorConfigData.p; p = p(end-aN+1:end);
    
    sp(ii,:) = p.^2./a;
end

sp = sp/4/pi;

yourFolder = append(pathname,'sp');
if ~exist(yourFolder, 'dir')
   mkdir(yourFolder)
end
data_name = append(pathname,'sp/',sim_name,'_sp.mat');
save(data_name,'sp');


end