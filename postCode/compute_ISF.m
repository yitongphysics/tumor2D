function compute_ISF(sim_name)

sim_name = string(sim_name*0.0002);
pathname = append(pwd,'/');
fstr = append(pathname,'pos/',sim_name,'.pos');

finfo = dir(fstr);
fprintf('-- Reading in %s\n',finfo.name);
fprintf('-- File size = %f MB\n',finfo.bytes/1e6);

% open file stream
fid = fopen(fstr);

fnum = 1;
FSTART = 1;
FSTEP = 1;
FEND = 1000;
ff = 1;

X=zeros(1000,1000);
Y=zeros(1000,1000);

for ii = FSTART:FEND

    % read in data
    tumorConfigData = readTumor2DInterface(fid);
    if tumorConfigData.NCELLS==0
        ii=ii-1;
        X(ii+1:end,:) = [];
        Y(ii+1:end,:) = [];
        break;
    end
    
    x = tumorConfigData.x;
    y = tumorConfigData.y;
    %fprintf('reading frame ii = %d/%d\n',ii,FEND);
    
    X(ii,:) = [x{1:1000}];
    Y(ii,:) = [y{1:1000}];
end
X(ii+1:end,:) = [];
Y(ii+1:end,:) = [];
Ly=7.9330033363271;
for jj=1:ii-1
    Y(jj+1:end,:) = Y(jj+1:end,:)-(Y(jj+1,:)-Y(jj,:)>Ly/2)*Ly;
    Y(jj+1:end,:) = Y(jj+1:end,:)+(Y(jj+1,:)-Y(jj,:)<-Ly/2)*Ly;
end
for jj=1:ii
    Y(jj,:) = Y(jj,:) - mean(Y(jj,:));
end
%%
Nframe = length([X(:,1)]);
dt=[0:Nframe-1];
isf = zeros(1,Nframe);
k=50;
for tt=0:Nframe-1
    X1 = X(1:Nframe-tt,:);X1 = reshape(X1,1,[]);
    X2 = X(1+tt:end,:);X2 = reshape(X2,1,[]);
    Y1 = Y(1:Nframe-tt,:);Y1 = reshape(Y1,1,[]);
    Y2 = Y(1+tt:end,:);Y2 = reshape(Y2,1,[]);
    isf(tt+1) = compute_ISF_2D(X1, Y1, X2, Y2, k, Ly);
end

%%
yourFolder = append(pathname,'ISF');
if ~exist(yourFolder, 'dir')
   mkdir(yourFolder)
end
data_name = append(pathname,'ISF/',sim_name,'_ISF.mat');
save(data_name,'dt','isf');

end
