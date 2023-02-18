function Vt_phi(sim_name)
%% voronoi tesselation for packing fraction
% assume all adipocytes have the same number of vertices

pathname = append(pwd,'/');
sim_name = string(sim_name);

FEND = 1000;

phi = zeros(1,FEND);

fstr = append(pathname,sim_name,'.pos');
fid = fopen(fstr);

for ii = 1:1:FEND
    vpos=[];
    % read in data
    tumorConfigData = readTumor2DInterface(fid);
    if tumorConfigData.NCELLS==0
        break;
    end
    NCELLS = tumorConfigData.NCELLS;
    tN = tumorConfigData.tN;
    aN = NCELLS - tN;
    L = tumorConfigData.L(1,:);
    Ly = L(2);
    Lx = L(1);
    wpos = tumorConfigData.wpos;
    a  = tumorConfigData.a;
    x = tumorConfigData.x;
    y = tumorConfigData.y;
    r = tumorConfigData.r;
    r0 = r{tN+1};
    aV = length(r0);
    r0=r0(1);
    x_v = [x{1:tN}];
    y_v = [y{1:tN}];

    for jj=tN+1:NCELLS
        x_v = [x_v,[x{jj}].'];
        y_v = [y_v,[y{jj}].'];
    end
    nV = length(x_v);
    for yy = -1:1
        vtmp = [x_v.', y_v.' + yy.'*Ly];
        vpos = [vpos;vtmp];
    end
    DT = delaunayTriangulation(vpos);
    [V_vor,r_vor] = voronoiDiagram(DT);

    A0 = 0;
    for aa = 1:nV-tN
        xtmp = V_vor(r_vor{aa+tN+nV},1);
        ytmp = V_vor(r_vor{aa+tN+nV},2);
        xtmp = xtmp(xtmp~=Inf);
        ytmp = ytmp(ytmp~=Inf);
        xtmp = [xtmp; xtmp(1)];
        ytmp = [ytmp; ytmp(1)];

        ytmp1 = 100;
        ytmp2 = 100;
        vv_list=find(xtmp>wpos & xtmp<Lx);
        if(length(find(xtmp<wpos)))
            for vv = 1:length(xtmp)-1
                if((xtmp(vv)-wpos)*(xtmp(vv+1)-wpos)<0)
                    if(ytmp1 == 100)
                        ytmp1 = (ytmp(vv+1)-ytmp(vv))/(xtmp(vv+1)-xtmp(vv))*(wpos-xtmp(vv)) + ytmp(vv);
                    elseif(ytmp2 ==100)
                        ytmp2 = (ytmp(vv+1)-ytmp(vv))/(xtmp(vv+1)-xtmp(vv))*(wpos-xtmp(vv)) + ytmp(vv);
                    else
                        fprintf("error");
                        %exit(1);
                    end
                end
            end
            xtmp = [xtmp(vv_list);wpos;wpos];
            ytmp = [ytmp(vv_list);ytmp1;ytmp2];
            if(vv_list(1) == 1)
                xtmp(1)=[];
                ytmp(1)=[];
            end
            [xtmp,ytmp] = poly2cw(xtmp,ytmp);
        end
        if(length(find(xtmp>Lx)))
            for vv = 1:length(xtmp)-1
                if((xtmp(vv)-Lx)*(xtmp(vv+1)-Lx)<0)
                    if(ytmp1 == 100)
                        ytmp1 = (ytmp(vv+1)-ytmp(vv))/(xtmp(vv+1)-xtmp(vv))*(Lx-xtmp(vv)) + ytmp(vv);
                    elseif(ytmp2 ==100)
                        ytmp2 = (ytmp(vv+1)-ytmp(vv))/(xtmp(vv+1)-xtmp(vv))*(Lx-xtmp(vv)) + ytmp(vv);
                    else
                        fprintf("error");
                        %exit(1);
                    end
                end
            end
            xtmp = [xtmp(vv_list);Lx;Lx];
            ytmp = [ytmp(vv_list);ytmp1;ytmp2];
            if(vv_list(1) == 1)
                xtmp(1)=[];
                ytmp(1)=[];
            end
            [xtmp,ytmp] = poly2cw(xtmp,ytmp);
        end
        [~,I] = sort(angle(complex(xtmp-mean(xtmp),ytmp-mean(ytmp))));
        xtmp = xtmp(I);
        ytmp = ytmp(I);
        A0 = A0 + polyarea(xtmp,ytmp);
    end
    A = sum(a(tN+1:NCELLS)) + aN*(aV/2+1)*pi*r0^2;
    phi(ii) = A/A0;
end


yourFolder = append(pathname,'phi');
if ~exist(yourFolder, 'dir')
   mkdir(yourFolder)
end
data_name = append(pathname,'phi/',sim_name,'_phi.mat');
save(data_name,'phi');


end
