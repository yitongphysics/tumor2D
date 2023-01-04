%% voronoi tesselation

%clear;close all;clc

%sim_name = '0.021_354';
%pathname = '/Users/yitongzheng/Documents/Corey/tumor2D/OverDamped/1214_long/';
%plotit = 1;
%videoit = 1;
%datait = 1;

%FEND = 1;

if(videoit)
    vid_name = append(pathname,sim_name,'_vor.mp4');
    vobj=VideoWriter(vid_name,'MPEG-4');
    vobj.FrameRate = 10;
    open(vobj)
end

fstr = append(pathname,sim_name,'.pos');
finfo = dir(fstr);
fid = fopen(fstr);


L_sur = zeros(1,FEND);
for ii = 1:1:FEND
    vpos=[];
    % read in data
    tumorConfigData = readTumor2DInterface(fid);
    if tumorConfigData.NCELLS==0
        break;
    end
    NCELLS = tumorConfigData.NCELLS;
    tN = tumorConfigData.tN;
    L = tumorConfigData.L(1,:);
    Ly = L(2);
    Lx = L(1);
    wpos = tumorConfigData.wpos;
    x = tumorConfigData.x;
    y = tumorConfigData.y;
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
    
    if(plotit)
        f=figure(1); clf, hold on, box on;
        f.Position = [100 100 1200 700];
        for jj=1:tN
            patch(V_vor(r_vor{jj+nV},1),V_vor(r_vor{jj+nV},2),'r','FaceAlpha',0.3)
            patch(V_vor(r_vor{jj+2*nV},1),V_vor(r_vor{jj+2*nV},2),'g','FaceAlpha',0.3)
            patch(V_vor(r_vor{jj},1),V_vor(r_vor{jj},2),'b','FaceAlpha',0.3)
            %text(vpos(jj+nV,1),vpos(jj+nV,2),string(jj))
        end
        for jj=tN+1:nV
            patch(V_vor(r_vor{jj+nV},1),V_vor(r_vor{jj+nV},2),'w','FaceAlpha',0.3)
            patch(V_vor(r_vor{jj+2*nV},1),V_vor(r_vor{jj+2*nV},2),'w','FaceAlpha',0.3)
            patch(V_vor(r_vor{jj},1),V_vor(r_vor{jj},2),'w','FaceAlpha',0.3)
            %text(vpos(jj+nV,1),vpos(jj+nV,2),string(jj))
        end
        hold on
        plot([wpos Lx],[0 0],'k-','linewidth',2);
        plot([wpos Lx],[Ly Ly],'k-','linewidth',2);
        plot([Lx Lx],[-Ly/4 5*Ly/4],'k-','linewidth',2);
        plot([wpos wpos],[-Ly/4 5*Ly/4],'k-','linewidth',2);
    end
    % find out shared side of polygons
    % loop over V_vor, find vertices that belongs to both T and A.
    % record cells and vertices
    nVT = size(V_vor);nVT=nVT(1);% number of voronoi vertices
    V_list=[];
    C_list={};
    C_tuple = [];
    for jj=2:nVT% exclude the first vertice (Inf,Inf).
        for kk=1:tN
            if(~isempty(find(r_vor{kk+nV}==jj,1)))
                for ll=tN+1:nV
                    if(~isempty(find(r_vor{ll+nV}==jj,1))||...
                        ~isempty(find(r_vor{ll}==jj,1))||...
                        ~isempty(find(r_vor{ll+2*nV}==jj,1)))
                        if(isempty(find(V_list==jj,1)))
                            V_list=[V_list;jj];
                            C_list{length(C_list)+1}=ll;
                            C_tuple = [C_tuple;ll];
                        else
                            C_list{end}=[C_list{end},ll];
                            C_tuple = [C_tuple;ll];
                        end
                    end
                end
                break;
            end
        end
    end
    %V_list(1)=[];
    %C_list(1)=[];
    [C_tuple,Index]=unique(C_tuple);
    C_tuple = [C_tuple;C_tuple + nV;C_tuple + 2*nV];
    
    %scatter(V_vor(V_list,1),V_vor(V_list,2),'b')
    %scatter(V_vor(V_list,1),V_vor(V_list,2)+Ly,'y')
    %scatter(V_vor(V_list,1),V_vor(V_list,2)-Ly,'g')
    %scatter(vpos(C_tuple,1),vpos(C_tuple,2))
    
    %compute length
    % 1st: tumor-adipocyte surface
    % 2nd: tumor-box surface
    L=0;
    L_right=[];
    L_left =[];
    %[UniXY,Index]=unique(C_tuple);
    for jj = 1:length(C_tuple)
        C_test = C_tuple(jj);
        V_test = r_vor{C_test};
        V_test = [V_test,V_test(1)];
        for kk = 1:length(V_test)-1
            if(~isempty(find(V_list==V_test(kk),1)) && ~isempty(find(V_list==V_test(kk+1),1)))
                x1 = V_vor(V_test(kk),1);
                y1 = V_vor(V_test(kk),2);
                x2 = V_vor(V_test(kk+1),1);
                y2 = V_vor(V_test(kk+1),2);
                
                if((x1<wpos && x2<wpos)||(x1>Lx && x2>Lx))
                    continue;
                end
                
                if(x1<wpos || x2 < wpos)
                    y0 = (y2-y1)/(x2-x1)*(wpos-x1) + y1;
                    if(x1<wpos)
                        L = L + sqrt((wpos-x2)^2+(y0-y2)^2);
                        if(plotit)
                            plot([wpos x2],[y0 y2],'r','LineWidth',2)
                            plot([wpos x2],[y0 y2]+Ly,'g','LineWidth',2)
                            plot([wpos x2],[y0 y2]-Ly,'b','LineWidth',2)
                        end
                    else
                        L = L + sqrt((x1-wpos)^2+(y1-y0)^2);
                        if(plotit)
                            plot([x1 wpos],[y1 y0],'r','LineWidth',2)
                            plot([x1 wpos],[y1 y0]+Ly,'g','LineWidth',2)
                            plot([x1 wpos],[y1 y0]-Ly,'b','LineWidth',2)
                        end
                    end
                    
                    if(y0<0)
                        y0=y0+Ly;
                    end
                    if(y0>Ly)
                        y0=y0-Ly;
                    end
                    L_left = [L_left; y0];
                    
                elseif(x1>Lx || x2>Lx)
                    y0 = (y2-y1)/(x2-x1)*(Lx-x1) + y1;
                    if(x1>Lx)
                        L = L + sqrt((Lx-x2)^2+(y0-y2)^2);
                        if(plotit)
                            plot([Lx x2],[y0 y2],'r','LineWidth',2)
                            plot([Lx x2],[y0 y2]+Ly,'g','LineWidth',2)
                            plot([Lx x2],[y0 y2]-Ly,'b','LineWidth',2)
                        end
                    else
                        L = L + sqrt((x1-Lx)^2+(y1-y0)^2);
                        if(plotit)
                            plot([x1 Lx],[y1 y0],'r','LineWidth',2)
                            plot([x1 Lx],[y1 y0]+Ly,'g','LineWidth',2)
                            plot([x1 Lx],[y1 y0]-Ly,'b','LineWidth',2)
                        end
                    end
                    
                    if(y0<0)
                        y0=y0+Ly;
                    end
                    if(y0>Ly)
                        y0=y0-Ly;
                    end
                    L_right = [L_right; y0];
                    
                    
                else
                    L = L + sqrt((x1-x2)^2+(y1-y2)^2);
                    if(plotit)
                        plot([x1 x2],[y1 y2],'r','LineWidth',2)
                        plot([x1 x2],[y1 y2]+Ly,'g','LineWidth',2)
                        plot([x1 x2],[y1 y2]-Ly,'b','LineWidth',2)
                    end
                end
            end
        end
    end
    
    % find out which cell (wpos,0) and (Lx,0) belongs to
    L_left = sort(L_left);
    if(mod(length(L_left),2)==1)
        for vi = 1:length(L_left)-1
            if(L_left(vi+1)-L_left(vi)<0.0001)
                L_left(vi)=[];
                break;
            end
        end
    end
    L_right= sort(L_right);
    if(mod(length(L_right),2)==1)
        for vi = 1:length(L_right)-1
            if(L_right(vi+1)-L_right(vi)<0.0001)
                L_right(vi)=[];
                break;
            end
        end
    end
    
    for jj=1:nV
        VVV = V_vor(r_vor{jj+nV},1);
        UUU = V_vor(r_vor{jj+nV},2);
        VVV = VVV(isfinite(VVV));
        UUU = UUU(isfinite(UUU));
        if(inpolygon(wpos,0,VVV,UUU)...
            ||inpolygon(wpos,-Ly,VVV,UUU)...
            ||inpolygon(wpos,Ly,VVV,UUU))
            if(jj>tN)
                L = L + sum(L_left(2:2:end)-L_left(1:2:end));
                if(plotit)
                    plot(ones(2,length(L_left)/2)*wpos,[L_left(2:2:end) L_left(1:2:end)].','r','LineWidth',2)
                    plot(ones(2,length(L_left)/2)*wpos,[L_left(2:2:end) L_left(1:2:end)].'+Ly,'g','LineWidth',2)
                    plot(ones(2,length(L_left)/2)*wpos,[L_left(2:2:end) L_left(1:2:end)].'-Ly,'b','LineWidth',2)
                end
            else
                L_left = sort([0;L_left;Ly]);
                
                L = L + sum(L_left(2:2:end)-L_left(1:2:end));
                if(plotit)
                    plot(ones(2,length(L_left)/2)*wpos,[L_left(2:2:end) L_left(1:2:end)].','r','LineWidth',2)
                    plot(ones(2,length(L_left)/2)*wpos,[L_left(2:2:end) L_left(1:2:end)].'+Ly,'g','LineWidth',2)
                    plot(ones(2,length(L_left)/2)*wpos,[L_left(2:2:end) L_left(1:2:end)].'-Ly,'b','LineWidth',2)
                end
            end
            break;
        end
    end
    
    for jj=1:nV
        VVV = V_vor(r_vor{jj+nV},1);
        UUU = V_vor(r_vor{jj+nV},2);
        VVV = VVV(isfinite(VVV));
        UUU = UUU(isfinite(UUU));
        if(inpolygon(Lx,0,VVV,UUU)...
            ||inpolygon(Lx,-Ly,VVV,UUU)...
            ||inpolygon(Lx,Ly,VVV,UUU))
            if(jj>tN)
                L = L + sum(L_right(2:2:end)-L_right(1:2:end));
                if(plotit)
                    plot(ones(2,length(L_right)/2)*Lx,[L_right(2:2:end) L_right(1:2:end)].','r','LineWidth',2)
                    plot(ones(2,length(L_right)/2)*Lx,[L_right(2:2:end) L_right(1:2:end)].'+Ly,'g','LineWidth',2)
                    plot(ones(2,length(L_right)/2)*Lx,[L_right(2:2:end) L_right(1:2:end)].'-Ly,'b','LineWidth',2)
                end
            else
                L_right = sort([0;L_right;Ly]);
                
                L = L + sum(L_right(2:2:end)-L_right(1:2:end));
                if(plotit)
                    plot(ones(2,length(L_right)/2)*Lx,[L_right(2:2:end) L_right(1:2:end)].','r','LineWidth',2)
                    plot(ones(2,length(L_right)/2)*Lx,[L_right(2:2:end) L_right(1:2:end)].'+Ly,'g','LineWidth',2)
                    plot(ones(2,length(L_right)/2)*Lx,[L_right(2:2:end) L_right(1:2:end)].'-Ly,'b','LineWidth',2)
                end
            end
            break;
        end
    end
    if(plotit)
        axis equal;
        ax = gca;
        ax.XTick = [];
        ax.YTick = [];
        ax.XLim = [-0.3 1.1]*Lx;
        ax.YLim = [-0.1 1.1]*Ly;
    end
    L_sur(ii) = L;
    
    if(videoit)
        writeVideo(vobj,getframe);
    end
end

if(videoit)
    close(vobj);
end
if(datait)
    data_name = append(pathname,'vT/',sim_name,'_VT.mat');
    save(data_name,'L_sur');
end