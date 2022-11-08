function [U1_tumor,U2_tumor,L_con, r_dis] = tumorEnergy(x,y,r,tN,NCELLS,l1,l2,Ly,kc,nv)
%calculate tumor-tumor energy

U1_tumor=0;
U2_tumor=0;
L_con=0;
K_tumor = 0;
cut_num = 939;
x_tumor = [x{1:tN}];
x_tumor = sort(x_tumor,'descend');
x_cut = x_tumor(cut_num);

kint = l1/(l2 - l1);
r_0 = [];
r_min=1;
xy_min=[1,1;1,1];
ind_min=[0,0];
counter=0;
for ii=1:tN-1
    xi=x{ii};
    yi=y{ii};
    
    if(xi>=x_cut && counter<=cut_num)
        counter = counter + 1;
        for jj=ii+1:tN
            sij = r{ii}+r{jj};
            xj=x{jj};
            yj=y{jj};
            dx=xi-xj;
            dy=yi-yj;
            dy = dy - Ly * round(dy / Ly);
            if(dy<sij*(1+l2))
                rij = sqrt(dx * dx + dy * dy);
                if(dy<sij*(1+l2))
                    xij=rij/sij;
                    if(xij<1+l2)
                        r_0 = [r_0, xij];
                        if(xij>1+l1)
                           U1_tumor = U1_tumor - 0.5*kint*power(1.0 + l2 - xij,2.0); 
                        else
                           if(xij>1)
                               U1_tumor = U1_tumor + 0.5*(power(1.0 - xij,2.0) - l1*l2); 
                           else
                               %{
                               if xij<r_min
                                   r_min=xij;
                                   xy_min = [xi,xj;yi,yj];
                                   ind_min = [ii,jj];
                               end
                               %}
                               U2_tumor = U2_tumor + 0.5*(power(1.0 - xij,2.0) - l1*l2); 
                           end
                        end
                    end
                end
            end
        end
    end
end

for ii=tN+1:NCELLS
    xii=x{ii};
    yii=y{ii};
    for vv=1:nv(ii)
        xi=xii(vv);
        yi=yii(vv);
        
        breaker=0;
        for jj=1:tN
            sij = r{ii}+r{jj};
            xj=x{jj};
            
            if(xj > x_cut)
                yj=y{jj};
                dx=xi-xj;
                dy=yi-yj;
                dy = dy - Ly * round(dy / Ly);
                if(dy<sij)
                    rij = sqrt(dx * dx + dy * dy);
                    if(rij<sij)
                        L_con = L_con + 1;
                        breaker=1;
                    end
                end
            end
            
            
       
            if(breaker==1)
                break;
            end
            
            
        end
    end
end

%{
x_bound = [];
y_bound = [];
for ii = 1:tN
    xii=x{ii};
    yii=y{ii};
    breaker=0;
    for jj = tN+1:NCELLS
        for vv=1:nv(jj)
            xjj=x{jj}(vv);
            yjj=y{jj}(vv);
            sij=2*r{ii}+r{jj}(1);
            
            dx=xii-xjj;
            dy=yii-yjj;
            dy = dy - Ly * round(dy / Ly);
            if(dy<sij)
                rij = sqrt(dx * dx + dy * dy);
                if(rij<sij)
                    x_bound = [x_bound,xii];
                    y_bound = [y_bound,yii];
                    breaker=1;
                    break;
                end
            end
            
        end
        
        if(breaker==1)
            break;
        end
    end
end
nB=length(x_bound);
List = zeros(nB,2);
Dist = [ones(nB,1)*10,ones(nB,1)*9];
for ii=1:nB
    xii=x_bound(ii);
    yii=y_bound(ii);
    for jj = 1:nB
        if(jj~=ii)
            xjj=x_bound(jj);
            yjj=y_bound(jj);
            dx=xii-xjj;
            dy=yii-yjj;
            dy = dy - Ly * round(dy / Ly);
            if(dy<sij)
                rij = sqrt(dx * dx + dy * dy);
                if(rij==0)
                    ii
                    jj
                end
                if(rij<Dist(ii,1))
                    if(rij<Dist(ii,2))
                        Dist(ii,2)=rij;
                        List(ii,2)=jj;
                    else
                        Dist(ii,1)=rij;
                        List(ii,1)=jj;
                    end
                end
            end
        end
    end
end

figure;hold on
for ii=1:nB
    plot([x_bound(ii),x_bound(List(ii,1))],[y_bound(ii),y_bound(List(ii,1))],'k')
    plot([x_bound(ii),x_bound(List(ii,2))],[y_bound(ii),y_bound(List(ii,2))],'k')
end
axis equal;
%}
U1_tumor = U1_tumor*kc;
U2_tumor = U2_tumor*kc;
r_dis=r_0;
end
