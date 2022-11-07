function [U1_tumor,U2_tumor,K_tumor,L_con] = tumorEnergy(x,y,vx,vy,r,tN,NCELLS,l1,l2,Ly,kc,nv)
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

counter=0;
for ii=1:tN-1
    xi=x{ii};
    yi=y{ii};
    
    if(xi>x_cut && counter<=cut_num)
        counter = counter + 1;
        K_tumor = K_tumor + vx{ii}^2+vy{ii}^2;
        for jj=ii:tN
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
                        if(xij>1+l1)
                           U1_tumor = U1_tumor - 0.5*kint*power(1.0 + l2 - xij,2.0); 
                        else
                           if(xij>1)
                               U1_tumor = U1_tumor + 0.5*(power(1.0 - xij,2.0) - l1*l2); 
                           else
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

U1_tumor = U1_tumor*kc;
U2_tumor = U2_tumor*kc;
K_tumor = 1/2 * K_tumor;
end