%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Bezier_Plotter
%
% Input:  p_1 = polynomial degree in direction 1
%         p_2 = polynomial degree in direction 2
%         p_3 = polynomial degree in direction 3
%         n_el = number of elements
%         P_b = array of Bezier control points
%         w_b = array of Bezier weights
%
% Output: N/A
%
% Purpose: Plot a given Bezier surface
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Bezier_Plotter(p_1,p_2,p_3,n_el,P_b,w_b)

n_loc = (p_1+1)*(p_2+1)*(p_3+1);

itr = 8;

x_plot = 0:1/(itr-1):1;
y_plot = 0:1/(itr-1):1;
z_plot = 0:1/(itr-1):1;

for i1 = 1:itr
    B1 = bernstein_poly(p_1,x_plot(i1));
    for i2 = 1:itr
        B2 = bernstein_poly(p_2,y_plot(i2));
        for i3 = 1:itr     
            B3 = bernstein_poly(p_3,z_plot(i3));
            
            for j1 = 1:(p_1+1)
                for j2 = 1:(p_2+1)
                    for j3 = 1:(p_3+1)
                        tensor_prod(j3,j2,j1) = B1(j1)*B2(j2)*B3(j3);
                    end
                end
            end
            
            temp = reshape(tensor_prod,n_loc,1);
        
            B(i1,i2,i3,:) = temp;
        end
    end    
end

figure
hold on
axis equal

for e = 1:n_el
    WX = zeros(itr,itr,itr);
    WY = zeros(itr,itr,itr);
    WZ = zeros(itr,itr,itr);
    W = zeros(itr,itr,itr);
    
    for a = 1:n_loc
        WX(:,:,:) = WX(:,:,:) + B(:,:,:,a)*P_b(a,1,e)*w_b(a,e);
        WY(:,:,:) = WY(:,:,:) + B(:,:,:,a)*P_b(a,2,e)*w_b(a,e);
        WZ(:,:,:) = WZ(:,:,:) + B(:,:,:,a)*P_b(a,3,e)*w_b(a,e);
        W(:,:,:) = W(:,:,:) + B(:,:,:,a)*w_b(a,e);
    end    
    
    X = WX./W;
    Y = WY./W;
    Z = WZ./W;
    
    X_plot = reshape(X,[itr^3,1]);
    Y_plot = reshape(Y,[itr^3,1]);
    Z_plot = reshape(Z,[itr^3,1]);

    scatter3(X_plot,Y_plot,Z_plot)
%     ids = abs(Z_plot - 1)<1e-1;
%     scatter3(X_plot(ids),Y_plot(ids),Z_plot(ids))
    scatter3(P_b(:,1,e),P_b(:,2,e),P_b(:,3,e),'MarkerEdgeColor','r','MarkerFaceColor','r');
end

view([45,20])
xlabel('x','FontName','Times New Roman','FontSize',20)
ylabel('y','FontName','Times New Roman','FontSize',20)
zlabel('z','FontName','Times New Roman','FontSize',20)
set(gca,'FontName','Times New Roman','FontSize',14)
camlight right; lighting phong


hold off