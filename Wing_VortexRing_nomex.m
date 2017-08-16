function Wing_VortexRing
    % Wing input Data
    b = .25*8; % Wingspan
    cr = .25; % Chord at root
    ct = .25; % Chord at tip
    s = 0; % Sweep angle
    alpha = 0.1; % Initial angle of attack
    CL_target = 0.85; % Required CL

    % Plot and output
    doplot = true;
    plotwake = false;
    
    % Discretization
    x_w = 100; % Distance of the wake as times c_r
    n_w = 1;
    n_y = 50;
    n_x = 30;
    max_iter = 15;
    
    % Pre calculationad and memory allocation
    S = b*(cr+ct)/2;
    Ar = b^2/S;
    
    coord = zeros(n_x,n_y,4,3);
    ring = zeros(n_x,n_y,4,3);
    ring_wk = zeros(n_w,n_y,4,3);
    collp = zeros(n_x,n_y,3);
    A = zeros(n_x*n_y);
    rhs = zeros(n_x*n_y,1);
    
    close all
    hold on
    disp('Assembling wing geometry...');
    for i=1:n_y
        y1 = (b/2)*(i-1)/n_y;
        y2 = (b/2)*i/n_y;
        y_m = (y1+y2)/2;
        for u=1:n_x
            x_ac_1 = cr/4+y1*tan(s);
            x_ac_2 = cr/4+y2*tan(s);
            x_ac_m = cr/4+y_m*tan(s);
            c1 = (ct-cr)*y1*2/b+cr;
            c2 = (ct-cr)*y2*2/b+cr;
            c_m = (ct-cr)*y_m*2/b+cr;
            x_max_1 = x_ac_1-c1*.25;
            x_max_2 = x_ac_2-c2*.25;
            x_max_m = x_ac_m-c_m*.25;
            x_m = x_max_m+c_m*(u-1)/n_x;
            
            % Panel coordinates
            coord(u,i,1,:) = [x_max_1+c1*(u-1)/n_x y1 0];
            coord(u,i,2,:) = [x_max_2+c2*(u-1)/n_x y2 0];
            coord(u,i,3,:) = [x_max_2+c2*(u)/n_x y2 0];
            coord(u,i,4,:) = [x_max_1+c1*(u)/n_x y1 0];
            
            % Ring coordinates
            ring(u,i,1,:) = [coord(u,i,1,1)+c1/(n_x*4) y1 0];
            ring(u,i,2,:) = [coord(u,i,2,1)+c2/(n_x*4) y2 0];
            ring(u,i,3,:) = [coord(u,i,3,1)+c2/(n_x*4) y2 0];
            ring(u,i,4,:) = [coord(u,i,4,1)+c1/(n_x*4) y1 0];
            
            % Collocation points
            collp(u,i,:) = [x_m+c_m*3/(n_x*4) y_m 0];
            
            if doplot
                patch(squeeze(coord(u,i,:,1)),squeeze(coord(u,i,:,2)),'w');
                patch(squeeze(ring(u,i,:,1)),squeeze(ring(u,i,:,2)),'g','EdgeColor','green','FaceColor','none','LineStyle','--');
                plot(collp(u,i,1),collp(u,i,2),'kx');
            end
        end
        for u=1:n_w
            % Wake Rings
            ring_wk(u,i,1,:) = [ring(n_x,i,4,1)+x_w*cr*(u-1)/n_w y1 0];
            ring_wk(u,i,2,:) = [ring(n_x,i,3,1)+x_w*cr*(u-1)/n_w y2 0];
            ring_wk(u,i,3,:) = [ring(n_x,i,3,1)+x_w*cr*u/n_w y2 0];
            ring_wk(u,i,4,:) = [ring(n_x,i,4,1)+x_w*cr*u/n_w y1 0];
            if doplot & plotwake
                patch(squeeze(ring_wk(u,i,:,1)),squeeze(ring_wk(u,i,:,2)),'b','EdgeColor','blue','FaceColor','none','LineStyle','--');
            end
        end
    end
    if doplot, axis equal; xlabel('x'); ylabel('y'); title('Discretization'); end
    
    disp('Assembling matrix...');
    tic;
    t = toc;
    % For each panel
    for ei=1:n_y
        if toc-t>1, fprintf('%.1f%%, remaining %.1f s\n', ei*100/n_y, toc*(n_y/ei-1)) ;t = toc; end
        for eu=1:n_x
            e_index = eu+n_x*(ei-1);
                % For each control point
                for ci=1:n_y
                    for cu=1:n_x
                        % Compute the collocation point index
                        c_index = cu+n_x*(ci-1);
                        % Compute the influence of the panel on the same
                        % semi.wing
                        vi1 = voring4(ring(eu,ei,:,:),collp(cu,ci,:));
                        % Add contribution of the other semi.wing
                        collp_aux(1,1,:) = [collp(cu,ci,1) -collp(cu,ci,2) collp(cu,ci,3)];
                        vi2 = voring4(ring(eu,ei,:,:),collp_aux(1,1,:));
                        vi2(2) = -vi2(2);
                        % Add to influence matrix
                        A(c_index,e_index) = (vi1+vi2)*[0 0 1]';
                        % Add to RHS
                        %rhs(c_index,1) = -u_inf*[0 0 1]';
                    end
                end
        end
        for eu=1:n_w
                % For each control point
                for ci=1:n_y
                    for cu=1:n_x
                        % Compute the collocation point index
                        c_index = cu+n_x*(ci-1);
                        % Compute the influence of the panel on the same
                        % semi.wing
                        vi1 = voring4(ring_wk(eu,ei,:,:),collp(cu,ci,:));
                        % Add contribution of the other semi.wing
                        collp_aux(1,1,:) = [collp(cu,ci,1) -collp(cu,ci,2) collp(cu,ci,3)];
                        vi2 = voring4(ring_wk(eu,ei,:,:),collp_aux(1,1,:));
                        vi2(2) = -vi2(2);
                        % Add to influence matrix
                        A(c_index,n_x) = A(c_index,n_x) + (vi1+vi2)*[0 0 1]';
                    end
                end
        end
    end
    
    % Solve the system of equations
    fprintf('Solving system of equations...');
    CL_ant = inf;
    for i=1:max_iter
        u_inf = [cos(alpha) 0 sin(alpha)];
        rhs(:) = -u_inf*[0 0 1]';
        Gamma = A\rhs;
        % Final computations and plots
        G = reshape(Gamma,n_x,n_y);
        Cl = 2*G(n_x,:)/cr;
        CL = 2*sum(G(n_x,:))*(b/(2*n_y))*(2/S);
        CLa = CL/alpha;
        alpha = CL_target/CLa;
        if abs(CL-CL_ant)<1e-3; break; end
        CL_ant = CL;
    end
    fprintf('done with %i iterations!\n',i);
    
    % Drag Estimation using Trefft plane
    G_vort = [0 -diff(G(n_x,:)) G(n_x,n_y)];
    cp_pos = [(ring(n_x,:,4,2)+ring(n_x,:,3,2))/2;...
              (ring(n_x,:,4,3)+ring(n_x,:,3,3))/2]';
    vort_pos = [ring(n_x,:,4,2) ring(n_x,n_y,3,2);...
                ring(n_x,:,4,3) ring(n_x,n_y,3,3)]';
    wi = zeros(1,n_y);
    for i1=1:n_y
        for i2=1:n_y+1
            % Contribution of trailing vortex of the same side
            Ay = vort_pos(i2,1)-cp_pos(i1,1);
            Az = vort_pos(i2,2)-cp_pos(i1,2);
            if i1~=i2
                wi(i1) = wi(i1) + (-G_vort(i2)/(2*pi))*Ay/(Ay^2+Az^2);
            end
            % Contribution of trailing vortex of the other side
            if i2>1
                % Trailing votex is on the other side
                Ay = -vort_pos(i2,1)-cp_pos(i1,1);
                % Gamma is the opposite
                wi(i1) = wi(i1) + (+G_vort(i2)/(2*pi))*Ay/(Ay^2+Az^2);
            end
        end
    end
    Cd = -(wi.*G(n_x,:))*(b/(2*n_y*S));
    CDi = -2*sum((wi.*G(n_x,:)))*(b/(2*n_y*S));
    e = CL^2/(pi*Ar*CDi);
    
    % Pressure coefficient estimation
    for ci=1:n_y
        cu = 1;
        c_mean = (coord(cu,ci,4,1)-coord(cu,ci,1,1)+coord(cu,ci,3,1)-coord(cu,ci,2,1))/2;
        cp(cu,ci) = 2*G(cu,ci)/c_mean;
        for cu=2:n_x
            c_mean = (coord(cu,ci,4,1)-coord(cu,ci,1,1)+coord(cu,ci,3,1)-coord(cu,ci,2,1))/2;
            cp(cu,ci) = 2*(G(cu,ci)-G(cu-1,ci))/c_mean;
        end
    end
    if doplot
        figure
        cmap = colormap;
        cp_max = max(cp(:));
        for ci=1:n_y
            for cu=1:n_x
                patch(squeeze(coord(cu,ci,:,1)),squeeze(coord(cu,ci,:,2)),cp(cu,ci),'EdgeColor','none');
            end
        end
        colorbar;
        xlabel('x'); ylabel('y'); title('\Delta Cp Distribution');
        axis equal;
    end
    
    % CD and CL plots
    figure
    yyaxis left
    plot(linspace(1/(2*n_y),1-1/(2*n_y),n_y),Cl/CL)
    xlabel('2y/b'); ylabel('Cl/CL'); axis([0 1 0 1.4])
    grid on
    yyaxis right
    plot(linspace(1/(2*n_y),1-1/(2*n_y),n_y),Cd/CDi)
    ylabel('Cd/CD');
    
    figure;
    plot(linspace(0,1,n_x),2*[G(1,1) diff(G(1:end,1))']*n_x/cr);
    grid on
    xlabel('x/c'); ylabel('\Delta Cp'); title('Cp at central section')
    
    fprintf(['Aspect Ratio = ' num2str(Ar) '\n']);
    fprintf(['Alpha = ' num2str(alpha) 'rad\n']);
    fprintf(['CL = ' num2str(CL) '\n']);
    fprintf(['CL_alpha = ' num2str(CLa) '\n']);
    fprintf(['CDi = ' num2str(CDi) '\n']);
    fprintf(['Efficiency Factor e = ' num2str(e) '\n']);
    fprintf(['CL/CDi = ' num2str(CL/CDi) '\n']);
    hold off
end

function [v,v2,v3] = voring4(r,p)
    persistent in1 in2 in3 in4
    if isempty(in1)
        in1 = sub2ind([1 1 4 3],[1 1 1],[1 1 1],[1 1 1],[1 2 3]);
        in2 = sub2ind([1 1 4 3],[1 1 1],[1 1 1],[2 2 2],[1 2 3]);
        in3 = sub2ind([1 1 4 3],[1 1 1],[1 1 1],[3 3 3],[1 2 3]);
        in4 = sub2ind([1 1 4 3],[1 1 1],[1 1 1],[4 4 4],[1 2 3]);
    end
    
    p2 = [p(1) p(2) p(3)]';
    v2 = vortexl(p2,r(in2),r(in3));
    v3 = vortexl(p2,r(in3),r(in4));
    v = vortexl(p2,r(in1),r(in2))+...
        v2+...
        v3+...
        vortexl(p2,r(in4),r(in1));
end

function v = vortexl(p,v1,v2)
    n1 =[p(1)-v1(1) p(2)-v1(2) p(3)-v1(3)];
    n2 =[p(1)-v2(1) p(2)-v2(2) p(3)-v2(3)];
    cp = [n1(2)*n2(3)-n1(3)*n2(2),...
          n1(3)*n2(1)-n1(1)*n2(3),...
          n1(1)*n2(2)-n1(2)*n2(1)];  
    r1 = sqrt(n1(1)^2+n1(2)^2+n1(3)^2);
    r2 = sqrt(n2(1)^2+n2(2)^2+n2(3)^2);
    cp2 = cp(1)^2+cp(2)^2+cp(3)^2;
    vn =[v2(1)-v1(1) v2(2)-v1(2) v2(3)-v1(3)];
    r0r1 = vn(1)*n1(1)+vn(2)*n1(2)+vn(3)*n1(3);
    r0r2 = vn(1)*n2(1)+vn(2)*n2(2)+vn(3)*n2(3);
    K = (r0r1/r1-r0r2/r2)/(4*pi*cp2);
    v = cp*K;
end