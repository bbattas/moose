function [eta,sumetasqu, grain_area] = anisotropic_grgr_incl_in_gamma(invoer,n,n_old,eta,save_freq,grain_area)

%input
%invoer = name of input file with all input variable values
%n = the simulation runs till time step n
%n_old : the simulations starts at time step n_old+1
%eta : cell of 2-d matrices with the values of the order parameters at time
%step n_old
%save_freq : intermediate results will be save every save_freq time steps
%grain_area : vector with grain_area values stored for the previous n_old
%time steps or an empty vector

%output
%eta : cell of 2-d matrices with the values of the order parameters at time
%step n
%sumetasqu : sum of the squared order parameter values at time step n
%grain_area : every save_freq time steps the number of grid points in grain
%represented by eta{1} is counted and stored in the vector grain_area

[kappa,g,L,m,deltax,deltat,order_par,bound_cond,delta_s,delta_s_kin,n_fold,n_fold_kin,ofset_angle,ofset_angle_kin,anis] = feval(invoer); % evaluation of the input file

%%coefficients for interpolation of g-function
p1 =      -3.0944 ;
p2 =    -1.8169;
p3 =       10.323;
p4 =      -8.1819;
p5 =       2.0033;


N = size(eta{1});
k1 = L*deltat;
k1_average = sum(sum(k1))/(0.5*order_par * (order_par - 1));

e = ones(N(1),1);
laplace_matrix = (1/deltax^2)*spdiags([e -2*e e], -1:1, N(1), N(1));
derivative_matrixr = (1/deltax)*spdiags([-1*e 1*e ], 0:1, N(1), N(1));
derivative_matrixl = (1/deltax)*spdiags([-1*e 1*e ], -1:0, N(1), N(1));
derivative_matrixc = 0.5*(1/deltax)*spdiags([-1*e 1*e ], -1:2:1, N(1), N(1));
e = ones(N(2),1);
laplace_matrix2 = (1/deltax^2)*spdiags([e -2*e e], -1:1, N(2), N(2));
derivative_matrixr2 = (1/deltax)*spdiags([-1*e 1*e ], 0:1, N(2), N(2));
derivative_matrixl2 = (1/deltax)*spdiags([-1*e 1*e ], -1:0, N(2), N(2));
derivative_matrixc2 = 0.5*(1/deltax)*spdiags([-1*e 1*e ], -1:2:1, N(2), N(2));


if bound_cond == 'pp' % only periodic boundary conditions were considered so far
    
    laplace_matrix(1,N(1))=1/deltax^2; laplace_matrix(N(1),1)=1/deltax^2;
    derivative_matrixr(N(1),1) = 1/deltax;
    derivative_matrixl(1,N(1)) = -1/deltax;
    derivative_matrixc(1,N(1)) = -0.5*(1/deltax);
    derivative_matrixc(N(1),1) = 0.5*(1/deltax);
    laplace_matrix2(1,N(2))=1/deltax^2; laplace_matrix2(N(2),1)=1/deltax^2;
    derivative_matrixr2(N(2),1) = 1/deltax;
    derivative_matrixl2(1,N(2)) = -1/deltax;
    derivative_matrixc2(1,N(2)) = -0.5*(1/deltax);
    derivative_matrixc2(N(2),1) = 0.5*(1/deltax);
    % elseif bound_cond == 'pn' % the program could still be extended to
    % situations with other types of boundary conditions
    % elseif bound_cond == 'nn'
    % elseif bound_cond == 'dd'
end
clear e

for k = 1 + n_old:n
    
    if mod(k,50) == 0 % every 50 times steps the time step is written on the screen
        fprintf('iterationstep %g \n',k)
    end
    
    eta_kwad = cell(order_par,1);
    eta_kwad{1} = eta{1}.^2;
    sumetasqu = eta_kwad{1};
    for ii = 2 : order_par
        eta_kwad{ii} = eta{ii}.^2;
        sumetasqu = sumetasqu + eta_kwad{ii};
    end
    
    sumetacrosssqu = zeros(N(1),N(2));
    sumetacrosssqu_k1ij = zeros(N(1),N(2));
    for ii = 1 : order_par
        for jj = ii + 1 : order_par
            
            extra = eta_kwad{ii}.*eta_kwad{jj};
            sumetacrosssqu = sumetacrosssqu + extra;
            phi_2 = (derivative_matrixc2*(eta{ii}-eta{jj})')';
            phi_1 = derivative_matrixc*(eta{ii}-eta{jj});
            theta_angle = (pi/2)*ones(N(1),N(2));
            aa = find((abs(phi_1) > 1e-6) & (abs(phi_2) > 1e-6));
            theta_angle(aa) = (atan(phi_2(aa)./phi_1(aa)));
            theta_angle_kin = mod((theta_angle + pi/2 + ofset_angle_kin(ii,jj)),pi) -pi/2;
            theta_angle = mod((theta_angle + pi/2 + ofset_angle(ii,jj)),pi) -pi/2;
            if anis(ii,jj) == 'w' %f_sigma = (1+delta_sigma*cos(n*(theta+theta_ofset)))
                g_g = 1 + delta_s(ii,jj)*(cos(n_fold(ii,jj)*theta_angle));
            elseif anis(ii,jj) == 'c' %cubes scaled : f_sigma = (1/1+delta_sigma)*(1+delta_sigma(abs(cos(theta+theta_ofset))+abs(sin(theta+theta_ofset))))
                g_g = (1/(1+delta_s(ii,jj)))*(1 + ...
                    delta_s(ii,jj)*((sign(cos(theta_angle))).*cos(theta_angle)...
                    + (sign(sin(theta_angle))).*sin(theta_angle)));
            elseif anis(ii,jj) == 's' % cubes not scaled : f_sigma = (1+delta_sigma(abs(cos(theta+theta_ofset))+abs(sin(theta+theta_ofset))))
                g_g =  (1 + ...
                    delta_s(ii,jj)*((sign(cos(theta_angle))).*cos(theta_angle)...
                    + (sign(sin(theta_angle))).*sin(theta_angle)));
            end
            
            kin_anis = 1 + delta_s_kin(ii,jj)*(cos(n_fold_kin(ii,jj)*theta_angle_kin)); % f_mu = (1+delta_kin*cos(n*(theta+theta_ofset_kin)))
            sumetacrosssqu_k1ij = sumetacrosssqu_k1ij + k1(ii,jj)*kin_anis.*g_g.*extra;
            
        end
    end
    
    
    interface = find(sumetacrosssqu > 1e-6);
    k1_ij = (k1_average)*ones(N(1),N(2));
    k1_ij(interface) = (sumetacrosssqu_k1ij(interface))./sumetacrosssqu(interface);
    
    
    for ii = 1 : order_par
        
        sumetasqu_gamma = zeros(N(1),N(2));
        sumetasqu_dgamma1 = zeros(N(1),N(2));
        sumetasqu_dgamma2 = zeros(N(1),N(2));
        
        for jj = 1 : ii - 1
            phi_2 = (derivative_matrixc2*(eta{jj}-eta{ii})')';
            phi_1 = derivative_matrixc*(eta{jj}-eta{ii});
            theta_angle = (pi/2)*ones(N(1),N(2));
            aa = find((abs(phi_1) > 1e-6) & (abs(phi_2) > 1e-6));
            theta_angle(aa) = (atan(phi_2(aa)./phi_1(aa)));
            theta_angle = mod((theta_angle + pi/2 + ofset_angle(jj,ii)),pi) -pi/2;
            if anis(jj,ii) == 'w' %f_sigma = (1+delta_sigma*cos(n*(theta+theta_ofset)))
                g_g = g(jj,ii)*(1 + delta_s(jj,ii)*(cos(n_fold(jj,ii)*theta_angle)));
                g_dg = - g(jj,ii)*(n_fold(jj,ii)*delta_s(jj,ii)*sin(n_fold(jj,ii)*theta_angle));
            elseif anis(jj,ii) == 'c' %cubes scaled : f_sigma = (1/1+delta_sigma)*(1+delta_sigma(abs(cos(theta+theta_ofset))+abs(sin(theta+theta_ofset))))
                g_g = g(jj,ii)*(1/(1+delta_s(jj,ii)))*(1 + ...
                    delta_s(jj,ii)*((sign(cos(theta_angle))).*cos(theta_angle) ...
                    + (sign(sin(theta_angle))).*sin(theta_angle)));
                g_dg = g(jj,ii)*(1/(1+delta_s(jj,ii)))*delta_s(jj,ii)*...
                    ((sign(cos(theta_angle))).*(-sin(theta_angle))...
                    + (sign(sin(theta_angle))).*cos(theta_angle));
            elseif anis(jj,ii) == 's' % cubes not scaled : f_sigma = (1+delta_sigma(abs(cos(theta+theta_ofset))+abs(sin(theta+theta_ofset))))
                g_g = g(jj,ii)*(1 + ...
                    delta_s(jj,ii)*((sign(cos(theta_angle))).*cos(theta_angle) ...
                    + (sign(sin(theta_angle))).*sin(theta_angle)));
                g_dg= g(jj,ii)*delta_s(jj,ii).*...
                    ((sign(cos(theta_angle))).*(-sin(theta_angle))...
                    + (sign(sin(theta_angle))).*cos(theta_angle));
            end
            g_g_sq = g_g.^2;
            derivative_polynomial = 4*p1*g_g_sq.^3 + 3*p2*g_g_sq.^2 + 2*p3*g_g_sq + p4;
            gamma_inverse = p1*g_g_sq.^4 + p2*g_g_sq.^3 + p3*g_g_sq.^2 + p4*g_g_sq + p5;
            gamma = 1./gamma_inverse;
            sumetasqu_gamma = sumetasqu_gamma + gamma.*eta_kwad{jj};
            theta_angle_1 = zeros(N(1),N(2));
            theta_angle_2 = zeros(N(1),N(2));
            unit = phi_1.^2 + phi_2.^2;
            extra = eta_kwad{ii}.*eta_kwad{jj};
            aa = find(unit > 1e-4); % only for unit > 1e-4 we divide by unit;
            %otherwise it is numerically not stable;
            
            theta_angle_1 = -phi_2; % this value is only kept for unit < 1e-4; this only occurs when well within the bulk,
            %where the evaluation of the interface energy will not affect the behavior
            theta_angle_2 = phi_1;
            theta_angle_1(aa) = -(phi_2(aa)./unit(aa)).*extra(aa); % at boundaries, where unit > 1e-4, we fully evaluate theta_angle_1 and _2
            theta_angle_2(aa) = (phi_1(aa)./unit(aa)).*extra(aa);
            
            
            sumetasqu_dgamma1 = sumetasqu_dgamma1 +  (-2*(gamma.^2).*derivative_polynomial.*g_g.*g_dg).*(-theta_angle_1);
            sumetasqu_dgamma2 = sumetasqu_dgamma2 +  (-2*(gamma.^2).*derivative_polynomial.*g_g.*g_dg).*(-theta_angle_2);
        end
        for jj = ii + 1 : order_par
            
            phi_2 = (derivative_matrixc2*(eta{ii}-eta{jj})')';
            phi_1 = derivative_matrixc*(eta{ii}-eta{jj});
            theta_angle = (pi/2)*ones(N(1),N(2));
            aa = find((abs(phi_1) > 1e-6) & (abs(phi_2) > 1e-6));
            theta_angle(aa) = (atan(phi_2(aa)./phi_1(aa)));
            theta_angle = mod((theta_angle + pi/2 + ofset_angle(ii,jj)),pi) -pi/2;
            if anis(ii,jj) == 'w' %f_sigma = (1+delta_sigma*cos(n*(theta+theta_ofset)))
                g_g = g(ii,jj)*(1 + delta_s(ii,jj)*(cos(n_fold(ii,jj)*theta_angle)));
                g_dg = - g(ii,jj)*(n_fold(ii,jj)*delta_s(ii,jj)*sin(n_fold(ii,jj)*theta_angle));
            elseif anis(ii,jj) == 'c' %cubes scaled : f_sigma = (1/1+delta_sigma)*(1+delta_sigma(abs(cos(theta+theta_ofset))+abs(sin(theta+theta_ofset))))
                g_g = g(ii,jj)*(1/(1+delta_s(ii,jj)))*(1 + ...
                    delta_s(ii,jj)*((sign(cos(theta_angle))).*cos(theta_angle)...
                    + (sign(sin(theta_angle))).*sin(theta_angle)));
                g_dg = g(ii,jj)*(1/(1+delta_s(ii,jj)))*delta_s(ii,jj)*...
                    ((sign(cos(theta_angle))).*(-sin(theta_angle)) ...
                    + (sign(sin(theta_angle))).*cos(theta_angle));
            elseif anis(ii,jj) == 's' % cubes not scaled : f_sigma = (1+delta_sigma(abs(cos(theta+theta_ofset))+abs(sin(theta+theta_ofset))))
                g_g = g(ii,jj)*(1 + ...
                    delta_s(ii,jj)*((sign(cos(theta_angle))).*cos(theta_angle)...
                    + (sign(sin(theta_angle))).*sin(theta_angle)));
                g_dg = g(ii,jj)*delta_s(ii,jj).*...
                    ((sign(cos(theta_angle))).*(-sin(theta_angle)) ...
                    + (sign(sin(theta_angle))).*cos(theta_angle));
            end
            g_g_sq = g_g.^2;
            derivative_polynomial = 4*p1*g_g_sq.^3 + 3*p2*g_g_sq.^2 + 2*p3*g_g_sq + p4;
            gamma_inverse = p1*g_g_sq.^4 + p2*g_g_sq.^3 + p3*g_g_sq.^2 + p4*g_g_sq + p5;
            gamma = 1./gamma_inverse;
            sumetasqu_gamma = sumetasqu_gamma + gamma.*eta_kwad{jj};
            theta_angle_1 = zeros(N(1),N(2));
            theta_angle_2 = zeros(N(1),N(2));
            unit = phi_1.^2 + phi_2.^2;
            extra = eta_kwad{ii}.*eta_kwad{jj};
            unit = (m/kappa)*g_g;
            aa = find(unit>1e-4);% only for unit > 1e-4 we divide by unit;
            %otherwise it is numerically not stable;
            
            theta_angle_1 = -phi_2; % this value is only kept for unit < 1e-4; this only occurs when well within the bulk,
            %where the evaluation of the interface energy will not affect the behavior
            theta_angle_2 = phi_1;
            theta_angle_1(aa) = -(phi_2(aa).*extra(aa))./unit(aa); % at boundaries, where unit > 1e-4, we fully evaluate theta_angle_1 and _2
            theta_angle_2(aa) = phi_1(aa).*extra(aa)./unit(aa);
            
            sumetasqu_dgamma1 = sumetasqu_dgamma1 +  (-2*(gamma.^2).*derivative_polynomial.*g_g.*g_dg).*theta_angle_1;
            sumetasqu_dgamma2 = sumetasqu_dgamma2 +  (-2*(gamma.^2).*derivative_polynomial.*g_g.*g_dg).*theta_angle_2;
            
        end
        
        gradient_anis = - m*(derivative_matrixc*sumetasqu_dgamma1 + ...
            (derivative_matrixc2*(sumetasqu_dgamma2'))');
        
        gradient_energy = -kappa.*(laplace_matrix*eta{ii} + (laplace_matrix2*eta{ii}')') ...
            + gradient_anis;
        bulk = m.*eta{ii}.*(eta_kwad{ii} - 1 + 2*sumetasqu_gamma);
        
        eta{ii} = eta{ii} - k1_ij.*(bulk + gradient_energy);
    end
    
    if mod(k,save_freq)==0
        microstructure = [];
        grain_area = [grain_area sum(sum(eta{1}))];
        for s = 1 : N(1)
            for t = 1 : N(2)
                max_val = abs(eta{1}(s,t));
                microstructure(s,t) = 1;
                for kk = 2 : order_par
                    if abs(eta{kk}(s,t)) > max_val;
                        max_val = abs(eta{kk}(s,t));
                        microstructure(s,t) = kk;
                    end
                end
            end
        end
        save(strcat(invoer,'nn',int2str(k)),'sumetasqu','eta','microstructure','grain_area')
    end
end

microstructure = [];
for s = 1 : N(1)
    for t = 1 : N(2)
        max_val = abs(eta{1}(s,t));
        microstructure(s,t) = 1;
        for kk = 2 : order_par
            if abs(eta{kk}(s,t)) > max_val;
                max_val = abs(eta{kk}(s,t));
                microstructure(s,t) = kk;
            end
        end
    end
end

sumetasqu = eta_kwad{1};
for ii = 2 : order_par
    sumetasqu = sumetasqu + eta_kwad{ii};
end

save(strcat(invoer,'einde',int2str(n)),'eta','sumetasqu','microstructure','grain_area')
%sumetasqu can be used to represent the grain structure; it is 1 within
%grains and smaller than 1 at grain boundaries.
%microstructure can be used as a sharp interface representation of the
%grain structures; microstructure gives at each grid point the number of
%the phase-field variable with the highest value.
