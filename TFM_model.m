%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TFM model
% version 1.0, 1 April 2020, by Pardon, G.
% Publication in progress
% Contact: Gaspard Pardon: gaspard@stanford.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TFM model was developed as a validation and benchmarking too for the 
% Streamlined TFM GUI published in Pardon, G. et al., Contra-X: a 
% streamlined and versatile pipeline for time-course analysis of the 
% contractile force of single hiPSC-cardiomyocytes, 2020, submitted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tfm_mains.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instructions:
% To analyse TFM video, execute tfm_main.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all
clear all

%load base image
base_im_frame = imread('beads_2.tif');

%===========================================
%Simulation parameters
%===========================================

%Virtual cell parameters
%Cell dimension and position parameters
cell_l = 90e-6;%cell length
cell_w = 18e-6;%cell width
%Cell main axis orientation (*****LOOPED OVER LIST*****)
cell_angle = [-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90];

%Simulated force time profile (*****CHOOSE ONE OR SEVERAL IN CELL ARRAY*****)
F_shape = {'Custom'};%Choose one of: 'Triangle';'Square';'Gaussian';'Custom';
n_peaks = 3;% number of peaks
peak_freq = 1;%frequency of peaks in Hz
peak_dr = 0.5;%Duty ratio of force: peak force/contraction cycle duration
m_force = 1e-8;%Force max amplitude
%List of force application gaussian area radius in micron (*****LOOPED OVER LIST*****)
fr = [8]*1e-6; 
F_ar = 1;%Force application guassian radius dimension aspect ratio

%Material parameters
E = 10e3;%Yongs modulus
nu = 0.45;%Poisson ratio

%Video parameters
dimx =size(base_im_frame,2);%Image dimension in pixels (*****MUST MATCH THAT OF THE MODEL IMAGE*****)
dimy = size(base_im_frame,1);
fps = 19; %frame per seconds
mu_pix = 0.27;%micron per pixel

%Image processing and analysis parameters 
%TO IMPLEMENT IN FUTURE RELEASE
%binning = 1;%1 is no binning, 2 is binning 2 by 2 px etc...
%spac_coef = 0;%enable matching the output of the piv step computed with a given spac_coef in the ContraX

%===========================================

%parameters conversion
A=(1+nu)/(pi*E);
D = mu_pix*1e-6;%conversion factor
Apx=D^2;%px area in micron
cell_x = dimx/2*D;%
cell_y = dimy/2*D;
dur_vid = n_peaks/peak_freq;% duration in seconds
dur_peak = 1/peak_freq;% duration in seconds
n_frame = dur_vid*fps;%number of frame per video
n_frame_peak = floor(fps*dur_peak);%number of frames for one contraction cycle
end_peaks = n_frame_peak;
for n = n_peaks+1
    end_peaks = [end_peaks, n_frame_peak*n_peaks];
end
n_frame_n0peak = floor(fps*dur_peak*peak_dr);%number of frames during non-zero active contraction
time_vec = linspace(0,dur_vid,n_frame);%time vector
time_vec_peak = linspace(0,dur_peak,n_frame_peak);%time vector for one contraction cycle
time_vec_n0peak = linspace(0,dur_peak*peak_dr,n_frame_n0peak);%time vector

% initialize some variables for speed
Trx = zeros(dimy,dimx, n_frame);
Try = zeros(dimy,dimx,n_frame);
fu = zeros(dimy,dimx,n_frame);
fv = zeros(dimy,dimx,n_frame);

F_tot = zeros(1,n_frame);
Fx_tot = zeros(1,n_frame);
Fy_tot = zeros(1,n_frame);

%Meashgrid
[X,Y] = meshgrid(1:dimx,1:dimy);
xdata(:,:,1) = X*D;
xdata(:,:,2) = Y*D;

%===========================================

%loop over Force time profiles
for m = 1:size(F_shape,2)
    F_type = F_shape{m};
    
    %Loop over the list of force application radius
    for n = 1:size(fr,2)
        
        %Force application guassian radius dimension aspect ratio
        fr_x = fr(n)*F_ar;
        fr_y = fr(n);
        %Norm of the force application area radiuses
        norm = 1/(2*pi*(fr_x^2+fr_y^2));
        
        %%Square force up and down in time
        sq_peak = [zeros(1,floor((n_frame_peak-n_frame_n0peak)/2)),...
            m_force*ones(1,n_frame_n0peak),...
            zeros(1,n_frame_peak-floor((n_frame_peak-n_frame_n0peak)/2)-n_frame_n0peak)];
        %%Triangle force up-down
        tri_peak = [zeros(1,floor((n_frame_peak-n_frame_n0peak)/2)),...
            (n_frame_n0peak/fps-abs(linspace(0,n_frame_n0peak/fps,n_frame_n0peak)-linspace(n_frame_n0peak/fps,0,n_frame_n0peak)))*m_force/(n_frame_n0peak/fps),...
            zeros(1,floor((n_frame_peak-n_frame_n0peak)/2))];
         %%Gaussian force peak in time
        gau_peak = m_force*exp(-((time_vec_peak-time_vec_n0peak(floor(n_frame_peak/2)))/(peak_dr/2)).^2);
        gau_peak(1) = 0; gau_peak(end) = 0;
        f_sh_vecs = {sq_peak,tri_peak,gau_peak};

            switch F_type
                
                case 'Gaussian'
                    force_vec = f_sh_vecs{1};
                    for i = 1:n_peaks-1
                        force_vec = [force_vec, f_sh_vecs{1}];
                    end
                case 'Square'
                    force_vec = f_sh_vecs{2};
                    for i = 1:n_peaks-1
                        force_vec = [force_vec, f_sh_vecs{2}];
                    end
                case 'Triangle'
                    force_vec = f_sh_vecs{3};
                    for i = 1:n_peaks-1
                        force_vec = [force_vec, f_sh_vecs{3}];
                    end
                case 'Custom'
                    force_vec = f_sh_vecs{1};
                    for i = 1:n_peaks-1
                        force_vec = [force_vec, f_sh_vecs{i+1}];
                    end
            end
            
        
        figure(1)
        plot(time_vec, force_vec) %force time profile
        %plot(linspace(0,peak_dr,length(force)), force) %force time profile
    
        %mask parameter
        phi = linspace(0,2*pi,180);
        cosphi = cos(phi);
        sinphi = sin(phi);
        
        %prepare movie figures
        
        % Plot displacement
        zlim = m_force;
        fig1 = figure(1);%'Name','Displacement magnitude');
        surf(xdata(:,:,1),xdata(:,:,2),sqrt(fu(:,:,1).^2+fv(:,:,1).^2),sqrt(fu(:,:,1).^2+fv(:,:,1).^2),'EdgeColor','none')
        %hold on
        %quiver3(xdata(:,:,1),xdata(:,:,2),-ones(dim,dim)*zlim,u,v,zeros(dim,dim),2,'k');
        axis tight manual
        ax = fig1.CurrentAxes;
        ax.XLim = [0 dimx*D];
        ax.YLim = [0 dimy*D];
        ax.ZLim = [-zlim zlim];
        ax.CLim = [0 zlim];
        ax.NextPlot = 'replaceChildren';
        
        %Plot tractions
        zlim2 = 2.5e2;
        fig2 = figure(2);%'Name','Traction magnitude');
        surf(xdata(:,:,1),xdata(:,:,2),Trx(:,:,1),'EdgeColor','none')
        %hold on
        %quiver3(xdata(:,:,1),xdata(:,:,2),-ones(dim,dim)*zlim,u,v,zeros(dim,dim),2,'k');
        axis tight manual
        ax2 = fig2.CurrentAxes;
        ax2.XLim = [0 dimx*D];
        ax2.YLim = [0 dimy*D];
        ax2.ZLim = [-zlim2 zlim2];
        ax2.CLim = [0 zlim2];
        ax2.NextPlot = 'replaceChildren';
        

        %loop over cell orientation angle
        for k = 1:size(cell_angle,2)

            
            % initialize some variables for speed
            u_fft = zeros(dimy,dimx);
            v_fft = zeros(dimy,dimx);
            Trx_save = zeros(dimy,dimx,n_frame);
            Try_save = zeros(dimy,dimx,n_frame);
            fu_save = zeros(dimy,dimx,n_frame);
            fv_save = zeros(dimy,dimx,n_frame);
            x = zeros(dimy,dimx,n_frame);
            y = zeros(dimy,dimx,n_frame);
            
            Mov1 = struct('cdata',cell(1,n_frame),'colormap',cell(1,n_frame));
            Mov2 = struct('cdata',cell(1,n_frame),'colormap',cell(1,n_frame));
            U = zeros(1,n_frame);
            av_disp = zeros(1,n_frame);
            disp2_prev = zeros(dimy,dimx,2);
            prev_frame = base_im_frame;

            
            
            %loop over force magnitude
            for l = 1:length(force_vec)
                
                %cell orientation angle
                cell_theta = cell_angle(k);
                
                            
                %calculate ellipse mask
                % Calculate cell ellipse
                R = [ cosd(cell_theta)   sind(cell_theta)
                    -sind(cell_theta)   cosd(cell_theta)];
                xy = [cell_l/2*cosphi; cell_w/2*sinphi];
                xy = R*xy;
                x_el = xy(1,:) + cell_x;
                y_el = xy(2,:) + cell_y;
                mask=poly2mask(x_el./D,y_el./D,dimy,dimx);

                
                %calculate blue ellipse
                %calculate based on predicted force propagation
                fu_pred=D;
                d=m_force*(1+nu)/(fu_pred*E*pi)*1e-6/D;
                an=0.75*sqrt(3)*cell_l/2;
                bn=1/0.75*sqrt(3)*cell_w/2;
                xyn = [an*cosphi; bn*sinphi];
                xyn = R*xyn;
                xn = xyn(1,:) + cell_x;
                yn = xyn(2,:) + cell_y;
                %generate blue ellipse mask
                mask_ellipse=poly2mask(xn/D,yn/D,dimy,dimx);
                
%                 figure('Name','Mask')
%                 imagesc(mask)
                
                %force application position
                fx_p_1 = cell_x-cell_l/2*cosd(-cell_theta);
                fy_p_1 = cell_y-cell_l/2*sind(-cell_theta);
                fx_p_2 = cell_x+cell_l/2*cosd(-cell_theta);
                fy_p_2 = cell_y+cell_l/2*sind(-cell_theta);
                
                %Distribute force as a 2D Gaussian
                para = [norm,fx_p_1,fr_x,fy_p_1,fr_y,cell_theta];
                %Rotate Gaussian axis
                xdatarot(:,:,1)= xdata(:,:,1)*cos(para(6)) - xdata(:,:,2)*sin(para(6));
                xdatarot(:,:,2)= xdata(:,:,1)*sin(para(6)) + xdata(:,:,2)*cos(para(6));
                x0rot = para(2)*cos(para(6)) - para(4)*sin(para(6));
                y0rot = para(2)*sin(para(6)) + para(4)*cos(para(6));
                %Generate normalized Gaussian stress distribution 1st end of the cell
                T = para(1)*exp(   -((xdatarot(:,:,1)-x0rot).^2/(2*para(3)^2) + (xdatarot(:,:,2)-y0rot).^2/(2*para(5)^2) )    );
                
                %Repeat for 2nd end of the cell
                para = [-norm,fx_p_2,fr_x,fy_p_2,fr_y,cell_theta];
                
                xdatarot(:,:,1)= xdata(:,:,1)*cos(para(6)) - xdata(:,:,2)*sin(para(6));
                xdatarot(:,:,2)= xdata(:,:,1)*sin(para(6)) + xdata(:,:,2)*cos(para(6));
                x0rot = para(2)*cos(para(6)) - para(4)*sin(para(6));
                y0rot = para(2)*sin(para(6)) + para(4)*cos(para(6));
                
                %Sum two stress fields
                T = T+(  para(1)*exp(   -((xdatarot(:,:,1)-x0rot).^2/(2*para(3)^2) + (xdatarot(:,:,2)-y0rot).^2/(2*para(5)^2) )    ));
                
%                 %Plot simulated stress field
%                 figure('Name','Simulated stress field')
%                 surf(T,'EdgeColor','none')
                
                %Project traction stress magnitude to x and y axis
                Trx = T*force_vec(l)*cosd(-cell_theta);
                Try = T*force_vec(l)*sind(-cell_theta);
                
                %Integrate force in the cellmask
                Fx = Apx*Trx.*mask_ellipse;
                Fy = Apx*Try.*mask_ellipse;
                Trt = sqrt( Trx.^2 + Try.^2 );
                Fr = Apx*Trt;%;.*mask;
                F_tot(l)=sum(sum(Fr))
                Fx_tot(l)=sum(sum(Fx));
                Fy_tot(l)=sum(sum(Fy));
                
                %Remove mean before fft2
                Tx1_0=(Trx)-nanmean(nanmean(Trx)).*ones(size(Trx,1),size(Trx,2));
                Ty1_0=(Try)-nanmean(nanmean(Try)).*ones(size(Try,1),size(Try,2));
                
                %Calculate the Fourier transform of the traction
                Tx_k = fft2(Tx1_0);
                Ty_k = fft2(Ty1_0);
                
                % Generate Fourier space coordinates vectors
                Nx=size(Tx_k,2);
                Ny=size(Ty_k,1);
                dkx = 1/(Nx*D);
                dky = 1/(Ny*D);
                kx = [0:fix(Nx/2)-1,-fix(Nx/2):-1]*dkx*2*pi;
                ky = [0:fix(Ny/2)-1,-fix(Ny/2):-1]*dky*2*pi;
                
                %Operate forward stress to displacement calculation in the Fourier
                %space
                ii=0;
                %loop over image pixels
                for i=ky(1:end)
                    ii=ii+1;
                    jj=0;
                    for j=kx(1:end)
                        jj=jj+1;
                        
                        kn=sqrt(i^2+j^2);
                        Tnx=Tx_k(ii,jj);
                        Tny=Ty_k(ii,jj);
                        Tn=[Tnx;Tny];
                        
                        r = sqrt(i^2+j^2);
                        
                        % Implement the Green tensor in the Fourier space
                        G=A*2*pi/(kn^3)*[(1-nu)*kn^2+nu*i^2,-nu*i*j; -nu*i*j,(1-nu)*kn^2+nu*j^2];
                        
                        %Take the product of the matrix in the Fourier domain
                        dn = G*Tn;
                        
                        u_fft(ii,jj)=dn(1);
                        v_fft(ii,jj)=dn(2);
                        
                        %Make sure no drift displacement is generated by nulling
                        %the zero-th order displacement in Fourier space
                        u_fft(1,1)=0;
                        v_fft(1,1)=0;
                    end
                end
                
                %Take the inverse Fourier transform of the displacement field
                fu(:,:,l) = real(ifft2(u_fft));
                fv(:,:,l) = real(ifft2(v_fft));
                % Apx=(conversion*1e-6)^2*(x1(1,2)-x1(1,1))*(y1(2,1)-y1(1,1));
                % Fx=Apx*Trx;
                
                
                % %         h = figure('Name','Simulated traction stress in x-direction');
                % %         subplot(1,2,1)
                % %         imagesc(Tx1_0)
                % %         subplot(1,2,2)
                % %         imagesc(Ty1_0)
                % %
                % %         %axis([0 max(xdata(:,:,1) 0 max(xdata(:,:,1) -1e-5 1e-5])
                % %         figure('Name','Simulated traction stress in y-direction')
                % %         plot(t,force)
                % %         axis([-1 1 0 1e-6])
                
                %Plot x-displacement
                fig1 = figure(fig1);
                surf(ax,xdata(:,:,1),xdata(:,:,2),sqrt(fu(:,:,l).^2+fv(:,:,l).^2),sqrt(fu(:,:,l).^2+fv(:,:,l).^2),'EdgeColor','none')
                %hold on
                %quiver3(ax,xdata(:,:,1),xdata(:,:,2),-ones(dim,dim)*zlim,u,v,zeros(dim,dim),2);
                drawnow
                Mov1(l) = getframe(ax);
                ax.NextPlot = 'replaceChildren';
                
                %Plot x-traction
                fig2 = figure(fig2);
                surf(ax2,xdata(:,:,1),xdata(:,:,2),sqrt(Trx.^2+Try.^2),sqrt(Trx.^2+Try.^2),'EdgeColor','none')
                %hold on
                %quiver3(ax,xdata(:,:,1),xdata(:,:,2),-ones(dim,dim)*zlim,u,v,zeros(dim,dim),2);
                drawnow
                Mov2(l) = getframe(ax2);
                ax2.NextPlot = 'replaceChildren';
                
                newdir=['Simulated videos/',F_type,'_signal/F_area_radius_',num2str(fr_x)];
                mkdir(newdir)
                %Write movie of displacement
                disps = sqrt(fu(:,:,l).^2+fv(:,:,l).^2);
                disps = (double(disps))./zlim*20;%(max(disps(:))-min(disps(:)));%Scale values for display
                imwrite(disps,parula,[newdir,'/',F_type(1:3),'_',num2str(m_force),'N_rad_',num2str(fr_x),'m_phi_',num2str(cell_theta),'_disp_tot.tif'],'WriteMode','append')
              
%                 newdir=['Simulated videos/',F_type,'_signal/F_area_radius_',num2str(fr_x)];
%                 mkdir(newdir)
                %Write movie of displacement
                imwrite(Trt./20,parula,[newdir,'/',F_type(1:3),'_',num2str(m_force),'N_rad_',num2str(fr_x),'m_phi_',num2str(cell_theta),'_Trt.tif'],'WriteMode','append')
                
                %Warp beads image and write movie of deforming image
                % Warp base image
                disp2(:,:,1)=fu(:,:,l);
                disp2(:,:,2)=fv(:,:,l);
                beads_d = imwarp(base_im_frame,-disp2./D);
                imwrite(beads_d,[newdir,'/',F_type(1:3),'_',num2str(m_force),'N_rad_',num2str(fr_x),'m_phi_',num2str(cell_theta),'.tif'],'WriteMode','append')

%                 %Warp beads image and write movie of deforming image
%                 %warp previous frame by another step from previous step
%                 disp2_inc(:,:,1)=fu(:,:,l) - disp2_prev(:,:,1);
%                 disp2_inc(:,:,2)=fv(:,:,l) - disp2_prev(:,:,2);
%                 beads_d2 = imwarp(prev_frame,-disp2_inc./D);
%                 imwrite(beads_d,[newdir,'/',F_type(1:3),'_',num2str(m_force),'N_rad_',num2str(fr_x),'m_phi_',num2str(cell_theta),'_sbs.tif'],'WriteMode','append')
%                 
%                 disp2_prev(:,:,1) = fu(:,:,l);
%                 disp2_prev(:,:,2) = fv(:,:,l);
%                 prev_frame =beads_d2;
%                 
%                 %reset image and displacement to initial between each peaks
%                 if any(l == end_peaks)
%                     prev_frame = base_im_frame;
%                     disp2_prev = zeros(dimx,dimy,2);
%                 end
                                
                %%Calculate Strain Energy
                Trx_vec = Trx(:);
                Try_vec = Try(:);
                fu_f = fu(:,:,l); 
                fv_f = fv(:,:,l); 
                u_vec = fu_f(:);
                v_vec = fv_f(:);
                dx = D;
                dy = D;
                U(l) = .5*sum((Trx_vec.*u_vec+Try_vec.*v_vec)*dx*dy);
                %%Calculate mean displacement
                av_disp(l) = mean(sqrt(fu_f(:).^2+fv_f(:).^2))
                
                %create x and y position array for saving
                x(:,:,l) = int16(xdata(1:end,:,1)/(mu_pix*1e-6));
                y(:,:,l) = int16(xdata(1:end,:,2)/(mu_pix*1e-6));
                fu(:,:,l) =fu(:,:,l)/(mu_pix*1e-6);
                fv(:,:,l) =fv(:,:,l)/(mu_pix*1e-6);
                Trx_save(:,:,l) = Trx;
                Try_save(:,:,l) = Try;
            end
            
            Trx = Trx_save;
            Try = Try_save;

            %         figure('Name','Displacement in y-direction')
            %         imagesc(v)
            
            %         %%Save data for Matlab GUI calculation
            newdir=['Simulated videos/',F_type,'_signal/F_area_radius_',num2str(fr_x),'/',F_type(1:3),'_',num2str(m_force),'N_rad_',num2str(fr_x),'m_phi_',num2str(cell_theta)];
            mkdir(newdir)
            mkdir([newdir,'/Datasets/Displacements'])
            mkdir([newdir,'/Datasets/Traction Forces'])
            mkdir([newdir,'/Datasets/Forces'])
            mkdir([newdir,'/Datasets/Strain Energy'])
            mkdir([newdir,'/Datasets/Average Displacement'])
            save([newdir,'/Datasets/Displacements','/tfm_piv_x.mat'],'x','-v7.3')
            save([newdir,'/Datasets/Displacements','/tfm_piv_y.mat'],'y','-v7.3')
            save([newdir,'/Datasets/Displacements','/tfm_piv_u.mat'],'fu','-v7.3')
            save([newdir,'/Datasets/Displacements','/tfm_piv_v.mat'],'fv','-v7.3')
            save([newdir,'/Datasets/Traction Forces','/sim_Trx.mat'],'Trx','-v7.3')
            save([newdir,'/Datasets/Traction Forces','/sim_Try.mat'],'Try','-v7.3')
            
            save([newdir,'/Datasets/Forces','/sim_mFt.mat'],'m_force','-v7.3')
            save([newdir,'/Datasets/Forces','/sim_F_tot.mat'],'F_tot','-v7.3')
            save([newdir,'/Datasets/Forces','/sim_Fx_tot.mat'],'Fx_tot','-v7.3')
            save([newdir,'/Datasets/Forces','/sim_Fy_tot.mat'],'Fy_tot','-v7.3')
            save([newdir,'/Datasets/Strain Energy','/sim_U.mat'],'U','-v7.3')
            save([newdir,'/Datasets/Average Displacement','/sim_av_disp.mat'],'av_disp','-v7.3')
            mkdir([newdir,'/Mask'])
            save([newdir,'/Mask/',F_type(1:3),'_',num2str(m_force),'N_rad_',num2str(fr_x),'m_phi_',num2str(cell_theta),'.mat'],'mask','-v7.3')
            
        end
    end
    
end


% %%====================
% %% Publication figure scripts
% %%====================
% 
% figure(5)
% %imagesc('XData',[0 xdata(end,end,1)],'Ydata',[0 xdata(end,end,2)],'CData',sqrt(u.^2+v.^2))
% contour(xdata(:,:,1), xdata(:,:,2),sqrt(u.^2+v.^2))
% hold on
% quiver(xdata(1:6:end,1:6:end,1),xdata(1:6:end,1:6:end,2),u(1:6:end,1:6:end),v(1:6:end,1:6:end),'AutoScaleFactor',0.5,'Color','r');
% view(-25,-30)
% axis off
% 
disp('Simulation finished')