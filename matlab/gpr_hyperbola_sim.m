% #####################################################################
% ########## Generate simulated GPR data from point targets ###########
% #####################################################################
% Point target locations are determined by p_scats matrix:
% format is [x0,y0,z0;x1,y1,z1;...]
% #####################################################################

xvec = [-5:.1:5];
rxstep = xvec(2)-xvec(1);

% target scatterer locations [x0,y0,z0;x1,y1,z1;...]
p_scats = [-2,3,0; 0,3,0; 2,3.5,0];

% permittivity of the ground
er = 1;

% antenna height above the ground (should really be y not z)
z_gnd = 2.5;

d1 = zeros(size(xvec));
d2=d1;
d3=d1;

dists = zeros(size(p_scats,1),numel(xvec));
for i=1:size(dists,2)
    for j = 1:size(dists,1)
        dists(j,i) = ray_range_2medium([xvec(i),0,0],p_scats(j,:),z_gnd,1,er);
    end
end

figure(101); hold on;
for i=1:size(dists,1)
    plot(xvec,dists(i,:)); 
end
plot(xvec,z_gnd*ones(size(xvec)));
scatter(p_scats(:,1),p_scats(:,2),'x');
title('Simulated target signature ranges');

t = linspace(0,1.5*max(max(dists))*2/3e8,500);

f = 6e8;
fs = 1/(t(2)-t(1));
w_rick = -(2*((pi*f)^2)*((t-(sqrt(2)/f)).^2)-1).*exp(-((pi*f)^2)*((t-(sqrt(2)/f)).^2));
w_rick = 64*w_rick;

Z = zeros(numel(t),numel(xvec));

for i=1:numel(xvec)
    Z(:,i) = (1/z_gnd^2)*shift(w_rick,z_gnd*2*fs/3e8);
    for j=1:size(dists,1)
        Z(:,i) = Z(:,i) + (1/dists(j,i)^2)*shift(w_rick,dists(j,i)*2*fs/3e8).';
    end
end

range = t*3e8/2;
[~,imax] = max(w_rick);
range = range - range(imax);

[X,Y] = meshgrid(xvec,range);

figure;

subplot(2,1,1);  surf(X,Y,abs(Z),'EdgeColor','none','FaceColor','interp'); view(90,90); axis tight;
xlabel('xpos'); ylabel('range'); title('mag Original');

subplot(2,1,2); surf(X,Y,Z,'EdgeColor','none','FaceColor','interp'); view(90,90); axis tight;
xlabel('xpos'); ylabel('range'); title('Original');

set(gcf, 'Units','centimeters', 'Position', [0 0 15 20])
movegui(gcf, 'west');


zpad1 = 0;
zpad2 = 0;
% field = Z;

field = [zeros(size(Z,1),zpad1/2),Z,zeros(size(Z,1),zpad1/2)];
field = [field;zeros(zpad2,size(field,2))];

rxlocs = xvec;

[~,imin]=min(abs(range-(z_gnd)));

range_gnd = range(imin:end);
field_gnd = field(imin:end,:);

% #####################################################################
% ########### Focus the image using k-omega algorithm #################
% #####################################################################
% K-omega algorithm works on spectrum of image to do focusing. Is
% computationally effecient, but gives poor results.
% #####################################################################

field_fft = fft2(ifftshift(field_gnd,2));

field_fft = field_fft./(max(max(abs(field_fft))));
freq = (linspace(0,fs,size(field_fft,1)));
k = 2*pi./(3e8./(sqrt(er)*freq));

kx = linspace(0,2*pi/(rxstep),size(field_fft,2));
Kz = zeros(size(field_fft));

[~,indr] = min(abs(4*k.^2-kx(end).^2));
if ((4*k(indr)^2)<(kx(end)^2))
    indr = indr+1;
end
    
kzlin = linspace(sqrt(4*k(indr)^2-kx(end)^2),sqrt(4*k(end-indr+1)^2-kx(end)^2),numel(k));

field_fft_focus = zeros(size(field_fft));
% figure; hold on;
for i=1:size(field_fft,2)
    [~,indr] = min(abs(4*k.^2-kx(i).^2));
    if ((4*k(indr)^2)<(kx(i)^2))
        indr = indr+1;
    end
    kztemp = sqrt(4*k(indr:end-indr+1).^2-kx(i).^2); 

    Kz(1:numel(kztemp),i) = kztemp;
%     plot(kztemp);
    field_fft_focus(:,i) = interp1(kztemp,field_fft(indr:end-indr+1,i),kzlin,'linear');

end
% figure; imagesc(real(Kz)); colorbar;
% figure; imagesc(isfinite(field_fft_focus)); colorbar;


figure; surf(abs(field_fft_focus),'EdgeColor','none','FaceColor','interp'); title('abs spectrum field fft focus');
figure; 
subplot(4,1,1:2); surf(real(field_fft_focus),'EdgeColor','none','FaceColor','interp'); title('real spectrum field fft focus');
subplot(4,1,3:4); surf(imag(field_fft_focus),'EdgeColor','none','FaceColor','interp'); title('imag spectrum field fft focus');

% field_foc_sq = fftshift(ifft2(field_fft),2);
field_fft_focus(~isfinite(field_fft_focus))=0;
field_foc = fftshift(ifft2(field_fft_focus),2);

field_foc = field_foc(1:end-zpad2,zpad1/2+1:end-zpad1/2);

figure; 

subplot(4,2,[1,3]); surf((abs(field_foc)),'EdgeColor','none','FaceColor','interp'); view(90,90); axis tight;
title('mag FFT focused image: sim');

subplot(4,2,[2,4]);  surf(20*log10(abs(field_foc)),'EdgeColor','none','FaceColor','interp'); view(90,90); axis tight;
title('db mag FFT focused image: sim');

subplot(4,2,[5,7]);  surf((real(field_foc)),'EdgeColor','none','FaceColor','interp'); view(90,90); axis tight;
title('real FFT focused image: sim');

subplot(4,2,[6,8]);  surf((imag(field_foc)),'EdgeColor','none','FaceColor','interp'); view(90,90); axis tight;
title('imag FFT focused image: sim');

set(gcf, 'Units','centimeters', 'Position', [0 0 40 20])
movegui(gcf, 'center');
    
% figure;
% imagesc(abs(field_foc)); title('mag fft focused image: sim');
% figure;
% imagesc(real(field_foc));title('real fft focused image: sim');
% figure;
% imagesc(imag(field_foc));title('imag fft focused image: sim');
% figure;
% imagesc(20*log10(abs(field_foc)));title('db mag fft focused image: sim');

% #####################################################################
% ########### Focus the image time domain backprojection ##############
% #####################################################################
% Time domain backprojection gives better results but is more
% computationally complex...
% #####################################################################

% ##### Run this to regenerated dist_mat: ######
% clear dist_mat_saved

% precompute lookup table for all possible p_ant-p_scat relative positions
% (only need first and last antenna location)
    
if (exist('dist_mat_saved'))
    dist_mat = dist_mat_saved;
else
    options = optimset('Display','off');
    x_offs = zeros(1,numel(rxlocs)*2);
    dist_mat = zeros(numel(rxlocs)*2,numel(range_gnd));
    off_ind = 0;
    for ii=[1,numel(rxlocs)]
        for j=1:size(field_gnd,2)

    % ######## takes er of ground into account: much slower ######## 
            dist = sqrt(range_gnd.^2-(rxlocs(ii)-rxlocs(j))^2);
            rind = find(imag(dist)==0,1);
    %         dist = zeros(1,numel(range_gnd));
            for jj=rind:numel(range_gnd)
                fun = @(x) range_gnd(jj) - ray_range_2medium([rxlocs(ii),0,0],[rxlocs(j),x,0],z_gnd,1,er); 
                dist(jj) = fzero(fun,dist(jj),options);
    %             w = warning('query','last')
            end
            if (mod(j,20)==0)
                figure(201); hold on; plot(real(dist)); title(sprintf('antloc %g, imloc %g',ii,j));
            end
            off_ind = off_ind+1;
            x_offs(off_ind) = rxlocs(j)-rxlocs(ii);
            dist_mat(off_ind,:) = dist;
        end
    end
    size(dist_mat)
    % remove repeated value at end
    if (x_offs(1)==x_offs(end))
        x_offs = x_offs(1:end-1);
        dist_mat = dist_mat(1:end-1,:);
    end

end

dist_mat_saved = dist_mat;

field_autofoc = zeros(size(field_gnd));
for ii=1:numel(rxlocs)
    ascan = field_gnd(:,ii);
    for j=1:size(field_autofoc,2)
% ######## Use precomputed matrix lookup ######## 
        
        x_off = rxlocs(j)-rxlocs(ii);
        dist_row = find(x_offs==x_off);
        dist = dist_mat(dist_row,:);
% ######## takes er of ground into account: much slower ######## 
%         dist = sqrt(range_gnd.^2-(rxlocs(ii)-rxlocs(j))^2);
%         rind = find(imag(dist)==0,1);
%         dist = zeros(1,numel(range_gnd));
%         for jj=rind:numel(range_gnd)
%             fun = @(x) range_gnd(jj) - ray_range_2medium([rxlocs(ii),0,0],[rxlocs(j),x,0],z_gnd,1,er); 
%             dist(jj) = fzero(fun,dist(jj));
%         end        
        % ######## doesn't take er of ground into account ######## 
%         dist = sqrt(range_gnd.^2-(rxlocs(ii)-rxlocs(j))^2);
        
        distprev = dist;
%         rinds = find(imag(dist)==0,1)
%         dist = dist(rinds:end);
%         finds = find(isfinite(dist),1)
%         dist = dist(finds:end);
%         dinds = rinds+finds-1;
        dinds = find((imag(dist)==0)&(isfinite(dist)));
        if (numel(dist)>0)
            field_temp =interp1(dist(dinds).',field_gnd(dinds,ii),range_gnd(:));
            field_autofoc(isfinite(field_temp),j)=field_autofoc(isfinite(field_temp),j)+field_temp(isfinite(field_temp));
%             field_autofoc(:,j)=field_autofoc(:,j)+interp1(dist(:),field_gnd(min(rinds):end,ii),range_gnd(:),'spline');
        end
    end
end

[X,Y] = meshgrid(rxlocs,range_gnd);

figure; 

subplot(4,2,[1,3]); surf(X,Y,(abs(field_autofoc)),'EdgeColor','none','FaceColor','interp'); view(90,90); axis tight;
title('mag Back-projection focused image: sim');

subplot(4,2,[2,4]);  surf(X,Y,20*log10(abs(field_autofoc)),'EdgeColor','none','FaceColor','interp'); view(90,90); axis tight;
title('db mag Back-projection focused image: sim');

subplot(4,2,[5,7]);  surf(X,Y,(real(field_autofoc)),'EdgeColor','none','FaceColor','interp'); view(90,90); axis tight;
title('real Back-projection focused image: sim');

subplot(4,2,[6,8]);  surf(X,Y,(imag(field_autofoc)),'EdgeColor','none','FaceColor','interp'); view(90,90); axis tight;
title('imag Back-projection focused image: sim');

set(gcf, 'Units','centimeters', 'Position', [0 0 40 20])
movegui(gcf, 'center');

