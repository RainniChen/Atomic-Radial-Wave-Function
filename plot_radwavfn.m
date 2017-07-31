function plot_radwavfn(n,l,h)
    rs = 2*n*(n+15);
    xs = ceil(log(rs));
    
    % get wave functions
    X = radial_wav_fn(n,l,h,0);
    R = radial_wav_fn(n,l,h,1);
    
    % plot
    figure('rend','painters','pos',[0 150 850 550])
    x_space = linspace(0,xs,xs/h+1);
    r_space = exp(x_space);
    subplot(2,1,1)
    plot(r_space, R, '.-')
    subplot(2,1,2)
    plot(x_space,X, '.-')
  
end