% Following program is used to calculate Complex Hurst, Version 2.0
% Contributed to github server, Jun Steed Huang on 12/24/2025
% Calculate the dementia vs normal dataset's Hurst value
% Clear the work space
  clear
  clf
  clc
% Input data ms & pT values
  abc=xlsread('control.xls');
  abd=xlsread('dementia.xls');
% Pick the interested column for this complex analysis
  ac=abc(:,1);     % this is column for time tick, we don't use it here
  bc=abc(:,2);     % this is MEG value in pT, can change to other
  ad=abd(:,1);     % sick data
  bd=abd(:,2);     % sick data
% Calculate the original length of data
  num_data=size(bc);    % this is the number of samples counted
% F=1.5 is the Fractional moment between Mean and Variance 
  F=1.5 
% Plot together
  figure(1)
% Calculate the Index of Dispersion for Count controlled
  k=1; %control round
  d=bc; 
% The level of grouping, minimum data points is 700
    L=[1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192];
% 1 is for 1 ms, 2 is for 2 ms, 128 ms = 7.83Hz
    m=length(L);   % number of dots on IDC curve
% Calculate each point of the graph one by one
  for r=1:m
      clear sum V E  num_batch  Y
      num_batch=fix(num_data/ L(r));
        % Add varaible within the group
        for i=1:num_batch
            Y(i)=0;
             for j=1:L(r)
                 Y(i)=Y(i)+ d((i-1)*L(r)+j);
             end
        end
        % Calculate the standard variance and normalize it
        n=length(Y);
        E=sum(Y)/n;
        sum=0;
        for j=1:n
            % This is where Complex Number started
             sum = sum + (Y(j)-E)^F;
        end
        % Normalize the Complex Number 
        V = sum/(n-1);
        % IDC of Complex Number
        I(r) = V/E;
  end 
  % Plot value part, ignore the sign
%  figure(1)
  loglog(L,abs(real(I)), '-ws','LineWidth',2,'MarkerSize',14,'MarkerFaceColor',[k/5,0.5,0.271]);
  % to add on next imag curve
  hold on;
  loglog(L,abs(imag(I)), '-ws','LineWidth',2,'MarkerSize',7,'MarkerFaceColor',[(k+1)/5,0.5,0.271]);
  % Use the last and first point to calculate the Complex Hurst
  slope=(log(I(m))-log(I(1)))/(log(L(m))-log(L(1)));
  Hurst_Control=(1+slope)/2    % this is the number we looking for!
  
% Calculate the Index of Dispersion for Count dementia
  k=3; %AD round
  d=bd; 
% The level of grouping, minimum data points is 700
    L=[1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192];
% 1 is for 1 ms, 2 is for 2 ms, 128 ms = 7.83Hz
    m=length(L);   % number of dots on IDC curve
% Calculate each point of the graph one by one
  for r=1:m
      clear sum V E  num_batch  Y
      num_batch=fix(num_data/ L(r));
        % Add varaible within the group
        for i=1:num_batch
            Y(i)=0;
             for j=1:L(r)
                 Y(i)=Y(i)+ d((i-1)*L(r)+j);
             end
        end
        % Calculate the standard variance and normalize it
        n=length(Y);
        E=sum(Y)/n;
        sum=0;
        for j=1:n
            % This is where Complex Number started
             sum = sum + (Y(j)-E)^F;
        end
        % Normalize the Complex Number 
        V = sum/(n-1);
        % IDC of Complex Number
        I(r) = V/E;
  end 
% Plot value part, ignore the sign
%  figure(1)
  loglog(L,abs(real(I)), '-ws','LineWidth',2,'MarkerSize',14,'MarkerFaceColor',[k/5,0.5,0.271]); 
  loglog(L,abs(imag(I)), '-ws','LineWidth',2,'MarkerSize',7,'MarkerFaceColor',[(k+1)/5,0.5,0.271]);
  % Use the last and first point to calculate the Complex Hurst
  slope=(log(I(m))-log(I(1)))/(log(L(m))-log(L(1)));
  Hurst_AD=(1+slope)/2    % this is the number we looking for!
  
% Sync the title with the F value, legend
legend('Normal Real','Normal Imag','Sick Real','Sick Imag');
  title('Index of Dispersion Plot for Semi Complex Hurst');
  xlabel('From 1 ms to 20 sec');
  ylabel('Index of Dispersion for Root of Pico Tesla');
  grid on;
hold off;
% end  