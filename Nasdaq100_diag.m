addpath('functions');
% data
imgs = readMNIST('../mnist', '../mnist_labels', 0);

N = length(imgs);
for i = 1:N
    X{i} = mat2cell(imgs(:, :,i) , 28, 28);
    X{i} = X{i}{1};
end



%Total 5421
Ntrain =5410;



% Setting the type of noise-----------------------------------------------

Type = struct('Q0','Full','Q','Full','R','Full');
I = size(X{1});
J = [28 28];
M = numel(I);
Ntest = N - Ntrain;

% %reexpressed tensors as vectors to fit LDS--------------------------------
% 
% X_vectorized = vec2ten(ten2vec(X), prod(I));
% 
% % LDS---------------------------------------------------------------------
% 
% disp('Fitting LDS...')
% sizes = zeros(prod(I), 1);
% for i = 1:prod(I)
%   sizes(i) = number_of_parameters(prod(I), i,Type);
% end
% clear i
% J_lds = find(sizes >= number_of_parameters(I, J,Type),1);
% clear i sizes
% model_lds = learn_mlds(subcell(X_vectorized, 1:Ntrain),Type, 'J',J_lds );
% result_lds = Generator(X, Ntrain, model_lds);
% err_lds=err(X, Ntrain, result_lds);
% Result_lds=ten_form(X,result_lds);
% 

% %preprocess the dataset to fit L-MLDS well--------------------------------
% for i=1:N
%      XX{i}=X{ i}'; 
% end
XX = X;
%compute J----------------------------------------------------------------
for i = 1:I(2)
  sizes(i) = L_number_of_parameters(XX,I(2), i,Type);
end
J_lmlds = find(sizes >= number_of_parameters(I, J,Type),1);
if isempty(J_lmlds)
    J_lmlds=I(2);
end




%dft-MLDS-----------------------------------------------------------------
result_dft= dct_mlds(XX,J_lmlds,Ntrain,Type);

for n = 1:10
    results(:, :, n) = result_dft{n};
end

T = [1:10]+Ntrain; 
figure(2)
montage(results(:,:, 1:10));
figure(3)
montage(imgs(:,:,T));


imwrite(results(:,:,1), 'forecasted_1.png');
imwrite(results(:,:,2), 'forecasted_2.png');
imwrite(results(:,:,3), 'forecasted_3.png');
imwrite(results(:,:,4), 'forecasted_4.png');
imwrite(results(:,:,5), 'forecasted_5.png');
imwrite(results(:,:,6), 'forecasted_6.png');
imwrite(results(:,:,7), 'forecasted_7.png');
imwrite(results(:,:,8), 'forecasted_8.png');
imwrite(results(:,:,9), 'forecasted_9.png');
imwrite(results(:,:,10), 'forecasted_0.png');

err_dft = Err_dct(XX,result_dft,Ntrain);


% MLDS--------------------------------------------------------------------


% plot results------------------------------------------------------------
disp('Plotting results...')

%% Error
% figure(4)
% subplot(1,1,1);
% hold on;
% 
% 
% plot(T, err_dft, 'Color', 'yellow');
% hold off;
% legend('dft-MLDS');
% xlim([1 Ntest] + Ntrain);
% xlabel('Time slice');
% ylabel('Error');
% 
%% Prediction value
% figure(5)
% 
% real_dft=zeros(1,Ntest);
% real=zeros(1,Ntest);
% 
% for i=1:Ntest
%     
%     real_dft(i)  = result_dft{i}(1,6);
%     real(i)      = X{i+Ntrain}(6,1); 
% end
% subplot(1,1,1);
% hold on;
% T = [1:Ntest]+Ntrain; 
% 
% plot(T, abs(real_dft), 'Color', 'yellow');
% 
% plot(T, real, 'Color', 'cyan');
% hold off;
% legend('dft-MLDS','real');
% xlim([1 Ntest] + Ntrain);
% xlabel('Time slice');
% ylabel('Price');
