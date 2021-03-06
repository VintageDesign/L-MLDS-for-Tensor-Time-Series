addpath('functions');
% data
load data/tesla
N = numel(X);
Ntrain =1100;

% Setting the type of noise-----------------------------------------------

Type = struct('Q0','Full','Q','Full','R','Full');
I = size(X{1});
J = [5 2];
M = numel(I);
Ntest = N - Ntrain;

%reexpressed tensors as vectors to fit LDS--------------------------------

X_vectorized = vec2ten(ten2vec(X), prod(I));

% LDS---------------------------------------------------------------------

disp('Fitting LDS...')
sizes = zeros(prod(I), 1);
for i = 1:prod(I)
  sizes(i) = number_of_parameters(prod(I), i,Type);
end
clear i
J_lds = find(sizes >= number_of_parameters(I, J,Type),1);
clear i sizes
model_lds = learn_mlds(subcell(X_vectorized, 1:Ntrain),Type, 'J',J_lds );
result_lds = Generator(X, Ntrain, model_lds);
err_lds=err(X, Ntrain, result_lds);
Result_lds=ten_form(X,result_lds);


%preprocess the dataset to fit L-MLDS well--------------------------------
for i=1:N
    XX{i}=X{i}'; 
end
%compute J----------------------------------------------------------------
for i = 1:I(2)
  sizes(i) = L_number_of_parameters(XX,I(2), i,Type);
end
J_lmlds = find(sizes >= number_of_parameters(I, J,Type),1);
if isempty(J_lmlds)
    J_lmlds=I(2);
end
%dct-MLDS-----------------------------------------------------------------
result_dct= dct_mlds(XX,J_lmlds,Ntrain,Type);
err_dct = Err_dct(XX,result_dct,Ntrain);

%dwt-MLDS-----------------------------------------------------------------

J_dwt=J_lmlds;
result_dwt =dwt_mlds( XX,J_dwt,Ntrain,Type);
err_dwt=Err_dct( XX,result_dwt,Ntrain);


%dft-MLDS-----------------------------------------------------------------

result_dft = dft_mlds(XX,J_lmlds,Ntrain,Type);
err_dft=Err_dct(XX,result_dft,Ntrain);

% MLDS--------------------------------------------------------------------

disp('Fitting MLDS with matching number of parameters...')
J_mlds = prod(J);
model_mlds = learn_mlds(subcell(X, 1:Ntrain),Type, 'J', J);
result_mlds = Generator(X, Ntrain, model_mlds);
err_mlds=err(X, Ntrain, result_mlds);
Result_mlds=ten_form(X,result_mlds);



% plot results------------------------------------------------------------
disp('Plotting results...')

%% Error
figure(1)
subplot(1,1,1);
hold on;
T = [1:Ntest]+Ntrain; 
plot(T, err_lds, 'Color', 'blue');
plot(T, err_mlds, 'Color', 'black');
plot(T, err_dft, 'Color', 'yellow');
plot(T, err_dct, 'Color', 'red');
plot(T, err_dwt, 'Color', 'green');
hold off;
legend('LDS','MLDS','dft-MLDS','dct-MLDS','dwt-MLDS');
xlim([1 Ntest] + Ntrain);
xlabel('Time slice');
ylabel('Error');

%% Prediction value
figure(2)
real_lds=zeros(1,Ntest);
real_mlds=zeros(1,Ntest);
real_dct=zeros(1,Ntest);
real_dwt=zeros(1,Ntest);
real_dft=zeros(1,Ntest);
real=zeros(1,Ntest);

for i=1:Ntest
    real_lds(i)  = Result_lds{i}(6,1);
    real_mlds(i) = Result_mlds{i}(6,1);
    real_dct(i)  = result_dct{i}(1,6);
    real_dwt(i)  = result_dwt{i}(1,6);
    real_dft(i)  = result_dft{i}(1,6);
    real(i)      = X{i+Ntrain}(6,1); 
end
subplot(1,1,1);
hold on;
T = [1:Ntest]+Ntrain; 
plot(T, real_lds, 'Color', 'blue');
plot(T, real_mlds, 'Color', 'black');
plot(T, abs(real_dft), 'Color', 'yellow');
plot(T, real_dct, 'Color', 'red');
plot(T, real_dwt, 'Color', 'green');
plot(T, real, 'Color', 'cyan');
hold off;
legend('LDS','MLDS','dft-MLDS','dct-MLDS','dwt-MLDS','real');
xlim([1 Ntest] + Ntrain);
xlabel('Time slice');
ylabel('Price');
