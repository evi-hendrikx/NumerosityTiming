function meanRho = NumTim_make_mean_correlation(corrs_1, corrs_2)
%% calculate mean correlation between two correlations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% corrs_1: first correlation coefficients; 
% corrs_2: second correlation coefficients 
%
% Output
% meanRho: mean correlation for each pair of coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

corrs_1 = reshape(corrs_1,1,length(corrs_1));
corrs_2 = reshape(corrs_2,1,length(corrs_2));

M = [corrs_1 ;corrs_2];

z = atanh(M);
meanRho = real(tanh(mean(z)));

% if one of the two is NaN, take the other one
meanRho(isnan(corrs_1)&~isnan(corrs_2))= corrs_2(isnan(corrs_1)&~isnan(corrs_2));
meanRho(isnan(corrs_2)&~isnan(corrs_1))= corrs_1(isnan(corrs_2)&~isnan(corrs_1));

end
