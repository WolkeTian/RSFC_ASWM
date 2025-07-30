%% GCA analysis; MVGC toolbox
% 参数设定
nvar = 2;
regmode ='OLS'; % 回归模型
icregmode = regmode; % 用于AIC/BIC
morder = 'AIC'; %使用AIC自动确定模型阶数
momax = 10; % 最大可能阶数
alpha = 0.05; % 显著性水平
% 开始
for i = 1:size(bold, 2)
    % testdata = [subi_LH_roi1, subi_roi2]';

    % z-score标准化
    testdata = zscore(testdata,0,2);
    
    % 模型阶数选择
    [AIC, BIC, moAIC, moBIC] = tsdata_to_infocrit(testdata, momax, icregmode);
    morder = 'AIC'; % 使用AIC
    if strcmp(morder, 'AIC')
        morder_used = moAIC;
    else
        morder_used = moBIC;
    end
    
    % 模型拟合
    [A, SIG, E] = tsdata_to_var(testdata, morder_used, regmode);
    
    % 检查稳定性
    %assert(isstable(A), 'VAR模型不稳定');
    
    % 计算GC
    F = var_to_autocov(A, SIG, [], [], true); 
    %assert(~info.error, '自协方差计算失败');
    GC = autocov_to_pwcgc(F);  % 结果为 2x2 的矩阵
    
    GC_1to2(i) = GC(1,2);
    GC_2to1(i) = GC(2,1);
end

