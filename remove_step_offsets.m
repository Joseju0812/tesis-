% function x_out = remove_step_offsets_v2(x, jumps_idx)
% 
%  N = numel(x);
% jumps_idx = unique(round(jumps_idx));
% jumps_idx = jumps_idx(jumps_idx>=2 & jumps_idx<=N);
% 
% 
% x_out = x;
% M = numel(jumps_idx);
% 
% for k = 1:M
%     j = jumps_idx(k);
%     mseanial_presalto = mean(x_out(1: j-1),'omitnan');
%     svmseanial_presalto = std(x_out(1: j-1),0,'omitnan');
%     coef_var = svmseanial_presalto / mseanial_presalto ;
% 
%      % w = round(interp1([0,100],[1,(0.1*j)], coef_var));
% 
%     pre_start = max(1, j - w); 
%     pre_end = j-1;
%     post_start = j+1 ;
%     post_end = min(N, j + w );
% 
%     pre_seg = x_out(pre_start:pre_end);
%     post_seg = x_out(post_start:post_end);
% 
%     m_pre  = mean(pre_seg, 'omitnan');
%     m_post = mean(post_seg, 'omitnan');
% 
%     % faltan valores para mean, con omitnan
%     d = m_post - m_pre;
%     % d =  m_pre - m_post;
%     % if abs(d) <= eps, continue; end
%     % if k < M, apply_until = jumps_idx(k+1)-1; 
%     %     else apply_until = N; 
%     % end
%     % if apply_until >= j
%     %     x_out(j:apply_until) = x_out(j:apply_until) - d;
%     % end
% 
% 
%     x_out(j:end) = x_out(j:end) - d;
% 
% end
% 
% end


% 
% 
% function x_out = remove_step_offsets_v2(x, jumps_idx)
% 
% N = numel(x);
% jumps_idx = unique(round(jumps_idx));
% jumps_idx = jumps_idx(jumps_idx>=2 & jumps_idx<=N);
% x_out = x;
% M = numel(jumps_idx);
% 
% for k = 1:M
%     j = jumps_idx(k);
% 
% 
%     pre_all = x(1:j-1);
%     if all(isnan(pre_all)), continue; end
% 
%     % MAD escalada ~ estimador robusto de std
%     med_pre_full = median(pre_all,'omitnan');
%     mad_pre = median(abs(pre_all - med_pre_full),'omitnan');
%     s_pre = 1.4826 * mad_pre;  
% 
%     % coeficiente de variación robusto: (1.4826*MAD) / median
%     if abs(med_pre_full) < eps
%         coef_var = 0;
%     else
%         coef_var = s_pre / med_pre_full;
%     end
%     coef_var = min(max(coef_var,0),100);  % clamp a [0,100] para interp1 seguro
% 
%     % interp1: mapear 0->1, 100->0.1*j
%     w = round(interp1([0,100],[1,(0.1*j)], coef_var, 'linear', 'extrap'));
%     w = max(1, min(w, max(1, round(0.1*j))));  % asegurar rango entero razonable
% 
%     % ventanas pre/post (usando la señal original para promedios)
%     pre_start = max(1, j - w); pre_end = j - 1;
%     post_start = j + 1;        post_end = min(N, j + w);
%     pre_seg = x(pre_start:pre_end); post_seg = x(post_start:post_end);
% 
%     % if isempty(pre_seg) || isempty(post_seg) || all(isnan(pre_seg)) || all(isnan(post_seg))
%     %     continue;
%     % end
% 
%     m_pre  = median(pre_seg,'omitnan');   % medianas locales
%     m_post = median(post_seg,'omitnan');
%     d = m_post - m_pre;
% 
%     % if k < M
%     %     apply_until = jumps_idx(k+1) - 1;
%     % else
%     %     apply_until = N;
%     % end
%     % if apply_until >= j
%         % x_out(j:apply_until) = x_out(j:apply_until) - d;
%          x_out(j:end) = x_out(j:end) - d;
%     % end
% end
% end


function x_out = remove_step_offsets(x, jumps_idx)
N = numel(x);
jumps_idx = unique(round(jumps_idx));
jumps_idx = jumps_idx(jumps_idx>=2 & jumps_idx<=N);
x_out = x;
M = numel(jumps_idx);

for k = 1:M
    j = jumps_idx(k);

    
    pre_all = x(1:j-1);
    if all(isnan(pre_all)), continue; end

    s_pre = std(pre_all,0,'omitnan');
    m_pre_full = median(pre_all,'omitnan');   

    
    if abs(m_pre_full) < eps
        coef_var = 0;
    else
        coef_var = s_pre / m_pre_full;
    end
    coef_var = min(max(coef_var,0),100); 

    
    w = round(interp1([0,100],[1,(0.1*j)], coef_var, 'linear', 'extrap'));
    w = max(1, min(w, max(1, round(0.1*j))));  

  
    pre_start = max(1, j - w); pre_end = j - 1;
    post_start = j + 1;        post_end = min(N, j + w);
    pre_seg = x(pre_start:pre_end); post_seg = x(post_start:post_end);

    if isempty(pre_seg) || isempty(post_seg) || all(isnan(pre_seg)) || all(isnan(post_seg))
        continue;
    end

    m_pre  = median(pre_seg,'omitnan');   % mediana local
    m_post = median(post_seg,'omitnan');  % mediana local
    d = m_post - m_pre;
 x_out(j:end) = x_out(j:end) - d;
end
end

