function [er] = err()

X_r = fscanfMat('wave_wv_10000/frame.x');
q_r = fscanfMat('wave_wv_10000/frame.q2');


X_1 = fscanfMat('wave_wv_100/frame.x');
q_1 = fscanfMat('wave_wv_100/frame.q2');
qq = interpln([X_r'; q_r(:,1)'], X_1);
nx = max(size(qq));
er_1 = norm(qq-q_1(:,1)')/nx;

X_1 = fscanfMat('wave_wv_200/frame.x');
q_1 = fscanfMat('wave_wv_200/frame.q2');
qq = interpln([X_r'; q_r(:,1)'], X_1);
nx = max(size(qq));
er_2 = norm(qq-q_1(:,1)')/nx;

X_1 = fscanfMat('wave_wv_400/frame.x');
q_1 = fscanfMat('wave_wv_400/frame.q2');
qq = interpln([X_r'; q_r(:,1)'], X_1);
nx = max(size(qq));
er_3 = norm(qq-q_1(:,1)')/nx;

X_1 = fscanfMat('wave_wv_800/frame.x');
q_1 = fscanfMat('wave_wv_800/frame.q2');
qq = interpln([X_r'; q_r(:,1)'], X_1);
nx = max(size(qq));
er_4 = norm(qq-q_1(:,1)')/nx;

X_1 = fscanfMat('wave_wv_1600/frame.x');
q_1 = fscanfMat('wave_wv_1600/frame.q2');
qq = interpln([X_r'; q_r(:,1)'], X_1);
nx = max(size(qq));
er_5 = norm(qq-q_1(:,1)')/nx;

er = [er_1, er_2, er_3, er_4, er_5];

endfunction
