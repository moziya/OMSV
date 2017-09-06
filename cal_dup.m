function cal_dup(ref_file,align_file,outFile,z_num)
if nargin == 0
	ref_file = 'ref_mat2.txt';
	align_file = 'cov_mat2.txt';
	outFile = 'C6661_dup.osv';
	z_num=2;
elseif nargin ~= 4
	fprintf('There are 4 parameters required: 1.input ref file, 2.input alignment file, 3.output file, 4.N\n');
	return;
end

outFile = [outFile(1:end-4) '_' num2str(z_num) outFile(end-3:end)];
for i = 1:24
[real_dup{i},real_del{i}] = red_smot(700,i,ref_file,align_file,z_num);
end
%save('Dup_Del_call.mat','real_dup','real_del');

fid = fopen(outFile,'w');
fprintf(fid,'#Called by Matlab code\n');
fprintf(fid,'#chr	start	stop	variant_call_type	variant_call_id	Reserved	sample_id	SV size	zygosity	likelihood	score	Coverage	Against\n');
cnt = 1;
for i = 1:length(real_dup)
for j = 1:size(real_dup{i},1)
siz = real_dup{i}(j,2)-real_dup{i}(j,1);
fprintf(fid,'%d\t%d\t%d\tDuplication\t%d\t0\tC666-1\tsizeChange=%d\tNone\t%d\t0\t0\t0\n',i,real_dup{i}(j,1),real_dup{i}(j,2),cnt,siz,real_dup{i}(j,3));
cnt = cnt + 1;
end
end
fclose(fid);

end


function [tot_dup,tot_del] = red_smot(res_size,chr,ref_file,align_file,z_stat)
%need to recover the segments lost in split map
%refmat = importdata('ref_mat.txt');
%z_stat = 2;
%smt_vec = [0.0228 0.9544 0.0228];
ref_pos = importdata(ref_file);
wd_size = res_size*z_stat; % std_res = resolution/2, then here z-statistic = 2 when wd_size = res_size*2 (i.e. X=wd_size/2=std_res*2)
%A = textread(['hg38/chr' num2str(chr) '.fa'],'%s',-1,'headerlines',1);
%seq = cell2mat(A(1:end-1));
%seq = seq';
%seq = seq(:);
%%generate marks of large repeats
% idx = seq>='a'&seq<='z';
% cnt = 0;m_v = 0;m_p = 0;
rep_mask = [];
% len_thre = 200000;
% for i = 1:length(idx)
%     if idx(i)
%         cnt=cnt+1;
%     else
%         if cnt >= 10000
%             rep_mask = [rep_mask;[floor((i-cnt)/wd_size) ceil(i/wd_size)]];
%         end
%         cnt = 0;
%     end
% end
% if size(rep_mask,1)>0
%     new_mask = rep_mask(1,:);
%     cnt = 1;
%     for i = 2:size(rep_mask,1)
%         if rep_mask(i,1)<=new_mask(cnt,2)+5 || mean(idx((new_mask(cnt,2)+1)*res_size:rep_mask(i,1)*res_size))>0.9
%             new_mask(cnt,2) = max(rep_mask(i,2),new_mask(cnt,2));
%         else
%             new_mask = [new_mask;rep_mask(i,:)];
%             cnt = cnt + 1;
%         end
%     end
%     rep_mask = new_mask(new_mask(:,2)-new_mask(:,1)>=200000/res_size,:);
% end

%%also need to recover the split maps
%seq = upper([seq' A{end}]);
clear new_mask idx

covmat = importdata(align_file);
covmat(:,3:4) = [];
rms = [];
for i = 2:size(covmat,1)
    if covmat(i-1,1) == covmat(i,1)&&covmat(i-1,2) == covmat(i,2)
        rms = [rms i-1];
        covmat(i,3) = min(covmat(i,3),covmat(i-1,3));
        covmat(i,4) = max(covmat(i,4),covmat(i-1,4));
    end
end
covmat(:,1) = [];
covmat(rms,:) = [];

[~,ids] = sort(covmat(:,2));
covmat = covmat(ids,:);
[~,ids] = sort(covmat(:,1));
covmat = covmat(ids,:);
covmat = covmat(covmat(:,1) == chr,:);
%for i = 1:size(refmat,1)
%    vec{refmat(i,1)} = zeros(1,round(refmat(i,2)/wd_size));
%end
%for i =1:size(covmat,1)
%    vec{covmat(i,1)}(max(round(covmat(i,2)/wd_size),1):round(covmat(i,3)/wd_size)) = vec{covmat(i,1)}(max(round(covmat(i,2)/wd_size),1):round(covmat(i,3)/wd_size)) + 1;
%end
%%add masks of large repeats in the reference
len_chr = ref_pos(find(ref_pos(:,1)==chr,1,'last'),2);
vec = zeros(1,round(len_chr/wd_size));
NumS = zeros(1,round(len_chr/wd_size));
NS=20;
clear seq
for i =1:size(covmat,1)
    vec(max(round(covmat(i,2)/wd_size),1):round(covmat(i,3)/wd_size)) = vec(max(round(covmat(i,2)/wd_size),1):round(covmat(i,3)/wd_size)) + 1;
end
% for i = 2:length(vec)-1
%     vec(i) = sum(vec(i-1:i+1).*smt_vec);
% end


ref_pos = ref_pos(ref_pos(:,1)==chr,2);

NumS(end) = eps*rand;
for i = 1:length(NumS)-NS
    NumS(i) = sum(ref_pos>(i-1)*wd_size&ref_pos<(i+NS-1)*wd_size);
end

NS_rat = ones(1,length(NumS));
%vec = vec{chr};

vec(vec < 0.1*mean(vec(vec>0))) = 0;
avail_reg = vec>0;

% avail_reg2 = vec>0;
% for i = 1:size(rep_mask,1)
%     avail_reg2(rep_mask(i,1):rep_mask(i,2)) = false;
% end
median_vec = median(vec(avail_reg));
for i = 1:length(avail_reg)
    if avail_reg(i)
        NS_rat(i) = median_vec/median(vec(avail_reg&NumS==NumS(i)));
        if isinf(NS_rat(i))
            disp('')
        end
    end
end
vec = vec.*NS_rat;

mean_vec = mean(vec(avail_reg));
std_vec = std(vec(avail_reg&vec<mean_vec*3));

pr_upper = 1 - normcdf((vec-mean_vec)/std_vec);
pr_lower = normcdf((vec-mean_vec)/std_vec);

OM_ava_len = 244075;

FPR = 0.05;
L = length(vec);
tot_dup = [];
tot_del = [];
for ts = 30:-1:2
    real_dup = [];
    real_del = [];
    cand_dup = [];
    cand_del = [];
    ts_score = (FPR/L*ts)^(1/ts);
    if ts_score >= 0.5
        continue;
    end
    disp(['Start test merge ' num2str(ts) ' blocks']);
    dup_yes = pr_upper < ts_score;
    del_yes = pr_lower < ts_score;
    dup_yes = dup_yes & avail_reg;
    del_yes = del_yes & avail_reg;
    for j = 1:length(dup_yes)-ts
        if (dup_yes(j)>0&&all(dup_yes(j:j+ts-1)>0))
            cand_dup = [cand_dup; [j j+ts-1]];
        end
    end
    dup_idx = [];
    for k = 1:size(cand_dup,1)
%        [mean(vec(cand_dup(k,1):cand_dup(k,2))/mean_vec) 1 - normcdf(mean(vec(cand_dup(k,1):cand_dup(k,2))),mean_vec,std_vec)]
        if (mean(vec(cand_dup(k,1):cand_dup(k,2))/mean_vec)>1.3 || mean(vec(cand_dup(k,1):cand_dup(k,2))/mean_vec)<0.7)% && 1 - normcdf(mean(vec(cand_dup(k,1):cand_dup(k,2))),mean_vec,std_vec)<FPR
            real_dup = [real_dup; cand_dup(k,:)];
            
        end
    end
    cadi_dup = [];
    cadi_del = [];
    if ~isempty(real_dup)
        [~,dup_idx] = sort(real_dup(:,1));
        real_dup = real_dup(dup_idx,:);
        new_dup = real_dup(1,:);
        ct = 1;
        for l = 2:size(real_dup,1)
            if real_dup(l,1) <= new_dup(ct,2)+1
                new_dup(ct,2) = max(real_dup(l,2),new_dup(ct,2));
            else
                new_dup = [new_dup;real_dup(l,:)];
                ct = ct + 1;
            end
        end
        real_dup = [];
        
        for xy = 1:size(new_dup,1)
            p_val = 1 - normcdf(mean(vec(new_dup(xy,1):new_dup(xy,2))),mean_vec,std_vec);
            if p_val < FPR/max((size(new_dup,1))/ts,1) && (new_dup(xy,2)-new_dup(xy,1))*res_size*2>OM_ava_len/2
                real_dup = [real_dup; [new_dup(xy,:)*wd_size -(log(p_val)/log(2))]];
                avail_reg(new_dup(xy,1):new_dup(xy,2)) = 0;
            elseif 1 - normcdf(mean(vec(new_dup(xy,1):new_dup(xy,2))),mean_vec,std_vec) < FPR
                cadi_dup = [cadi_dup; new_dup(xy,:)];
    %            cct = cct+1;
    %            [mean(vec(new_dup(xy,1):new_dup(xy,2))/mean_vec) 1 - normcdf(mean(vec(new_dup(xy,1):new_dup(xy,2))),mean_vec,std_vec)]
            end
        end
    end
    disp(['There are ' num2str(size(cadi_dup,1)) ' candidate dups']);
%     real_dup = new_dup;
        
    for j = 1:length(del_yes)-ts
        try
        if (del_yes(j)>0&&all(del_yes(j:j+ts-1)>0))
            cand_del = [cand_del; [j j+ts-1]];
        end
        catch
            disp('')
        end
    end
    for k = 1:size(cand_del,1)
 %       [mean(vec(cand_del(k,1):cand_del(k,2))/mean_vec) normcdf(mean(vec(cand_del(k,1):cand_del(k,2))),mean_vec,std_vec)]
        if (mean(vec(cand_del(k,1):cand_del(k,2))/mean_vec)>1.3 || mean(vec(cand_del(k,1):cand_del(k,2))/mean_vec)<0.7)% && normcdf(mean(vec(cand_del(k,1):cand_del(k,2))),mean_vec,std_vec)<FPR
            real_del = [real_del; cand_del(k,:)];
            %avail_reg(cand_del(k,1):cand_del(k,2)) = 0;
        end
    end
    if ~isempty(real_del)
        [~,del_idx] = sort(real_del(:,1));
        real_del = real_del(del_idx,:);
        new_del = real_del(1,:);
        ct = 1;
        for l = 2:size(real_del,1)
            if real_del(l,1) <= new_del(ct,2)+1
                new_del(ct,2) = max(real_del(l,2),new_del(ct,2));
            else
                new_del = [new_del;real_del(l,:)];
                ct = ct + 1;
            end
        end
        real_del = [];
        for xy = 1:size(new_del,1)
            p_val = normcdf(mean(vec(new_del(xy,1):new_del(xy,2))),mean_vec,std_vec);
            if p_val < FPR/max((size(new_del,1))/ts,1) && (new_del(xy,2)-new_del(xy,1))*res_size*2>OM_ava_len/2
                real_del = [real_del; [new_del(xy,:)*wd_size -(log(p_val)/log(2))]];
                avail_reg(new_del(xy,1):new_del(xy,2)) = 0;
            elseif normcdf(mean(vec(new_del(xy,1):new_del(xy,2))),mean_vec,std_vec) < FPR
                cadi_del = [cadi_del; new_del(xy,:)];
    %            cct = cct+1;
    %            [mean(vec(new_dup(xy,1):new_dup(xy,2))/mean_vec) 1 - normcdf(mean(vec(new_dup(xy,1):new_dup(xy,2))),mean_vec,std_vec)]
            end
        end
    end
    tot_dup = [tot_dup;real_dup];
    tot_del = [tot_del;real_del];
    disp(['There are ' num2str(size(cadi_del,1)) ' candidate dels']);
    disp(['There are ' num2str(size(tot_dup,1)) ' duplications, and ' num2str(size(tot_del,1)) ' deletions']);
    tot_dup
    tot_del
%     real_del = new_del;
end
if ~isempty(tot_dup)
    rem_idx = (tot_dup(:,2)-tot_dup(:,1))>OM_ava_len/2;
    tot_dup = tot_dup(rem_idx,:);
end

if ~isempty(tot_dup)
    [~,id1] = sort(tot_dup(:,1));
    tot_dup = tot_dup(id1,:);
%     for i = 1:size(tot_dup,1)
%         close all
%         figure('Visible','off');
%         d_sz = min([tot_dup(i,2)-tot_dup(i,1),tot_dup(i,1)*2,(length(vec)-tot_dup(i,2))*2]);
%         if mod(d_sz,2)==1
%             d_sz = d_sz-1;
%         end
%         tp_dat = (vec(tot_dup(i,1)/res_size/2-d_sz/2:tot_dup(i,2)/res_size/2+d_sz/2)/mean_vec-1)*2;
%         plot(tp_dat);
%         set(gca,'xtick',[d_sz/2+1 length(tp_dat)-d_sz/2],'xticklabel',tot_dup(i,1:2)*res_size,'ytick',-2:2,'yticklabel',[{'-2'} {'-1'} {'0'} {'+1'} {'+2'} ]);
%         saveas(gca,['fig_res/Dup_chr' num2str(chr) '_' num2str(tot_dup(i,1)) '_' num2str(tot_dup(i,2)) '.jpg'],'jpg');
%     end
end
if ~isempty(tot_del)
    [~,id2] = sort(tot_del(:,1));
    tot_del = tot_del(id2,:);
end

end
% for i = 1:size(tot_del,1)
%     d_sz = min(tot_del(i,2)-tot_del(i,1),tot_del(i,1)*2,length(vec)-tot_del(i,2)*2);
%     if mod(d_sz,2)==1
%         d_sz = d_sz-1;
%     end
%     tp_dat = (vec(tot_del(i,1)-d_sz/2:tot_del(i,2)+d_sz/2)/mean_vec-1)*2;
%     plot(tp_dat);
%     set(gca,'xtick',[d_sz/2+1 length(tp_dat)-d_sz/2],'xticklabel',tot_del(i,:)*res_size,'ytick',-2:2,'yticklabel',[{'-2'} {'-1'} {'0'} {'+1'} {'+2'} ]);
%     saveas(gca,['fig_res/Del_chr' num2str(chr) '_' num2str(tot_del(i,1)*res_size) '_' num2str(tot_del(i,2)*res_size) '.jpg'],'jpg');
% end
% vec = vec/mean_vec;
% vec_1st = vec;
% dec = 0.25;
% cnt = 0;
% for i = 1:iter
%     p_vec = vec;
%     ngb = round(log(iter-i+1)/log(2));
%     %ngb = 10;
%     plot(vec-1,'.');
%     vec_o = vec;
%     pause(0.1)
%     %saveas(gca,['fig_' num2str(i-1) '.jpg'],'jpg');
%     idx_1 = vec_o>=1-dec&vec_o<=1+dec;
%     idx_0 = vec_o<=dec;
%     st = 0;
%     for j = 1:length(idx_1)
%         if st == 0 && idx_1(j)==1
%             st = j;
%         end
%         if st~=0 && idx_1(j)==0
%             vec_o(st:j-1) = mean(vec_o(st:j-1));
%             st = 0;
%         end
%     end
%     st = 0;
%     for j = 1:length(idx_0)
%         if st == 0 && idx_0(j)==1
%             st = j;
%         end
%         if st~=0 && idx_0(j)==0
%             vec_o(st:j-1) = mean(vec_o(st:j-1));
%             st = 0;
%         end
%     end
%     for j =11:length(vec)-10
%         vec(j) = mean(vec_o(j-ngb:j+ngb));
%     end
%     if sum(abs(vec-p_vec))/length(vec) < 0.001
%         cnt=cnt+1;
%     else
%         cnt = 0;
%     end
%     if cnt == 5
%         disp(['Converged after ' num2str(i) ' iterations!']);
%         break;
%     end
% end
% idx_1 = vec_1st>=1-dec&vec_1st<=1+dec;
% vec_1st = vec_1st-1;
% vec = vec - 1;
% for i = 2:length(vec)-1
%     if (vec(i)-vec(i-1))*(vec(i)-vec(i+1))>=0
%         tp_vec = vec(i);
% %         
% plot(vec,'.');
% saveas(gca,['fig_' num2str(iter) '.jpg'],'jpg');
