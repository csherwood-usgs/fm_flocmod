% fm_kernal_stats
knames = {'f_coll_prob_sh',...
'f_coll_prob_ds',...
'f_g1_sh',...
'f_g1_ds',...
'f_g3',...
'f_l3'};

for i=1:length(knames);
fprintf('%s ',knames{i})
%eval(['size( ',char(knames{i}),')'])
eval(['val(1)=min( ',char(knames{i}),'(:));'])
eval(['val(2)=max( ',char(knames{i}),'(:));'])
eval(['val(3)=sum( ',char(knames{i}),'(:));'])
eval(['val(4)=sum(abs( ',char(knames{i}),'(:)));'])
%fprintf(1,'Min: %g Max: %g Sum: %g SumAbs: %g\n',val)
fprintf('%g\n',val(3))
end