function [] = txt_export_all(Params,Gr,mom,flag)
% In this function we export all model results to txt files. We do not
% display on the screen
% screen

if ~isstruct(Params)
    error('Input Params must be a structure')
end
if ~isstruct(Gr)
    error('Input Gr must be a structure')
end
if ~isstruct(mom)
    error('Input mom must be a structure')
end
if ~isstruct(flag)
    error('Input flag must be a structure')
end

% Unpack variables from struct "mom"
% ave = mom.ave;
% ave_young = mom.ave_young;
% ave_health = mom.ave_health;
% ave_health_young = mom.ave_health_young;
% gini = mom.gini;
% gini_health = mom.gini_health;
% corr = mom.corr;

dir_name = flag.dir_name;
myfolder = fullfile(flag.results,dir_name);
if ~isfolder(myfolder)
    mkdir(myfolder)
end

% Get names ave, ave_young, etc
Names_mom = fieldnames(mom); 
for ii = 1:numel(Names_mom)
    Name_ii = Names_mom{ii};
    Names_subfield = fieldnames(mom.(Name_ii));
    for jj = 1:numel(Names_subfield)
        Name_jj = Names_subfield{jj};
        name_var = [Name_ii,'_',Name_jj,'.txt'];
        file_name = fullfile(flag.results,dir_name,name_var);
        myexport(mom.(Name_ii).(Name_jj),file_name);
    end
end

end %end function


function [] = myexport(X,file_name)
% Export variable X to file file_name (it is a txt file)

X_vec = reshape(X,[numel(X),1]);
writematrix(X_vec, file_name);
% dlmwrite is obsolete but it may be faster
%dlmwrite('filename.txt',a)

end