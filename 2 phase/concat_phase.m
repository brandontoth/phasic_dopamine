function ph_struct = concat_phase(save_loc)
%% initialize struct
ph_struct = struct;

%% navigate to desired dir
filedir = save_loc;
cd(filedir);

%% identify files in dir
file_pat  = fullfile(filedir, '*.mat');
the_files = dir(file_pat);

%% loop through files
j = 1;
for i = 1:length(the_files)
    cur_file   = load(the_files(i).name);
    cur_struct = getVarName(cur_file);

    if j > 1
        [store, cur_struct] = checkFields(store, cur_struct);
    end
    
    cur_struct = cell2struct(cellfun(@(x)(cur_struct.(x)), ...
        fieldnames(cur_struct), 'uni', false), fieldnames(cur_struct), 1);
    store(i) = cur_struct;

    j = j + 1;
end

fields = string(fieldnames(store));

for i = 1:length(fields)
    if size(store(1).(fields(i)), 1) == 1
        for j = 1:length(store)
            if isnan(store(j).(fields(i)))
                store(j).(fields(i)) = [];
            end
            store(j).(fields(i)) = store(j).(fields(i))';
        end
    end

    fet_field = fields(i);
    
    ph_struct.(fet_field) = vertcat(store.(fields(i)))';
    ph_struct.(fet_field) = rmmissing(ph_struct.(fet_field), 2);
end

end