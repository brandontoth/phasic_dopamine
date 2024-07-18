function final = append_data

selpath = uigetdir();
filedir = selpath;
cd(filedir);

filePattern = fullfile(filedir, '*.mat');
theFiles = dir(filePattern);

j = 1;
for i = 1:length(theFiles)
    curFile = load(theFiles(i).name);
    thisStruct = getVarName(curFile);

    if j > 1
        [store, thisStruct] = checkFields(store, thisStruct);
    end
    
    thisStruct = cell2struct(cellfun(@(x)(thisStruct.(x)), ...
        fieldnames(thisStruct), 'uni', false), fieldnames(thisStruct), 1);
    store(i) = thisStruct;

    j = j + 1;
end

fields = string(fieldnames(store));

for i = 1:length(fields)
    if contains(fields(i), 'loc', 'IgnoreCase', true)
        store = rmfield(store, fields(i));
    end
end

fields(contains(fields, 'lick', 'IgnoreCase', true)) = [];

for i=1:length(fields)
    if size(store(1).(fields(i)), 1)>size(store(1).(fields(i)),2)
        for j=1:length(store)
            store(j).(fields(i))=store(j).(fields(i))'; 
        end 
    end
    getField = fields(i);
    
    final.(getField) = horzcat(store.(fields(i)))';
end

end