function [photo_signal, raw, labels, food, lick, rw] = load_data()

    fprintf('Loading data. \n');
    
    files = dir('*.mat');
    for i = 1:length(files)
        if contains(files(i).name, 'Signal') || contains(lower(files(i).name), 'labels', 'IgnoreCase', true) ...
            || contains(files(i).name, 'RW') || contains(files(i).name, 'Food') || contains(files(i).name, 'Lick') ...
            || contains(files(i).name, 'raw')
            load(files(i).name)
        end
    end
     
    if exist('zPhotoSyncRight', 'var')
        photoSignal = zPhotoSyncRight;
    elseif exist('zPhotoSyncLeft', 'var')
        photoSignal = zPhotoSyncLeft;
    end
    
    if exist('adjustedLabels', 'var')
        labels = adjustedLabels;
    end

    photo_signal = photoSignal;
    rw = RW; food = Food; lick = Lick; 
    if exist('sync_array', 'var'); raw = sync_array; else; raw = []; end
end