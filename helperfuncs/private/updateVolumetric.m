function updateVolumetric(h_f,M)

h_f.UserData(1:end-1,1:end-1,1:end-1) = single(M);
h_f.UserData(end,:,:) = h_f.UserData(end-1,:,:);
h_f.UserData(:,end,:) = h_f.UserData(:,end-1,:);
h_f.UserData(:,:,end) = h_f.UserData(:,:,end-1);

i = 1;
while i<length(h_f.Children) % Find out which graphics element is the handle for the checkbox
    if isprop(h_f.Children(i),'Style') && strcmp(h_f.Children(i).Style,'checkbox')
        h_f.Children(i).Callback{1}(h_f.Children(i),[],h_f.Children(i).Callback{2}); % invoke callback function to update plot
        break;
    else
        i = i+1;
    end
end