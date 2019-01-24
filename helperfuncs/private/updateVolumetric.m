function updateVolumetric(h_f,M)
h_f.UserData = padarray(M,[1 1 1],'replicate','post');
i = 1;
while i<length(h_f.Children) % Find out which graphics element is the handle for the checkbox
    if isprop(h_f.Children(i),'Style') && strcmp(h_f.Children(i).Style,'checkbox')
        h_f.Children(i).Callback{1}(h_f.Children(i),[],h_f.Children(i).Callback{2}); % invoke callback function to update plot
        break;
    else
        i = i+1;
    end
end