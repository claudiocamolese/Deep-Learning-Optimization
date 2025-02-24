function [f,choice]=menu()

datasets = {'Extended Freudenstein and Roth function', 'Problem 75', 'Problem 76'};

% Pop up the menu with the choice of the function
[choice, ok] = listdlg('PromptString', 'Select a functions:', ...
                          'SelectionMode', 'single', ...
                          'ListString', datasets, ...
                          'ListSize', [250, 250], ...
                          'Name', 'Function selection', ...
                          'OKString', 'Select', ...
                          'CancelString', 'Cancel');

% Check if the user correctly select a function
if ok
    % Choice of the function
    switch choice
        case 1
            f = @(x) extendedFreudensteinRoth(x);      
        case 2
          f = @(x) problem75(x);   
        case 3
          f=@(x) problem76(x); 
    end
    disp('Function loaded!');
else
    error('Abort operation');
end
