%Routine for saving a matrix to a file
function saveThisMatrixToFile(a, filename)

[Cols, Rows] = ndgrid( 1 : size(a, 2),   1: size(a, 1) );
M = [Rows(:), Cols(:), a'(:) ];
dlmwrite(strcat(filename,'.dat'), M, 'delimiter', ' ');

endfunction