currentDirectory = pwd;

if(ispc)
	breakCharactor = ';';
else
	breakCharactor = ':';
end

paths = [...
	[currentDirectory, '/matlab/', breakCharactor], ...
	[currentDirectory, '/matlab/mex', breakCharactor], ...
];

addpath(paths);
