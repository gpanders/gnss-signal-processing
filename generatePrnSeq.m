function [ prn ] = generatePrnSeq( txId )
%GENERATEPRNSEQ Generate GPS L1 C/A PRN sequences
%   GENERATEPRNSEQ(I) produces the 1023-by-1 +/- 1 valued spreading code for the
%   satellite with PRN I

% Define parameters for GPS L1 C/A Gold codes
ciVec1 = [3, 10]';
ciVec2 = [2, 3, 6, 8, 9, 10]';
a0Vec = ones(10, 1);

F1 = generateLfsrSequence(10, ciVec1, a0Vec);
F2 = generateLfsrSequence(10, ciVec2, a0Vec);

% TODO: refactor this to a switch statement to avoid overhead of loading
% this array every time this function is called
DELAY_CHIPS = [5, 6, 7, 8, 17, 18, 139, 140, 141, 251, 252, 254, 255, ...
              256, 257, 258, 469, 470, 471, 472, 473, 474, 509, 512, ...
              513, 514, 515, 516, 859, 860, 861, 862, 863];
              
prn = 2*mod(F1 + circshift(F2, DELAY_CHIPS(txId)), 2) - 1;

end

