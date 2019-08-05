
%Location for scrabble dictionary. 
url = 'https://www.wordgamedictionary.com/sowpods/download/sowpods.txt';
disp(['Loading data from: ' url ]);
options = weboptions('ContentType','text');
data = webread(url,options);
inputWordArray = splitlines(data);

%Fist 6 lines are preamble,  
inputWordArray = string(inputWordArray(7:end));