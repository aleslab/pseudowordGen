%% Load a dictionary file.

%These functions download and load in a couple of open access list of
%words.
%Drugbank open access data rom: https://www.drugbank.ca
loadDrugnameDict;

%The European SOWPODS scrabble dictionary:
%loadScrabbleDict;

%% Generate substring list that describe each markov "state"

disp('Generating substring list');
subStrLen = 4; %length of substrings to consider

%The following is a crazy matlab way of doing things.  Most people would
%use nested for loops.  But it makes it hard to extend to arbitrary
%substring length sequences.  Instead I'm using a hard to understand
%combination of ndgrid and cell arrays.  ndgrid is a great way to create
%all possible combinations of a list of things.  

%Build list of strings, use asci character values first for ease:
%charList = ([97:(97+25)]); %Ascii lower case a-z
charList = ([97:(97+25) 32 36]); %Ascii lower case a-z plus space ' ' and $

repCharList = repmat({charList}, 1, subStrLen);
%Use a cell array to grab different number of output arguments changeable using single subStrLen parameter
allSubStr = cell(1,subStrLen);
[allSubStr{:}] = ndgrid(repCharList{:});
%Funky way to reshape the matrices above.  Use implicit indexing:  "(:)" to
%reshape matrices into vecotrs.  Use cellfun() to apply the vector reshape to each cell. 
%Could also have used reshape(), but that requires more thinking;
allSubStr = cellfun(@(x) x(:),allSubStr,'UniformOutput',false);

%Finaly concatenate the ascii values together (this crazy matlab syntax [ myCell{:} ]
%and change it back to a chartacter array.
subStringChar = char( [allSubStr{:}]);

%Now let's make a lookuptable to quickly find this sustrings
strLookupTbl = containers.Map(string(subStringChar),num2cell(1:length(subStringChar)));


%% Splice up input word dictionary into subStrLen tuplet sequences +1 for next letter:
%Again this uses crazy matlab syntax to accomplish the goal. 


%First prepend with spaces to indicate start of word and $ to indicate end of
%word.  We use spaces for the start so we can make the start of words have
%a different probabilities than the middle.  
inputWordArray = strcat(string(char(repmat(32,1,subStrLen))),inputWordArray,'$');
%Now convert to numeric (no reason to use double(inefficient), but lazy
%matlab default because matlab doesn't always play nice when numeric values
%are different) Usinging numeric so I can use efficient numeric
%calculations below. I'm taking the input words concatenating them all
%toghether, changing them to lowecase, then changing them to be interpreted
%as numbers instead of letters.
inputWordArray = double(lower([inputWordArray{:}]));

%Ok.  So this line was stupid, broke the nesting out below to comment  
%inputTuplets = char(fliplr(toeplitz(inputWordArray,inputWordArray(1:(subStrLen+1)))));

%To get all the substrings from the input word we're using the toeplitz
%function.  Often used to create sliding timewindows in timeseries data.
%This is a weird use to create sliding string windows in the word input. 
%We make our sliding window 1 greater than the subtring length we want so
%we can keep track of the input substring and the letter that follows the
%substring. 
subStringMatrix = toeplitz(inputWordArray,inputWordArray(1:(subStrLen+1)));
%By default the toeplitz represents things shiftinng going backwards(blame
%matrix algebra), to make these easier to read for humans I'm flipping.
subStringMatrix = fliplr(subStringMatrix);
%Now take the numbers and turn them back into string representation.
% (Note: this ultimately slows things down, but makes it much easier for people to
% read)
inputTuplets = char(subStringMatrix);

%% Create markov state space. 

disp('Analyzing probabilities from input word list');

%Now generate a state table: Number of input substrings x number of
%possible output characters. 

%This is the heart of the "learning" part.  Basicaly we just go through the
%input word list substrings see what letter follows the substring, we find
%the corresponding markov "state"  and increment a counter for the output
%letter.

%Predfine things here 
nInputStates = size(subStringChar,1);
nInputTuplets = size(inputTuplets,1);

%This is important to pre-allocate because it's big and the loop below is
%long. 
stateTable = zeros(nInputStates,length(charList));

%Generate a simple lookuptable. 
%Doing this for efficiency reasons/understanding reasons.  
%Techinically a markov state transition table should be nStates x nStates,  
%But in this specific case we are generating a single letter so we just
%need 26 outputs not nState outputs
%Technically this is a new state: e.g. if the state is ABC and the next
%char is D then the new state is BCD
%
%Could use a sparse matrix, to represent the full nstate x nstate but it
%gets large and unwieldy for this use case.


%I want to easily index into the output characters, but Ascii characters
%start at 97, so build a simple lookup table that takes an ascii character
%value and remaps it into the the allowed output chars. 
%Don't really need a sparse matrix here, using it for example.
charIdxLookup = sparse(1,300);
charIdxLookup(charList) = 1:length(charList);

%Now loop over each of spliced input tuplets and build up the observed transition
%matrix
waitH = waitbar(0,'Analyzing substring probabilities in input');
tic
for iInTup = 1:nInputTuplets,
    
    
    if mod(iInTup,1000)==0,
        elapsedTime = toc;
        timePerTuplet = elapsedTime/iInTup;
        timeLeft =(nInputTuplets-iInTup)*timePerTuplet;
       
        waitbar(iInTup/nInputTuplets,waitH,...
            ['Analyzing ' num2str(nInputTuplets) ' substring probabilities, time left: ' num2str(timeLeft,3) ' seconds']);

    end
    
    %Grab the first substrlen characters that describe the state
    thisTuplet = inputTuplets(iInTup,1:subStrLen);
    
    %if this tuplet doesn't exist in our table move on
    %This catches wierd inputs in the training set that we don't want.
    if ~strLookupTbl.isKey(thisTuplet)
        continue;
    end
    
    %The last char of the tuplut the character that follows thisTuplet
    %is the output letter
    nextChar = double(inputTuplets(iInTup,end));
   
    %Use the lookuptable to index into the table.
    nextCharIdx = charIdxLookup(nextChar);
    
    %Now we need to find wich row of the table the current tuplet
    %correponds with.
    stateIdx = strLookupTbl(thisTuplet);
    
    %If we cant find the tuplet in our state list, or can't find the
    %character in our output character list skip this tuplet.
    if isempty(stateIdx) || nextCharIdx==0,
        continue;
    end
    
    %Finaly, increment our count of frequency of what the next character is
    %given the current substring value:
    stateTable(stateIdx,nextCharIdx) = stateTable(stateIdx,nextCharIdx)+1;
end

close(waitH);


%The above code just counts raw frequencies.  To use this we need to
%convert to probabilities.  That is, for any given state what is the
%probabiility to generate the next char.  
%Normalize the output probabilities.
for iState = 1:size(stateTable,1),
    
 stateTable(iState,:) =    stateTable(iState,:)./sum(stateTable(iState,:));
end


%% Generate pseudowords
%Given a state table created above we can now generate strings that follow
%the markov probability. NOTE: These are not guaranteed to be unique.  It
%is quite probably that entries from the training set show up. 


%Initialize the state of the markov sequence generation. We mark out word
%beginnings by being precceded by only ' '(spaces). So we initialize the
%sequence with an the all space character state.
markovSeq =  repmat(' ',1,subStrLen) ;

finished = false;

while length(markovSeq)<2000;
    
    %Check the current string in the sequence and use the lookup table to
    %find out what row the the state table that corresponds too. 
    stateIdx = strLookupTbl(markovSeq((end-subStrLen+1):end));
    
    %Now in order to generate a markov sequence we need to draw a random
    %number that follows the probability distribution given for this state.
    %We do that by converting the probability density to a cummulative
    %probabilty.  We then generate a random number between 0 and 1 and

    %add a 0 because cumulative probability function starts at 0.
    thisTransitionProb = [0 stateTable(stateIdx,:)];
    
    %Now cummulative summate all probabilitiies to create cumulative distribution
    cumProb  = cumsum( thisTransitionProb')';

    %Generatye a random number between 0 and 1 (rand()).  Then compare it
    %with the cumulative distribution and find the index in the table it
    %corresponds with.  
    %This is an inefficient way to generate draws. But works.   
    newChar = find(cumProb<rand(),1,'last');
    
    %Take the output index value and turn it into the correct character
    %value using the lookup table..
    newChar = char(charList(newChar));
    
    %Finaly, grow the string sequence. 
    markovSeq(end+1) = newChar;
end


