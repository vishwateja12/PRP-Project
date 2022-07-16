function [msg] = huffmandec(symb,code,bitstream)

n=length(symb);

% finding maximum length of code[1111010101010101010101010][101,1,0,10]
lengths=[];
for i=1:n
    len=length(char(code(i)));
    lengths=[lengths len];
end

maxlen=max(lengths);

% initialise msg
msg='';

%bitstream length is denoted by streamlen
streamlen=length(bitstream);

% initialising i
i=1;
while i<=streamlen
    % initialising j
    j=0;
    while j<maxlen
        % selecting the bits from i to i+j,where j goes to maxlen of bit
        c=bitstream(i:i+j);
        
        %initialising index
        ind=1;
        
        while (ind<=n && ~isequal(char(code(ind)),c))
            ind=ind+1;
        end
        % if it cant find any equal it will give n+1 so to avoid error we
        % use below if else
        if ind<=n
            msg=[msg char(symb(ind))];
            break
        else
            % increase j to check for the increased bit length
            j=j+1;
        end
    end
    % if we get a message for some i,j then next we should start from i+j+1
    % th bit so the below expression
    i=i+j+1;
end

end


