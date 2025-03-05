function colorStr=setColor(type)
    switch type
        case 2
            colorStr="\color{green}";
        case 1
            colorStr="\color{blue}";
        case 0
            colorStr="\color{black}";
        case -1
            colorStr="\color{red}";
        otherwise
            colorStr="\color{magenta}";
    end
end