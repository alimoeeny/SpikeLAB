function n = SelectByMonkey(OldList, Monkey)
    n = {};
    for i = 1:length(OldList)
        if (upper(OldList{i}(1)) == upper(Monkey(1)))
            n{length(n)+1} = OldList{i};
        end
    end