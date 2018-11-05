function btdstruct = subMeanFromStim(btdstruct)



for i = 1:length(btdstruct.btd)
    btdstruct.btd(1, i).glt(1, 1).yData = btdstruct.btd(1, i).glt(1, 1).yData - 128;
end



end

