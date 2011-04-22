for iN= 1:length(AllNeurons),
   if (0<size(EyeTrack{iN},1))
        SaccPC(iN) = size(EyeTrack{iN}{26},1);
        SaccNC(iN) = size(EyeTrack{iN}{27},1);
        SacRPC(iN) = size(EyeTrack{iN}{28},1);
        SacRNC(iN) = size(EyeTrack{iN}{29},1);

        SaccP(iN,:) = mean(EyeTrack{iN}{26});
        SaccN(iN,:) = mean(EyeTrack{iN}{27});
        SacRP(iN,:) = mean(EyeTrack{iN}{28});
        SacRN(iN,:) = mean(EyeTrack{iN}{29});
   end
end
