function [CellValid, ValidityScore] = CellValidity(Expt, testType, varargin)
% for TWO we check to see if the performance had been acceptable, meaning
% the effect of background cylinder disparity on the dx=0 had been larger
% than th eflip effect also there is an additional clause that could be
% enabled indicating background and forground can not be in the same hemifield!

CellValid =  1; % All cells are valid unless proven otherwise 
ValidityScore = [];

if isempty(Expt)
    MonkeyName = varargin{2};
    NeuronNumber = varargin{3};
    ClusterName = varargin{4};
    FileType = varargin{5};
    DataPath = GetDataPath();
    StimulusType = 'cylinder';
    %filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');
    %filename = MakeFileName(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    Neuron = load(filepath);
    Expt = Neuron.Expt;
end


switch testType
    case 'TWO' % (Expt, 'TWO', pD, )
        pD = varargin{1};
        dxValues = unique([Expt.Trials(:).dx]);
        
        p1 = sum([Expt.Trials(:).dx] == dxValues(1)   & [Expt.Trials(:).bd] == dxValues(1)   & [Expt.Trials(:).RespDir] == 1) / sum([Expt.Trials(:).dx] == dxValues(1)   & [Expt.Trials(:).bd] == dxValues(1)   & [Expt.Trials(:).RespDir] ~= 0);
        p2 = sum([Expt.Trials(:).dx] == dxValues(end) & [Expt.Trials(:).bd] == dxValues(end) & [Expt.Trials(:).RespDir] == 1) / sum([Expt.Trials(:).dx] == dxValues(end) & [Expt.Trials(:).bd] == dxValues(end) & [Expt.Trials(:).RespDir] ~= 0);
        q1 = sum([Expt.Trials(:).dx] == dxValues(1)   & [Expt.Trials(:).bd] == -dxValues(1)   & [Expt.Trials(:).RespDir] == 1) / sum([Expt.Trials(:).dx] == dxValues(1)   & [Expt.Trials(:).bd] == -dxValues(1)   & [Expt.Trials(:).RespDir] ~= 0);
        q2 = sum([Expt.Trials(:).dx] == dxValues(end) & [Expt.Trials(:).bd] == -dxValues(end) & [Expt.Trials(:).RespDir] == 1) / sum([Expt.Trials(:).dx] == dxValues(end) & [Expt.Trials(:).bd] == -dxValues(end) & [Expt.Trials(:).RespDir] ~= 0);
        z1 = sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd] == dxValues(1)   & [Expt.Trials(:).RespDir] == 1) / sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd] == dxValues(1)   & [Expt.Trials(:).RespDir] ~= 0);
        z2 = sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd] == dxValues(end) & [Expt.Trials(:).RespDir] == 1) / sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd] == dxValues(end) & [Expt.Trials(:).RespDir] ~= 0);
       
%         if abs(z1-z2) > ((abs(p1-q1)+ abs(p2-q2))/2) 
%             CellValid = 1;
%         else 
%             CellValid = 0;
%         end
%         ValidityScore = abs(z1-z2) - ((abs(p1-q1)+ abs(p2-q2))/2);

        if ~isfield(Expt.Stimvals, 'backxo'), Expt.Stimvals.backxo = - sign(Expt.Stimvals.xo); end
        if abs(z1-z2) > max(abs(p1-q1), abs(p2-q2)) ...
                && (sign(Expt.Stimvals.xo)~=sign(Expt.Stimvals.backxo)) % This is the additional clause indicating background and forground can not be in the same hemifield!
            CellValid = 1;
        else 
            CellValid = 0;
        end
        ValidityScore = abs(z1-z2) - abs(z1-z2) > max(abs(p1-q1), abs(p2-q2));

end

