% declare 
imagedir = 'y:\michael\';
experimentdir = [imagedir,'TF dynamics data\Figure S3\Fig 3A\'];
datadir = [experimentdir,'Data\'];
% declare row, col, and sites to analyze
rows = 2:7;numrows = numel(rows);
cols = [2:11]; numcols = numel(cols);
sites = 1:4; numsites = numel(sites);
plates = [1:2]; numplates = numel(plates);
posmat = 1:numrows*numcols;
posmat = reshape(posmat,numcols,numrows)';

yfpmat = zeros(numrows,numcols,numplates);
mchymat = zeros(numrows,numcols,numplates);
cy5mat = zeros(numrows,numcols,numplates);


%% Loading Data
for plate = 1:numplates
% figure
% set(gcf,'color','white')

    for col = 1:numcols

        for row = 1:numrows
          
            cy5well = [];
            mchywell = [];
            yfpwell = [];
            dnawell = [];

            for site = 1:numsites

                shot = [num2str(rows(row)),'_',num2str(cols(col)),'_',num2str(sites(site))];
                
                try
                    load([datadir,'fixdata_Plate_',num2str(plates(plate)),'_',shot,'.mat']); 
                catch
                    continue
                end
                
                yfp = fixdata(:,8);
                yfpwell = [yfpwell;yfp];
                cy5 = fixdata(:,12);
                cy5well = [cy5well;cy5];
                mchy = fixdata(:,11);
                mchywell = [mchywell; mchy];
                dnawell = [dnawell;fixdata(:,4)];
                    
            end
            badyfp = yfpwell<0 | isnan(yfpwell) | yfpwell > 4000;
            badcy5  = cy5well <0 | cy5well>10000;
            badmchy  = mchywell <0 | isnan(mchywell) | mchywell> 10000;            
            baddna = dnawell<5e5 | dnawell > 2e6;
            badcells = badyfp | badmchy | badcy5 | baddna; 
            yfpwell(badcells) = [];
            mchywell(badcells) = [];
            cy5well(badcells)  =[];
            dnawell(badcells) = []; 

            yfpmat(row,col,plate) = mean(yfpwell);
            mchymat(row,col,plate) = mean(mchywell);
            cy5mat(row,col,plate) = mean(cy5well);

        end

    end

end

%% Rosi Plotting 
rosiconcs = [0.019073486	0.038146973	0.076293945	0.152587891	0.305175781	0.610351563	1.220703125	2.44140625	4.8828125	9.765625	19.53125	39.0625	78.125	156.25	312.5	625	1250	2500	5000	10000
];
rosiconcs = round(rosiconcs*100)/100;

adiponectin = [mchymat(1:2,:,2),mchymat(3:4,:,2)];
glut4 = [yfpmat(1:2,:,2),yfpmat(3:4,:,2)];
cebpa = [cy5mat(1:2,:,2),cy5mat(3:4,:,2)];

dmiadipo = mean(mchymat(6,6:10,1));
dmiglut4 = mean(yfpmat(6,6:10,1));
dmicebpa = mean(cy5mat(6,6:10,1));

figure,notBoxPlot(adiponectin)
set(gca,'XLim',[0 21],'XTick',1:20,'XTickLabel',rosiconcs,'XTickLabelRotation',45,...
    'FontName','Arial','FontSize',15)
xlabel('Rosiglitazone \muM');ylabel('Adiponectin');
line(get(gca,'XLim'),[dmiadipo dmiadipo],'LineWidth',2)

figure,notBoxPlot(glut4)
set(gca,'XLim',[0 21],'XTick',1:20,'XTickLabel',rosiconcs,'XTickLabelRotation',45,...
    'FontName','Arial','FontSize',15)
xlabel('Rosiglitazone \muM');ylabel('glut4');
line(get(gca,'XLim'),[dmiglut4 dmiglut4],'LineWidth',2)

figure,notBoxPlot(cebpa)
set(gca,'XLim',[0 21],'XTick',1:20,'XTickLabel',rosiconcs,'XTickLabelRotation',45,...
    'FontName','Arial','FontSize',15)
xlabel('Rosiglitazone \muM');ylabel('cebpa');
line(get(gca,'XLim'),[dmicebpa dmicebpa],'LineWidth',2)


