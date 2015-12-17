classdef pairedRDS < StimulusGenerator
    
    % Default properties
    properties
        dotsize=3; % dot width in pixels
        density=0.24; % dot density; proportion of stimulus coverage assuming no overlap
        dx=0; % horizontal disparity in pixels
        dy=0; % vertical disparity in pixels
        correlation=1; %-1 is anti-correlated; 0 is uncorrelated; 1 is correlated;
        
        Nx=292; % stimulus width in pixels
        Ny=292; % stimulus height in pixels
                
        annulus_dx;
        annulus_correlation;
        
        secondary_dx=0;
        secondary_dx_fraction=0;
        secondary_dx_dotmatch=0;
        
        rho1=0;
        rho2=0;
        
        
    end
    
    properties (Hidden)
        % None
    end
    
    methods
        
        % Constructor
        function rds = pairedRDS(varargin)
            % Constructs an RDS object which can be used to generate
            % binocular stimuli (left and right images).             
            % Usage: rds = RDS(varargin);
            % rds = RDS() will initialise an RDS object with default
            % stimulus parameters
            % 
            % ===optional arguments===
            % dotsize: dot width in pixels
            % density: dot density (fraction of occupied space assuming no
            % dot occlusion; note: dots are allowed to occlude)
            % correlation: binocular correlation of stimulus (-1 to 1)
            % dx: horizontal disparity of stimulus in pixels
            % rds.dy: vertical disparity of stimulus in pixels            
            % Nx: stimulus width in pixels
            % Ny: stimulus height in pixels
                                                
            
            for j = 1:2:length(varargin);
                assert(isprop(bem,varargin{j}),'Property not in generator structure');
                rds.(varargin{j}) = varargin{j+1};
            end
            
            
        end
        
        function [LeftImage,RightImage] = generate(rds)
            % The structure of this function:
            % Section 1: Defines the parameters we need to create our
            % paired RDS. It also divides the image into subsets so that we can
            % induce a lateral shift of the central disk and paint the surrounding
            % area appropriately. There are three key subsets of the stimulus: the
            % annulus (left and right will differ if there is disparity), the
            % central disk (will be identical, but with positional shift) and the
            % new area left after the positional shift of the disk. The dots in the
            % new area in the left and right images are uncorrelated.

            % Section 2: Specifies the fill of the dots that we are
            % painting (i.e. 1 or -1). What exactly we're doing here will depend on
            % dotmatch, rho1 and rho2.

            % Section 3: defines the location of the dots.

            % The final for loop paints the dots.

            assert((rds.secondary_dx_fraction <= 1) && (rds.secondary_dx_fraction >= 0),'Proportion disparity2 dots must be in range [0,1]');
            assert((rds.secondary_dx_dotmatch <= 1) && (rds.secondary_dx_dotmatch >= 0),'Dot match for disparity2 dots must be in range [0,1]');
            
            if isempty(rds.annulus_correlation);
                rds.annulus_correlation = rds.correlation;
            end
            
            if isempty(rds.annulus_dx);
                rds.annulus_dx = 0;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Section 1: Parameters and subsets %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            

            % If we want the annulus to have a different binocular sign than the
            % disk, we can specify that here.
            ActualLeftDot = new_dot(rds.dotsize);
            [leftrow,leftcol] = find(abs(ActualLeftDot) > 10e-2);
            LeftDot(:,1) = leftrow;
            LeftDot(:,2) = leftcol;
            LeftPixels = ActualLeftDot(abs(ActualLeftDot) > 10e-2);
            RightDot = LeftDot; RightPixels = LeftPixels;
            AnnulusDot = RightDot; AnnulusPixels = RightPixels;

            % move from correlation to dotmatch
            AnnulusDotmatch = (rds.annulus_correlation+1)*0.5;
            DiskDotmatch = (rds.correlation+1)*0.5;


            DiskSize = 2.5/3.5; % This is the ratio used by Doi et al. (2013)


            PairFlip = 2*rds.rho2 -1; % -1 if rho2 = 0, 1 if rho2 = 1;
            AnnulusRadius = 0.8*rds.Nx/2;
            DiskRadius = AnnulusRadius*DiskSize;

            DotDiameter = max(LeftDot(:,1))-min(LeftDot(:,1)) + 1; % Get the diameter of the dot
            DotArea = length(LeftDot);

            voffset = round(1.1*DotDiameter); % Vertical distance between the two dots in a pair

            Outcomes = [1,1;-1,-1;1,-1;-1,1];

            % Define the coordinates for the annulus and the disk
            HalfNy = round(rds.Ny/2); HalfNx = round(rds.Nx/2);
            [Y,X] = meshgrid(-linspace(-HalfNy,HalfNy,rds.Ny),linspace(-HalfNx,HalfNx,rds.Nx));
            K = sqrt(Y.^2 + X.^2);
            [Disk_Y,Disk_X] = find(K <= DiskRadius);
            Annulus = find(K <= AnnulusRadius);

            LeftDiskY = Disk_Y - ceil(rds.dy/2); LeftDiskX = Disk_X - ceil(rds.dx/2);
            RightDiskY = Disk_Y + floor(rds.dy/2); RightDiskX = Disk_X + floor(rds.dx/2);

            Disk = sub2ind([rds.Ny,rds.Nx],Disk_Y,Disk_X);
            LeftDisk = sub2ind([rds.Ny,rds.Nx],LeftDiskY,LeftDiskX);
            RightDisk = sub2ind([rds.Ny,rds.Nx],RightDiskY,RightDiskX);

            %Create left and right annuli
            LeftAnnulus = setdiff(Annulus,LeftDisk);
            RightAnnulus = setdiff(Annulus,RightDisk);

            IntersectAnnulus = intersect(LeftAnnulus,RightAnnulus);
            [IntersectAnnulusY,IntersectAnnulusX] = ind2sub([rds.Ny,rds.Nx],IntersectAnnulus);

            DiffDiskL = setdiff(RightDisk,LeftDisk); % These are the areas that are left "clear" after shifting
            DiffDiskR = setdiff(LeftDisk,RightDisk);

            [DiffDiskLY,DiffDiskLX] = ind2sub([rds.Ny,rds.Nx],DiffDiskL);
            [DiffDiskRY,DiffDiskRX] = ind2sub([rds.Ny,rds.Nx],DiffDiskR);

            AreaCleared = round((length(DiffDiskL) + length(DiffDiskR))/2);

            % Number of dots to paint; painting annulus and disk separately

            nDotsDisk = round(length(Disk) * rds.density / DotArea); % For the central disk
            nDotsAnnulus = round(length(IntersectAnnulus) * rds.density / DotArea); % For the surrounding annulus
            nDotsNewArea = round(AreaCleared * rds.density / DotArea);

            % Number of vertical pairs formed
            nPairedDotsDisk = round(rds.rho1/2*nDotsDisk);
            nPairedDotsAnnulus = round(rds.rho1/2*nDotsAnnulus);
            nPairedDotsNewArea = round(rds.rho1/2*nDotsNewArea);

            % Number of paired dots that are mismatched
            % We don't need this for NewArea because dots in this area aren't
            % binocularly paired
            nDotsDiskMismatch = nPairedDotsDisk - round(nPairedDotsDisk*DiskDotmatch);
            nDotsAnnulusMismatch = nPairedDotsAnnulus - round(nPairedDotsAnnulus*AnnulusDotmatch);

            % Number of singleton dots
            nSingletonDotsDisk = round((1-rds.rho1)*nDotsDisk);
            nSingletonDotsAnnulus = round((1-rds.rho1)*nDotsAnnulus);
            nSingletonDotsNewArea = round((1-rds.rho1)*nDotsNewArea);

            % Number of singleton dots that are mismatched
            nSingletonDotsDiskMismatch = nSingletonDotsDisk - round(nSingletonDotsDisk*DiskDotmatch);
            nSingletonDotsAnnulusMismatch = nSingletonDotsAnnulus - round(nSingletonDotsAnnulus*AnnulusDotmatch);



            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Section 2: Generate fill %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            LeftImage = zeros(rds.Ny,rds.Nx); RightImage = zeros(rds.Ny,rds.Nx);

            % First we paint the disk
            % Assign the fill for the initial dots

            LeftDiskInitialFill = [zeros(1,ceil(nPairedDotsDisk/2))-1, ones(1,floor(nPairedDotsDisk/2))];
            RightDiskInitialFill = LeftDiskInitialFill;

            DiskMismatchIndices = randperm(nPairedDotsDisk,nDotsDiskMismatch);
            RightDiskInitialFill(DiskMismatchIndices) = RightDiskInitialFill(DiskMismatchIndices) * -1;

            % Assign the fill for the complementary dots -- This is
            % a bit involved...
            DiskNewContrastIndices = abs(LeftDiskInitialFill .* RightDiskInitialFill * PairFlip - 1);
            DiskWhichOutcome = randi(2,1,length(DiskNewContrastIndices));
            LeftDiskComplementaryFill = Outcomes(DiskNewContrastIndices+DiskWhichOutcome,1)';
            RightDiskComplementaryFill = Outcomes(DiskNewContrastIndices+DiskWhichOutcome,2)';

            % Assign the fill for the remaining (singleton) dots
            LeftDiskRemainingFill = [zeros(1,ceil(nSingletonDotsDisk/2))-1, ones(1,floor(nSingletonDotsDisk/2))];
            RightDiskRemainingFill = LeftDiskRemainingFill;
            DiskMismatchIndices = randperm(nSingletonDotsDisk,nSingletonDotsDiskMismatch);
            RightDiskRemainingFill(DiskMismatchIndices) = RightDiskRemainingFill(DiskMismatchIndices) * -1;

            % All fills together now...
            LeftDiskFill = [LeftDiskInitialFill, LeftDiskComplementaryFill, LeftDiskRemainingFill];
            RightDiskFill = [RightDiskInitialFill, RightDiskComplementaryFill, RightDiskRemainingFill];

            % And now the fill for the dots painted in the area cleared after the
            % disk shift
            LeftClearedFill = randi([0,1],1,nDotsNewArea)*2 -1;
            RightClearedFill = randi([0,1],1,nDotsNewArea)*2 -1;

            LeftAnnulusInitialFill = [zeros(1,ceil(nPairedDotsAnnulus/2))-1, ones(1,floor(nPairedDotsAnnulus/2))];
            RightAnnulusInitialFill = LeftAnnulusInitialFill;

            AnnulusMismatchIndices = randperm(nPairedDotsAnnulus,nDotsAnnulusMismatch);
            RightAnnulusInitialFill(AnnulusMismatchIndices) = RightAnnulusInitialFill(AnnulusMismatchIndices) * -1;

            % Assign the fill for the complementary dots -- This is
            % a bit involved...
            AnnulusNewContrastIndices = abs(LeftAnnulusInitialFill .* RightAnnulusInitialFill * PairFlip - 1);
            AnnulusWhichOutcome = randi(2,1,length(AnnulusNewContrastIndices));
            LeftAnnulusComplementaryFill = Outcomes(AnnulusNewContrastIndices+AnnulusWhichOutcome,1)';
            RightAnnulusComplementaryFill = Outcomes(AnnulusNewContrastIndices+AnnulusWhichOutcome,2)';

            % Assign the fill for the remaining (singleton) dots
            LeftAnnulusRemainingFill = [zeros(1,floor(nSingletonDotsAnnulus/2))-1, ones(1,ceil(nSingletonDotsAnnulus/2))];
            RightAnnulusRemainingFill = LeftAnnulusRemainingFill;
            AnnulusMismatchIndices = randperm(nSingletonDotsAnnulus,nSingletonDotsAnnulusMismatch);
            RightAnnulusRemainingFill(AnnulusMismatchIndices) = RightAnnulusRemainingFill(AnnulusMismatchIndices) * -1;



            % All fills together now...
            LeftAnnulusFill = [LeftAnnulusInitialFill, LeftAnnulusComplementaryFill, LeftAnnulusRemainingFill];
            RightAnnulusFill = [RightAnnulusInitialFill, RightAnnulusComplementaryFill, RightAnnulusRemainingFill];


            LeftFill = [LeftAnnulusFill, LeftDiskFill, LeftClearedFill];
            RightFill = [RightAnnulusFill, RightDiskFill, RightClearedFill];



            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Section 3: Generating dot centres %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Disk centres
            DiskInitialIndices = randperm(length(Disk),nPairedDotsDisk);
            DiskRemainingIndices = randperm(length(Disk),nSingletonDotsDisk);

            DiskInitialCenters = [Disk_Y(DiskInitialIndices),Disk_X(DiskInitialIndices)];
            DiskComplementaryCenters = [DiskInitialCenters(:,1)+voffset,DiskInitialCenters(:,2)];
            DiskRemainingCenters = [Disk_Y(DiskRemainingIndices),Disk_X(DiskRemainingIndices)];

            DiskCenters = [DiskInitialCenters; DiskComplementaryCenters; DiskRemainingCenters];

            % Annulus centres
            AnnulusInitialIndices = randperm(length(IntersectAnnulus),nPairedDotsAnnulus);
            AnnulusRemainingIndices = randperm(length(IntersectAnnulus),nSingletonDotsAnnulus);

            AnnulusInitialCenters = [IntersectAnnulusY(AnnulusInitialIndices),IntersectAnnulusX(AnnulusInitialIndices)];
            AnnulusComplementaryCenters = [AnnulusInitialCenters(:,1)+voffset,AnnulusInitialCenters(:,2)];
            AnnulusRemainingCenters = [IntersectAnnulusY(AnnulusRemainingIndices), IntersectAnnulusX(AnnulusRemainingIndices)];

            AnnulusCenters = [AnnulusInitialCenters; AnnulusComplementaryCenters; AnnulusRemainingCenters];

            % NewArea centres, these are specified for the left and right eyes
            % separately as the dots in these areas are not binocularly paired (that is, they are monocular singletons)

            % First the left image's singletons
            NewAreaLeftInitialIndices = randperm(length(DiffDiskL),nPairedDotsNewArea);
            NewAreaLeftRemainingIndices = randperm(length(DiffDiskL),nSingletonDotsNewArea);

            NewAreaLeftInitialCenters = [DiffDiskLY(NewAreaLeftInitialIndices),DiffDiskLX(NewAreaLeftInitialIndices)];
            NewAreaLeftComplementaryCenters = [NewAreaLeftInitialCenters(:,1)+voffset,NewAreaLeftInitialCenters(:,2)];
            NewAreaLeftRemainingCenters = [DiffDiskLY(NewAreaLeftRemainingIndices),DiffDiskLX(NewAreaLeftRemainingIndices)];

            NewAreaLeftCenters = [NewAreaLeftInitialCenters;NewAreaLeftComplementaryCenters;NewAreaLeftRemainingCenters];

            % Now the right image's singletons
            NewAreaRightInitialIndices = randperm(length(DiffDiskR),nPairedDotsNewArea);
            NewAreaRightRemainingIndices = randperm(length(DiffDiskR),nSingletonDotsNewArea);

            NewAreaRightInitialCenters = [DiffDiskRY(NewAreaRightInitialIndices),DiffDiskRX(NewAreaRightInitialIndices)];
            NewAreaRightComplementaryCenters = [NewAreaRightInitialCenters(:,1)+voffset,NewAreaRightInitialCenters(:,2)];
            NewAreaRightRemainingCenters = [DiffDiskRY(NewAreaRightRemainingIndices),DiffDiskRX(NewAreaRightRemainingIndices)];

            NewAreaRightCenters = [NewAreaRightInitialCenters;NewAreaRightComplementaryCenters;NewAreaRightRemainingCenters];

            % Now we put all the dots together in a single matrix so that we can
            % paint the dots sequentially
            AllLeftCenters = [[AnnulusCenters(:,1)-ceil(rds.annulus_dx/2),AnnulusCenters(:,2)];[DiskCenters(:,1) - ceil(rds.dy/2),DiskCenters(:,2) - ceil(rds.dx/2)]; NewAreaLeftCenters];
            AllRightCenters = [[AnnulusCenters(:,1)+floor(rds.annulus_dx/2),AnnulusCenters(:,2)];[DiskCenters(:,1) + floor(rds.dy/2),DiskCenters(:,2) + floor(rds.dx/2)];NewAreaRightCenters];

            nDotsToPaint = min(length(LeftFill),length(AllLeftCenters));

            % This is for the disparity2 parameter: We simply choose a subset of
            % all of these dots and set their centers to be rds.secondary_dx. This is
            % a hack, but the option isn't oging to be very common so whatevs
            %%% =========
            nDotsDisparity2 = round(nDotsToPaint*rds.secondary_dx_fraction);
            disparity2Indices = randperm(nDotsToPaint,nDotsDisparity2);
            AllRightCenters(disparity2Indices,2) = AllLeftCenters(disparity2Indices,2)+rds.secondary_dx;

            nDisparity2Revert = round(nDotsDisparity2*(1-rds.secondary_dx_fraction)); % Number of dots to reverse contrast

            RightFill(disparity2Indices(1:nDisparity2Revert)) = RightFill(disparity2Indices(1:nDisparity2Revert))*-1;
            RightFill(disparity2Indices((nDisparity2Revert+1):nDotsDisparity2)) = RightFill(disparity2Indices((nDisparity2Revert+1):nDotsDisparity2));
            %%% ==========

            RandomisedIndices = randperm(nDotsToPaint);

            % Let's play around with this to see what gives us the most consistent
            % placement of the disk
            RandomisedLeftCenters = AllLeftCenters(RandomisedIndices,:) - 2*DotDiameter; 
            RandomisedRightCenters = AllRightCenters(RandomisedIndices,:) - 2*DotDiameter;
            LeftFill = LeftFill(RandomisedIndices);
            RightFill = RightFill(RandomisedIndices);


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Section 4: Paint the dots %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            for dots = 1:nDotsToPaint;
                LeftYCenter = RandomisedLeftCenters(dots,1);
                LeftXCenter = RandomisedLeftCenters(dots,2);

                RightYCenter = RandomisedRightCenters(dots,1);
                RightXCenter = RandomisedRightCenters(dots,2);
                if RandomisedIndices(dots) > nDotsAnnulus
                    LeftImage = paint_dot(LeftImage,LeftYCenter,LeftXCenter,LeftDot,LeftFill(dots),LeftPixels);
                    RightImage = paint_dot(RightImage,RightYCenter,RightXCenter,RightDot,RightFill(dots),RightPixels);
                else
                    LeftImage = paint_dot(LeftImage,LeftYCenter,LeftXCenter,AnnulusDot,LeftFill(dots),AnnulusPixels);
                    RightImage = paint_dot(RightImage,RightYCenter,RightXCenter,AnnulusDot,RightFill(dots),AnnulusPixels);
                end
            end            
        end
    end
    
end

function [LeftDot,RightDot] = new_dot(r)

    n = r*10;

    MotherDot = zeros(n, n);

    x = 1:n;
    y = 1:n;

    [X, Y] = meshgrid(x, y);

    cx = median(x);
    cy = median(y);

    d = sqrt(abs(X-cx).^2 + abs(Y-cy).^2);
    actual_dot = (d < r);
    anti_aliasing_index = logical((d < r*1.125).*(d>r));
    anti_aliasing = d(anti_aliasing_index) / max(d(anti_aliasing_index)) * 0.5;

    MotherDot(actual_dot) = 1;
    MotherDot(anti_aliasing_index) = anti_aliasing;
    
    LeftDot = MotherDot; RightDot = MotherDot;
    
end

function A = paint_dot(A,yc,xc,mydot,fill,pixels);
    % Usage: paintdot(A,yc,xc,mydot,fill)
    % A is your input matrix
    % yc and xc are your dot "centres"
    % mydot is the dot you want to paint
    % fill is any value you want to fill the dot with (e.g. 1 or -1)
    
    
    for px = 1:length(mydot);
        currentpx = pixels(px);
        ydc = yc+mydot(px,1)-1; xdc = xc+mydot(px,2)-1;
        A(ydc,xdc) = fill*currentpx;
    end
end
