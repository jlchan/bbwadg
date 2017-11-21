%% Surface Fitting With Custom Equations to Biopharmaceutical Data
%
% This example shows how to use Curve Fitting Toolbox(TM) to fit response
% surfaces to some anesthesia data to analyze drug interaction effects. Response
% surface models provide a good method for understanding the pharmacodynamic
% interaction behavior of drug combinations.
%
% This data is based on the results in this paper: Kern SE, Xie G, White JL, Egan
% TD. Opioid-hypnotic synergy: A response surface analysis of
% propofol-remifentanil pharmacodynamic interaction in volunteers. Anesthesiology
% 2004; 100: 1373-81.
%
% Anesthesia is typically at least a two-drug process, consisting of an opioid
% and a sedative hypnotic. This example uses Propofol and Reminfentanil as drug
% class prototypes. Their interaction is measured by four different measures of
% the analgesic and sedative response to the drug combination. Algometry, Tetany,
% Sedation, and Laryingoscopy comprise the four measures of surrogate drug
% effects at various concentration combinations of Propofol and Reminfentanil.
%
% The following code, using Curve Fitting Toolbox methods, reproduces the
% interactive surface building with the Curve Fitting Tool described in
% "Biopharmaceutical Drug Interaction Surface Fitting".

% Copyright 2012 The MathWorks, Inc.

%% Load Data
% Load the data from file.
data = importdata( 'OpioidHypnoticSynergy.txt' );
Propofol      = data.data(:,1);
Remifentanil  = data.data(:,2);
Algometry     = data.data(:,3);
Tetany        = data.data(:,4);
Sedation      = data.data(:,5);
Laryingoscopy = data.data(:,6);

%% Create the Model Fit Type
% You can use the |fittype| function to define the model from the paper, where
% |CA| and |CB| are the drug concentrations, and |IC50A|, |IC50B|, |alpha|, and
% |n| are the coefficients to be estimated. Create the model fit type.
ft = fittype( 'Emax*( CA/IC50A + CB/IC50B + alpha*( CA/IC50A ) * ( CB/IC50B ) )^n /(( CA/IC50A + CB/IC50B + alpha*( CA/IC50A ) * ( CB/IC50B ) )^n  + 1 )', ...
    'independent', {'CA', 'CB'}, 'dependent', 'z', 'problem', 'Emax' )
%%
% Assume |Emax = 1| because the effect output is normalized.
Emax = 1;

%% Set Fit Options
% Set fit options for robust fitting, bounds, and start points.
opts = fitoptions( ft );
opts.Lower = [0, 0, -5, -0];
opts.Robust = 'LAR';
opts.StartPoint = [0.0089, 0.706, 1.0, 0.746];

%% Fit and Plot a Surface for Algometry
[f, gof] = fit( [Propofol, Remifentanil], Algometry, ft,...
    opts, 'problem', Emax )
plot( f, [Propofol, Remifentanil], Algometry );


%% Fit a Surface to Tetany
% Reuse the same |fittype| to create a response surface for tetany.
[f, gof] = fit( [Propofol, Remifentanil], Tetany, ft, opts, 'problem', Emax )
plot( f, [Propofol, Remifentanil], Tetany );


%% Fit a Surface to Sedation
[f, gof] = fit( [Propofol, Remifentanil], Sedation, ft, opts, 'problem', Emax )
plot( f, [Propofol, Remifentanil], Sedation );


%% Fit a Surface to Laryingoscopy
[f, gof] = fit( [Propofol, Remifentanil], Laryingoscopy, ft, opts, 'problem', Emax )
plot( f, [Propofol, Remifentanil], Laryingoscopy );

displayEndOfDemoMessage(mfilename)
