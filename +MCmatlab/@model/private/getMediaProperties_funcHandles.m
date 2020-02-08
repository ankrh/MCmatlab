function model = getMediaProperties_funcHandles(model,simType)

switch simType
  case 1 % Monte Carlo
    wavelength = model.MC.wavelength;
  case 2 % Fluorescence Monte Carlo
    wavelength = model.FMC.wavelength;
  case 3 % Heat Simulations
    wavelength = NaN;
end

G = model.G;

mP_raw = G.mediaPropertiesFunc(wavelength,G.mediaPropParams); % The raw media properties, as defined by the user in the model file
mP_fH(length(mP_raw)) = struct();

% Loop over all the media and throw error if the char arrays did not contain FR or T
anyFRdependence = false;
anyTdependence = false;
anyFDdependence = false;
for i=1:length(mP_raw)
  mP_fH(i).name = mP_raw(i).name;
  if simType <= 2
    if ~isfield(mP_raw,'mua') || isempty(mP_raw(i).mua)
      error('Error: Medium %s has no mua.',mP_raw(i).name);
    elseif ischar(mP_raw(i).mua)
      FRdependence = contains(mP_raw(i).mua,'FR');
      anyFRdependence = anyFRdependence || FRdependence;
      Tdependence = contains(mP_raw(i).mua,'T');
      anyTdependence = anyTdependence || Tdependence;
      FDdependence = contains(mP_raw(i).mua,'FD');
      anyFDdependence = anyFDdependence || FDdependence;
      if ~FRdependence && ~Tdependence && ~FDdependence
        error('Error: mua of %s is a char array but depends on neither FR (Fluence Rate), T (Temperature) nor FD (Fractional Damage)',mP_raw(i).name);
      end
      mP_fH(i).mua = str2func(['@(FR,T,FD)(' formulaToElementwiseFormula(mP_raw(i).mua) ')']);
    else
      mP_fH(i).mua = mP_raw(i).mua;
    end
    if ~isfield(mP_raw,'mus') || isempty(mP_raw(i).mus)
      error('Error: Medium %s has no mus.',mP_raw(i).name);
    elseif ischar(mP_raw(i).mus)
      FRdependence = contains(mP_raw(i).mus,'FR');
      anyFRdependence = anyFRdependence || FRdependence;
      Tdependence = contains(mP_raw(i).mus,'T');
      anyTdependence = anyTdependence || Tdependence;
      FDdependence = contains(mP_raw(i).mus,'FD');
      anyFDdependence = anyFDdependence || FDdependence;
      if ~FRdependence && ~Tdependence && ~FDdependence
        error('Error: mus of %s is a char array but depends on neither FR (Fluence Rate), T (Temperature) nor FD (Fractional Damage)',mP_raw(i).name);
      end
      mP_fH(i).mus = str2func(['@(FR,T,FD)(' formulaToElementwiseFormula(mP_raw(i).mus) ')']);
    else
      mP_fH(i).mus = mP_raw(i).mus;
    end
    if ~isfield(mP_raw,'g') || isempty(mP_raw(i).g)
      error('Error: Medium %s has no g.',mP_raw(i).name);
    elseif ischar(mP_raw(i).g)
      FRdependence = contains(mP_raw(i).g,'FR');
      anyFRdependence = anyFRdependence || FRdependence;
      Tdependence = contains(mP_raw(i).g,'T');
      anyTdependence = anyTdependence || Tdependence;
      FDdependence = contains(mP_raw(i).g,'FD');
      anyFDdependence = anyFDdependence || FDdependence;
      if ~FRdependence && ~Tdependence && ~FDdependence
        error('Error: g of %s is a char array but depends on neither FR (Fluence Rate), T (Temperature) nor FD (Fractional Damage)',mP_raw(i).name);
      end
      mP_fH(i).g = str2func(['@(FR,T,FD)(' formulaToElementwiseFormula(mP_raw(i).g) ')']);
    else
      mP_fH(i).g = mP_raw(i).g;
    end
    
    if ~isfield(mP_raw,'Y') || isempty(mP_raw(i).Y)
      mP_fH(i).Y = NaN;
    elseif ischar(mP_raw(i).Y)
      FRdependence = contains(mP_raw(i).Y,'FR');
      Tdependence = contains(mP_raw(i).Y,'T');
      FDdependence = contains(mP_raw(i).Y,'FD');
      if ~FRdependence && ~Tdependence && ~FDdependence
        error('Error: Y of %s is a char array but depends on neither FR (Fluence Rate), T (Temperature) nor FD (Fractional Damage)',mP_raw(i).name);
      end
      mP_fH(i).Y = str2func(['@(FR,T,FD)(' formulaToElementwiseFormula(mP_raw(i).Y) ')']);
    else
      mP_fH(i).Y = mP_raw(i).Y;
    end
    
    if ~isfield(mP_raw,'n') || isempty(mP_raw(i).n)
      mP_fH(i).n = NaN;
    else
      mP_fH(i).n = mP_raw(i).n;
    end
  else
    if ~isfield(mP_raw,'VHC') || isempty(mP_raw(i).VHC)
      error('Error: Medium %s has no VHC.',mP_raw(i).name);
    elseif ischar(mP_raw(i).VHC)
      FRdependence = contains(mP_raw(i).VHC,'FR');
      Tdependence = contains(mP_raw(i).VHC,'T');
      anyTdependence = anyTdependence || Tdependence;
      FDdependence = contains(mP_raw(i).VHC,'FD');
      anyFDdependence = anyFDdependence || FDdependence;
      if FRdependence
        error('Error: Thermal properties cannot depend on FR (Fluence Rate)');
      end
      if ~Tdependence && ~FDdependence
        error('Error: VHC of %s is a char array but depends on neither T (Temperature) nor FD (Fractional Damage)',mP_raw(i).name);
      end
      mP_fH(i).VHC = str2func(['@(T,FD)(' formulaToElementwiseFormula(mP_raw(i).VHC) ')']);
    else
      mP_fH(i).VHC = mP_raw(i).VHC;
    end
    if ~isfield(mP_raw,'TC') || isempty(mP_raw(i).TC)
      error('Error: Medium %s has no TC.',mP_raw(i).name);
    elseif ischar(mP_raw(i).TC)
      FRdependence = contains(mP_raw(i).TC,'FR');
      Tdependence = contains(mP_raw(i).TC,'T');
      anyTdependence = anyTdependence || Tdependence;
      FDdependence = contains(mP_raw(i).TC,'FD');
      anyFDdependence = anyFDdependence || FDdependence;
      if FRdependence
        error('Error: Thermal properties cannot depend on FR (Fluence Rate)');
      end
      if ~Tdependence && ~FDdependence
        error('Error: TC of %s is a char array but depends on neither T (Temperature) nor FD (Fractional Damage)',mP_raw(i).name);
      end
      mP_fH(i).TC = str2func(['@(T,FD)(' formulaToElementwiseFormula(mP_raw(i).TC) ')']);
    else
      mP_fH(i).TC = mP_raw(i).TC;
    end
    
    if ~isfield(mP_raw,'E') || isempty(mP_raw(i).E)
      mP_fH(i).E = NaN;
    else
      mP_fH(i).E = mP_raw(i).E;
    end
    if ~isfield(mP_raw,'A') || isempty(mP_raw(i).A)
      mP_fH(i).A = NaN;
    else
      mP_fH(i).A = mP_raw(i).A;
    end
  end
  
  if ~isfield(mP_raw,'nBins') || isempty(mP_raw(i).nBins) % nBins is number of bins used for fluence rate, temperature or damaged fraction dependences
    mP_fH(i).nBins = NaN;
  else
    mP_fH(i).nBins = mP_raw(i).nBins;
  end
end

switch simType
  case 1
    model.MC.mediaProperties_funcHandles = mP_fH;
    model.MC.FRdependent = anyFRdependence;
    model.MC.Tdependent  = anyTdependence;
    model.MC.FDdependent = anyFDdependence;
  case 2
    model.FMC.mediaProperties_funcHandles = mP_fH;
    model.FMC.FRdependent = anyFRdependence;
    model.FMC.Tdependent  = anyTdependence;
    model.FMC.FDdependent = anyFDdependence;
  case 3
    model.HS.mediaProperties_funcHandles = mP_fH;
    model.HS.Tdependent  = anyTdependence;
    model.HS.FDdependent = anyFDdependence;
end
end

function str = formulaToElementwiseFormula(str)
str = strrep(str,'/','./');
str = strrep(str,'*','.*');
str = strrep(str,'^','.^');
str = strrep(str,'..','.');
end