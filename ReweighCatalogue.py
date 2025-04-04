import numpy as np
import pandas as pd
import warnings
from tqdm import tqdm
import glob
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
import gc

""" Python Script that matches isochrones to the population catalogue and weighs them according to the popolation mass."""

abbrev = 'SanityChecks/'+ sys.argv[1] #"NucSynthOnfall"

number_workers = 15

min_logL = 1 #None #log10(L/L_sun)


isochrone_path = './Resources/Isochrones/NewIsochroneChabrier/'



def find_closest_isochrone(row, padova, age_epsilon=1e-9, met_epsilon=1e-9):
    # Calculate age differences
    age_diff = (padova['logAge'] - row['LogAge']).abs()
    age_min = age_diff.min()


    # Keep all entries matching the minimal age difference (within a small epsilon if needed)
    age_subset = padova[age_diff <= age_min + age_epsilon]

    # Among those, find the minimal difference in metallicity
    met_diff = (age_subset['Zini'] - row['Metallicity']).abs()
    met_min = met_diff.min()

    # Keep all entries matching the minimal metallicity difference
    subset = age_subset[met_diff <= met_min + met_epsilon]

    # if (age_min < age_epsilon):
    if ((subset['Zini'].min()- subset['Zini'].max()) > 0 or (subset['Age'].min()- subset['Age'].max() > 0)):
        warnings.warn("Selected more than one isochrone")


    return subset


def weigh_isochrones(row, padova):
    
    isochrone = find_closest_isochrone(row, padova).copy()

    isochrone['weight'] = isochrone['IMF'] *row['PopulationMass']*1e9
    # isochrone['Nstars'] = isochrone['weight']/isochrone['Mini']

    isochrone['Radius'] = row['Radius']
    isochrone['BirthRadius'] = row['BirthRadius']

    isochrone['PopulationMass'] = row['PopulationMass']
    isochrone['Metallicity'] = row['Metallicity']
    isochrone['HeH'] = row['HeH']
    isochrone['ZH'] = row['ZH']
    isochrone['FeH'] = row['FeH']
    isochrone['OH'] = row['OH']
    isochrone['MgH'] = row['MgH']
    isochrone['CH'] = row['CH']
    isochrone['SiH'] = row['SiH']
    isochrone['CaH'] = row['CaH']
    isochrone['MnH'] = row['MnH']
    isochrone['CrH'] = row['CrH']
    isochrone['CoH'] = row['CoH']
    isochrone['EuH'] = row['EuH']
    isochrone['R2Age'] = row['TrueAge']
  

    return isochrone


def make_isochrone_catalogue(min_logL=None):
    cols =["Zini","MH","logAge","Mini","int_IMF","Mass","logL","logTe","logg","label","McoreTP","C_O","period0","period1","period2","period3","period4","pmode","Mloss","tau1m","X","Y","Xc","Xn","Xo","Cexcess","Z","mbolmag","Umag","Bmag","Vmag","Rmag","Imag","Jmag","Hmag","Kmag"]

    # Gather all isochorne .dat files
    files_newpadova = sorted(glob.glob(isochrone_path + '/*.dat'))

    # Load each file into a DataFrame and collect them in a list
    dfs = []
    for filename in files_newpadova:
        data = np.loadtxt(filename)
        df_temp = pd.DataFrame(data, columns=cols)
        dfs.append(df_temp)

    # Concatenate all DataFrames into one
    padova = pd.concat(dfs, ignore_index=True)

    padova = padova.sort_values(by=['logAge', 'Zini']).reset_index(drop=True)

    padova['Age'] = 10**padova['logAge'] / 1e9

    # get the IMF from the integrated IMF
    padova['IMF_diff'] = padova['int_IMF'].diff()
    # Take care of the first entry in a new isochrone
    padova['IMF'] = padova['IMF_diff'].where(padova['IMF_diff'] >= 0, padova['int_IMF'])

    # For the very first row, diff() yields NaN, so also set it to its own int_IMF:
    first_idx = padova.index[0]
    padova.loc[first_idx, 'IMF'] = padova.loc[first_idx, 'int_IMF']

    zero_mask = padova['IMF'] == 0

    ## padove isochrones have several points where the IMF has the same value. 
    # This happens at quickly evolving states, where we want to keep all points.
    # As our precision is bigger than isochrone one, we can distribute the IMF value over the zero points.
     
    # Identify contiguous blocks among the zeros using a vectorized group marker.
    group_ids = (zero_mask != zero_mask.shift(1)).cumsum()

    # For each contiguous group of zeros, look for the row immediately preceding the block,
    # then distribute its IMF value over (1 + number of rows in the zero block).
    for _, group in padova[zero_mask].groupby(group_ids[zero_mask]):
        first_zero_idx = group.index[0]
        # Only proceed if there is a preceding row and it is nonzero.
        if first_zero_idx == 0:
            warnings.warn("The very first entry had IMF = 0 - something went wrong.")
            continue
        if padova.at[first_zero_idx - 1, 'IMF'] == 0:
            warnings.warn("The value before the group was 0 too - something went wrong.")
            continue
        # Include the preceding row in the distribution.
        group_indices = [first_zero_idx - 1] + list(group.index)
        total_count = len(group_indices)
        preceding_value = padova.at[first_zero_idx - 1, 'IMF']
        new_value = preceding_value / total_count
        padova.loc[group_indices, 'IMF'] = new_value


        # padova.drop(padova[padova['IMF'] == 0].index, inplace=True)


    padova.drop([ 'McoreTP', 'C_O', 'period0', 'period1', 'period2',
        'period3', 'period4', 'pmode', 'Mloss', 'tau1m', 'X', 'Y', 'Xc', 'Xn',
        'Xo', 'Cexcess', 'IMF_diff'], inplace = True, axis = 1)

    #drop post AGB isochrone points
    padova.drop(padova[padova['label'] == 9].index, inplace=True)

    # drop all entries definitely not in the selction function
    if min_logL is not None:
        padova.drop(padova[padova['logL'] <min_logL].index, inplace=True)

    return padova




ddf = pd.read_csv('./Output/' +abbrev+ '/StellarCatalogue.dat', sep=', ', engine ='python')

ddf.drop(ddf.PopulationMass[ddf.PopulationMass == 0].index, inplace=True)

ddf.loc[ddf['TrueAge'] == 0, 'TrueAge'] = 0.001
ddf['LogAge']  = np.log10(ddf['TrueAge']*1e9)

padova = make_isochrone_catalogue(min_logL=min_logL)


## actually run the isochrone functions
rows = ddf.to_dict(orient='records')

def process_row(row_dict):
    row = pd.Series(row_dict)
    weighted_isochrones = weigh_isochrones(row, padova)
    return weighted_isochrones


with ProcessPoolExecutor(max_workers=number_workers) as executor:
    futures = [executor.submit(process_row, r) for r in rows]
    results = []
    for future in tqdm(as_completed(futures), total=len(futures), desc="Processing Populations"):
        results.append(future.result())


## writing the weighted isochrones to one file
for i, chunk in enumerate(tqdm(results, desc="Writing chunks to CSV", total=len(results))):
    df_chunk = chunk if isinstance(chunk, pd.DataFrame) else pd.DataFrame(chunk)
    if i == 0:
        df_chunk.to_csv("./Output/"+abbrev+"/WeightedPopulation.csv", mode="w", index=False, header=True)
    else:
        df_chunk.to_csv("./Output/"+abbrev+"/WeightedPopulation.csv", mode="a", index=False, header=False)
 

del results
del padova

gc.collect()