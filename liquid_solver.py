#!/usr/local/bin/python
# -*- coding: utf-8 -*-

import sys
import getopt
import numpy as np
import pandas as pd
from scipy.optimize import brentq

verbose = False
input_filename = 'input_data_petrolog.csv'
output_filename = 'output_data.csv'
eps = 0.1
fugacity_delta = 0
exp_T_flg = False
exp_fug_flg=False

options, remainder = getopt.getopt(sys.argv[1:], 'i:o:vt:f:', ['verbose', 'oxid=', 'oliv-melt=', 'tol=', 'fugacity=', 'fugacity-delta=', 'exp-temp=', 'exp-fug='])

for opt, arg in options:
    if opt in ('-o', '--output'):
        output_filename = arg
    elif opt in ('-i', '--input'):
        input_filename = arg
    elif opt in ('-v', '--verbose'):
        verbose = True
    elif opt == '--oxid':
        oxidation_model = arg
    elif opt == '--oliv-melt':
        olivine_melt_model = arg
#        oxid_list = oxid_list_dict[olivine_melt_model]
    elif opt in ('-t', '--tol'):
        eps = float(arg)
    elif opt in ('-f', '--fugacity'):
        fugacity_method = arg
    elif opt in ('--fugacity-delta'):
        fugacity_delta = float(arg)
    elif opt == '--exp-temp':
        exp_T_flg = (int(arg) == 1)
    elif opt == '--exp-fug':
        exp_fug_flg = (int(arg) == 1)

def manipulate_liquid(row, oxid_w, oxidation_model, olivine_melt_model, full_flag=False):
    # leaving required oxids
#    oxid_list = oxid_list_dict[olivine_melt_model][:]
    oxid_list = oxid_w.index.values.tolist()
    if not full_flag:
        oxid_list.remove('Fe2O3')

    row = row[oxid_list]

    # dividing by oxid weights
    row = row / np.array(oxid_w[oxid_list])

    # some must be halfed
    if ((not full_flag) and oxidation_model == 'Borisov-Shapkin') or \
            (full_flag and olivine_melt_model in('Roeder-Emslie', 'Toplis')):
        two_kat_list = ['Li2O', 'Na2O', 'Al2O3', 'P2O5', 'K2O', 'Sc2O3', 'V2O5', 'Cr2O3', 'As2O5', 'Nb2O3',
                            'Sb2O3', 'Fe2O3']
        two_kat_filtered_list = list(set(two_kat_list) & set(row.keys().tolist()))
        row[two_kat_filtered_list] = row[two_kat_filtered_list] / 2.0
    # normalizing to sum = 1 again
    row = row / np.array(row).sum()
    return row

def compute_fugacity(T_K, fugacity_delta):
    fug_K1_list = {'NNO': -57403.5, 'QFM': -56280.0, 'WM': -85083.0, 'FeW': -61790.0}
    fug_K2_list = {'NNO': 21.55, 'QFM': 19.09, 'WM': 37.04, 'FeW': 14.9}

    if fugacity_method in fug_K1_list.keys():
        fugacity = (fug_K1_list[fugacity_method] / prev_T_K + fug_K2_list[fugacity_method]) * np.log10(np.e)
    else:
        sys.exit('Error: fugacity computation method {} not supported'.format(fugacity_method))
    return fugacity + fugacity_delta


# reading the input liquid data
input_data = pd.read_csv(input_filename, sep = ';')

# column names can be unclean — cleaning them (leave only letters and numbers)
input_data.columns = [''.join([letter for letter in col if letter.isalpha() or letter.isdigit() or letter == '_'])
                      for col in input_data.columns.tolist()]

if exp_T_flg:
    exp_T_list = (input_data['exp_T_C'] + 273.15).tolist()
if exp_fug_flg:
    exp_fug_list = input_data['log_fO2'].tolist()

if 'exp_T_C' in input_data.columns:
    input_data = input_data.drop('exp_T_C', axis=1)
if 'log_fO2' in input_data.columns:
    input_data = input_data.drop('log_fO2', axis=1)

# normalizing to sum = 100
input_data = pd.DataFrame(np.array(input_data) * 100.0 / np.array(input_data).sum(axis=1).reshape(-1, 1), columns=input_data.columns)

#also reading the oxid weights
oxid_w = pd.read_csv('elements_weights.csv', sep = ';').iloc[0, :]
oxid_list = oxid_w.index.values.tolist()
oxid_list.remove('Fe2O3')

for col in set(oxid_list) - set(input_data.columns.tolist()):
    input_data[col] = 0
input_data = input_data[oxid_list]
#oxid_w = oxid_w[oxid_list]

molec_data = pd.DataFrame(np.hstack(
    (input_data, input_data.apply(lambda x: manipulate_liquid(x, oxid_w, oxidation_model, olivine_melt_model), axis=1))),
                                                columns = input_data.columns.tolist() + \
                                      [col+'_molec_proc' for col in input_data.columns.tolist()])

# table for final results
final_df = input_data.copy()
opt_T_list = []
Mg_volume_list = []
Fe_volume_list = []
for i, row in molec_data.iterrows():
    # repeating 2 steps until convergence or current liquid
    cur_T_K = 1520
    prev_T_K = 2700
    if verbose:
        print 'liquid {}'.format(str(i))

    iter = 0
    while np.abs(cur_T_K - prev_T_K) > eps:
        if olivine_melt_model == 'Toplis':
            prev_T_K = cur_T_K
            full_liquid_info = manipulate_liquid(row, oxid_w, oxidation_model, olivine_melt_model)
            if exp_T_flg:
                cur_T_K = exp_T_list[i]
            else:
                cur_T_K = 1446 - 144 * full_liquid_info['SiO2'] - 50 * full_liquid_info['FeO'] + 1232 * full_liquid_info['MgO'] - 389.9 * full_liquid_info['CaO']
            if exp_fug_flg:
                fugacity = exp_fug_list[i]
            else:
                fugacity = compute_fugacity(cur_T_K, fugacity_delta)
            Fe_frac = np.exp(0.2185 * fugacity * np.log(10) + 12670.0 / cur_T_K - 7.54 - 2.24 * full_liquid_info['Al2O3'] + 1.55 * full_liquid_info['FeO'] +
                             2.96 * full_liquid_info['CaO'] + 8.42 * full_liquid_info['Na2O'] + 9.59 * full_liquid_info['K2O'])
            old_FeO = full_liquid_info['FeO']
            full_liquid_info['FeO'] = old_FeO / (1.0 + 2.0 * Fe_frac)
            full_liquid_info['Fe2O3'] = old_FeO - full_liquid_info['FeO']
            full_liquid_info = full_liquid_info * 100
            if full_liquid_info['SiO2'] < 60:
                coef = (46.0 / (100.0 - full_liquid_info['SiO2']) - 0.93) * (full_liquid_info['Na2O'] + full_liquid_info['K2O']) + (-533.0 / (100.0 - full_liquid_info['SiO2']) + 9.69)
                full_liquid_info['SiO2'] = full_liquid_info['SiO2'] + coef * (full_liquid_info['Na2O'] + full_liquid_info['K2O'])
            else:
                full_liquid_info['SiO2'] = full_liquid_info['SiO2'] + (5.5 * np.exp(-0.13 * (full_liquid_info['Na2O'] + full_liquid_info['K2O'])) * (2.0 - 100.0 / (100.0 - full_liquid_info['SiO2'])) - 1) * \
                                                                      (full_liquid_info['Na2O'] + full_liquid_info['K2O'])
            Fe_volume = brentq(lambda x: (0.036 * full_liquid_info['SiO2'] - 0.22) * np.exp(-(6766 + 7.34 * cur_T_K) / (8.314 * cur_T_K) + 3000.0 * (3 * x - 1) / (8.314 * cur_T_K)) - \
                               full_liquid_info['MgO'] * x / (full_liquid_info['FeO'] * (2.0 / 3 - x)), 0.0001, 2.0/3 - 0.0001, xtol=1e-16)
            Mg_volume = 2.0 / 3 - Fe_volume
        else:
            if verbose:
                print 'iter {}'.format(iter)
            iter += 1
            # ferrum quantity computation
            prev_T_K = cur_T_K
            if verbose:
                print 'current temperature is {}'.format(str(prev_T_K))
            if exp_fug_flg:
                fugacity = exp_fug_list[i]
            else:
                fugacity = compute_fugacity(cur_T_K, fugacity_delta)
            if verbose:
                print 'fugacity = {}'.format(str(fugacity))

            if oxidation_model == 'Borisov-Shapkin':
                a = 0.1735
                b = 771.7
                c = 1.914
                oxid_const = pd.read_csv('ferrum_const.csv', sep = ';').T
                oxid_const.columns = ['d', 'd1', 'd2']
                oxid_const['D'] = oxid_const['d'] + oxid_const['d1'] / prev_T_K + fugacity * oxid_const['d2']
                oxid_const['D*X'] = np.array(oxid_const['D']) * np.array(row[[col + '_molec_proc' for col in
                                                          ['SiO2', 'TiO2', 'Al2O3', 'FeO', 'MgO', 'CaO', 'Na2O',  'K2O']]] * 100.0)
                ferrum_frac= np.power(10.0, np.sum(oxid_const['D*X']) + a * fugacity + b / prev_T_K + c)
                print ferrum_frac
                if verbose:
                    print 'Fe_3+ / Fe_2+ = {}'.format(str(ferrum_frac))
                full_liquid_info = pd.Series({oxid: row[oxid + '_molec_proc'] for oxid in oxid_list})
                full_liquid_info['FeO'] = full_liquid_info['FeO'] / (1.0 + ferrum_frac)
                full_liquid_info['Fe2O3'] = full_liquid_info['FeO'] * ferrum_frac / 2
                full_liquid_info = full_liquid_info / np.sum(full_liquid_info)
                print 'final fractions'
                print full_liquid_info
                if verbose:
                    print 'final liquid fractions:'
                    print full_liquid_info
            else:
                sys.exit('Error: oxidation model {} not supported'.format(oxidation_model))

            # optimal temperature computation
            if olivine_melt_model == 'Ford':
                Fe_C_list = [-4.325, 6515.1, 0.035572, -0.94621, -0.19957, 1.4172, 0.64833, 0.18455, 0.51707, 3.2151, 1.1978]
                Mg_C_list = [-2.2008, 4896.7, 0.013804, -1.4993, -1.2674, 2.2394, 4.3926, -1.5124, 0.2268, 2.4709, 4.5012]
                # exp(Fe_K1 + Fe_K2 / t) + exp(Mg_K1 + Mg_K2 / T) = 2/3
                # computing consts

                Fe_K1 = np.log(full_liquid_info['FeO']) + Fe_C_list[0] + Fe_C_list[3] * \
                np.log(1.5 * (full_liquid_info['MgO'] + full_liquid_info['FeO'] + full_liquid_info['CaO'] + \
                             full_liquid_info['MnO'] + full_liquid_info['Cr2O3'])) +\
                Fe_C_list[4] * np.log(3 * full_liquid_info['SiO2'])  + Fe_C_list[5] * np.log(1 - full_liquid_info['Al2O3']) + \
                Fe_C_list[6] * np.log(1 - full_liquid_info['Fe2O3']) + Fe_C_list[7] * np.log(1 - full_liquid_info['CaO']) + \
                Fe_C_list[8] * np.log(1 - full_liquid_info['Na2O'] - full_liquid_info['K2O']) +\
                Fe_C_list[9] * np.log(1 - full_liquid_info['TiO2']) + Fe_C_list[10] * np.log(1 - full_liquid_info['P2O5'])
                Fe_K2 = Fe_C_list[1]

                Mg_K1 = np.log(full_liquid_info['MgO']) + Mg_C_list[0] + Mg_C_list[3] * \
                np.log(1.5 * (full_liquid_info['MgO'] + full_liquid_info['FeO'] + full_liquid_info['CaO'] + \
                             full_liquid_info['MnO'] + full_liquid_info['Cr2O3'])) +\
                Mg_C_list[4] * np.log(3 * full_liquid_info['SiO2'])  + Mg_C_list[5] * np.log(1 - full_liquid_info['Al2O3']) + \
                Mg_C_list[6] * np.log(1 - full_liquid_info['Fe2O3']) + Mg_C_list[7] * np.log(1 - full_liquid_info['CaO']) + \
                Mg_C_list[8] * np.log(1 - full_liquid_info['Na2O'] - full_liquid_info['K2O']) +\
                Mg_C_list[9] * np.log(1 - full_liquid_info['TiO2']) + Mg_C_list[10] * np.log(1 - full_liquid_info['P2O5'])
                Mg_K2 = Mg_C_list[1]
                cur_T_K = brentq(lambda x: np.exp(Fe_K1 + Fe_K2 / x) +
                    np.exp(Mg_K1 + Mg_K2 / x) - 2.0 / 3,
                                1, 3000, xtol=1e-16)
                if verbose:
                    print 'new optimal T = {}'.format(str(cur_T_K))
            elif olivine_melt_model == 'Roeder-Emslie':
                Mg_A, Mg_B = 3740.0, 1.87
                Fe_A, Fe_B = 3911.0, 2.5
                cur_T_K = brentq(lambda x: full_liquid_info['MgO'] * np.power(10.0, Mg_A / x - Mg_B) +
                                            full_liquid_info['FeO'] * np.power(10.0, Fe_A / x - Fe_B) - 2.0 /3 , 1, 3000, xtol=1e-16)
            elif olivine_melt_model == 'Herzberg':
                A, B, C = 0.381, 0.79, 1.039
                Fe_volume = 2.0 * (A * full_liquid_info['FeO'] - B * full_liquid_info['FeO'] / row['MgO'] + \
                             C * full_liquid_info['FeO'] / (row['MgO'] * row['MgO'])) / \
                            (3.0 * (full_liquid_info['MgO'] + A * full_liquid_info['FeO'] - \
                                    B * full_liquid_info['FeO'] / row['MgO'] + C * full_liquid_info['FeO'] / (row['MgO'] * row['MgO'])))
                Mg_volume = 2.0 / 3 - Fe_volume
                H, S, R = 113100, 52.05, 8.314
                NF = 7.0 * (np.log(1 - full_liquid_info['Al2O3']) / 2.0 + np.log(1 - full_liquid_info['TiO2']))
                cur_T_K = (H / R) /(S / R + 2 * np.log(Mg_volume / full_liquid_info['MgO']) + \
                                    2 * np.log(3.0 * (full_liquid_info['FeO'] + full_liquid_info['MgO'] + \
                                                      full_liquid_info['MnO'] + full_liquid_info['CaO']) / 2) +\
                                    2 * np.log(3.0 * full_liquid_info['SiO2']) - NF)
            elif olivine_melt_model == 'Beattie':
                A, B = 0.299, 0.027
                Fe_volume = (2.0 * A * full_liquid_info['FeO'] / 3.0 + B * full_liquid_info['MgO'] * full_liquid_info['FeO']) / \
                            (full_liquid_info['MgO'] + A * full_liquid_info['FeO'])
                Mg_volume = 2.0 / 3 - Fe_volume
                H, S, R = 113100, 52.05, 8.314
                NF = 7.0 * (np.log(1 - full_liquid_info['Al2O3']) / 2.0 + np.log(1 - full_liquid_info['TiO2']))
                cur_T_K = (H / R) /(S / R + 2 * np.log(Mg_volume / full_liquid_info['MgO']) + \
                                    2 * np.log(3.0 * (full_liquid_info['FeO'] + full_liquid_info['MgO'] + \
                                                      full_liquid_info['MnO'] + full_liquid_info['CaO']) / 2) +\
                                    2 * np.log(3.0 * full_liquid_info['SiO2']) - NF)
            elif olivine_melt_model == 'Putirka_AB':
                A_Mg, B_Mg = 4490.5, 2.02
                A_Fe, B_Fe = 3793.3, 2.66
                cur_T_K = 273.15 + brentq(lambda x: full_liquid_info['MgO'] * np.exp(A_Mg / x - B_Mg) +
                                                    full_liquid_info['FeO'] * np.exp(A_Fe / x - B_Fe) - 2.0 / 3,
                                          1, 3000, xtol=1e-16)
                Mg_volume = full_liquid_info['MgO'] * np.exp(A_Mg / (cur_T_K - 273.15) - B_Mg)
                Fe_volume = full_liquid_info['FeO'] * np.exp(A_Fe / (cur_T_K - 273.15) - B_Fe)
            elif olivine_melt_model == 'Putirka_CD':
                A_Mg, A_Fe = 3063.2, 2556.4
                B_Mg, B_Fe = 2.106193, 3.25
                Si_Mg, Si_Fe = 0.019, 0.028
                NaK_Mg, NaK_Fe = 0.08, 0.052
                cur_T_K = 273.15 + brentq(lambda x: full_liquid_info['MgO'] * np.exp(A_Mg / x - B_Mg + Si_Mg * row['SiO2'] +
                                                                                     NaK_Mg * (row['Na2O'] + row['K2O'])) +
                                                    full_liquid_info['FeO'] * np.exp(A_Fe / x - B_Fe + Si_Fe * row['SiO2'] +
                                                                                     NaK_Fe * (row['Na2O'] + row['K2O'])) - 2.0 / 3,
                                          1, 3000, xtol=1e-16)
                Mg_volume = full_liquid_info['MgO'] * np.exp(A_Mg / (cur_T_K - 273.15) - B_Mg + Si_Mg * row['SiO2'] +
                                                                                     NaK_Mg * (row['Na2O'] + row['K2O']))
                Fe_volume = full_liquid_info['FeO'] * np.exp(A_Fe / (cur_T_K - 273.15) - B_Fe + Si_Fe * row['SiO2'] +                                                                                 NaK_Fe * (row['Na2O'] + row['K2O']))
            else:
                sys.exit('Error: olivine-melt model {} not supported'.format(olivine_melt_model))
    # writing result for this liquid to table
    opt_T_list.append(cur_T_K - 273.15)
    if olivine_melt_model == 'Ford':
        Mg_volume_list.append(np.exp(Mg_K1 + Mg_K2 / cur_T_K))
        Fe_volume_list.append(np.exp(Fe_K1 + Fe_K2 / cur_T_K)) # 2.0 / 3 - full_liquid_info['Mg_volume']
    elif olivine_melt_model == 'Roeder-Emslie':
    	Mg_volume_list.append(full_liquid_info['MgO'] * np.power(10.0, Mg_A / cur_T_K - Mg_B))
    	Fe_volume_list.append(full_liquid_info['FeO'] * np.power(10.0, Fe_A / cur_T_K - Fe_B))
    elif olivine_melt_model in ('Herzberg', 'Beattie', 'Toplis', 'Putirka_AB', 'Putirka_CD'):
        Mg_volume_list.append(Mg_volume)
        Fe_volume_list.append(Fe_volume)
final_df['Temp'] = opt_T_list
final_df['Mg_volume'] = Mg_volume_list
final_df['Fe_volume'] = Fe_volume_list
final_df['sum_volume_check'] = final_df['Mg_volume'] + final_df['Fe_volume']
olivine_sum = final_df['Mg_volume'] + final_df['Fe_volume']
final_df['Mg_part'] = final_df['Mg_volume'] * 100.0 / olivine_sum
final_df['Fe_part'] = final_df['Fe_volume'] * 100.0 / olivine_sum
final_df.to_csv(output_filename, index=False, sep=';')