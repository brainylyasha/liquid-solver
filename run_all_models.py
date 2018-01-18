import subprocess
import pandas as pd


oxid_model_list = ['Borisov-Shapkin', 'Kilinc', 'Sack']
oliv_model_list = ['Ford', 'Roeder-Emslie', 'Herzberg', 'Beattie', 'Putirka_AB', 'Putirka_CD', 'Toplis']

fugacity_list = ['NNO', 'QFM', 'WM', 'FeW']

print 'Would you like to use experimental fugacity (recorded in input data file) for calculations? (0 - No, 1 - Yes)'
exp_fug_flg = raw_input()

if int(exp_fug_flg) != 1:
    print 'Choose fugacity computation method ({}):'.format(', '.join([str(i) + ' - ' + method for i, method in enumerate(fugacity_list)]))

    fugacity_str = raw_input()
    if '+' in fugacity_str:
        fugacity_method = fugacity_list[int(fugacity_str.split('+')[0])]
        fugacity_delta = float(fugacity_str.split('+')[1])
    elif '-' in fugacity_str:
        fugacity_method = fugacity_list[int(fugacity_str.split('-')[0])]
        fugacity_delta = -float(fugacity_str.split('-')[1])
    else:
        fugacity_method = fugacity_list[int(fugacity_str)]
        fugacity_delta = 0

print 'Choose oxidation model ({}):'.\
    format(', '.join([str(i) + ' - ' + method for i, method in enumerate(oxid_model_list)]))
oxid_model = oxid_model_list[int(raw_input().strip())]

print 'Choose olivine-melt models (divide them with whitespace; {} or -1 for all of them):'.\
    format(', '.join([str(i) + ' - ' + method for i, method in enumerate(oliv_model_list)]))
user_oliv_model_list = map(int, raw_input().split())
if user_oliv_model_list != [-1]:
    oliv_model_list = [model for i, model in enumerate(oliv_model_list) if i in map(int, user_oliv_model_list)]
print oliv_model_list

if 'Toplis' in oliv_model_list:
    print 'Would you like to use experimental temperature (recorded in input data file) for calculations in Toplis&Sugawara olivine-melt model? (0 - No, 1 - Yes)'
    exp_T_flg = raw_input()


tech_fields = ['Temp', 'Mg_volume', 'Fe_volume', 'sum_volume_check', 'Mg_part', 'Fe_part']

for i, oliv_model in enumerate(oliv_model_list):
    cmd = 'python liquid_solver.py --oxid {} --oliv-melt {} --tol 0.00001 -o output_{}.csv'.\
        format(oxid_model, oliv_model, oliv_model)
    if oliv_model=='Toplis':
        cmd = cmd + ' --exp-temp {}'.format(exp_T_flg)
    if int(exp_fug_flg) != 1:
        cmd = cmd + ' -f {} --fugacity-delta {}'.format(fugacity_method, str(fugacity_delta))
    else:
        cmd = cmd + ' --exp-fug {}'.format(exp_fug_flg)
    print cmd
    subprocess.call(cmd, shell=True)

    cur_output = pd.read_csv('output_{}.csv'.format(oliv_model), sep=';')
    if i == 0:
        full_data = cur_output[[col for col in cur_output.columns if col not in tech_fields]].copy()
    full_data['T_{}'.format(oliv_model)] = cur_output['Temp']
    full_data['Olv_Fo {}'.format(oliv_model)] = cur_output['Mg_part']

full_data.to_csv('output_data.csv', sep=';', index=False)