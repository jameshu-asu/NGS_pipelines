import os
import sys
import shutil
import subprocess
import pandas as pd

#______INPUTS______#

'''By James C. Hu
This section stores input variables for the run.
Variables are pulled from the RVP_main_input.csv file.
sequencer input field in the input file is case sensitive. Must be "INSTRUMENT_NAME1" or "INSTRUMENT_NAME2" to reflect sequencer directory path.
'''

infile = pd.read_csv('RVP_main_input.csv').copy()
run = infile['run'].to_list()[0]  # alternatively, variables can be set using the following format: run = infile.at[0, 'run'] which would be faster.
sequencer = infile['sequencer'].to_list()[0]
instrament_run_number = infile['instrament_run_number'].to_list()[0]


#____STATIC_VARIABLES____#
'''Notes
The list and dictionaries in this section are used to modify dataframes.
add_col_list is used in the add_col funciton.
virus_dict is used in the sort_groups funciton.
norm_counts_groups is used in the sum_groups function.
'''

# for adding indexing columns
add_col_list = ['Site', 'Date', 'Sample_ID']

# for formatting.
# This dictionary replaces spacs in virus names with underscores
virus_dict = {'Human adenovirus B1': 'Human_adenovirus_B1', 'Human adenovirus C': 'Human_adenovirus_C', 'Human adenovirus E': 'Human_adenovirus_E',
              'Human bocavirus 2c PK isolate PK-5510': 'Human_bocavirus_2c_PK_isolate_PK-5510', 'Human bocavirus 3': 'Human_bocavirus_3', 'Human bocavirus 4 NI strain HBoV4-NI-385': 'Human_bocavirus_4_NI_strain_HBoV4-NI-385', 'Primate bocaparvovirus 1 isolate st2': 'Primate_bocaparvovirus_1_isolate_st2',
              'Human coronavirus 229E': 'Human_coronavirus_229E', 'Human coronavirus HKU1': 'Human_coronavirus_HKU1', 'Human Coronavirus NL63': 'Human_Coronavirus_NL63', 'Human coronavirus OC43': 'Human_coronavirus_OC43',
              'Enterovirus C104': 'Enterovirus_C104', 'Human enterovirus 109': 'Human_enterovirus_109',
              'Human metapneumovirus': 'Human_metapneumovirus',
              'Human parainfluenza virus 1': 'Human_parainfluenza_virus_1', 'Human parainfluenza virus 3': 'Human_parainfluenza_virus_3', 'Human parainfluenza virus 4a': 'Human_parainfluenza_virus_4a',
              'Human parechovirus 6': 'Human_parechovirus_6', 'Human parechovirus type 1 PicoBank/HPeV1/a': 'Human_parechovirus_type_1_PicoBank/HPeV1/a',
              'Human Respiratory syncytial virus 9320': 'Human_Respiratory_syncytial_virus_9320', 'Respiratory syncytial virus': 'Respiratory_syncytial_virus',
              'Human rhinovirus 89': 'Human_rhinovirus_89', 'Human rhinovirus C': 'Human_rhinovirus_C', 'Rhinovirus B14': 'Rhinovirus_B14',
              'Human rubulavirus 2': 'Human_rubulavirus_2',
              'Influenza A': 'Influenza_A', 'Influenza B': 'Influenza_B',
              'KI polyomavirus Stockholm 60': 'KI_polyomavirus_Stockholm_60', 'WU Polyomavirus': 'WU_Polyomavirus',
              'SARS-CoV-2': 'SARS-CoV-2'
              }

# for sum by groups
# this dictionary is used to group viruses.
norm_counts_groups = {'Adenovirus': ['Human_adenovirus_B1', 'Human_adenovirus_C', 'Human_adenovirus_E'],
                      'Bocavirus': ['Human_bocavirus_2c_PK_isolate_PK-5510', 'Human_bocavirus_3', 'Human_bocavirus_4_NI_strain_HBoV4-NI-385', 'Primate_bocaparvovirus_1_isolate_st2'],
                      'Coronavirus': ['Human_coronavirus_229E', 'Human_coronavirus_HKU1', 'Human_Coronavirus_NL63', 'Human_coronavirus_OC43'],
                      'Enterovirus': ['Enterovirus_C104', 'Human_enterovirus_109'],
                      'Metapneumovirus': ['Human_metapneumovirus'],
                      'Parainfluenza': ['Human_parainfluenza_virus_1', 'Human_parainfluenza_virus_3', 'Human_parainfluenza_virus_4a'],
                      'Parechovirus': ['Human_parechovirus_6', 'Human_parechovirus_type_1_PicoBank/HPeV1/a'],
                      'RSV': ['Human_Respiratory_syncytial_virus_9320', 'Respiratory_syncytial_virus'],
                      'Rhinovirus': ['Human_rhinovirus_89', 'Human_rhinovirus_C', 'Rhinovirus_B14'],
                      'Rubulavirus': ['Human_rubulavirus_2'],
                      'Influenza': ['Influenza_A', 'Influenza_B'],
                      'Polyomavirus': ['KI_polyomavirus_Stockholm_60', 'WU_Polyomavirus'],
                      'SARS-CoV-2': ['SARS-CoV-2']
                      }


#______FUNCTIONS______#
'''Notes
This section contains function for running the RVP_main_v2.py script.
The main scrip will do the following:
1) Set up directory structure for specified run based on input file.
    a) A log file will be generated to track the workflow, both success and errors are tracked.
    b) Script will stop executing after encoutering an error using sys.exit() and an error message will be printed to the terminal.
2) Execute the RVP monitoring script and all other associated scripts.
3) Perform R analysis/normaization/data merging to database files in the metadata directory.
    a) Most common errors in this section is caused by delimiter missmatch for .txt files when performing R analysis.
       Try adjusting the delimiter (sep="," to sep="\t", sep="," is default : line 17 of the R script) as needed to fix R analysis related errors.
4) Generate outputs for Streamlit dashboard in the metadata directory.
'''


def make_RVP_directory(run: str) -> None:
    '''
    This function generates the directory structure for RVP analysis
    1) Creates directory structure.
    2) Populates directory with necessary scripts.
    3) Runs monitoring script.
    '''
    os.chdir('../')
    os.makedirs(f'NextSeq1000_{run}', mode=0o777)
    shutil.copy('RVP_scripts/wastewater_respiratory_Monitor_v3.py',
                f'NextSeq1000_{run}')
    shutil.copy('RVP_scripts/ReadsInOut_2.4.py', f'NextSeq1000_{run}')
    shutil.copy('RVP_scripts/QCWorkflow_WW_Agave_v3.2.py', f'NextSeq1000_{run}')
    shutil.copy('../metadata/ABCTL Samples.xlsx', f'../respiratory_virus/NextSeq1000_{run}')
    return None


def get_fastq_path(sequencer: str, instrament_run_number: str) -> str:
    '''
    This function parses the sequencer directory for the run specific file name.
    Sequencer: Ganymede or Callisto, case sensitive
    Instrament_run: instrament run number as a string.
    '''
    path = f'/mnt/storage/data/{sequencer}/outputDirectory'
    sequencer_path = ''
    os.chdir(path)
    instrament_files = [file.split('_') for file in os.listdir() if file.startswith('23')]
    for file in instrament_files:
        if instrament_run_number == file[2]:
            sequencer_path = '_'.join(file)
    return os.path.join(path, sequencer_path)


def start_RVP_monitoring(run: str, sequencer: str, sequencer_path: str) -> None:
    '''
    This funciton cds into the run directory and initiates the monitoring script on the provided path.
    '''
    os.chdir(f'/mnt/storage/wastewater_analysis/respiratory_virus/NextSeq1000_{run}')
    subprocess.call(['python', 'wastewater_respiratory_Monitor_v3.py', f'{sequencer_path}'])
    return None


def setup_r_analysis_directory(run: str) -> None:
    '''
    This function populates the run specific r_analysis folder with relevant files.
    '''
    os.chdir('/mnt/storage/wastewater_analysis/metadata')
    # make sure to move back after updates.
    shutil.copy(f'RVP_metadata.txt', f'../respiratory_virus/NextSeq1000_{run}/r_analysis')
    # make sure to move back after updates.
    shutil.copy(f'RVP_counts_virusOnly.txt', f'../respiratory_virus/NextSeq1000_{run}/r_analysis')
    os.chdir('/mnt/storage/wastewater_analysis/respiratory_virus/RVP_scripts')
    shutil.copy('AccessionToName_RVP.txt', f'../NextSeq1000_{run}/r_analysis')
    shutil.copy('RVP_analysis_v0.2.R', f'../NextSeq1000_{run}/r_analysis')
    shutil.copy(f'../NextSeq1000_{run}/read_metrics.txt', f'../NextSeq1000_{run}/r_analysis')
    shutil.copy(f'../NextSeq1000_{run}/counts_2022/countsTable.txt',
                f'../NextSeq1000_{run}/r_analysis')
    return None


def add_QC_metrics(run: str) -> None:
    '''
    This function updates the the RVP_metadata.txt file with QC read metrics.
    '''
    os.chdir(f'/mnt/storage/wastewater_analysis/respiratory_virus/NextSeq1000_{run}/r_analysis')
    df = pd.read_csv('countsTable.txt', sep='\t', index_col=0)
    df.columns = df.columns.str.strip('.txt')
    df.to_csv('countsTable.txt', sep='\t')
    df = pd.read_csv('read_metrics.txt', sep='\t', index_col=False)
    df1 = df.copy()
    df1 = df1.rename(columns={'Sample ID': 'Respiratory Seq ID',
                              'DeDupe_Merged_Dedupe_Filtered_Output': 'QC'}).set_index('Respiratory Seq ID')
    df = pd.read_csv('RVP_metadata.txt', sep='\t')
    df2 = df.copy()
    df2 = df2.set_index('Respiratory Seq ID')
    df2.update(df1)
    df2.to_csv('RVP_metadata.txt', sep='\t')
    return None


def run_RVP_r_script(run: str) -> None:
    '''
    This function executes the R script and moves the output files to the metadata directory.
    '''
    os.chdir(f'/mnt/storage/wastewater_analysis/respiratory_virus/NextSeq1000_{run}/r_analysis')
    subprocess.call(['Rscript', 'RVP_analysis_v0.2.R'])
    shutil.copy('RVP_metadata.txt', '../../../metadata/RVP_metadata.txt')
    shutil.copy('RVP_counts_virusOnly_updated.txt', '../../../metadata/RVP_counts_virusOnly.txt')
    # shutil needs full path to destination to overwrite existing file.
    shutil.copy('RVP_counts_virusOnly_normalized.txt',
                '/mnt/storage/wastewater_analysis/metadata/RVP_counts_virusOnly_normalized.txt')
    return None


def sort_groups(infile: str) -> pd.DataFrame:
    '''
    This funciton formats the infile using virus_dict.
    '''
    df = pd.read_csv(infile)
    df = df.rename(columns={'Unnamed: 0': 'Virus'})
    df = df.set_index('Virus')
    df = df.T
    df = df.rename(columns=virus_dict)
    df = df[virus_dict.values()]
    return df


def sum_groups(df: pd.DataFrame) -> pd.DataFrame:
    '''
    This function sums the reads by species.
    '''
    df0 = pd.DataFrame(index=df.index)
    for i, j in zip(norm_counts_groups.keys(), norm_counts_groups.values()):
        df0[i] = df[j].sum(axis=1)
    df0.index.name = 'Respiratory_Seq_ID'
    return df0


def add_col(df: pd.DataFrame) -> pd.DataFrame:
    '''
    This function adds empty columns to the df based on the add_col_list.
    columns are inserted starting at index = 0
    '''
    for i in add_col_list:
        df.insert(0, i, '')
    return df


def fetch_data(df_in: pd.DataFrame) -> pd.DataFrame:
    '''
    This function fills in the metadata for the df based on the columns created using add_col.
    '''
    df = pd.read_csv('combined_keys.csv')
    # Tempe
    df_tempe = df[df['Site'].str.startswith('ASUT')]
    df_tempe_RVP = df_tempe[df_tempe['Sample_ID'].str.endswith('(2)')]
    df_tempe_RVP = df_tempe_RVP.rename(
        columns={'Seq_ID': 'Respiratory_Seq_ID'}).set_index('Respiratory_Seq_ID')

    # Poly
    df_polytechnic = df[df['Site'].str.startswith('ASUP')]
    df_polytechnic_RVP = df_polytechnic[df_polytechnic['Sample_ID'].str.endswith('(2)')]
    df_polytechnic_RVP = df_polytechnic_RVP.rename(
        columns={'Seq_ID': 'Respiratory_Seq_ID'}).set_index('Respiratory_Seq_ID')

    # West
    df_west = df[df['Site'].str.startswith('ASUW')]
    df_west_RVP = df_west[df_west['Sample_ID'].str.endswith('(2)')]
    df_west_RVP = df_west_RVP.rename(
        columns={'Seq_ID': 'Respiratory_Seq_ID'}).set_index('Respiratory_Seq_ID')

    # Downtown
    df_downtown = df[df['Site'].str.startswith('ASUDT')]
    df_downtown_RVP = df_downtown[df_downtown['Sample_ID'].str.endswith('(2)')]
    df_downtown_RVP = df_downtown_RVP.rename(
        columns={'Seq_ID': 'Respiratory_Seq_ID'}).set_index('Respiratory_Seq_ID')

    # Preparing updated grouped file
    df_out = df_in
    df_out.update(df_tempe_RVP)
    df_out.update(df_polytechnic_RVP)
    df_out.update(df_west_RVP)
    df_out.update(df_downtown_RVP)
    df_out['Date'] = df_out['Date'].astype(str).str.split('.').str[0]
    return df_out


def sort_df_data(df: pd.DataFrame) -> pd.DataFrame:
    '''
    This function converts dates in the Datetime column.
    There is a much easier way to do this.
    '''
    df['Date'] = df['Date'].astype(str)
    new_dates = []
    for i in df['Date']:
        j = i[:2] + '-' + i[2:4] + '-' + i[4:]
        k = j.split('-')
        l = k[1] + '-' + k[2] + '-' + k[0]
        new_dates.append(l)
    df['Date'] = new_dates
    df['Date'] = df['Date'].sort_values()
    df['Date'] = pd.to_datetime(df['Date'])
    df['Date'] = df['Date'].dt.strftime('%b-%d-%Y')
    df = df.set_index('Date')
    return df


def clean_data(df: pd.DataFrame) -> pd.DataFrame:
    '''
    This function removes select non RVP samples from the df.
    Update the drop_these_i_numbers.csv file as needed.
    '''
    df1 = pd.read_csv('drop_these_i_numbers.csv')
    df1 = df1.set_index('Respiratory_Seq_ID')
    df = df.drop([i for i in df1.index])
    return df


def group_RVP_outputs(infile: str) -> None:
    '''
    This function executes other helper functions sequentially to generate the RVP streamlit input file.
    '''
    os.chdir('/mnt/storage/wastewater_analysis/metadata')
    sorted_infile = sort_groups(infile)
    output_df = sum_groups(sorted_infile)
    idx_output_df = add_col(output_df)
    update_grouped_df = fetch_data(idx_output_df)
    clean_grouped_df = clean_data(update_grouped_df)
    # This is the RVP input file for the streamlit app.
    clean_grouped_df.to_csv(f'output_for_streamlit/RVP_counts_virusOnly_normalized_filled.csv')
    return None


def group_by_campus(df: pd.DataFrame) -> pd.DataFrame:
    '''
    This function is called during RVP_location_outputs to group campus data by dates.
    '''
    df_out = pd.DataFrame(columns=df.columns)
    for sub_df in df.groupby('Date'):
        df_len = len(sub_df[1])
        adjusted_df = sub_df[1].groupby('Date').sum(numeric_only=True) / df_len
        df_out = pd.concat([df_out, adjusted_df])
    return df_out


def modify_grouped_output(df: pd.DataFrame) -> pd.DataFrame:
    '''
    This function is called during RVP_location_outputs to fix specified dataframes.
    '''
    df = df.reset_index().drop(columns=['Respiratory_Seq_ID', 'Sample_ID', 'Date', 'Site'])
    df = df.rename(columns={'index': 'Date'})
    df = sort_df_data(df)
    return df


def RVP_location_outputs(infile: str) -> pd.DataFrame:
    '''
    This function calls prior functions to splits up the final output into campus specific datasets.
    '''
    os.chdir('/mnt/storage/wastewater_analysis/metadata/output_for_streamlit')
    df = pd.read_csv(infile).copy()
    df_dropna = df.dropna(axis=0)

    if len(df) - len(df_dropna) > 0:
        print('\nWARNING!\n')
        print(f'Total samples dropped : {len(df) - len(df_dropna)}\n')
        df_na = df[df.isnull().any(axis=1)]
        print(df_na)
        df_na.to_csv(f'log_files/Run{run}_seq_ids_missing_data.csv')
        print('The following file is missing data: /mnt/storage/wastewater_analysis/metadata/output_for_streamlit/RVP_counts_virusOnly_normalized_filled.csv')
        print(
            f'Please review the output file: "{run}_seq_id_missing_data.csv".\n Exiting script.\n')
        sys.exit()
    else:
        df_tempe = df_dropna[df_dropna['Site'].str.startswith('ASUT')]
        df_polytechnic = df_dropna[df_dropna['Site'].str.startswith('ASUP')]
        df_west = df_dropna[df_dropna['Site'].str.startswith('ASUW')]
        df_downtown = df_dropna[df_dropna['Site'].str.startswith('ASUDT')]
        df_tempe_out = modify_grouped_output(group_by_campus(df_tempe))
        df_polytechnic_out = modify_grouped_output(group_by_campus(df_polytechnic))
        df_west_out = modify_grouped_output(group_by_campus(df_west))
        df_downtown_out = modify_grouped_output(group_by_campus(df_downtown))
        df_tempe_out.to_csv('RVP_counts_Tempe_normalized_avg_groupby_date.csv')
        df_polytechnic_out.to_csv('RVP_counts_Polytechnic_normalized_avg_groupby_date.csv')
        df_west_out.to_csv('RVP_counts_West_normalized_avg_groupby_date.csv')
        df_downtown_out.to_csv('RVP_counts_Downtown_normalized_avg_groupby_date.csv')
    return None


#_____SCRIPT_EXECUTION____#
'''Notes
can add main function if needed

______Workflow______

Step1 : make_RVP_directory
    Directory setup
Step2 : get_fastq_path
    Generate path to sequencer data directory based on input
Step3 : start_RVP_monitoring
    Trigger monitoring script and other associated scripts
Step4 : setup_r_analysis_directory
    Prepares r analysis directory for processing.
Step5 : add_QC_metrics
    Fills metadata files (eg. RVP_metadata.txt)
Step6 : run_RVP_r_script
    Executes R analysis scripts
Step7 : group_RVP_outputs
    Prepares the output: "RVP_counts_virusOnly_normalized_filled.csv"
Step8 : RVP_location_outputs
    Uses the output from step 7 to generate files for streamlit.
'''

#____TROUBLESHOOTING____#
'''Notes
Uncomment block below for troubleshooting
'''
# make_RVP_directory(run)  # STEP_1
# sequencer_path = get_fastq_path(sequencer, instrament_run_number)  # STEP_2
# start_RVP_monitoring(run, sequencer, sequencer_path)  # STEP_3
# setup_r_analysis_directory(run)  # STEP_4
# add_QC_metrics(run)  # STEP_5
# run_RVP_r_script(run)  # STEP_6
# group_RVP_outputs('RVP_counts_virusOnly_normalized.txt')  # STEP_7
# RVP_location_outputs('RVP_counts_virusOnly_normalized_filled.csv')  # STEP_8

#_________________________#
# Comment out code below during troubleshooting

#_________________________#

sequencer_path = ''

with open(f'RVP_main_logs/Run{run}_RVP_main_log.txt', 'w') as f:
    f.write(f'{run}\n')

with open(f'RVP_main_logs/Run{run}_RVP_main_log.txt', 'w') as f:
    # STEP_1
    try:
        print('====================================================================\n')
        print('Creating RVP directory structure.\n')
        make_RVP_directory(run)
        f.write('Step 1: make_RVP_directory was processed successfully without error.\n')
    except Exception as e:
        f.write('An error occured during Step 1: make_RVP_directory.\n')
        f.write(f'{e}\n')
        print(e)
        print('Now exiting script.')
        sys.exit()
    # STEP_2
    try:
        sequencer_path = get_fastq_path(sequencer, instrament_run_number)
        print('====================================================================\n')
        print(f'Curent sequencer path: {sequencer_path}\n')
        f.write('Step 2: get_fastq_path was processed successfully without error.\n')
    except Exception as e:
        f.write('An error occured during Step 2: get_fastq_path.\n')
        f.write(f'{e}\n')
        print(e)
        print('Now exiting script.')
        sys.exit()
    # STEP_3
    try:
        print('====================================================================\n')
        print('Begin running monitoring script.')
        start_RVP_monitoring(run, sequencer, sequencer_path)
        f.write('Step 3: start_RVP_monitoring was processed successfully without error.\n')
    except Exception as e:
        f.write('An error occured during Step 3: get_fastq_path.\n')
        f.write(f'{e}\n')
        print(e)
        print('Now exiting script.')
        sys.exit()
    # STEP_4
    try:
        print('====================================================================\n')
        print('Populating r_analysis directory.\n')
        setup_r_analysis_directory(run)
        f.write('Step 4: setup_r_analysis_directory was processed successfully without error.\n')
    except Exception as e:
        f.write('An error occured during Step 4: setup_r_analysis_directory.\n')
        f.write(f'{e}\n')
        print(e)
        print('Now exiting script.')
        sys.exit()
    # STEP_5
    try:
        print('====================================================================\n')
        print('Updating RVP_metadta.txt file.\n')
        add_QC_metrics(run)
        f.write('Step 5: add_QC_metrics was processed successfully without error.\n')
    except Exception as e:
        f.write('An error occured during Step 5: add_QC_metrics.\n')
        f.write(f'{e}\n')
        print(e)
        print('Now exiting script.')
        sys.exit()
    # STEP_6
    try:
        print('====================================================================\n')
        print('Executing Rscript.\n')
        run_RVP_r_script(run)
        f.write('Step 6: run_RVP_r_script was processed successfully without error.\n')
    except Exception as e:
        f.write('An error occured during Step 6: run_RVP_r_script.\n')
        f.write(f'{e}\n')
        print(e)
        print('Now exiting script.')
        sys.exit()
    # STEP_7
    try:
        print('====================================================================\n')
        print('Merging normalized output data to database file.\n')
        group_RVP_outputs('RVP_counts_virusOnly_normalized.txt')
        f.write('Step 7: group_RVP_outputs was processed successfully without error.\n')
    except Exception as e:
        f.write('An error occured during Step 7: group_RVP_outputs.\n')
        f.write(f'{e}\n')
        print(e)
        print('Now exiting script.')
        sys.exit()
    # STEP_8
    try:
        print('====================================================================\n')
        print('Updating output files for Streamlit.\n')
        RVP_location_outputs('RVP_counts_virusOnly_normalized_filled.csv')
        f.write('Step 8: RVP_location_outputs was processed successfully without error.\n')
    except Exception as e:
        f.write('An error occured during Step 8: RVP_location_outputs.\n')
        f.write(f'{e}\n')
        print(e)
        print('Now exiting script.')
        sys.exit()
