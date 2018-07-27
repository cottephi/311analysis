#Created by Philippe Cotte philippe.cotte@cea.fr on 26/02/2018 for WA105 collaboration
#Requiers python3 with numpy and pandas

#This script sorts the 311 runs by "scan" type : drift, amplification, extraction and induction. It will look for runs that differ by only one field and store them in a list.
#The output is printed directly on the terminal, and can be either a raw information or formatted to C++ if/else conditions. Set output_format = "c++" to have the later.
#For a given scan type, it is possible to have several different configurations of the other 3 fields. So more than one list of runs is produced by scan type.
#The script assumes that a variation of Drift field inferior to 5% is not significant

#The input must be runs.csv, formatted this way: "run;Drift;Extraction;Amplification;Induction" (the field order does not matter but "run" must come first). The programm will load that in a pandas.DataFrame using the first row as header (columns names) and the first column (should be run numbers) as key. All fields from runs.csv are divided by 1000!

#The sorting loops basically work as follows:
  #First loop is over all the runs. Stores the fields values.
    #Second loop is over all the runs following the current run of loop 1. It will compare the field of the 2 runs. If one and only one field is different, it will store the field name in the variable "field_to_scan" and do the third loop
      #Third loop is over all the remaining runs. It will store all the runs whose only field difference with the run of loop 1 is "field_to_scan"

#The script is able to handle runs with exact same field configuration. If a run is strictly equal to the run of loop 1 is found, it is stored in a separate list, and later on added alongside its corresponding run of loop 1 in every list where the later appears.

import pandas
import os
import csv
import numpy
import sys

def check_diff(df_run_ref, df_run_to_compare):#Count how many fields are different between the two input DataFrames. A difference inferior to 5% in Drift will be considered as an equality
  diff = 0
  field_to_scan = ""
  drift_limit = .0
  for field in df_run_ref:
    if (float(df_run_ref[field]) != float(df_run_to_compare[field]) and field != "Drift") or ( (abs(float(df_run_ref[field]) - float(df_run_to_compare[field]))/float(df_run_ref[field]) > drift_limit or abs(float(df_run_ref[field]) - float(df_run_to_compare[field]))/float(df_run_to_compare[field]) > drift_limit) and field == "Drift"):
      diff = diff + 1
      field_to_scan = field
  return diff, field_to_scan
  

runs_file = "all_runs.csv"
if len(sys.argv) > 1:
  runs_file = sys.argv[1]
if not os.path.isfile(runs_file):
  print("ERROR: can't find",runs_file)
  exit()

df_raw_runs = pandas.read_csv(runs_file,sep=";",header=0,index_col=0)
runs_sorted_by_scan = { field:[] for field in df_raw_runs }
found_field_to_scan = { field:[] for field in df_raw_runs }
runs_to_add_later = []
#output_format = ""#
output_format = "c++"

for run_ref in df_raw_runs.index.values:#First loop
  sub_df_run_ref = df_raw_runs.loc[[run_ref]]
  l_df_runs_list = [sub_df_run_ref]
  
  for run in df_raw_runs.index.values:#Second loop
    if run <= run_ref:
      continue
    sub_df_run = df_raw_runs.loc[[run]]
    diff, field_to_scan = check_diff(sub_df_run_ref, sub_df_run)
    if diff != 1 or run in found_field_to_scan[field_to_scan]:
      if diff == 0:
        runs_to_add_later.append([sub_df_run_ref,sub_df_run])
        for ifield in found_field_to_scan:
          found_field_to_scan[ifield].append(run)
      continue
    l_df_runs_list.append(sub_df_run)
    found_field_to_scan[field_to_scan].append(run)
    if run_ref not in found_field_to_scan[field_to_scan]:
      found_field_to_scan[field_to_scan].append(run_ref)
    
    for third_run in df_raw_runs.index.values:#Third loop
      if third_run <= run:
        continue
      third_sub_df_run = df_raw_runs.loc[[third_run]]
      third_diff, third_field_to_scan = check_diff(sub_df_run_ref, third_sub_df_run)
      if third_diff != 1 or third_field_to_scan != field_to_scan or third_run in found_field_to_scan[field_to_scan]:
        if third_diff == 0:
          runs_to_add_later.append([sub_df_run_ref,third_sub_df_run])
        continue
      l_df_runs_list.append(third_sub_df_run)
      found_field_to_scan[field_to_scan].append(third_run)

    if len(l_df_runs_list) > 1:
      runs_sorted_by_scan[field_to_scan].append(l_df_runs_list)
      l_df_runs_list = [sub_df_run_ref]

#Output result (and add runs with all same fields)
if output_format == "c++":
  first = True
  for scan_field in runs_sorted_by_scan:
    i = 0
    if len(runs_sorted_by_scan[scan_field]) == 0:
      continue
    print("\n  //",scan_field,":")
    fields_to_print = [ field for field in runs_sorted_by_scan[scan_field][0][0] if field != scan_field ]
    for runs in runs_sorted_by_scan[scan_field]:
      i += 1
      if first:
        print('  if( scan_type == "' + scan_field + '" && scan_num ==',i,'){')
        first = False
      else:
        print('  else if( scan_type == "' + scan_field + '" && scan_num ==',i,"){")
      [ print("   //",field,":",float(runs[0][field])/1000.,' ', end='') for field in fields_to_print ]
      outstring = '\n    run_list = {'
      for run in runs:
        for run_to_add in runs_to_add_later:
          if run.index.values[0] == run_to_add[0].index.values[0] and run_to_add[1].index.values[0] not in [ irun.index.values[0] for irun in runs ]:
            runs.append(run_to_add[1])

      for run in runs:
        outstring += str(run.index.values[0]) + ", "
      outstring = outstring[:-2]
      outstring += '};\n  }'
      print(outstring)
  print('\n')
  print("  else{\n    cout<< \"Error: unknown combination of scan type and scan number (\" << scan_type << \", \" << scan_num << \").\" << endl;\n  continue;\n  }")


else:
  for scan_field in runs_sorted_by_scan:
    i = 0
    print("\n",scan_field,":")
    if len(runs_sorted_by_scan[scan_field]) == 0:
      print("No runs found")
      continue
    fields_to_print = [ field for field in runs_sorted_by_scan[scan_field][0][0] if field != scan_field ]
    for runs in runs_sorted_by_scan[scan_field]:
      [ print(" ",field,":",float(runs[0][field])/1000.,' ', end='') for field in fields_to_print ]
      outstring = '\n  runs : '
      for run in runs:
        for run_to_add in runs_to_add_later:
          if run.index.values[0] == run_to_add[0].index.values[0] and run_to_add[1].index.values[0] not in [ irun.index.values[0] for irun in runs ]:
            runs.append(run_to_add[1])

      for run in runs:
        outstring += str(run.index.values[0]) + ", "
      outstring = outstring[:-2]
      outstring += '\n  fields : '
      for run in runs:
        outstring += str(run.iloc[0][scan_field]/1000.) + ", "
      outstring = outstring[:-2]
      outstring += '\n'
      print(outstring)
  print('\n')
