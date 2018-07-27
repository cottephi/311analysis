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
  drift_limit = .5
  for field in df_run_ref:
    if (float(df_run_ref[field]) != float(df_run_to_compare[field]) and field != "Drift") or ( (abs(float(df_run_ref[field]) - float(df_run_to_compare[field]))/float(df_run_ref[field]) > drift_limit or abs(float(df_run_ref[field]) - float(df_run_to_compare[field]))/float(df_run_to_compare[field]) > drift_limit) and field == "Drift"):
      diff = 1
  return diff
  

runs_file = "all_runs.csv"
if len(sys.argv) > 1:
  runs_file = sys.argv[1]
if not os.path.isfile(runs_file):
  print("ERROR: can't find",runs_file)
  exit()

df_raw_runs = pandas.read_csv(runs_file,sep=";",header=0,index_col=0)
sorted_runs = []
runs_seen = []
#output_format = ""#
output_format = "c++"

for run_ref in df_raw_runs.index.values:#First loop
  sub_df_run_ref = df_raw_runs.loc[[run_ref]]
  if run_ref in runs_seen:
    continue
  runs_seen.append(run_ref)
  fields = [ [field,float(sub_df_run_ref[field])] for field in sub_df_run_ref ]
  sorted_runs.append([[run_ref,fields]])
  
  for run in df_raw_runs.index.values:#Second loop
    if run <= run_ref or run in runs_seen:
      continue
    sub_df_run = df_raw_runs.loc[[run]]
    diff = check_diff(sub_df_run_ref, sub_df_run)
    if diff == 0:
      sorted_runs[-1].append([run])
      runs_seen.append(run)
     
  if len(sorted_runs[-1]) == 1:
    sorted_runs = sorted_runs[:-1]

#Output result (and add runs with all same fields)
if output_format == "c++":
  first = True
  print("\n  //Identical runs:")
  i = 0
  for runs in sorted_runs:
    i = i + 1
    fields = runs[0][1]
    if first:
      print('  if( scan_type == "Identical" && scan_num ==',i,'){')
      first = False
    else:
      print('  else if( scan_type == "Identical" && scan_num ==',i,"){")
    [ print("   //",field[0],":",field[1]/1000.,' ', end='') for field in fields ]
    outstring = '\n    run_list = {'
    for run in runs:
      outstring += str(run[0]) + ", "
    outstring = outstring[:-2]
    outstring += '};\n  }'
    print(outstring)
  print('\n')
