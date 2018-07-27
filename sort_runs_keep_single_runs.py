import pandas
import os
import csv
import numpy
import sys

runs_file = "runs.csv"
if len(sys.argv) > 1:
  runs_file = sys.argv[1]
if not os.path.isfile(runs_file):
  print("ERROR: can't find",runs_file)
  exit()

df_raw_runs = pandas.read_csv(runs_file,sep=";",header=0,index_col=0)
runs_sorted_by_fields = {}

#output_format = ""#
output_format = "c++"

for run_ref in df_raw_runs.index.values:#First loop
  df_run_ref = df_raw_runs.loc[[run_ref]]
  tmp_field_list = tuple([ df_run_ref[field].values[0] for field in df_run_ref ])
  drift_current = tmp_field_list[0]
  other_fields_current = [ tmp_field_list[i] for i in range(1,len(tmp_field_list)) if i != 2 ]
  recorded_run = False
  if len(runs_sorted_by_fields) == 0:
    runs_sorted_by_fields[tmp_field_list] = [run_ref]
  else:
    for field in runs_sorted_by_fields:
      drift = field[0]
      other_fields = [ field[i] for i in range(1,len(field)) if i != 2 ]
      if other_fields_current == other_fields and float(drift_current) < float(drift)*1.12 and float(drift_current) > float(drift)*(1.-0.12):
        runs_sorted_by_fields[field].append(run_ref)
        recorded_run = True
        break
  if not recorded_run:
    runs_sorted_by_fields[tmp_field_list] = [run_ref]
    
i = 0
for stuff in runs_sorted_by_fields:
  i = i+1
  to_print = ""
  if i == 1:
    to_print = "  if( scan_type == \"All\" && scan_num == " + str(i) + " ){\n    run_list = {"
  else:
    to_print = "  else if( scan_type == \"All\" && scan_num == " + str(i) + " ){\n    run_list = {"
  for run in runs_sorted_by_fields[stuff]:
    to_print = to_print + str(run) + ", "
  to_print = to_print[:-2] + "};\n  }"
  print(to_print)
