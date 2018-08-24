import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import sys
import os
import getopt
import math

pedestal_runs = [729,1168, 1169, 1170, 1171]
pedestal_run_to_use = [729,1168]
pedestal_data_by_strip = {}
strips_to_plot = [0,25,320,319,345,1280]
strips_to_ignore = [200]
string_outdir = ""

def usage():
  print("Usage is:")
  print("-r run")
  print("-f first event")
  print("-l last event")
  print("-h")
  print("-p")
  
def DrawEvent(eventfilename, isped):
  data_view0 = {}
  data_view1 = {}
  print("  Opening file " + eventfilename)
  if len(pedestal_data_by_strip) != 1280 and not isped:
    print("You need a pedestal run! Please run this script on a pedestal run first.")
    exit(1)
  
  with open(eventfilename,"r") as eventfile:
    lines = eventfile.readlines()
    strip = 0
    for line in lines:
      ADCs = [ float(adc) for adc in line.split("\n")[0].split(" ")[:-1] ]
      if strip < 320:
        if not isped:
          data_view0[strip] = np.array([ adc - pedestal_data_by_strip[strip][0] for adc in ADCs ])
        else:
          data_view0[strip] = np.array(ADCs)
      elif strip < 1280:
        if not isped:
          data_view1[strip] = np.array([ adc - pedestal_data_by_strip[strip][0] for adc in ADCs ])
        else:
          data_view1[strip] = np.array(ADCs)
      else:
        print("    ERROR: unexpected strip number " + str(strip))
      strip += 1
  print("  Closed file " + eventfilename)

  outfile = string_outdir
  if isped:
    outfile = outfile + "/events_plots"
    if not os.path.isdir(outfile):
      print("  Creating directory " + outfile)
      os.mkdir(outfile)

  plt.ylabel("ADC")
  plt.xlabel("time bin")
  for strip in data_view0:
    if strip in strips_to_ignore:
      continue
    plt.plot(data_view0[strip])
  plt.savefig(outfile + "/" + eventfilename.split("/")[-1].split(".")[0] + ".pdf")
  plt.clf()
    
  
def AnalysePed(int_run, int_first_event, int_last_event, list_eventfiles, list_event_numbers):

  int_nevents = 0
  if not os.path.isdir(string_outdir + "/strips_plots"):
    print("  Creating directory " + string_outdir + "/strips_plots")
    os.mkdir(string_outdir + "/strips_plots")
  data_by_strips = {}
  for eventfilename,eventnumber in zip(list_eventfiles,list_event_numbers):
    if eventnumber < int_first_event or eventnumber > int_last_event:
      continue
    int_nevents += 1
    print("  Opening file " + eventfilename)
    with open(eventfilename,"r") as eventfile:
      lines = eventfile.readlines()
      strip = 0
      for line in lines:
        ADCs = np.array([ int(adc) for adc in line.split("\n")[0].split(" ")[:-1] ])
        if strip < 1280:
          if not strip in data_by_strips:
            data_by_strips[strip] = [ADCs]
          else:
            data_by_strips[strip].append(ADCs)
        else:
          print("    ERROR: unexpected strip number " + str(strip))
        strip += 1
    print("  Closed file " + eventfilename)

  means_by_strips = {}
  rmss_by_strips = {}
  mean_means_by_strips = {}
  mean_rmss_by_strips = {}
  error_for_means = {}
  error_for_rmss = {}
  for strip in data_by_strips:
    means_by_strips[strip] = np.array([ np.mean(data_by_strips[strip][event]) for event in range(0,len(data_by_strips[strip])) ])
    mean_means_by_strips[strip] = np.mean(means_by_strips[strip])
    error_for_means[strip] = math.sqrt(np.var(means_by_strips[strip])/int_nevents)
    rmss_by_strips[strip]  = np.array([ math.sqrt(np.var(data_by_strips[strip][event]))  for event in range(0,len(data_by_strips[strip])) ])
    mean_rmss_by_strips[strip] = np.mean(rmss_by_strips[strip])
    error_for_rmss[strip] = math.sqrt(np.var(rmss_by_strips[strip])/int_nevents)
  
#  with open(string_outdir + "/mean_means_and_rmss.txt","w") as outfile:
#    for strip in mean_means_by_strips:
#      outfile.write(str(strip) + " " + str(mean_means_by_strips[strip]) + " " + str(mean_rmss_by_strips[strip]) + "\n")
  
  for strip in means_by_strips:
    if strip in strips_to_plot:
      ax0 = plt.subplot2grid((2,2),(0,0), colspan = 2)
      plt.xlabel("time bin")
      plt.ylabel("ADC")
      for event in data_by_strips[strip]:
        ax0.plot(event)
      
      ax1 = plt.subplot2grid((2,2),(1,0))
      plt.xlabel("mean adc in events")
      counts_mean, binned_means, _ = plt.hist(means_by_strips[strip])
      maxy_mean = max(counts_mean)
      ax1.set_ylim([0,maxy_mean*1.1])
      ax2 = plt.subplot2grid((2,2),(1,1))
      plt.xlabel("rms adc in events")
      counts_rms, binned_rms, _ = plt.hist(rmss_by_strips[strip])
      maxy_rms = max(counts_rms)
      ax2.set_ylim([0,maxy_rms*1.1])
      plt.savefig(string_outdir + "/strips_plots/" + str(strip) + ".pdf")
      plt.clf()
      
  for strip in strips_to_ignore:
    mean_means_by_strips.pop(strip,None)
    mean_rmss_by_strips.pop(strip,None)
    error_for_means.pop(strip,None)
    error_for_rmss.pop(strip,None)
  
  ax3 = plt.subplot(211)
  plt.xlabel("strip number")
  plt.ylabel("mean ADC")
  ax3.errorbar(list(mean_means_by_strips.keys()), list(mean_means_by_strips.values()), yerr=list(error_for_means.values()), marker="o", ms=2, mfc='red', mec='green', ls='None')
  ax3 = plt.subplot(212)
  plt.xlabel("strip number")
  plt.ylabel("rms ADC")
  ax3.errorbar(mean_rmss_by_strips.keys(),mean_rmss_by_strips.values(), yerr=error_for_rmss.values(), marker="o", ms=2, mfc='green', mec='red', ls='None')
  plt.savefig(string_outdir + "/mean_over_all_run.pdf")

def LoadPed(int_run):
  if not int_run in pedestal_runs:
    print("ERROR: run " + str(int_run) + " not a calibration run")
    exit(1)
  if not os.path.isfile("plots/" + str(int_run) + "/mean_means_and_rmss.txt"):
     print("ERROR: Calibration file not found. Please run this script for run " + str(int_run))
     exit(1)
     
  global pedestal_data_by_strip
  with open("plots/" + str(int_run) + "/mean_means_and_rmss.txt","r") as ifile:
    lines = ifile.readlines()
    for line in lines:
      strip = int(line.split(" ")[0])
      mean = float(line.split(" ")[1])
      rms = float(line.split("\n")[0].split(" ")[2])
      pedestal_data_by_strip[strip] = [mean, rms]

def main(argv):

  int_run = -1; 
  int_first_event = -1;
  int_last_event = -1;
  compute_ped = False
  
  try:
    opts, args = getopt.getopt(sys.argv[1:], "r:f:l:hp", [])
  except getopt.GetoptError as err:
    print(str(err))
    usage()
    exit(0)
  for opt, arg in opts:
    if opt == "-h":
      usage()
      exit(0)
    elif opt == "-r":
      int_run = int(arg)
    elif opt == "-f":
      int_first_event = int(arg)
    elif opt == "-l":
      int_last_event = int(arg)
    elif opt == "-p":
      compute_ped = True
    else:
      usage()
  
  string_indir = "textfiles/" + str(int_run)
  if not os.path.isdir(string_indir):
    print("ERROR: " + string_indir + " not found.")
    exit(1)
  global string_outdir
  string_outdir = "plots/"
  if not os.path.isdir(string_outdir):
    print("Creating directory " + string_outdir)
    os.mkdir(string_outdir)
  string_outdir = string_outdir + str(int_run)
  if not os.path.isdir(string_outdir):
    print("Creating directory " + string_outdir)
    os.mkdir(string_outdir)
    
  list_eventfiles = sorted(glob(string_indir + '/*'))
  list_event_numbers = [ int(evt.split("/")[-1].split(".")[0]) for evt in list_eventfiles ]
  if len(list_eventfiles) == 0:
    print("ERROR: no files found in " + str(string_indir))
    exit(1)
  if int_first_event < 0 and int_last_event < 0:
    if int_run in pedestal_runs:
      if compute_ped:
        print("Will plot pedestal histograms of all events of run " + str(int_run))
        AnalysePed(int_run, min(list_event_numbers), max(list_event_numbers), list_eventfiles, list_event_numbers)
      else:
        print("Will process all extracted events in run " + str(int_run))
        for eventfile in list_eventfiles:
          DrawEvent(eventfile, True)
    else:
      if int_run <= 840:
        LoadPed(pedestal_run_to_use[0])
      if int_run > 840:
        LoadPed(pedestal_run_to_use[1])
      print("Will process all extracted events in run " + str(int_run))
      for eventfile in list_eventfiles:
        DrawEvent(eventfile, False)
  else:
    if int_first_event < 0:
      int_first_event = min(list_event_numbers)
    if int_last_event < 0:
      int_last_event = max(list_event_numbers)
    if min(list_event_numbers) > int_first_event:
      print("Event " + str(int_first_event) + " not found. Will start at first found event: " + str(min(list_event_numbers)))
      int_first_event = min(list_event_numbers)
    if max(list_event_numbers) < int_last_event:
      print("Event " + str(int_last_event) + " not found. Will end at last found event: " + str(max(list_event_numbers)))
      int_last_event = max(list_event_numbers)
    if int_run in pedestal_runs:
      if compute_ped:
        print("Will plot pedestal histograms for event " + str(int_first_event) + " to " + str(int_last_event) + " in run " + str(int_run))
        AnalysePed(int_run, int_first_event, int_last_event, list_eventfiles, list_event_numbers)
      else:
        print("Will process events from " + str(int_first_event) + " to " + str(int_last_event) + " in run " + str(int_run))
        for eventfile,eventnumber in zip(list_eventfiles,list_event_numbers):
          if eventnumber < int_first_event or eventnumber > int_last_event:
            continue
          DrawEvent(eventfile, True)
    else:
      if int_run <= 840:
        LoadPed(pedestal_run_to_use[0])
      if int_run > 840:
        LoadPed(pedestal_run_to_use[1])
      print("Will process events from " + str(int_first_event) + " to " + str(int_last_event) + " in run " + str(int_run))
      for eventfile,eventnumber in zip(list_eventfiles,list_event_numbers):
        if eventnumber < int_first_event or eventnumber > int_last_event:
          continue
        DrawEvent(eventfile, False)
    
  
  
  
  
if __name__ == '__main__':
  if len(sys.argv) == 1:
    usage()
    exit(0)
  main(sys.argv[1:])
