import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import gridspec
import ROOT
from datetime import datetime
from dateutil import tz
import mpl_toolkits.mplot3d as mp3d
from scipy.optimize import curve_fit
import warnings

def file_exists(filename): #returns True if file exists, False if it does not
    return bool(glob(filename))

def get_number_of_subruns(run): #returns number of subruns for a given run
    number_of_subruns = 0
    subrun = 0
    while file_exists('/eos/experiment/wa105/data/311/rawdata/' + str(run) + '/' + str(run) + '-' + str(subrun) + '.dat'):
        number_of_subruns += 1
        subrun += 1
    return int(number_of_subruns)

def get_number_of_events_subrun(run, subrun): #returns number of events for a given run and subrun
    filename = '/eos/experiment/wa105/data/311/rawdata/' + str(run) + '/' + str(run) + '-' + str(subrun) + '.dat'
    file_size = os.path.getsize(filename) #get file size in bytes
    total_event_number = file_size - 5.0 - 4.0 #First 5 bytes contain run header, last 4 bytes contain run footer.
    total_event_number /= 35.0 + 1280.0 * 1667.0 * 1.5 #For each event, the first 35 bytes are the event header, then come the 1667 ADC counts for the 1280 channels stored in 12 bit (1.5 bytes) format.
    if not total_event_number.is_integer():
        print('Number of events is not an integer. Something is fucking wrong!')
    else:
        return int(total_event_number)

def get_number_of_events_run(run): #returns number of events for a given run for all its subruns
    total_event_number_run = 0
    subrun = 0
    while file_exists('/eos/experiment/wa105/data/311/rawdata/' + str(run) + '/' + str(run) + '-' + str(subrun) + '.dat'):
        total_event_number_run += get_number_of_events_subrun(run=run, subrun=subrun)
        subrun += 1
    return int(total_event_number_run)

def read_one_event(run, subrun, event): # Reads 3x1x1 binary files and returns ADC counts for each channel in a 1280(channels) x 1667(ticks) array
    filename = '/eos/experiment/wa105/data/311/rawdata/' + str(run) + '/' + str(run) + '-' + str(subrun) + '.dat'
    with open(filename, 'rb') as binary_file:
        binary_file.seek(5 + event * int(35.0 + 1280.0 * 1667.0 * 1.5) + 35)
        ADC_values = []
        for i in range(int(1280.0 * 1667.0 / 2.0)):
            byte_1 = binary_file.read(1)
            byte_2 = binary_file.read(1)
            byte_3 = binary_file.read(1)
            ADC_values.append((ord(byte_1) << 4) + (ord(byte_2) >> 4))
            ADC_values.append(((ord(byte_2) & 15) << 8) + ord(byte_3))
    ADC_values_reshaped = np.array(ADC_values).reshape(1280, 1667)
    
    # Remap channel from DAQ mapping (binary file mapping) to view mapping
    channel_map = []
    channel_arrangement = [295, 135, 263, 103, 231, 71, 199, 39, 167, 7, 487, 327, 519, 359, 551, 391, 583, 423, 615, 455, 807, 647, 839, 679, 871, 711, 903, 743, 935, 775, 1127, 967, 1159, 999, 1191, 1031, 1223, 1063, 1255, 1095]
    for i in channel_arrangement:
        for j in range(4):
            channel_map.extend(range(i+8*j, i+8*j-8, -1))
    ADC_values_viewchannel = np.zeros((1280, 1667))
    for i in range(1280):
        ADC_values_viewchannel[channel_map[i]] = ADC_values_reshaped[i]
    return ADC_values_viewchannel

def subtract_pedestal(ADC_values, method='median'): # subtracts the base level of the readout, setting the true 0:
    ADC_values_minped = np.zeros((1280, 1667))
    if method == 'median':
        for i in range(1280):
            ADC_values_minped[i] = ADC_values[i] - np.median(ADC_values[i])
    elif method == 'mean':
        for i in range(1280):
            ADC_values_minped[i] = ADC_values[i] - np.mean(ADC_values[i])
    else:
        print 'Method not recognized. Fucking implement it yourself or check what you typed!'
    return ADC_values_minped

#def event_display(run, subrun, event, clim):
#    ADC_values = read_one_event(run=run, subrun=subrun, event=event)
#    ADC_values_minped = subtract_pedestal(ADC_values)
#    fig = plt.figure()
#    fig.suptitle('Run ' + str(run) + ', SubRun ' + str(subrun) + ', Event ' + str(event))
#    gs = gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[1, 3])
    
#    plt.subplot(gs[0])
#    plt.title('View 0')
#    ax0 = plt.gca()
#    ax0.set_xlabel('channel number')
#    ax0.set_ylabel('tick')
#    plt.xticks(np.arange(0, 321, 160))
#    img0 = ax0.imshow(np.rot90(ADC_values_minped[:320]), extent=(0, 320, 0, 1667), interpolation='nearest', cmap=plt.cm.jet, origin='upper', clim=clim, aspect=320.0/1667.0)
#    ax0.invert_yaxis()
    
#    plt.subplot(gs[1])
#    plt.title('View 1')
#    ax1 = plt.gca()
#    ax1.set_xlabel('channel number')
#    ax1.set_ylabel('tick')
#    plt.xticks(np.arange(0, 961, 160))
#    img1 = ax1.imshow(np.rot90(ADC_values_minped[320:]), extent=(0, 960, 0, 1667), interpolation='nearest', cmap=plt.cm.jet, origin='upper', clim=clim, aspect=960.0/1667.0/3.0)
#    ax1.invert_yaxis()
    
#    p0 = ax0.get_position().get_points().flatten()
#    p1 = ax1.get_position().get_points().flatten()
#    ax_cbar = fig.add_axes([p0[0], 0.45 * p0[2], p1[2]-p0[0], 0.02])
#    fig.colorbar(img0, cax=ax_cbar, orientation='horizontal', label='ADC counts')
#    plt.show()

def plot_waveform(run, subrun, event, channel, view=None):
    ADC_values = read_one_event(run=run, subrun=subrun, event=event)
    ADC_values_minped = subtract_pedestal(ADC_values)
    if view == 0:
        ADC_values_one_channel = ADC_values_minped[channel]
        title = 'Run ' + str(run) + ', SubRun ' + str(subrun) + ', Event ' + str(event) + ', Channel ' + str(channel)
    elif view == 1:
        ADC_values_one_channel = ADC_values_minped[320 + channel]
        title = 'Run ' + str(run) + ', SubRun ' + str(subrun) + ', Event ' + str(event) + ', Channel ' + str(320 + channel)
    else:
        ADC_values_one_channel = ADC_values_minped[channel]
        title = 'Run ' + str(run) + ', SubRun ' + str(subrun) + ', Event ' + str(event) + ', Channel ' + str(channel)
    
    fig = plt.figure()
    plt.title(title)
    ax = plt.gca()
    ax.grid(True)
    ax.set_xlabel('tick')
    ax.set_ylabel('ADC')
    ax.plot(np.arange(0.0, 1667.0, 1.0), ADC_values_one_channel, color='navy')
    plt.show()

def get_reconstruction_variables(run, subrun, event): # Reads all reconstruction variables from the root files and returns a dictionary with all of them.
    dictionary = {}
    filename = '/eos/experiment/wa105/offline/LArSoft/Data/Parser/2018_Feb_05/' + str(run) + '/' + str(run) + '-' + str(subrun) + '-Parser.root'
    root_file = ROOT.TFile.Open(filename)
    root_file_tree = root_file.Get('analysistree/anatree')
    root_file_tree.GetEntry(event)
    
    event_time_seconds = root_file_tree.EventTimeSeconds
    dictionary['event_time_seconds'] = event_time_seconds
    event_time_utc_datetime = datetime.fromtimestamp(event_time_seconds)
    event_time_utc_datetime = event_time_utc_datetime.replace(tzinfo=tz.gettz('UTC'))
    event_time_geneva_datetime = event_time_utc_datetime.astimezone(tz.gettz('Europe/Zurich'))
    dictionary['event_time_geneva_datetime'] = event_time_geneva_datetime
    dictionary['event_time_nanoseconds'] = root_file_tree.EventTimeNanoseconds
    number_of_hits = root_file_tree.NumberOfHits
    dictionary['number_of_hits'] = number_of_hits
    hit_tpc = []
    hit_view = []
    hit_channel = []
    hit_peak_time = []
    hit_charge_summed_adc = []
    hit_charge_integral = []
    hit_peak_height = []
    hit_start_time = []
    hit_end_time = []
    hit_width = []
    hit_goodness_of_fit = []
    hit_multiplicity = []
    hit_track_id = []
    hit_cluster_id = []
    for i in range(0, number_of_hits):
        hit_tpc.append(root_file_tree.Hit_TPC[i])
        hit_view.append(root_file_tree.Hit_View[i])
        hit_channel.append(root_file_tree.Hit_Channel[i])
        hit_peak_time.append(root_file_tree.Hit_PeakTime[i])
        hit_start_time.append(root_file_tree.Hit_StartTime[i])
        hit_end_time.append(root_file_tree.Hit_EndTime[i])
        hit_width.append(root_file_tree.Hit_Width[i])
        hit_goodness_of_fit.append(root_file_tree.Hit_GoodnessOfFit[i])
        hit_multiplicity.append(root_file_tree.Hit_Multiplicity[i])
        hit_track_id.append(root_file_tree.Hit_TrackID[i])
        hit_cluster_id.append(root_file_tree.Hit_ClusterID[i])
    dictionary['hit_tpc'] = np.array(hit_tpc)
    dictionary['hit_view'] = np.array(hit_view)
    dictionary['hit_channel'] = np.array(hit_channel)
    dictionary['hit_peak_time'] = np.array(hit_peak_time)
    dictionary['hit_start_time'] = np.array(hit_start_time)
    dictionary['hit_end_time'] = np.array(hit_end_time)
    dictionary['hit_width'] = np.array(hit_width)
    dictionary['hit_goodness_of_fit'] = np.array(hit_goodness_of_fit)
    dictionary['hit_multiplicity'] = np.array(hit_multiplicity)
    dictionary['hit_track_id'] = np.array(hit_track_id)
    dictionary['hit_cluster_id'] = np.array(hit_cluster_id)
    
    number_of_clusters = root_file_tree.NumberOfClusters
    dictionary['number_of_clusters'] = number_of_clusters
    cluster_id = []
    cluster_number_of_hits = []
    cluster_view = []
    cluster_charge_integral = []
    cluster_charge_integral_average_per_hit = []
    cluster_start_channel = []
    cluster_start_tick = []
    cluster_end_channel = []
    cluster_end_tick = []
    cluster_start_charge = []
    cluster_start_angle = []
    cluster_end_charge = []
    cluster_end_angle = []
    for i in range(0, number_of_clusters):
        cluster_id.append(root_file_tree.ClusterID[i])
        cluster_number_of_hits.append(root_file_tree.Cluster_NumberOfHits[i])
        cluster_view.append(root_file_tree.Cluster_View[i])
        cluster_charge_integral.append(root_file_tree.Cluster_ChargeIntegral[i])
        cluster_charge_integral_average_per_hit.append(root_file_tree.Cluster_ChargeIntegralAveragePerHit[i])
        cluster_start_channel.append(root_file_tree.Cluster_StartChannel[i])
        cluster_start_tick.append(root_file_tree.Cluster_StartTick[i])
        cluster_end_channel.append(root_file_tree.Cluster_EndChannel[i])
        cluster_end_tick.append(root_file_tree.Cluster_EndTick[i])
        cluster_start_charge.append(root_file_tree.Cluster_StartCharge[i])
        cluster_start_angle.append(root_file_tree.Cluster_StartAngle[i])
        cluster_end_charge.append(root_file_tree.Cluster_EndCharge[i])
        cluster_end_angle.append(root_file_tree.Cluster_EndAngle[i])
    dictionary['cluster_id'] = np.array(cluster_id)
    dictionary['cluster_number_of_hits'] = np.array(cluster_number_of_hits)
    dictionary['cluster_view'] = np.array(cluster_view)
    dictionary['cluster_charge_integral'] = np.array(cluster_charge_integral)
    dictionary['cluster_charge_integral_average_per_hit'] = np.array(cluster_charge_integral_average_per_hit)
    dictionary['cluster_start_channel'] = np.array(cluster_start_channel)
    dictionary['cluster_start_tick'] = np.array(cluster_start_tick)
    dictionary['cluster_end_channel'] = np.array(cluster_end_channel)
    dictionary['cluster_end_tick'] = np.array(cluster_end_tick)
    dictionary['cluster_start_charge'] = np.array(cluster_start_charge)
    dictionary['cluster_start_angle'] = np.array(cluster_start_angle)
    dictionary['cluster_end_charge'] = np.array(cluster_end_charge)
    dictionary['cluster_end_angle'] = np.array(cluster_end_angle)
    
    number_of_tracks = root_file_tree.NumberOfTracks_pmtrack
    dictionary['number_of_tracks'] = number_of_tracks
    track_id = []
    track_number_of_hits = []
    track_length = []
    track_start_point_x = []
    track_start_point_y = []
    track_start_point_z = []
    track_start_point_distance_to_boundary = []
    track_end_point_x = []
    track_end_point_y = []
    track_end_point_z = []
    track_end_point_distance_to_boundary = []
    track_start_direction_theta = []
    track_start_direction_phi = []
    track_start_direction_x = []
    track_start_direction_y = []
    track_start_direction_z = []
    track_end_direction_theta = []
    track_end_direction_phi = []
    track_end_direction_x = []
    track_end_direction_y = []
    track_end_direction_z = []
    track_pitch_in_views = []
    track_number_of_hits_per_view = []
    index_hit_inside_track = 0
    track_hit_x = []
    track_hit_y = []
    track_hit_z = []
    track_dx_local_track_direction = []
    track_dx_3d_position = []
    track_hit_tpc = []
    track_hit_view = []
    track_hit_channel = []
    track_hit_peak_time =[]
    track_hit_charge_summed_adc = []
    track_hit_charge_integral = []
    track_hit_peak_height = []
    track_hit_start_time = []
    track_hit_end_time = []
    track_hit_width = []
    track_hit_goodness_of_fit = []
    track_hit_multiplicity = []
    for i in range(0, number_of_tracks):
        track_id.append(root_file_tree.TrackID_pmtrack[i])
        current_track_number_of_hits = root_file_tree.Track_NumberOfHits_pmtrack[i]
        track_number_of_hits.append(current_track_number_of_hits)
        track_length.append(root_file_tree.Track_Length_pmtrack[i])
        track_start_point_x.append(root_file_tree.Track_StartPoint_X_pmtrack[i])
        track_start_point_y.append(root_file_tree.Track_StartPoint_Y_pmtrack[i])
        track_start_point_z.append(root_file_tree.Track_StartPoint_Z_pmtrack[i])
        track_start_point_distance_to_boundary.append(root_file_tree.Track_StartPoint_DistanceToBoundary_pmtrack[i])
        track_end_point_x.append(root_file_tree.Track_EndPoint_X_pmtrack[i])
        track_end_point_y.append(root_file_tree.Track_EndPoint_Y_pmtrack[i])
        track_end_point_z.append(root_file_tree.Track_EndPoint_Z_pmtrack[i])
        track_end_point_distance_to_boundary.append(root_file_tree.Track_EndPoint_DistanceToBoundary_pmtrack[i])
        track_start_direction_theta.append(root_file_tree.Track_StartDirection_Theta_pmtrack[i])
        track_start_direction_phi.append(root_file_tree.Track_StartDirection_Phi_pmtrack[i])
        track_start_direction_x.append(root_file_tree.Track_StartDirection_X_pmtrack[i])
        track_start_direction_y.append(root_file_tree.Track_StartDirection_Y_pmtrack[i])
        track_start_direction_z.append(root_file_tree.Track_StartDirection_Z_pmtrack[i])
        track_end_direction_theta.append(root_file_tree.Track_EndDirection_Theta_pmtrack[i])
        track_end_direction_phi.append(root_file_tree.Track_EndDirection_Phi_pmtrack[i])
        track_end_direction_x.append(root_file_tree.Track_EndDirection_X_pmtrack[i])
        track_end_direction_y.append(root_file_tree.Track_EndDirection_Y_pmtrack[i])
        track_end_direction_z.append(root_file_tree.Track_EndDirection_Z_pmtrack[i])
        track_pitch_in_views.append(np.array([root_file_tree.Track_PitchInViews_pmtrack[2*i], root_file_tree.Track_PitchInViews_pmtrack[2*i+1]]))
        track_number_of_hits_per_view.append(np.array([root_file_tree.Track_NumberOfHitsPerView_pmtrack[2*i], root_file_tree.Track_NumberOfHitsPerView_pmtrack[2*i+1]]))
        current_track_hit_x = []
        current_track_hit_y = []
        current_track_hit_z = []
        current_track_dx_local_track_direction = []
        current_track_dx_3d_position = []
        current_track_hit_tpc = []
        current_track_hit_view = []
        current_track_hit_channel = []
        current_track_hit_peak_time =[]
        current_track_hit_charge_summed_adc = []
        current_track_hit_charge_integral = []
        current_track_hit_peak_height = []
        current_track_hit_start_time = []
        current_track_hit_end_time = []
        current_track_hit_width = []
        current_track_hit_goodness_of_fit = []
        current_track_hit_multiplicity = []
        for j in range(0, current_track_number_of_hits):
            current_track_hit_x.append(root_file_tree.Track_Hit_X_pmtrack[index_hit_inside_track + j])
            current_track_hit_y.append(root_file_tree.Track_Hit_Y_pmtrack[index_hit_inside_track + j])
            current_track_hit_z.append(root_file_tree.Track_Hit_Z_pmtrack[index_hit_inside_track + j])
            current_track_dx_local_track_direction.append(root_file_tree.Track_Hit_dx_LocalTrackDirection_pmtrack[index_hit_inside_track + j])
            current_track_dx_3d_position.append(root_file_tree.Track_Hit_dx_3DPosition_pmtrack[index_hit_inside_track + j])
            current_track_hit_tpc.append(root_file_tree.Track_Hit_TPC_pmtrack[index_hit_inside_track + j])
            current_track_hit_view.append(root_file_tree.Track_Hit_View_pmtrack[index_hit_inside_track + j])
            current_track_hit_channel.append(root_file_tree.Track_Hit_Channel_pmtrack[index_hit_inside_track + j])
            current_track_hit_peak_time.append(root_file_tree.Track_Hit_PeakTime_pmtrack[index_hit_inside_track + j])
            current_track_hit_charge_summed_adc.append(root_file_tree.Track_Hit_ChargeSummedADC_pmtrack[index_hit_inside_track + j])
            current_track_hit_charge_integral.append(root_file_tree.Track_Hit_ChargeIntegral_pmtrack[index_hit_inside_track + j])
            current_track_hit_peak_height.append(root_file_tree.Track_Hit_PeakHeight_pmtrack[index_hit_inside_track + j])
            current_track_hit_start_time.append(root_file_tree.Track_Hit_StartTime_pmtrack[index_hit_inside_track + j])
            current_track_hit_end_time.append(root_file_tree.Track_Hit_EndTime_pmtrack[index_hit_inside_track + j])
            current_track_hit_width.append(root_file_tree.Track_Hit_Width_pmtrack[index_hit_inside_track + j])
            current_track_hit_goodness_of_fit.append(root_file_tree.Track_Hit_GoodnessOfFit_pmtrack[index_hit_inside_track + j])
            current_track_hit_multiplicity.append(root_file_tree.Track_Hit_Multiplicity_pmtrack[index_hit_inside_track + j])
        index_hit_inside_track += current_track_number_of_hits
        track_hit_x.append(np.array(current_track_hit_x))
        track_hit_y.append(np.array(current_track_hit_y))
        track_hit_z.append(np.array(current_track_hit_z))
        track_dx_local_track_direction.append(np.array(current_track_dx_local_track_direction))
        track_dx_3d_position.append(np.array(current_track_dx_3d_position))
        track_hit_tpc.append(np.array(current_track_hit_tpc))
        track_hit_view.append(np.array(current_track_hit_view))
        track_hit_channel.append(np.array(current_track_hit_channel))
        track_hit_peak_time.append(np.array(current_track_hit_peak_time))
        track_hit_charge_summed_adc.append(np.array(current_track_hit_charge_summed_adc))
        track_hit_charge_integral.append(np.array(current_track_hit_charge_integral))
        track_hit_peak_height.append(np.array(current_track_hit_peak_height))
        track_hit_start_time.append(np.array(current_track_hit_start_time))
        track_hit_end_time.append(np.array(current_track_hit_end_time))
        track_hit_width.append(np.array(current_track_hit_width))
        track_hit_goodness_of_fit.append(np.array(current_track_hit_goodness_of_fit))
        track_hit_multiplicity.append(np.array(current_track_hit_multiplicity))
    dictionary['track_id'] = np.array(track_id)
    dictionary['track_number_of_hits'] = np.array(track_number_of_hits)
    dictionary['track_length'] = np.array(track_length)
    dictionary['track_start_point_x'] = np.array(track_start_point_x)
    dictionary['track_start_point_y'] = np.array(track_start_point_y)
    dictionary['track_start_point_z'] = np.array(track_start_point_z)
    dictionary['track_start_point_distance_to_boundary'] = np.array(track_start_point_distance_to_boundary)
    dictionary['track_end_point_x'] = np.array(track_end_point_x)
    dictionary['track_end_point_y'] = np.array(track_end_point_y)
    dictionary['track_end_point_z'] = np.array(track_end_point_z)
    dictionary['track_end_point_distance_to_boundary'] = np.array(track_end_point_distance_to_boundary)
    dictionary['track_start_direction_theta'] = np.array(track_start_direction_theta)
    dictionary['track_start_direction_phi'] = np.array(track_start_direction_phi)
    dictionary['track_start_direction_x'] = np.array(track_start_direction_x)
    dictionary['track_start_direction_y'] = np.array(track_start_direction_y)
    dictionary['track_start_direction_z'] = np.array(track_start_direction_z)
    dictionary['track_end_direction_theta'] = np.array(track_end_direction_theta)
    dictionary['track_end_direction_phi'] = np.array(track_end_direction_phi)
    dictionary['track_end_direction_x'] = np.array(track_end_direction_x)
    dictionary['track_end_direction_y'] = np.array(track_end_direction_y)
    dictionary['track_end_direction_z'] = np.array(track_end_direction_z)
    dictionary['track_pitch_in_views'] = np.array(track_pitch_in_views)
    dictionary['track_number_of_hits_per_view'] = np.array(track_number_of_hits_per_view)
    dictionary['track_hit_x'] = np.array(track_hit_x)
    dictionary['track_hit_y'] = np.array(track_hit_y)
    dictionary['track_hit_z'] = np.array(track_hit_z)
    dictionary['track_dx_local_track_direction'] = np.array(track_dx_local_track_direction)
    dictionary['track_dx_3d_position'] = np.array(track_dx_3d_position)
    dictionary['track_hit_tpc'] = np.array(track_hit_tpc)
    dictionary['track_hit_view'] = np.array(track_hit_view)
    dictionary['track_hit_channel'] = np.array(track_hit_channel)
    dictionary['track_hit_peak_time'] = np.array(track_hit_peak_time)
    dictionary['track_hit_charge_summed_adc'] = np.array(track_hit_charge_summed_adc)
    dictionary['track_hit_charge_integral'] = np.array(track_hit_charge_integral)
    dictionary['track_hit_peak_height'] = np.array(track_hit_peak_height)
    dictionary['track_hit_start_time'] = np.array(track_hit_start_time)
    dictionary['track_hit_end_time'] = np.array(track_hit_end_time)
    dictionary['track_hit_width'] = np.array(track_hit_width)
    dictionary['track_hit_goodness_of_fit'] = np.array(track_hit_goodness_of_fit)
    dictionary['track_hit_multiplicity'] = np.array(track_hit_multiplicity)
    return dictionary

def reco_display_3d(run, subrun, event):
    reco_info = get_reconstruction_variables(run=run, subrun=subrun, event=event)
    bottom = [(0.0, -50.0, -50.0), (0.0, 50.0, -50.0), (300.0, 50.0, -50.0), (300.0, -50.0, -50.0)]
    top = [(0.0, -50.0, 50.0), (0.0, 50.0, 50.0), (300.0, 50.0, 50.0), (300.0, -50.0, 50.0)]
    toptop = [(0.0, -50.0, 60.0), (0.0, 50.0, 60.0), (300.0, 50.0, 60.0), (300.0, -50.0, 60.0)]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('Z')
    ax.set_ylabel('Y')
    ax.set_zlabel('X')
    opacity = 0.3
    face1 = mp3d.art3d.Poly3DCollection([bottom], alpha=opacity, linewidth=1)
    face2 = mp3d.art3d.Poly3DCollection([top], alpha=opacity, linewidth=1)
    face3 = mp3d.art3d.Poly3DCollection([toptop], alpha=opacity, linewidth=1)
    face1.set_facecolor((0.5, 0.5, 0.5, opacity))
    face2.set_facecolor((0, 0, 1, opacity))
    face3.set_facecolor((0, 0, 1, opacity))
    ax.add_collection3d(face1)
    ax.add_collection3d(face2)
    ax.add_collection3d(face3)
    ax.set_zlim(-155.0, 155.0)
    ax.set_ylim(-155.0, 155.0)
    ax.set_xlim(-5.0, 305.0)
    for i in range(len(reco_info['track_hit_x'])):
        for j in range(len(reco_info['track_hit_x'][i])):
            if reco_info['track_hit_view'][i][j] == 0:
                ax.scatter(reco_info['track_hit_z'][i], reco_info['track_hit_y'][i], reco_info['track_hit_x'][i], depthshade=False, marker='o')
            if reco_info['track_hit_view'][i][j] == 1:
                ax.scatter(reco_info['track_hit_z'][i], reco_info['track_hit_y'][i], reco_info['track_hit_x'][i], depthshade=False, marker='v')
    plt.show()

def point_to_line_dist(m, b, x0, y0):
    #Formulae derived analytically. 
    # minimum distance satisfies (slope <dot> (x1-x0, y1-y0) = 0)
    #     where (x0,y0) is the point, and (x1,y1) is the point on the line that is closest the point in question
    linePointX = (x0 + y0*m -b*m)/(m*m + 1)
    linePointY = (m*m*y0 + m*x0 +b)/(m*m + 1)
    distance = ((linePointX-x0)**2 + (linePointY-y0)**2)**0.5
    return distance

#uses the numpy polyfit to make a fit for a list of [x] points and a list of [y] points
# returns a dictionary of the fit parameters
def get_line_of_best_fit(x, y):
    fit_parameters = np.polyfit(x, y, 1)  # polyfit returns highest-degree parameters first. y=mx+b -> [m,b]
    angle_radians = np.arctan(fit_parameters[0])  
    angle_degrees = angle_radians * 180.0 / np.pi
    # note: this angle is measured from the 'horizontal' to the slope. 
    #       /|
    #      / |
    #     /  | m    a = arctan(m/1)
    #    /   |
    #   /a___|
    #     1
    return {'slope': fit_parameters[0], 'y_intercept': fit_parameters[1], 'angle_radians': angle_radians, 'angle_degrees': angle_degrees}

#takes a list of [x] points and [y] points, plus the width of a rectangle to put around the line
#         change- feed the fit into this function since it's calculated right before calculated again here
def get_points_rectangle_around_line(x, y, width, line_fit):
    start_point_x = np.min(x)
    end_point_x = np.max(x)
    
    slope = line_fit['slope']
    y_intercept = line_fit['y_intercept']
    angle_radians = line_fit['angle_radians']
    
    #find the y-intercepts of top and bottom lines of the box, given the width of the box
    #       to visualize, expand a narrow rectangle around an angled line, draw a line straight up from the middle of a narrow end and extend the top of the box to the drawn line
    #       the extra distance up is what we're adding to the y_intercept (the height of the box-side-midpoint
    parallel_line_1_y_intercept = y_intercept + width / 2.0 / np.cos(angle_radians)
    parallel_line_2_y_intercept = y_intercept - width / 2.0 / np.cos(angle_radians)
    
    #finding the slope of the other sides of the rectangle. SINCE these are perpendicular, thie is the slope
    perpendicular_lines_slope = - 1.0 / slope

    # these y-intercepts are calculated by knowing the slope of the perpendicular lines and that they intersect the OG fit at the smallest x-value of the fit.
    #       basicaly::   mx_{0} + b = (-1/m)x_{0} + b_{perp}, and solve for b_{perp}
    perpendicular_line_1_y_intercept = (slope * start_point_x + y_intercept) - perpendicular_lines_slope * start_point_x
    perpendicular_line_2_y_intercept = (slope * end_point_x + y_intercept) - perpendicular_lines_slope * end_point_x
    # these are similar. A parallel line and a perpendicular line meet at these points, and so you solve for (x) when the parametrizations are equal, then use that (x) to find (y). 
    point_1 = np.array([(perpendicular_line_1_y_intercept - parallel_line_1_y_intercept) / (slope - perpendicular_lines_slope), slope * (perpendicular_line_1_y_intercept - parallel_line_1_y_intercept) / (slope - perpendicular_lines_slope) + parallel_line_1_y_intercept])
    point_2 = np.array([(perpendicular_line_2_y_intercept - parallel_line_1_y_intercept) / (slope - perpendicular_lines_slope), slope * (perpendicular_line_2_y_intercept - parallel_line_1_y_intercept) / (slope - perpendicular_lines_slope) + parallel_line_1_y_intercept])
    point_3 = np.array([(perpendicular_line_2_y_intercept - parallel_line_2_y_intercept) / (slope - perpendicular_lines_slope), slope * (perpendicular_line_2_y_intercept - parallel_line_2_y_intercept) / (slope - perpendicular_lines_slope) + parallel_line_2_y_intercept])
    point_4 = np.array([(perpendicular_line_1_y_intercept - parallel_line_2_y_intercept) / (slope - perpendicular_lines_slope), slope * (perpendicular_line_1_y_intercept - parallel_line_2_y_intercept) / (slope - perpendicular_lines_slope) + parallel_line_2_y_intercept])
    area_rectangle = np.sqrt((point_1[0] - point_2[0])**2.0 + (point_1[1] - point_2[1])**2.0) * np.sqrt((point_2[0] - point_3[0])**2.0 + (point_2[1] - point_3[1])**2.0)
    return {'rectangle_points': np.array([point_1, point_2, point_3, point_4]), 'fit_slope': slope, 'fit_y_intercept': y_intercept, 'parallel_y_intercepts': np.array([parallel_line_1_y_intercept, parallel_line_2_y_intercept]), 'perpendicular_slope': perpendicular_lines_slope, 'perpendicular_y_intercepts': np.array([perpendicular_line_1_y_intercept, perpendicular_line_2_y_intercept]), 'area_rectangle': area_rectangle }

#calculates the area of a triangle given these three points. **deprecated**
#no longer used since now using a different method for checking if a point is in a rectangle
def area_triangle(point_1, point_2, point_3):
    area = 0.5 * np.abs(((point_2[0] * point_3[1] - point_3[0] * point_2[1]) - (point_1[0] * point_3[1] - point_3[0] * point_1[1]) + (point_1[0] * point_2[1] - point_2[0] * point_1[1])))
    return area

# vector product, used in is_point_in_rectangle
def dotProd(p0, p1):
    return(p0[0]*p1[0] + p0[1]*p1[1])

#vector difference 
def vSub(p0, p1):
    p = [p1[0]-p0[0], p1[1]-p0[1]]
    return p

def sort_rectangle_points(points_rectangle):
    # sort the rectangle points so they go in order around. 
    # this is important, becase the analytical way requires three rect-points ABC, such that AB is perp to BC
    xmax = 0
    xmin = 3
    for i in range(3):
        if points_rectangle[i+1][0]>points_rectangle[xmax][0]:
            xmax=i+1
        if points_rectangle[2-i][0]<points_rectangle[xmin][0]:
            xmin=2-i
    # this give points on extremes of the rectangle, which will be on opposite sides of the rectangle.
    ordRectangle = [points_rectangle[xmin], [0,0],points_rectangle[xmax], [0,0]]
    
    #found the minimum and maximum ones, so now we just need to choose one of the remainders to go in the middle
    #can go either clockwise around (ABCD) or (ADBC). These don't make a difference! 
    
    # I'm sure there's a better way...
    added = False
    for i in range(4):
        if (i!=xmax and i!=xmin):
            if (added):
                ordRectangle[3]=points_rectangle[i]
            else:
                ordRectangle[1]=points_rectangle[i]
                added = True
    
    return(ordRectangle)

def get_slope_between(p1, p2):
    return( (p2[1]-p1[1])/(p2[0]-p1[0]))

#calculates the cross product of vectors v1 and v2
def xProd(v1, v2):
    if len(v1) != len(v2):
        print("error. Vectors need be of the same dimensionality")
        return([0,0,0])
    if len(v1)==2:
        #assuming 3rd component is zero. 
        return([0.0, 0.0, 	v1[0]*v2[1] - v2[0]*v1[1]])
    if len(v1)==3:
        print("3d vectors not yet supported")
        return([0,0,0])
    else:
        return([0,0,0])

# given a an array of three points defining a triangle, determines whether a provided point lies within the triangle. Assumes 2D euclidean geometry.
def is_point_in_triangle( triangle, point):
    # tri points ABC, point P
    AB = vSub(triangle[0], triangle[1])
    BC = vSub(triangle[1], triangle[2])
    CA = vSub(triangle[2], triangle[0])
    AP = vSub(triangle[0], point)
    BP = vSub(triangle[1], point)
    CP = vSub(triangle[2], point)
    
    x1 = xProd(AB, AP)[2]
    x2 = xProd(BC, BP)[2]
    x3 = xProd(CA, CP)[2]
    
    x1x2 = x1*x2
    x2x3 = x2*x3
    
    if (x1x2>0 and x2x3>0):
        return True
    else:
        return False

#formula a fun little byproduct of analytical geometry!
def is_point_in_rectangle(points_rectangle, point, arg):

    #Starting from the min-x one, going sequentially, label the points ABC and D. Our point of curiosity is M
	# not sure if I can explain the logic behind this in the comments, sorry. 
	# 
	# I verified this by generating random rectangles, ordering the points with the sort_rectangle_points function
	# 	   plotting the points, and using this function to fill in the space contained by the points.
    
    AB = vSub(points_rectangle[0], points_rectangle[1])
    AM = vSub(points_rectangle[0], point)
    BC = vSub(points_rectangle[1], points_rectangle[2])
    BM = vSub(points_rectangle[1], point)
    dotABAM = dotProd(AB, AM)
    dotABAB = dotProd(AB, AB)
    dotBCBM = dotProd(BC, BM)
    dotBCBC = dotProd(BC, BC)
    
    if( (0 < dotABAM) and (dotABAM < dotABAB) and (0 < dotBCBM) and (dotBCBM < dotBCBC) ):
        return True
    else:
        slope = arg['slope']
        # it is possible that track-hits can be inside these little triangles at the edge of the boxes, though they aren't inside the boxes
        # as a result, even if a point is not in a box, we may still want to add it in the charge sum. 
        # there are a lot of cases depending on the box orientation and order of points in the rectangle array
        if slope > 0:
            # up-going box.
            if get_slope_between(points_rectangle[0], points_rectangle[1]) > 0:
                #clockwise from min ABCD
                tri1Point1 = points_rectangle[3]
                tri1Point2 = [0.5*points_rectangle[0][0]+0.5*points_rectangle[3][0], 0.5*points_rectangle[0][1]+0.5*points_rectangle[3][1]]
                tri1Point3 = [tri1Point2[0], tri1Point1[1]]
                
                tri2Point1 = points_rectangle[1]
                tri2Point2 = [0.5*points_rectangle[1][0]+0.5*points_rectangle[2][0], 0.5*points_rectangle[1][1]+0.5*points_rectangle[2][1]]
                tri2Point3 = [tri2Point2[0], tri2Point1[1]]
            else:
                #counterclockwise from min ABCD
                tri1Point1 = points_rectangle[1]
                tri1Point2 = [0.5*points_rectangle[0][0]+0.5*points_rectangle[1][0], 0.5*points_rectangle[0][1]+0.5*points_rectangle[1][1]]
                tri1Point3 = [tri1Point2[0], tri1Point1[1]]
                
                tri2Point1 = points_rectangle[3]
                tri2Point2 = [0.5*points_rectangle[3][0]+0.5*points_rectangle[2][0], 0.5*points_rectangle[3][1]+0.5*points_rectangle[2][1]]
                tri2Point3 = [tri2Point2[0], tri2Point1[1]]
        elif (slope < 0):
            #down going box
            if get_slope_between(points_rectangle[0], points_rectangle[1]) > 0:
                tri1Point1 = points_rectangle[1]
                tri1Point2 = [0.5*points_rectangle[0][0]+0.5*points_rectangle[1][0], 0.5*points_rectangle[0][1]+0.5*points_rectangle[1][1]]
                tri1Point3 = [tri1Point2[0], tri1Point1[1]]
                
                tri2Point1 = points_rectangle[3]
                tri2Point2 = [0.5*points_rectangle[3][0]+0.5*points_rectangle[2][0], 0.5*points_rectangle[3][1]+0.5*points_rectangle[2][1]]
                tri2Point3 = [tri2Point2[0], tri2Point1[1]]
                
            else:
                tri1Point1 = points_rectangle[3]
                tri1Point2 = [0.5*points_rectangle[0][0]+0.5*points_rectangle[3][0], 0.5*points_rectangle[0][1]+0.5*points_rectangle[3][1]]
                tri1Point3 = [tri1Point2[0], tri1Point1[1]]
                
                tri2Point1 = points_rectangle[1]
                tri2Point2 = [0.5*points_rectangle[1][0]+0.5*points_rectangle[2][0], 0.5*points_rectangle[1][1]+0.5*points_rectangle[2][1]]
                tri2Point3 = [tri2Point2[0], tri2Point1[1]]
        else:
            return False
        return( is_point_in_triangle([tri1Point1, tri1Point2, tri1Point3], point) or is_point_in_triangle([tri2Point1, tri2Point2, tri2Point3], point))



# ---- begin unit conversions ----

def convert_tick_to_x_position(tick, drift_velocity=0.00160562889065):
#    drift_velocity_cm_per_us = ((-0.0464*(87.0 - 105.749) + 1.0) * (1.88125 * 0.5 * np.log(1.0 + 0.99408 / 0.5) + 0.01172 * 0.5**4.20214) + 0.01712 * (87.0 - 105.749)) / 10.0 #in cm/us
#    drift_velocity_m_per_us = drift_velocity_cm_per_us * 0.01
    x_position = - float(tick) * 0.4 * drift_velocity + 0.5
    return x_position

def convert_x_position_to_tick(x_position, drift_velocity=0.00160562889065):
    tick = int(np.floor((0.5-x_position)/(0.4*drift_velocity)))
    return tick

def convert_x_position_to_us(x_position, drift_velocity=0.00160562889065):
#    drift_velocity_cm_per_us = ((-0.0464*(87.0 - 105.749) + 1.0) * (1.88125 * 0.5 * np.log(1.0 + 0.99408 / 0.5) + 0.01172 * 0.5**4.20214) + 0.01712 * (87.0 - 105.749)) / 10.0 #in cm/us
#    drift_velocity_m_per_us = drift_velocity_cm_per_us * 0.01
    time_us = (0.5 - x_position / 100.0) / drift_velocity
    return time_us

def convert_y_position_to_channel_number_view0(y_position):
    channel_number = int(np.floor((y_position+0.5)/0.003125))
    return channel_number

def convert_z_position_to_channel_number_view1(z_position):
    channel_number = int(np.floor(z_position/ 0.003125))
    return channel_number

def convert_channel_number_view0_to_y_position(channel_number):
    y_position = float(channel_number) * 0.003125 - 0.5
    return y_position

def convert_channel_number_view1_to_z_position(channel_number):
    z_position = float(channel_number) * 0.003125
    return z_position

# ---- end unit conversions ----





def get_reconstruction_track_hit_positions(run, subrun, event):
    dictionary = {}
    filename = '/eos/experiment/wa105/offline/LArSoft/Data/Parser/2018_Feb_05/' + str(run) + '/' + str(run) + '-' + str(subrun) + '-Parser.root'
    root_file = ROOT.TFile.Open(filename)
    root_file_tree = root_file.Get('analysistree/anatree')
    root_file_tree.GetEntry(event)
    
    number_of_hits = root_file_tree.NumberOfHits
    dictionary['number_of_hits'] = number_of_hits
    
    number_of_tracks = root_file_tree.NumberOfTracks_pmtrack
    dictionary['number_of_tracks'] = number_of_tracks
    index_hit_inside_track = 0
    track_hit_x = []
    track_hit_y = []
    track_hit_z = []
    track_hit_view = []
    for i in range(0, number_of_tracks):
        current_track_number_of_hits = root_file_tree.Track_NumberOfHits_pmtrack[i]
        current_track_hit_x = []
        current_track_hit_y = []
        current_track_hit_z = []
        current_track_hit_view = []
        for j in range(0, current_track_number_of_hits):
            current_track_hit_x.append(root_file_tree.Track_Hit_X_pmtrack[index_hit_inside_track + j])
            current_track_hit_y.append(root_file_tree.Track_Hit_Y_pmtrack[index_hit_inside_track + j])
            current_track_hit_z.append(root_file_tree.Track_Hit_Z_pmtrack[index_hit_inside_track + j])
            current_track_hit_view.append(root_file_tree.Track_Hit_View_pmtrack[index_hit_inside_track + j])
        index_hit_inside_track += current_track_number_of_hits
        track_hit_x.append(np.array(current_track_hit_x))
        track_hit_y.append(np.array(current_track_hit_y))
        track_hit_z.append(np.array(current_track_hit_z))
        track_hit_view.append(np.array(current_track_hit_view))
    dictionary['track_hit_x'] = np.array(track_hit_x)
    dictionary['track_hit_y'] = np.array(track_hit_y)
    dictionary['track_hit_z'] = np.array(track_hit_z)
    dictionary['track_hit_view'] = np.array(track_hit_view)
    return dictionary

def get_track_2D_positions_in_each_view(run, subrun, event, track):
    track_variables = get_reconstruction_track_hit_positions(run, subrun, event)
    reco_hits_x_view0 = []
    reco_hits_y_view0 = []
    reco_hits_x_view1 = []
    reco_hits_z_view1 = []
    for i in range(len(track_variables['track_hit_x'][track])):
        if track_variables['track_hit_view'][track][i] == 0:
            reco_hits_x_view0.append(track_variables['track_hit_x'][track][i])
            reco_hits_y_view0.append(track_variables['track_hit_y'][track][i])
        elif track_variables['track_hit_view'][track][i] == 1:
            reco_hits_x_view1.append(track_variables['track_hit_x'][track][i])
            reco_hits_z_view1.append(track_variables['track_hit_z'][track][i])
    reco_hits_x_view0 = np.array(reco_hits_x_view0) / 100.0
    reco_hits_y_view0 = np.array(reco_hits_y_view0) / 100.0
    reco_hits_x_view1 = np.array(reco_hits_x_view1) / 100.0
    reco_hits_z_view1 = np.array(reco_hits_z_view1) / 100.0
    return {'reco_hits_x_view0': reco_hits_x_view0, 'reco_hits_y_view0': reco_hits_y_view0, 'reco_hits_x_view1': reco_hits_x_view1, 'reco_hits_z_view1': reco_hits_z_view1}

def get_rectangle_parameters_for_track(run, subrun, event, track, rectangle_width):
    reco_hits_track = get_track_2D_positions_in_each_view(run, subrun, event, track)
    lin_fit_view0 = get_line_of_best_fit(reco_hits_track['reco_hits_y_view0'], reco_hits_track['reco_hits_x_view0'])
    lin_fit_view1 = get_line_of_best_fit(reco_hits_track['reco_hits_z_view1'], reco_hits_track['reco_hits_x_view1'])
    rectangle_points_view0 = get_points_rectangle_around_line(reco_hits_track['reco_hits_y_view0'], reco_hits_track['reco_hits_x_view0'], rectangle_width, lin_fit_view0)
    rectangle_points_view1 = get_points_rectangle_around_line(reco_hits_track['reco_hits_z_view1'], reco_hits_track['reco_hits_x_view1'], rectangle_width, lin_fit_view1)
    return {'lin_fit_view0': lin_fit_view0, 'lin_fit_view1': lin_fit_view1, 'rectangle_points_view0': rectangle_points_view0, 'rectangle_points_view1': rectangle_points_view1, 'reco_hits_x_view0': reco_hits_track['reco_hits_x_view0'], 'reco_hits_y_view0': reco_hits_track['reco_hits_y_view0'], 'reco_hits_x_view1': reco_hits_track['reco_hits_x_view1'], 'reco_hits_z_view1': reco_hits_track['reco_hits_z_view1']}

def segment_track_in_y(x, y, number_of_segments):
    max_y = np.max(y)
    min_y = np.min(y)
    segments_y = []
    for i in range(number_of_segments):
        segments_y.append([min_y + i * (max_y - min_y) / number_of_segments, min_y + (i + 1) * (max_y - min_y) / number_of_segments])
    x_segmented = []
    y_segmented = []
    for i in segments_y:
        x_1segment = []
        y_1segment = []
        for j in range(np.size(x)):
            if y[j] + 1.0e-6 >= i[0] and y[j] < i[1] + 1.0e-6:
                x_1segment.append(x[j])
                y_1segment.append(y[j])
        x_segmented.append(np.array(x_1segment))
        y_segmented.append(np.array(y_1segment))
    return np.array(x_segmented), np.array(y_segmented)

def get_rectangle_parameters_for_segment_of_track(run, subrun, event, track, rectangle_width, reco_hits_track):
    lin_fit_view0 = get_line_of_best_fit(reco_hits_track['reco_hits_y_view0'], reco_hits_track['reco_hits_x_view0'])
    lin_fit_view1 = get_line_of_best_fit(reco_hits_track['reco_hits_z_view1'], reco_hits_track['reco_hits_x_view1'])
    rectangle_points_view0 = get_points_rectangle_around_line(reco_hits_track['reco_hits_y_view0'], reco_hits_track['reco_hits_x_view0'], rectangle_width, lin_fit_view0)
    rectangle_points_view1 = get_points_rectangle_around_line(reco_hits_track['reco_hits_z_view1'], reco_hits_track['reco_hits_x_view1'], rectangle_width, lin_fit_view1)
    return {'lin_fit_view0': lin_fit_view0, 'lin_fit_view1': lin_fit_view1, 'rectangle_points_view0': rectangle_points_view0, 'rectangle_points_view1': rectangle_points_view1, 'reco_hits_x_view0': reco_hits_track['reco_hits_x_view0'], 'reco_hits_y_view0': reco_hits_track['reco_hits_y_view0'], 'reco_hits_x_view1': reco_hits_track['reco_hits_x_view1'], 'reco_hits_z_view1': reco_hits_track['reco_hits_z_view1']}

def sum_up_charge_in_rectangle(run, subrun, event, track, rectangle_width, include_all_points_in_rectangle=False):
    rectangle_parameters_track = get_rectangle_parameters_for_track(run, subrun, event, track, rectangle_width)
    raw_info = read_one_event(run, subrun, event)
    raw_info_minped = subtract_pedestal(raw_info)
    max_chan_view0_rectangle = np.max(rectangle_parameters_track['rectangle_points_view0']['rectangle_points'][:, 0])
    min_chan_view0_rectangle = np.min(rectangle_parameters_track['rectangle_points_view0']['rectangle_points'][:, 0])
    max_tick_view0_rectangle = np.max(rectangle_parameters_track['rectangle_points_view0']['rectangle_points'][:, 1])
    min_tick_view0_rectangle = np.min(rectangle_parameters_track['rectangle_points_view0']['rectangle_points'][:, 1])
    max_chan_view1_rectangle = np.max(rectangle_parameters_track['rectangle_points_view1']['rectangle_points'][:, 0])
    min_chan_view1_rectangle = np.min(rectangle_parameters_track['rectangle_points_view1']['rectangle_points'][:, 0])
    max_tick_view1_rectangle = np.max(rectangle_parameters_track['rectangle_points_view1']['rectangle_points'][:, 1])
    min_tick_view1_rectangle = np.min(rectangle_parameters_track['rectangle_points_view1']['rectangle_points'][:, 1])
    charge_rectangle_view0 = 0.0
    charge_rectangle_view1 = 0.0
    if include_all_points_in_rectangle:
        points_in_rectangle_view0_y = []
        points_in_rectangle_view0_x = []
        points_in_rectangle_view1_z = []
        points_in_rectangle_view1_x = []
    for channel_view0 in range(320):
        y_position_channel = convert_channel_number_view0_to_y_position(channel_view0)
        if min_chan_view0_rectangle < y_position_channel and y_position_channel < max_chan_view0_rectangle:
            for tick_view0 in range(1667):
                x_position_tick = convert_tick_to_x_position(tick_view0)
                if min_tick_view0_rectangle < x_position_tick and x_position_tick < max_tick_view0_rectangle:
                    if is_point_in_rectangle(rectangle_parameters_track['rectangle_points_view0']['rectangle_points'], np.array([y_position_channel, x_position_tick]), rectangle_parameters_track['rectangle_points_view0']['area_rectangle']):
                        charge_rectangle_view0 += raw_info_minped[channel_view0][tick_view0]
                        if include_all_points_in_rectangle:
                            points_in_rectangle_view0_y.append(y_position_channel)
                            points_in_rectangle_view0_x.append(x_position_tick)
    for channel_view1 in range(960):
        z_position_channel = convert_channel_number_view1_to_z_position(channel_view1)
        if min_chan_view1_rectangle < z_position_channel and z_position_channel < max_chan_view1_rectangle:
            for tick_view1 in range(1667):
                x_position_tick = convert_tick_to_x_position(tick_view1)
                if min_tick_view1_rectangle < x_position_tick and x_position_tick < max_tick_view1_rectangle:
                    if is_point_in_rectangle(rectangle_parameters_track['rectangle_points_view1']['rectangle_points'], np.array([z_position_channel, x_position_tick]), rectangle_parameters_track['rectangle_points_view1']['area_rectangle']):
                        charge_rectangle_view1 += raw_info_minped[channel_view1+320][tick_view1]
                        if include_all_points_in_rectangle:
                            points_in_rectangle_view1_z.append(z_position_channel)
                            points_in_rectangle_view1_x.append(x_position_tick)
    charge_fC_rectangle_view0 = charge_rectangle_view0 / 54.77942
    charge_fC_rectangle_view1 = charge_rectangle_view1 / 67.08538
    if include_all_points_in_rectangle:
        points_in_rectangle_view0_y = np.array(points_in_rectangle_view0_y)
        points_in_rectangle_view0_x = np.array(points_in_rectangle_view0_x)
        points_in_rectangle_view1_z = np.array(points_in_rectangle_view1_z)
        points_in_rectangle_view1_x = np.array(points_in_rectangle_view1_x)
    if include_all_points_in_rectangle:
        dictionary = {'ADC_counts_minped': raw_info_minped, 'charge_fC_rectangle_view0': charge_fC_rectangle_view0, 'charge_fC_rectangle_view1': charge_fC_rectangle_view1, 'points_in_rectangle_view0_y': points_in_rectangle_view0_y, 'points_in_rectangle_view0_x': points_in_rectangle_view0_x, 'points_in_rectangle_view1_z': points_in_rectangle_view1_z, 'points_in_rectangle_view1_x': points_in_rectangle_view1_x}
    else:
        dictionary = {'ADC_counts_minped': raw_info_minped, 'charge_fC_rectangle_view0': charge_fC_rectangle_view0, 'charge_fC_rectangle_view1': charge_fC_rectangle_view1}
    dictionary_tot = dict(dictionary.items() + rectangle_parameters_track.items())
    return dictionary_tot

#using our defined rectangles around the tracks, we find the charge enclosed.
def sum_up_charge_in_rectangle_for_track_segment(run, subrun, event, track, rectangle_width, rectangle_parameters_track, include_all_points_in_rectangle=False):
    raw_info = read_one_event(run, subrun, event)
    raw_info_minped = subtract_pedestal(raw_info)
    max_chan_view0_rectangle = np.max(rectangle_parameters_track['rectangle_points_view0']['rectangle_points'][:, 0])
    min_chan_view0_rectangle = np.min(rectangle_parameters_track['rectangle_points_view0']['rectangle_points'][:, 0])
    max_tick_view0_rectangle = np.max(rectangle_parameters_track['rectangle_points_view0']['rectangle_points'][:, 1])
    min_tick_view0_rectangle = np.min(rectangle_parameters_track['rectangle_points_view0']['rectangle_points'][:, 1])
    max_chan_view1_rectangle = np.max(rectangle_parameters_track['rectangle_points_view1']['rectangle_points'][:, 0])
    min_chan_view1_rectangle = np.min(rectangle_parameters_track['rectangle_points_view1']['rectangle_points'][:, 0])
    max_tick_view1_rectangle = np.max(rectangle_parameters_track['rectangle_points_view1']['rectangle_points'][:, 1])
    min_tick_view1_rectangle = np.min(rectangle_parameters_track['rectangle_points_view1']['rectangle_points'][:, 1])
    charge_rectangle_view0 = 0.0
    charge_rectangle_view1 = 0.0
    if include_all_points_in_rectangle:
        points_in_rectangle_view0_y = []
        points_in_rectangle_view0_x = []
        points_in_rectangle_view1_z = []
        points_in_rectangle_view1_x = []

    # get the channel range around the rectangle
    # I add the +/- 1 as a buffer. 
    # start:
    channel_view0 = convert_y_position_to_channel_number_view0(min_chan_view0_rectangle) - 1
    #channel_view0 = 0
    # end:
    last_channel =  convert_y_position_to_channel_number_view0(max_chan_view0_rectangle) + 1
    
    # here I call a function to sort the rectangle points into a proper sequential format needed by my 
    # is_point_in_rectangle function
    sorted_rectangle_points_0= sort_rectangle_points( rectangle_parameters_track['rectangle_points_view0']['rectangle_points'] )
    sorted_rectangle_points_1= sort_rectangle_points( rectangle_parameters_track['rectangle_points_view1']['rectangle_points'] )
    
    while(channel_view0<last_channel and channel_view0 < 320 ):
        y_position_channel = convert_channel_number_view0_to_y_position(channel_view0)
        #now, for these, we're looking for the range in the box. need fit
        # I'm using the fit, which is used to define the box as well, to find the bottom and top ticks of the box at this y-address
        # please note the weird notation. The "y_intercept" of this line is actually an x-intercept
        # ALSO, the conversion between tick and x-position has a sign change, meaning small tick correspond to high position and vice versa
        #      what thsi means is that the max position is the smallest tick and the min position is the biggest tick! 
        # similar to before, I add the  +/-1 as a buffer
        tick_view0=convert_x_position_to_tick(y_position_channel*rectangle_parameters_track['lin_fit_view0']['slope'] + rectangle_parameters_track['lin_fit_view0']['y_intercept'] + 0.5*rectangle_width/np.cos(rectangle_parameters_track['lin_fit_view0']['angle_radians']) ) -1
        last_tick= convert_x_position_to_tick(y_position_channel*rectangle_parameters_track['lin_fit_view0']['slope'] + rectangle_parameters_track['lin_fit_view0']['y_intercept'] - 0.5*rectangle_width/np.cos(rectangle_parameters_track['lin_fit_view0']['angle_radians']) ) +1

        # the buffer and the fit can give nonsense negative ticks (since I parametrizing the above conversion for every contingency is not worth it),
        # so I instead just make sure it doesn't go too low or too high
        if tick_view0<0:
            tick_view0 = 0
 
        while(tick_view0<1667 and tick_view0<last_tick):
            # scanning across all ticks in the tick range
            # convert tick to position
            x_position_tick = convert_tick_to_x_position(tick_view0)
                
            if is_point_in_rectangle( sorted_rectangle_points_0 , [y_position_channel, x_position_tick], rectangle_parameters_track['lin_fit_view0']  ):
                #print("returned true")
                charge_rectangle_view0 += raw_info_minped[int(channel_view0)][int(tick_view0)]
                if include_all_points_in_rectangle:
                    points_in_rectangle_view0_y.append(y_position_channel)
                    points_in_rectangle_view0_x.append(x_position_tick)
            tick_view0 += 1
        channel_view0+=1
    
    # doing the same for channel 1! Note there is a different dimensionality here! 
    # start:
    channel_view1 = convert_z_position_to_channel_number_view1(min_chan_view1_rectangle) - 1
    #channel_view1 =  0
    # end:
    last_channel =  convert_z_position_to_channel_number_view1(max_chan_view1_rectangle) + 1 
    while(channel_view1<last_channel and channel_view1 < 960 ):
        z_position_channel = convert_channel_number_view1_to_z_position(channel_view1)
        #again, we use the fits (now in view1) to get the top and bottom tick at this channel.
        # I know, the y-intercept part is weird. It's more like, the fit generator is general purpose
        # and in most cases your vertical axis is the y-axis. So the fit generator returns a "y_intercept" in its dictionary
        # and our geometry is weird, so that "y-axis" actually points along the x-axis: the channel axis. 
            
           
        tick_view1=convert_x_position_to_tick(z_position_channel*rectangle_parameters_track['lin_fit_view1']['slope'] + rectangle_parameters_track['lin_fit_view1']['y_intercept'] + 0.5*rectangle_width/np.cos(rectangle_parameters_track['lin_fit_view1']['angle_radians']) ) -1
        last_tick= convert_x_position_to_tick(z_position_channel*rectangle_parameters_track['lin_fit_view1']['slope'] + rectangle_parameters_track['lin_fit_view1']['y_intercept'] - 0.5*rectangle_width/np.cos(rectangle_parameters_track['lin_fit_view1']['angle_radians']) ) -1

        if tick_view1<0:
            tick_view1 = 0
 
 
        while(tick_view1<1667 and tick_view1<last_tick):
            x_position_tick = convert_tick_to_x_position(tick_view1)
            #print("outer if state")
            if is_point_in_rectangle(sorted_rectangle_points_1, [z_position_channel, x_position_tick], rectangle_parameters_track['lin_fit_view0'] ):
                #print("returned true")
                charge_rectangle_view1 += raw_info_minped[int(channel_view1+320)][int(tick_view1)]
                if include_all_points_in_rectangle:
                    points_in_rectangle_view1_z.append(z_position_channel)
                    points_in_rectangle_view1_x.append(x_position_tick)
            tick_view1 += 1
        channel_view1 += 1
    
    #scale the charge properly using the calibration constants:
    charge_fC_rectangle_view0 = charge_rectangle_view0 / 54.77942
    charge_fC_rectangle_view1 = charge_rectangle_view1 / 67.08538
    if include_all_points_in_rectangle:
        points_in_rectangle_view0_y = np.array(points_in_rectangle_view0_y)
        points_in_rectangle_view0_x = np.array(points_in_rectangle_view0_x)
        points_in_rectangle_view1_z = np.array(points_in_rectangle_view1_z)
        points_in_rectangle_view1_x = np.array(points_in_rectangle_view1_x)
    if include_all_points_in_rectangle:
        dictionary = {'charge_fC_rectangle_view0': charge_fC_rectangle_view0, 'charge_fC_rectangle_view1': charge_fC_rectangle_view1, 'points_in_rectangle_view0_y': points_in_rectangle_view0_y, 'points_in_rectangle_view0_x': points_in_rectangle_view0_x, 'points_in_rectangle_view1_z': points_in_rectangle_view1_z, 'points_in_rectangle_view1_x': points_in_rectangle_view1_x, 'run': run, 'subrun': subrun, 'event': event, 'track': track, 'rectangle_width': rectangle_width}
    else:
        dictionary = {'charge_fC_rectangle_view0': charge_fC_rectangle_view0, 'charge_fC_rectangle_view1': charge_fC_rectangle_view1, 'run': run, 'subrun': subrun, 'event': event, 'track': track, 'rectangle_width': rectangle_width}
    dictionary_tot = dict(dictionary.items() + rectangle_parameters_track.items())
    return dictionary_tot

def get_charge_in_multiple_rectangles(run, subrun, event, track, number_of_segments, small_rectangle_width, large_rectangle_width, include_all_points_in_rectangle=False):
    dictionaries_charge_sum_small_rectangle = []
    dictionaries_charge_sum_large_rectangle = []
    reco_hits_tot_track = get_track_2D_positions_in_each_view(run=run, subrun=subrun, event=event, track=track)
    reco_hits_segmented_track_view0_y, reco_hits_segmented_track_view0_x = segment_track_in_y(x=reco_hits_tot_track['reco_hits_y_view0'], y=reco_hits_tot_track['reco_hits_x_view0'], number_of_segments=number_of_segments)
    reco_hits_segmented_track_view1_z, reco_hits_segmented_track_view1_x = segment_track_in_y(x=reco_hits_tot_track['reco_hits_z_view1'], y=reco_hits_tot_track['reco_hits_x_view1'], number_of_segments=number_of_segments)
    for i in range(number_of_segments):
        reco_hits_track_segment = {'reco_hits_x_view0': reco_hits_segmented_track_view0_x[i], 'reco_hits_y_view0': reco_hits_segmented_track_view0_y[i], 'reco_hits_x_view1': reco_hits_segmented_track_view1_x[i], 'reco_hits_z_view1': reco_hits_segmented_track_view1_z[i]}
        small_rectangle_parameters_track_segment = get_rectangle_parameters_for_segment_of_track(run, subrun, event, track, small_rectangle_width, reco_hits_track_segment)
        large_rectangle_parameters_track_segment = get_rectangle_parameters_for_segment_of_track(run, subrun, event, track, large_rectangle_width, reco_hits_track_segment)
        charge_in_small_rectangle = sum_up_charge_in_rectangle_for_track_segment(run, subrun, event, track, small_rectangle_width, small_rectangle_parameters_track_segment, include_all_points_in_rectangle)
        charge_in_large_rectangle = sum_up_charge_in_rectangle_for_track_segment(run, subrun, event, track, large_rectangle_width, large_rectangle_parameters_track_segment, include_all_points_in_rectangle)
        dictionaries_charge_sum_small_rectangle.append(charge_in_small_rectangle)
        dictionaries_charge_sum_large_rectangle.append(charge_in_large_rectangle)
    return dictionaries_charge_sum_small_rectangle, dictionaries_charge_sum_large_rectangle

def get_charge_in_multiple_rectangles_handle_warnings(run, subrun, event, track, number_of_segments, small_rectangle_width, large_rectangle_width, include_all_points_in_rectangle=False):
    for i in np.arange(number_of_segments, 0, -1):
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                dictionaries_charge_sum_small_rectangle, dictionaries_charge_sum_large_rectangle = get_charge_in_multiple_rectangles(run, subrun, event, track, i, small_rectangle_width, large_rectangle_width, include_all_points_in_rectangle)
                break
            except (TypeError, np.RankWarning):
                print 'Fuck! Not enough hits per box for ' + str(i) + ' segments (run: ' + str(run) + ', subrun: ' + str(subrun) + ', event: ' + str(event) + ', track: ' + str(track) + '). Trying again with fewer boxes.'
    return dictionaries_charge_sum_small_rectangle, dictionaries_charge_sum_large_rectangle

def get_reconstruction_track_variables_for_purity_track_selection(run, subrun, event):
    dictionary = {}
    filename = '/eos/experiment/wa105/offline/LArSoft/Data/Parser/2018_Feb_05/' + str(run) + '/' + str(run) + '-' + str(subrun) + '-Parser.root'
    root_file = ROOT.TFile.Open(filename)
    root_file_tree = root_file.Get('analysistree/anatree')
    root_file_tree.GetEntry(event)
    
    event_time_seconds = root_file_tree.EventTimeSeconds
    dictionary['event_time_seconds'] = event_time_seconds
    event_time_utc_datetime = datetime.fromtimestamp(event_time_seconds)
    event_time_utc_datetime = event_time_utc_datetime.replace(tzinfo=tz.gettz('UTC'))
    event_time_geneva_datetime = event_time_utc_datetime.astimezone(tz.gettz('Europe/Zurich'))
    dictionary['event_time_geneva_datetime'] = event_time_geneva_datetime
    dictionary['event_time_nanoseconds'] = root_file_tree.EventTimeNanoseconds
    number_of_tracks = root_file_tree.NumberOfTracks_pmtrack
    dictionary['number_of_tracks'] = number_of_tracks
    track_id = []
    track_number_of_hits = []
    track_start_point_x = []
    track_start_point_y = []
    track_start_point_z = []
    track_end_point_x = []
    track_end_point_y = []
    track_end_point_z = []
    index_hit_inside_track = 0
    for i in range(0, number_of_tracks):
        track_id.append(root_file_tree.TrackID_pmtrack[i])
        current_track_number_of_hits = root_file_tree.Track_NumberOfHits_pmtrack[i]
        track_start_point_x.append(root_file_tree.Track_StartPoint_X_pmtrack[i])
        track_start_point_y.append(root_file_tree.Track_StartPoint_Y_pmtrack[i])
        track_start_point_z.append(root_file_tree.Track_StartPoint_Z_pmtrack[i])
        track_end_point_x.append(root_file_tree.Track_EndPoint_X_pmtrack[i])
        track_end_point_y.append(root_file_tree.Track_EndPoint_Y_pmtrack[i])
        track_end_point_z.append(root_file_tree.Track_EndPoint_Z_pmtrack[i])
        index_hit_inside_track += current_track_number_of_hits
    dictionary['track_id'] = np.array(track_id)
    dictionary['track_start_point_x'] = np.array(track_start_point_x)
    dictionary['track_start_point_y'] = np.array(track_start_point_y)
    dictionary['track_start_point_z'] = np.array(track_start_point_z)
    dictionary['track_end_point_x'] = np.array(track_end_point_x)
    dictionary['track_end_point_y'] = np.array(track_end_point_y)
    dictionary['track_end_point_z'] = np.array(track_end_point_z)
    return dictionary

def select_tracks_based_on_length_position(run, max_number_of_tracks=10, min_track_length=90.0):
    number_of_subruns = get_number_of_subruns(run=run)
    chosen_tracks = []
    for subrun in range(number_of_subruns):
        number_of_events_in_subrun = get_number_of_events_subrun(run=run, subrun=subrun)
        for event in range(number_of_events_in_subrun):
            reco_info = get_reconstruction_track_variables_for_purity_track_selection(run=run, subrun=subrun, event=event)
            if reco_info['number_of_tracks'] <= max_number_of_tracks:
                for track in range(reco_info['number_of_tracks']):
                    if np.abs(reco_info['track_start_point_x'][track] - reco_info['track_end_point_x'][track]) >= min_track_length:
                        if reco_info['track_start_point_z'][track] > 52.0 and reco_info['track_start_point_z'][track] < 248.0:
                            if reco_info['track_end_point_z'][track] > 52.0 and reco_info['track_end_point_z'][track] < 248.0:
                                chosen_tracks.append(np.array([run, subrun, event, track]))
                                print 'Possible track for Run: ' + str(run) + ', SubRun: ' + str(subrun) + ', Event: ' + str(event) + ', Track: ' + str(track)
    return np.array(chosen_tracks)

def get_box_charge_ratios_multiple_boxes(run, initial_tracks_file, small_rectangle_width, large_rectangle_width, number_of_boxes):
    initial_track_selection = np.load(initial_tracks_file)
    charge_ratios_tracks = []
    for i in initial_track_selection:
        dictionary_small_box, dictionary_large_box = get_charge_in_multiple_rectangles_handle_warnings(i[0], i[1], i[2], i[3], number_of_boxes, small_rectangle_width, large_rectangle_width)
        charge_ratios_view0 = []
        charge_ratios_view1 = []
        for j in range(len(dictionary_small_box)):
            charge_small_box_view0 = dictionary_small_box[j]['charge_fC_rectangle_view0']
            charge_small_box_view1 = dictionary_small_box[j]['charge_fC_rectangle_view1']
            charge_large_box_view0 = dictionary_large_box[j]['charge_fC_rectangle_view0']
            charge_large_box_view1 = dictionary_large_box[j]['charge_fC_rectangle_view1']
            charge_ratio_view0 = (charge_large_box_view0 - charge_small_box_view0) / charge_small_box_view0
            charge_ratio_view1 = (charge_large_box_view1 - charge_small_box_view1) / charge_small_box_view1
            charge_ratios_view0.append(charge_ratio_view0)
            charge_ratios_view1.append(charge_ratio_view1)
        charge_ratios_tracks.append({'run': i[0], 'subrun': i[1], 'event': i[2], 'track': i[3], 'charge_ratios_view0': charge_ratios_view0, 'charge_ratios_view1': charge_ratios_view1})
        print 'Charge ratios for Run ' + str(i[0]) + ', SubRun ' + str(i[1]) + ', Event ' + str(i[2]) + ', Track ' + str(i[3]) + ':   View 0 ' + str(charge_ratios_view0) + ', View 1 ' + str(charge_ratios_view1)
    return charge_ratios_tracks

#def event_display_fit_boxes(run, subrun, event, track, number_of_segments, rectangle_width_small, rectangle_width_large, clim, include_all_points_in_rectangle=False):
    
#    dict_list_small, dict_list_large = get_charge_in_multiple_rectangles_handle_warnings(run, subrun, event, track, number_of_segments, rectangle_width_small, rectangle_width_large, include_all_points_in_rectangle)
#    ADC_values = read_one_event(run, subrun, event)
#    ADC_values_minped = subtract_pedestal(ADC_values, method='median')
    
#    fig = plt.figure()
#    fig.suptitle('Run ' + str(run) + ', SubRun ' + str(subrun) + ', Event ' + str(event) + ', Track ' + str(track))
#    gs = gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[1, 3])
#    colors = ['magenta', 'black', 'cyan', 'lime', 'grey', 'mediumpurple', 'olive']
#    colors = colors = colors + colors + colors + colors + colors
    
#    plt.subplot(gs[0])
#    plt.title('View 0')
#    ax0 = plt.gca()
#    ax0.set_xlabel('Y')
#    ax0.set_ylabel('X')
#    ax0.set_xlim((-0.51, 0.51))
#    ax0.set_ylim((-0.58, 0.51))

#    img0 = ax0.imshow(np.flipud(np.rot90(ADC_values_minped[:320])), extent=(-0.5, 0.5, convert_tick_to_x_position(1666), 0.5), interpolation='nearest', cmap=plt.cm.jet, origin='upper', clim=clim)
    
#    for i in range(len(dict_list_small)):
#        ax0.plot(dict_list_small[i]['reco_hits_y_view0'], dict_list_small[i]['reco_hits_x_view0'], marker='o', linestyle='None', color=colors[i])
#        y_coordinate_extent = np.arange(np.min(dict_list_small[i]['reco_hits_y_view0']), np.max(dict_list_small[i]['reco_hits_y_view0']), 0.001)
#        ax0.plot(y_coordinate_extent, dict_list_small[i]['lin_fit_view0']['slope'] * y_coordinate_extent + dict_list_small[i]['lin_fit_view0']['y_intercept'], marker='None', color=colors[i+1])
        
#        for j in range(4):
#            ax0.plot([dict_list_small[i]['rectangle_points_view0']['rectangle_points'][j % 4][0], dict_list_small[i]['rectangle_points_view0']['rectangle_points'][(j+1) % 4][0]], [dict_list_small[i]['rectangle_points_view0']['rectangle_points'][j % 4][1], dict_list_small[i]['rectangle_points_view0']['rectangle_points'][(j+1) % 4][1]], color=colors[i+2])
#            ax0.plot([dict_list_large[i]['rectangle_points_view0']['rectangle_points'][j % 4][0], dict_list_large[i]['rectangle_points_view0']['rectangle_points'][(j+1) % 4][0]], [dict_list_large[i]['rectangle_points_view0']['rectangle_points'][j % 4][1], dict_list_large[i]['rectangle_points_view0']['rectangle_points'][(j+1) % 4][1]], color=colors[i+3])
        
#        if include_all_points_in_rectangle:
#            ax0.plot(dict_list_small[i]['points_in_rectangle_view0_y'], dict_list_small[i]['points_in_rectangle_view0_x'], marker='.', linestyle='None', color='#aec7e8')
    
#    plt.subplot(gs[1])
#    plt.title('View 1')
#    ax1 = plt.gca()
#    ax1.set_xlabel('Z')
#    ax1.set_ylabel('X')
#    ax1.set_xlim((-0.01, 3.01))
#    ax1.set_ylim((-0.58, 0.51))
    
#    img1 = ax1.imshow(np.flipud(np.rot90(ADC_values_minped[320:])), extent=(0.0, 3.0, convert_tick_to_x_position(1666), 0.5), interpolation='nearest', cmap=plt.cm.jet, origin='upper', clim=clim)
    
#    for i in range(len(dict_list_small)):
#        ax1.plot(dict_list_small[i]['reco_hits_z_view1'], dict_list_small[i]['reco_hits_x_view1'], marker='o', linestyle='None', color=colors[i])
#        y_coordinate_extent = np.arange(np.min(dict_list_small[i]['reco_hits_z_view1']), np.max(dict_list_small[i]['reco_hits_z_view1']), 0.001)
#        ax1.plot(y_coordinate_extent, dict_list_small[i]['lin_fit_view1']['slope'] * y_coordinate_extent + dict_list_small[i]['lin_fit_view1']['y_intercept'], marker='None', color=colors[i+1])
        
#        for j in range(4):
#            ax1.plot([dict_list_small[i]['rectangle_points_view1']['rectangle_points'][j % 4][0], dict_list_small[i]['rectangle_points_view1']['rectangle_points'][(j+1) % 4][0]], [dict_list_small[i]['rectangle_points_view1']['rectangle_points'][j % 4][1], dict_list_small[i]['rectangle_points_view1']['rectangle_points'][(j+1) % 4][1]], color=colors[i+2])
#            ax1.plot([dict_list_large[i]['rectangle_points_view1']['rectangle_points'][j % 4][0], dict_list_large[i]['rectangle_points_view1']['rectangle_points'][(j+1) % 4][0]], [dict_list_large[i]['rectangle_points_view1']['rectangle_points'][j % 4][1], dict_list_large[i]['rectangle_points_view1']['rectangle_points'][(j+1) % 4][1]], color=colors[i+3])
        
#        if include_all_points_in_rectangle:
#            ax1.plot(dict_list_small[i]['points_in_rectangle_view1_z'], dict_list_small[i]['points_in_rectangle_view1_x'], marker='.', linestyle='None', color='#aec7e8')
    
#    p0 = ax0.get_position().get_points().flatten()
#    p1 = ax1.get_position().get_points().flatten()
#    ax_cbar = fig.add_axes([p0[0], 0.45 * p0[2], p1[2]-p0[0], 0.02])
#    fig.colorbar(img0, cax=ax_cbar, orientation='horizontal', label='ADC counts')
    
#    charge_ratios_view0 = []
#    max_x_coordinate_boxes_view0 = []
#    charge_ratios_view1 = []
#    max_x_coordinate_boxes_view1 = []
#    for i in range(len(dict_list_small)):
#        charge_ratios_view0.append((dict_list_large[i]['charge_fC_rectangle_view0'] - dict_list_small[i]['charge_fC_rectangle_view0']) / dict_list_small[i]['charge_fC_rectangle_view0'])
#        max_x_coordinate_boxes_view0.append(np.max(dict_list_large[i]['rectangle_points_view0']['rectangle_points'][:, 1]))
#        charge_ratios_view1.append((dict_list_large[i]['charge_fC_rectangle_view1'] - dict_list_small[i]['charge_fC_rectangle_view1']) / dict_list_small[i]['charge_fC_rectangle_view1'])
#        max_x_coordinate_boxes_view1.append(np.max(dict_list_large[i]['rectangle_points_view1']['rectangle_points'][:, 1]))
    
#    inds_view0 = np.array(max_x_coordinate_boxes_view0).argsort()
#    sorted_charge_ratios_view0 = np.array(charge_ratios_view0)[inds_view0]
#    inds_view1 = np.array(max_x_coordinate_boxes_view1).argsort()
#    sorted_charge_ratios_view1 = np.array(charge_ratios_view1)[inds_view1]
    
#    print 'Charge ratios View0: ' + str(np.round(np.array(sorted_charge_ratios_view0), 3))
#    print 'Charge ratios View1: ' + str(np.round(np.array(sorted_charge_ratios_view1), 3))
    
#    plt.show()
