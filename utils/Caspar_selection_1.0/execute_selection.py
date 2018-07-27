from mip_track_selection_functions import *

run = 840
max_number_of_tracks_in_event = 10
min_track_length_x_direction = 90.0
initial_track_selection = select_tracks_based_on_length_position(run = run, max_number_of_tracks = max_number_of_tracks_in_event, min_track_length = min_track_length_x_direction)
np.save('selected_tracks/initial_track_selection_run' + str(run) + '.npy', initial_track_selection)
