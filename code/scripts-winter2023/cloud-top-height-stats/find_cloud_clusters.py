from scipy.signal import find_peaks
import numpy as np
import xarray as xr
import pandas as pd

# this algorithm checks if the cloud height from the return power at the next
# time / distance step is within the cluster_threshold (defaults to 250m).
# if so, include that next point in the cloud cluster. If it's too high or low,
# start a new cluster!
def clusters_within_250m( H, xaxis, cluster_threshold= 250, surface_threshold=50):
    # save total cloud clusters, and the current cluster, in individual lists
    # do the same for x axis values
    H_clusters = []
    xaxis_clusters = []
    H_cluster_current = []
    xaxis_cluster_current = []

    # cycle through every cloud height value looking for cloud clusters
    for i in range( len( H) - 1):

        # lidar reaches down to the surface for at least 3 data points in a row: no more cloud cluster
        if i > len( H) - 3:
            pass
        else:
            surf_points = 0
            for j in range( 3):
                if abs( H[ i + j]) < surface_threshold:
                    surf_points += 1

        # same cluster: heights between individual CRL scans are within 250m (or a user chosen value) of each other
        if abs( H[ i] - H[ i+1]) <= cluster_threshold and surf_points < 3:

            # append values to the current cluster list
            H_cluster_current.append( H[ i])
            xaxis_cluster_current.append( float( xaxis[ i].values))

        # else case: start of a new cluster
        else:
            # append current cloud height to this cluster: it's the last valid value!
            H_cluster_current.append( H[ i])
            xaxis_cluster_current.append( float( xaxis[ i].values))

            # append current cloud clusters to total list
            H_clusters.append(  H_cluster_current)
            xaxis_clusters.append(  xaxis_cluster_current)
            # clear current cluster lists
            H_cluster_current = []
            xaxis_cluster_current = []

        # check if there are any last clusters to append at the end of the for loop!
        if i == len( H) - 1:
            # make sure the current cluster has values!
            if len( H_cluster_current) > 0:
                H_clusters.append(  H_cluster_current)
                xaxis_clusters.append(  xaxis_cluster_current)
    return H_clusters, xaxis_clusters




# this algorithm is a slight improvement on the previous one. It's still based on clustering using a 
# cluster_threshold height difference, but clear air can now break up a cluster, and some smaller clusters 
# are combined if their mean (or adjacent?) cloud heights are close enough?
def clusters_within_threshold_updated( H, xaxis, cluster_threshold= 250, surface_threshold=50, min_size=False, combine_clusters=False):
    # save total cloud clusters, and the current cluster, in individual lists
    # do the same for x axis values
    H_clusters = []
    xaxis_clusters = []
    H_cluster_current = []
    xaxis_cluster_current = []

    # cycle through every cloud height value looking for cloud clusters
    for i in range( len( H) -1 ):
        # the lidar doesn't reach down to the surface- only add a nan to the dataset!
        if abs( H[ i]) < surface_threshold:
            # append current cloud clusters to total list, if present
            if len(H_cluster_current) > 0:
                H_clusters.append(  H_cluster_current)
                xaxis_clusters.append(  xaxis_cluster_current)
                # clear current cluster lists
                H_cluster_current = []
                xaxis_cluster_current = []

            # append a nan to this cluster and end it
            H_cluster_current.append( np.nan)
            xaxis_cluster_current.append( float( xaxis[ i-1].values))
            H_clusters.append(  H_cluster_current)
            xaxis_clusters.append(  xaxis_cluster_current)
            H_cluster_current = []
            xaxis_cluster_current = []

        # same cluster: heights between individual CRL scans are within 250m (or a user chosen value) of each other
        elif abs( H[ i] - H[ i+1]) <= cluster_threshold:
            # append values to the current cluster list
            H_cluster_current.append( H[ i])
            xaxis_cluster_current.append( float( xaxis[ i].values))

        # else case: start of a new cluster
        else:
            # append current cloud height to this cluster: it's the last valid value!
            H_cluster_current.append( H[ i])
            xaxis_cluster_current.append( float( xaxis[ i].values))

            # append current cloud clusters to total list
            H_clusters.append(  H_cluster_current)
            xaxis_clusters.append(  xaxis_cluster_current)
            # clear current cluster lists
            H_cluster_current = []
            xaxis_cluster_current = []

        # check if there are any last clusters to append at the end of the for loop!
        if i == len( H) - 2:
            # make sure the current cluster has values!
            if len( H_cluster_current) > 0:
                H_clusters.append(  H_cluster_current)
                xaxis_clusters.append(  xaxis_cluster_current)

    # after running the first algorithm to create basic cloud clusters, use the flags below to 
    # add additional complexity to improve results!

    # optional code to romove clusters that are too thin... must be 2 data points ~= 520m to count?
    # or 4 data points ~= 1040m? Can specify the width using min_size
    if min_size:
        clusters_final, xaxis_final = [], []
        for clusteri, cluster in enumerate( H_clusters):
            if len(cluster) >= min_size:
                clusters_final.append( cluster)
                xaxis_final.append( xaxis_clusters[clusteri])
        
        # just use the min_size sorting case
        if not combine_clusters:
            return clusters_final, xaxis_final
        # otherwise, do the following...
        # part two: try to combine larger clusters with small gaps between them? optional.
        # do this after removing small clusters: these will be factored into the cloud gaps.
        else:
            # save final return clusters here
            new_clusters = []
            new_xaxis_vals = []
            # if add flag is true, the current cluster has already been added. false -> not added
            add_flag = False

            # first, find where clusters are and aren't present
            # make a list like [ 0 0 1 1 1 1 1 0 0 2 0 3 3 3 ] representing where clusters are and aren't present
            # aka a list with 0s at non clusters, higher numbers for cluster order
            cluster_order_list = np.zeros(len(xaxis)) # initialize
            # search through each x axis cluster
            for timei, timeval in enumerate( xaxis_final):
                # indices where cloud clusters are present on x axis!
                cluster_inds =  np.where(np.in1d( xaxis.values, np.array(timeval)))[0] 
                cluster_order_list[ cluster_inds] = timei + 1

            # testing
            # print( len( xaxis.values))
            # print(len(xaxis_flat))
            # print( len(np.where(np.in1d( xaxis.values, np.array(xaxis_flat)))[0]))    
            # print( np.where(np.in1d( xaxis.values, np.array(xaxis_flat)))[0])
            print( cluster_order_list)

            # then, get the indices for the start of each cluster ( like [2 9 11] for example above)
            # AND the end of each cluster ( [6 9 13] above)
            # add flag starts false
            cluster_start, cluster_end = [], [] # store beginning and end points here
            current_start, current_end = 1, 1 # current starting and ending indices to search for
            for xi in range(len(cluster_order_list)-1):
                xval = cluster_order_list[xi]
                xnext = cluster_order_list[xi + 1]
                # found a new cluster
                if xval == current_start:
                    cluster_start.append(xi)
                    current_start += 1
                # found end of a cluster- next value is different
                # make sure current value isn't 0 either! this will cause an infinite cluster break loop lol
                if xval != 0 and xnext != current_end:
                    cluster_end.append(xi)
                    current_end += 1
                # otherwise, just continue!
            # testing
            # print(cluster_start)
            # print(cluster_end)

            # cycle through all of the clusters except the last
            for i in range( len(cluster_start) - 1):

                # set up cluster indices here!
                # + 1 includes last element in cluster... kinda annoying for other index cases, though :/
                starti, endi = cluster_start[i], cluster_end[i] + 1 
                nextstarti, nextendi = cluster_start[i+1], cluster_end[i+1]+1
                # slice the cluster count array here for the current cluster, the clear air / cloud gap, and the next cluster
                current_cluster = cluster_order_list[starti : endi]
                if endi == nextstarti:
                    gap = [] # no gap
                else:
                    gap = cluster_order_list[ endi : nextstarti] # find gap
                next_cluster = cluster_order_list[nextstarti : nextendi]
                
                # height values for the current and next clusters!!
                current_vals = H[ starti:endi]
                next_vals = H[nextstarti:nextendi]

                # testing
                # print(current_cluster)
                # print(gap)
                # print(next_cluster)
                # print("\n")

                # Next, go to the current cluster and check if the cluster is long enough... 
                # and there's a small gap between them (< min_size)...
                # and the next cloud heights have similar means. 
                if len( current_cluster) > min_size and len( gap) < min_size and (np.abs( np.nanmean(current_vals) - np.nanmean(next_vals)) <= cluster_threshold):              
                    
                    print( np.abs( np.nanmean(current_vals) - np.nanmean(next_vals)))
                    # If so, check add flag
                    # true -> current cluster added already, so just add small space and next cluster to last created cluster in total cluster list.
                    # false -> add current cluster, small space, and next cluster to new cluster list.
                    # change add flag to true. this way, we know that this cluster has already been added.
                    # gapnans signify the outlier clouds / clear air being skipped over. Maybe put the actual heights here instead? 
                    if add_flag:
                        gapnans = np.empty(len(gap))
                        gapnans[:] = np.nan
                        # gapnans = gapnans.tolist()
                        new_clusters[-1] = np.concatenate( [ new_clusters[-1], np.concatenate( [gapnans, H[nextstarti:nextendi] ]) ])
                        new_xaxis_vals[-1] = np.concatenate( [ new_xaxis_vals[-1], xaxis[starti : nextendi].values ])
                    else:
                        gapnans = np.empty(len(gap))
                        gapnans[:] = np.nan
                        # gapnans = gapnans.tolist()
                        new_clusters.append( np.concatenate( [ H[starti : endi], gapnans, H[nextstarti:nextendi] ] ))
                        new_xaxis_vals.append( xaxis[starti : nextendi].values)
                    add_flag = True

                # if not, just add the current cluster to the list and set add flag to false
                # the gap space will be skipped over here, and we'll just look at the next cluster
                else:
                    new_clusters.append( H[starti : endi])
                    new_xaxis_vals.append( xaxis[starti : endi].values)
                    add_flag = False

            # Old ideas: don't do this
            # if the next cluster is the last cluster (starting at 11 in example above):
            # if gap found sucessfully, just break the loop
            # if not, add the last cluster separately and break
            # not last cluster case
            # reset i to next cluster and continue.
            return new_clusters, new_xaxis_vals

    # original original case: just return the simplest algorithm output!
    else:
        return H_clusters, xaxis_clusters

    # 5/31 new code:
    # example: with the clouds shown below, combine the lower 3 clusters into 1 big cluster!
    #       --       --            
    #                         -- 
    #                       --  -------  ----------
    # ------  -------  -----           --
    # 
    # the algorithm searches cluster by cluster. If there's a pass case, distance representing all 3
    # clusters is added to the final list. Then, 


    '''
    H_clusters_len = []
    for cluster in H_clusters:
        H_clusters_len += [len(cluster)]
    # print(H_clusters)
    # print(H_clusters_len)
    # print(np.where(np.array(H_clusters_len) > 2)[0])

    clusters_final, xaxis_final = [], []

    # just iterate using a while loop! keep track of indices manually. Nice bc we can iterate by more than one
    clusteri = 0 
    while True:
        # end of indices case! 
        if clusteri >= len(H_clusters) - 2:
            break

        # nice local definitions
        thiscluster = H_clusters[clusteri]
        middlecluster = H_clusters[clusteri+1]
        nextcluster = H_clusters[clusteri+2]

        # small cloud case: just append a new cluster and continue
        if H_clusters_len[clusteri] < 3:
            clusters_final.append(H_clusters[clusteri])
            xaxis_final.append(xaxis_clusters[clusteri])

        # large cloud case: multiple routes
        else:
            # make sure there's a small cluster in between and that the nextcluster is large!
            if len(middlecluster) <= 3 and len(nextcluster) >=4:
                # finally, make sure that the heights are similar!
                    if np.abs( np.nanmean(thiscluster) - np.nanmean(nextcluster)) <= cluster_threshold:
                    # if all these tests pass, add it to new list

            else:
    '''


    '''
    # 5/31 new code: check how much space exists between clusters. If there's not that much space (1 or 2 spaces), 
    # and the mean heights of the two clusters are close, combine the clusters! pad the distance between them with nans
    for clusteri in range(len(H_clusters)-2):
        thiscluster = H_clusters[clusteri]
        middlecluster = H_clusters[clusteri+1]
        nextcluster = H_clusters[clusteri+2]

        # do the test here! make sure clusters are wide enough
        if len(thiscluster) >= 4 and len(nextcluster) >= 4:
            # make sure there's a small cluster in between
            if len(middlecluster) <= 3:
                # make sure the heights are similar!
                if np.abs( np.nanmean(thiscluster) = np.nanmean(nextcluster)) <= cluster_threshold:
                    # if all these tests pass, add it to new list
    

        # one of the cases above fail: just add the current cluster to the list and continue

        else:
            clusters_final.append(thiscluster)
            xaxis_final.append(xaxis_clusters[clusteri])
            continue
    '''



# new algorithm: apply a low pass filter to the heights to find the dominant clusters without the high
# frequency noise. 
def clusters_low_pass( H, xaxis, surface_threshold=50, prominence=500, width=4):
    pass



# this code works but not very well lol.
# this script is attempting to preserve larger cloud clusters by implementing
# a few less restrictive cases to sort out non clusters
def large_clusters( H, xaxis, surface_threshold=50):
    # save total cloud clusters, and the current cluster, in individual lists
    # do the same for x axis values
    H_clusters = []
    xaxis_clusters = []
    H_cluster_current = []
    xaxis_cluster_current = []

    # cycle through every cloud height value looking for cloud clusters
    for i in range( len( H) - 1):
        pass_case = True

        # run multiple tests to see if the next cloud height should be included in the cluster
        # if just one of these cases fail, set pass = False and start a new cluster

        # make sure we aren't at the last 5 indices for the next test... things would break :/
        if i > len( H) - 5:
            pass
        else:
            # first, check if the next cloud top is close to the current one... if so,
            # skip this test completely
            # within 250 m of each other?
            if abs( H[ i] - H[ i+1]) < 250:
                pass_case = True
                pass

            # next test:
            # iterate over the next five heights: if that average is past a certain threshold,
            # call it a new cloud
            height_avg = 0
            for j in range( 5):
                height_avg += H[ i + j]
            height_avg = height_avg / 5

            if abs( H[ i] - height_avg) > 1000:
                pass_case=False
                # print('i = ' + str( i) + ': break case from > 1 km averages at x = ' + str( xaxis[i].values))


        # lidar reaches down to the surface for at least 3 data points in a row: no more cloud cluster
        if i > len( H) - 3:
            pass
        else:
            surf_points = 0
            for j in range( 3):
                if abs( H[ i + j]) < surface_threshold:
                    surf_points += 1

            # all three points are under threshold
            if surf_points == 3:
                    pass_case = False
                    # print('i = ' + str( i) + ': break case from 3 surf points at x = ' + str( xaxis[i].values))
                    pass

        # case where all the clusters pass the test!!
        if pass_case:
            # append values to the current cluster list
            H_cluster_current.append( H[ i])
            xaxis_cluster_current.append( float( xaxis[ i].values))
            # print('i = ' + str( i) + ': pass case at x = ' + str( xaxis[i].values))


        # else case: start of a new cluster
        else:
            # append current cloud height to this cluster: it's the last valid value!
            H_cluster_current.append( H[ i])
            xaxis_cluster_current.append( float( xaxis[ i].values))

            # append current cloud clusters to total list
            H_clusters.append(  H_cluster_current)
            xaxis_clusters.append(  xaxis_cluster_current)
            # clear current cluster lists
            H_cluster_current = []
            xaxis_cluster_current = []

        # check if there are any last clusters to append at the end of the for loop!
        if i == len( H) - 2:
            # make sure the current cluster has values!
            if len( H_cluster_current) > 0:
                H_clusters.append(  H_cluster_current)
                xaxis_clusters.append(  xaxis_cluster_current)

    return H_clusters, xaxis_clusters


# code below isn't working :/
# this algorithm finds peaks using multiple cloud height data! 
# This will hopefully allow for finding clouds beneath initial layers
def multi_clusters_within_250m( H, xaxis, surface_threshold=50, cluster_threshold=250, prominence=500):
    # This code uses find_multi_cloud_heights() to get cloud heights.
    # currently, there are len(cloudtime) lists stored in H. Each list has at least one
    # cloud height, or a nan if too close to the surface. There can be multiple cloud heights per list.
    # there is one time value in xaxis for each height list: len(H) = len(xaxis)

    # goal: initialize a cluster matrix to hold cloud clusters from multiple levels! There are a max of 8 clusters with this method.
    H_clusters = np.empty( (len(xaxis), 8))   
    H_clusters[:,:] = np.nan

    # cycle through every x axis position looking for cloud clusters
    for i in range( len( H) - 1):

        # base case: establish new clusters
        if i == 0:
            # cycle through every initial height and add to empty matrix
            for hti in range(len( H[i])):
                current_H = H[i][hti]
                H_clusters[i, hti] = current_H

        # recursive case: find the nearest cluster for this height!
        else:
            # the last values in each current cluster
            # values are removed from here after a new height is added! no duplicate cases
            cluster_ends = H_clusters[i-1, :] 
        
            # cycle through every height found at this x axis position (0 to 8 heights possible)
            for hti in range(len( H[i])):
                # define the current cloud height

                # no more clusters case: make a new one! 
                if len( np.where(np.isnan(cluster_ends))[0]) == 0:
                    print("No valid remaining clusters! Adding a new one")
                    # search previous cluster ends for the first nan -> cluster break. add height next to this one
                    for j in range(len(8)):
                        if H_clusters[i-1, j] == np.nan:
                            H_clusters[i, j] = current_H
                            continue

                current_H = H[i][hti]
                nearest_ht_i = np.nanargmin( np.abs(cluster_ends - current_H))
                nearest_ht_difference = np.nanmin( np.abs(cluster_ends - current_H)) # cluster_ends[ nearest_ht_i]

                # see if the nearest cluster end is within the specified cluster threshold
                # if so, add to cluster array
                if nearest_ht_difference <= cluster_threshold:
                    print('pass!')
                    H_clusters[i, nearest_ht_i] = current_H

                # testing
                if i < 10:
                    print("time " + str(i) + ", height " + str(hti))
                    print(H_clusters[0:20, :])
                    # print("\n" + str(i))
                    # print(cluster_ends)
                    # print(current_H)
                    # print(nearest_ht_difference)


    return H_clusters


# maybe try checking if for the next five or ten data points, if four or eight of them
# are all above / below a threshold, and those four or eight are all similar heights,
# maybe call that a new cloud cluster?
# this kinda combines the ideas from the 250m and larger scripts! It makes sure that the
# next potential cluster has similar heights (with wiggle room for a couple outliers) like
# the 250m script, but it also looks at more than just the next point, like the larger cluster
# script
# or, make a script that ignores smaller jagged clouds in cluster calculations altogether??
# in some ways, they aren't true clusters at all- just noisy
# or, should I just filter these clusters out at the end?
# initial idea:
# maybe reverse the previous script? currently im searching for cases that break a cluster up,
# but should I instead be looking for cases that continue the current cluster?
# aka switch from finding false cases to true cases