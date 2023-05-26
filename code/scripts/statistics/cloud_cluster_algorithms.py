from scipy.signal import find_peaks


# this algorithm checks if the cloud height from the return power at the next
# time / distance step is within the cluster_threshold (defaults to 250m).
# if so, include that next point in the cloud cluster. If it's too high or low,
# start a new cluster!
def clusters_within_250m( H, xaxis, cluster_threshold= .250 ):
    # save total cloud clusters, and the current cluster, in individual lists
    # do the same for x axis values
    H_clusters = []
    xaxis_clusters = []
    H_cluster_current = []
    xaxis_cluster_current = []

    # cycle through every cloud height value looking for cloud clusters
    for i in range( len( H) - 1):

        # same cluster: heights between individual CRL scans are within 250m (or a user chosen value) of each other
        if abs( H[ i] - H[ i+1]) <= cluster_threshold:

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

    return H_clusters, xaxis_clusters



# this script is attempting to preserve larger cloud clusters by implementing
# a few less restrictive cases to sort out non clusters
def large_clusters( H, xaxis, surface_threshold= .050 ):
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
            if abs( H[ i] - H[ i+1]) < .250:
                pass_case = True
                pass

            # next test:
            # iterate over the next five heights: if that average is past a certain threshold,
            # call it a new cloud
            height_avg = 0
            for j in range( 5):
                height_avg += H[ i + j]
            height_avg = height_avg / 5

            if abs( H[ i] - height_avg) > 1.0:
                pass_case=False
                print('i = ' + str( i) + ': break case from > 1 km averages at x = ' + str( xaxis[i].values))


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
                    print('i = ' + str( i) + ': break case from 3 surf points at x = ' + str( xaxis[i].values))
                    pass

        # print( " pass case: " + str( pass_case))

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


'''

# initial idea:
# maybe reverse the previous script? currently im searching for cases that break a cluster up,
# but should I instead be looking for cases that continue the current cluster?
# aka switch from finding false cases to true cases
def large_clusters_additive( H, xaxis, surface_threshold= .050, cluster_points =5):
    # save total cloud clusters, and the current cluster, in individual lists
    # do the same for x axis values
    H_clusters = []
    xaxis_clusters = []
    H_cluster_current = []
    xaxis_cluster_current = []

    # cycle through every cloud height value looking for cloud clusters
    for i in range( len( H) - 1):
        # run multiple tests to see if the next cloud height should be included in the cluster
        # the cluster needs to pass all the following tests for pass_case=True and to be included!
        pass_case = False

        # lidar reaches down to the surface for at least 3 data points in a row: no more cloud cluster.
        # make sure the last 3 indices aren't used...
        if i > len( H) - 3:
            pass
        else:
            surf_points = 0
            for j in range( 3):
                if abs( H[ i + j]) < surface_threshold:
                    surf_points += 1
            # all three points are under threshold
            if surf_points == 3:
                pass


        # make sure we aren't at the last 5 indices for the next test... things would break :/
        if i > len( H) - 5:
            pass
        else:



# check if the next cluster points are all close together



'''



# use find_peaks() algorithm to find unique cloud clusters and individual cells from
# crl cloud height data! It worked before, so why shouldn't it work now?
# Update: this function works pretty poorly lol, not all cloud clusters have unique peaks!
def find_peaks_clusters( H, xaxis, surface_threshold= .050):
    peaks = find_peaks( H, prominence= 1.0, width=3)
    return H[ peaks[0]], xaxis[ peaks[0]]

