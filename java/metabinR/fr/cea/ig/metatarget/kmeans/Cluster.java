package fr.cea.ig.metatarget.kmeans;;

public class Cluster {

    // Indices of the member reads.
    private int[] mMemberIndexes;
    // The cluster center.
    private ReadkMeansCentroid mCenter;
    
    public Cluster(int[] memberIndexes, ReadkMeansCentroid center) {
        mMemberIndexes = memberIndexes;
        mCenter = center;
    }
    
    public int[] getMemberIndexes() {
        return mMemberIndexes;
    }
    
    public ReadkMeansCentroid getCenter() {
        return mCenter;
    }
    
}
