/*
 *
 * MetaTarget ClusterVectorCB
 *
 * Copyright (C) 2022 Anestis Gkanogiannis <anestis@gkanogiannis.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */
package fr.cea.ig.metatarget.kmeans;;

public class ClusterVectorCB {

    // Indices of the member reads.
    private int[] mMemberIndexes;
    // The cluster center.
    private SequencekMeansCentroid mCenter;
    
    public ClusterVectorCB(int[] memberIndexes, SequencekMeansCentroid center) {
        mMemberIndexes = memberIndexes;
        mCenter = center;
    }
    
    public int[] getMemberIndexes() {
        return mMemberIndexes;
    }
    
    public SequencekMeansCentroid getCenter() {
        return mCenter;
    }
    
}
