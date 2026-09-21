package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Annotation;
import edu.ucsd.msjava.msutil.Matter;
import edu.ucsd.msjava.msutil.Peptide;

import java.util.ArrayList;
import java.util.HashMap;

public abstract class DeNovoGraph<T extends Matter> {
    protected T source;
    protected T pmNode;
    protected ArrayList<T> sinkNodes;
    protected ArrayList<T> intermediateNodes;
    private HashMap<T, Integer> nodeIndexMap;

    public T getSource() {
        return source;
    }

    public T getPMNode() {
        return pmNode;
    }

    public ArrayList<T> getSinkList() {
        return sinkNodes;
    }

    public ArrayList<T> getIntermediateNodeList() {
        return intermediateNodes;
    }

    /**
     * A dense index of the graph's nodes, for tables indexed by node instead of keyed by node: an index in
     * 0..getNumNodeIndices()-1 for every node of the graph, or -1 for a node that cannot be in it. An index may
     * belong to no node (FlexAminoAcidGraph indexes by nominal mass, and not every mass is a node), so a table
     * indexed this way must allow for empty entries.
     */
    public int getNodeIndex(T node) {
        Integer index = getNodeIndexMap().get(node);
        return index == null ? -1 : index;
    }

    public int getNumNodeIndices() {
        return getNodeIndexMap().size();
    }

    private HashMap<T, Integer> getNodeIndexMap() {
        if (nodeIndexMap == null) {
            HashMap<T, Integer> map = new HashMap<T, Integer>();
            map.put(source, map.size());
            for (T n : intermediateNodes)
                if (!map.containsKey(n))
                    map.put(n, map.size());
            for (T n : sinkNodes)
                if (!map.containsKey(n))
                    map.put(n, map.size());
            nodeIndexMap = map;
        }
        return nodeIndexMap;
    }

    public abstract boolean isReverse();

    public abstract int getScore(Peptide pep);

    public abstract int getScore(Annotation annotation);

    public abstract int getNodeScore(T node);

    //	public abstract int getEdgeScore(T curNode, T prevNode);
    public abstract ArrayList<Edge<T>> getEdges(T curNode);

    public abstract T getComplementNode(T node);

    public abstract AminoAcidSet getAASet();

    public static class Edge<T extends Matter> {
        private T prevNode;
        private float probability;
        private int index;
        private float mass;

        // scores
        private int cleavageScore;
        private int errorScore;

        public Edge(T prevNode, float probability, int index, float mass) {
            this.prevNode = prevNode;
            this.probability = probability;
            this.index = index;
            this.mass = mass;
        }

        public T getPrevNode() {
            return prevNode;
        }

        public void setCleavageScore(int cleavageScore) {
            this.cleavageScore = cleavageScore;
        }

        public void setErrorScore(int errorScore) {
            this.errorScore = errorScore;
        }

        public void setEdgeMass(float mass) {
            this.mass = mass;
        }

        public int getEdgeScore() {
            return cleavageScore + errorScore;
        }

        public int getErrorScore() {
            return errorScore;
        }

        public float getEdgeProbability() {
            return probability;
        }

        public int getEdgeIndex() {
            return index;
        }

        public float getEdgeMass() {
            return mass;
        }
    }
}
