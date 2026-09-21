package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msutil.*;
import edu.ucsd.msjava.msutil.Modification.Location;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

public class FlexAminoAcidGraph extends DeNovoGraph<NominalMass> {
    public static final int MODIFIED_EDGE_PENALTY = 0;
    private ScoredSpectrum<NominalMass> scoredSpec;
    private Enzyme enzyme;
    private boolean direction;    // true: forward (e.g. Lys-C), false: reverse (e.g. Trypsin)
    private AminoAcidSet aaSet;
    private boolean useProtNTerm;
    private boolean useProtCTerm;

    // indexed by nominal mass + nodeIndexOffset (nodes can have negative mass with negative modifications)
    private int nodeIndexOffset;
    private ArrayList<DeNovoGraph.Edge<NominalMass>>[] edgeMap;
    private int[] nodeScore;

    private static AtomicInteger negativeCompNodeMassWarnCount;
    private static AtomicInteger negativeNodeMassWarnCount;

    private static AtomicInteger nullNodeCountGetNodeScore;
    private static AtomicInteger exceptionCountGetNodeScore;

    public FlexAminoAcidGraph(
            AminoAcidSet aaSet,
            int peptideMass,
            Enzyme enzyme,
            ScoredSpectrum<NominalMass> scoredSpec
    ) {
        this(aaSet, peptideMass, enzyme, scoredSpec, false, false);
    }

    public FlexAminoAcidGraph(
            AminoAcidSet aaSet,
            int peptideMass,
            Enzyme enzyme,
            ScoredSpectrum<NominalMass> scoredSpec,
            boolean useProteinNTerm,
            boolean useProteinCTerm
    ) {
        this.enzyme = enzyme;
        this.direction = scoredSpec.getMainIonDirection();
        this.scoredSpec = scoredSpec;
        this.aaSet = aaSet;
        this.useProtNTerm = useProteinNTerm;
        this.useProtCTerm = useProteinCTerm;

        super.source = new NominalMass(0);

        super.pmNode = new NominalMass(peptideMass);

        if (negativeNodeMassWarnCount == null) {
            negativeNodeMassWarnCount = new AtomicInteger();
        }

        if (negativeCompNodeMassWarnCount == null) {
            negativeCompNodeMassWarnCount = new AtomicInteger();
        }

        if (nullNodeCountGetNodeScore == null) {
            nullNodeCountGetNodeScore = new AtomicInteger();
        }

        if (exceptionCountGetNodeScore == null) {
            exceptionCountGetNodeScore = new AtomicInteger();
        }

        int minAANominalMass = 0;
        for (Location location : Location.values()) {
            ArrayList<AminoAcid> aaList = aaSet.getAAList(location);
            if (aaList == null)
                continue;
            for (AminoAcid aa : aaList)
                if (aa.getNominalMass() < minAANominalMass)
                    minAANominalMass = aa.getNominalMass();
        }
        // A wide precursor tolerance can put the sink below zero; the arrays must still include both it and the source.
        nodeIndexOffset = -Math.min(minAANominalMass, peptideMass);
        int numNodeIndices = Math.max(peptideMass, 0) + nodeIndexOffset + 1;
        @SuppressWarnings("unchecked")
        ArrayList<DeNovoGraph.Edge<NominalMass>>[] edges = new ArrayList[numNodeIndices];
        edgeMap = edges;
        edgeMap[getNodeIndex(source)] = new ArrayList<DeNovoGraph.Edge<NominalMass>>();
        setForwardEdgesFromSource();
        setForwardEdgesFromIntermediateNodes();
        super.intermediateNodes = new ArrayList<NominalMass>();
        for (int i = 0; i < edgeMap.length; i++)    // in ascending mass order
            if (edgeMap[i] != null)
                super.intermediateNodes.add(new NominalMass(i - nodeIndexOffset));
        this.setBackwardEdgesFromSink();

        super.sinkNodes = new ArrayList<NominalMass>();
        sinkNodes.add(pmNode);

        computeNodeScores();
    }

    @Override
    public NominalMass getComplementNode(NominalMass node) {
        return new NominalMass(pmNode.getNominalMass() - node.getNominalMass());
    }

    @Override
    public int getNodeIndex(NominalMass node) {
        int index = node.getNominalMass() + nodeIndexOffset;
        if (index < 0 || index >= edgeMap.length)
            return -1;
        return index;
    }

    @Override
    public int getNumNodeIndices() {
        return edgeMap.length;
    }

    @Override
    public ArrayList<DeNovoGraph.Edge<NominalMass>> getEdges(NominalMass curNode) {
        int index = getNodeIndex(curNode);
        if (index < 0)
            return null;
        return edgeMap[index];
    }

    @Override
    public int getNodeScore(NominalMass node) {

        if (node == null) {
            int errorCount = nullNodeCountGetNodeScore.addAndGet(1);
            if (notifyError(errorCount)) {
                System.out.println("Note: null node encountered in getNodeScore");
            }
            return 0;
        }

        try {
            return nodeScore[getNodeIndex(node)];
        } catch (Exception ex) {
            int errorCount = exceptionCountGetNodeScore.addAndGet(1);
            if (notifyError(errorCount)) {
                System.out.println("Note: Exception in getNodeScore retrieving node at nominal mass " +
                        node.getNominalMass() + ": " + ex.getMessage());
            }
            return 0;
        }

    }

    @Override
    public int getScore(Peptide pep) {
        int score = 0;

        NominalMass prevNode = source;
        int nominalMass = 0;
        for (int i = 0; i < pep.size() - 1; i++) {
            AminoAcid aa;
            if (direction == true)
                aa = pep.get(i);
            else
                aa = pep.get(pep.size() - 1 - i);

            nominalMass += aa.getNominalMass();
            NominalMass curNode = new NominalMass(nominalMass);
            int nodeScore = getNodeScore(curNode);
            int edgeScore = scoredSpec.getEdgeScore(curNode, prevNode, aa.getMass());
            if (prevNode == source && direction == false && enzyme != null) {
                if (enzyme.isCleavable(aa))
                    edgeScore += aaSet.getPeptideCleavageCredit();
                else
                    edgeScore += aaSet.getPeptideCleavagePenalty();
            }
            prevNode = curNode;
            score += nodeScore + edgeScore;
        }
        if (direction == true && enzyme != null) {
            if (enzyme.isCleavable(pep.get(pep.size() - 1)))
                score += aaSet.getPeptideCleavageCredit();
            else
                score += aaSet.getPeptideCleavagePenalty();
        }
        if (direction == true)
            nominalMass += pep.get(pep.size() - 1).getNominalMass();
        else
            nominalMass += pep.get(0).getNominalMass();

        if (nominalMass != pmNode.getNominalMass())
            return Integer.MIN_VALUE;
        else
            return score;
    }

    @Override
    public int getScore(Annotation annotation) {
        int score = getScore(annotation.getPeptide());
        if (enzyme != null) {
            AminoAcid neighboringAA;
            if (enzyme.isCTerm())
                neighboringAA = annotation.getPrevAA();
            else
                neighboringAA = annotation.getNextAA();
            if (neighboringAA == null || enzyme.isCleavable(neighboringAA))
                score += aaSet.getNeighboringAACleavageCredit();
            else
                score += aaSet.getNeighboringAACleavagePenalty();
        }
        return score;
    }

    @Override
    public boolean isReverse() {
        return !direction;
    }

    @Override
    public AminoAcidSet getAASet() {
        return aaSet;
    }

    private void computeNodeScores() {
        nodeScore = new int[edgeMap.length];
        nodeScore[getNodeIndex(source)] = 0;

        boolean warnNegativeNodeMass = false;
        boolean warnNegativeCompNodeMass = false;

        for (int i = 1; i < intermediateNodes.size(); i++) {
            NominalMass node = intermediateNodes.get(i);
            NominalMass compNode = this.getComplementNode(node);
            if (node.getNominalMass() < 0 && !warnNegativeNodeMass) {
                warnNegativeNodeMass = true;
                // Mass of the node is negative
                // This can happen if we have a negative dynamic mod at the C-terminus, for example Lys-Loss
                int warnCount = negativeNodeMassWarnCount.addAndGet(1);
                if (notifyError(warnCount)) {
                    System.out.println("Note: negative node mass in computeNodeScores " +
                            "(count = " + Integer.toString(warnCount) + ")");
                }
            }
            if (compNode.getNominalMass() < 0 && !warnNegativeCompNodeMass) {
                warnNegativeCompNodeMass = true;
                int warnCount = negativeCompNodeMassWarnCount.addAndGet(1);
                if (notifyError(warnCount)) {
                    System.out.println("Note: negative compnode mass in computeNodeScores " +
                            "(count = " + Integer.toString(warnCount) + ")");
                }
            }
            int score;
            if (isReverse())
                score = scoredSpec.getNodeScore(compNode, node);
            else
                score = scoredSpec.getNodeScore(node, compNode);
            nodeScore[getNodeIndex(node)] = score;
        }
        for (NominalMass node : this.sinkNodes)
            nodeScore[getNodeIndex(node)] = 0;
    }

    private boolean notifyError(int errorCount) {
        if (errorCount < 5 || errorCount == 100 || errorCount == 1000 || errorCount % 10000 == 0) {
            return true;
        } else {
            return false;
        }
    }

    private void setForwardEdgesFromSource() {
        Location location;
        if (direction) {
            if (!this.useProtNTerm)
                location = Location.N_Term;
            else
                location = Location.Protein_N_Term;
        } else {
            if (!this.useProtCTerm)
                location = Location.C_Term;
            else
                location = Location.Protein_C_Term;
        }

        ArrayList<AminoAcid> aaList = aaSet.getAAList(location);
        makeForwardEdges(source, aaList, getAAIndices(aaList), enzyme != null && direction == enzyme.isNTerm());
    }

    private void setForwardEdgesFromIntermediateNodes() {
        ArrayList<AminoAcid> aaList = aaSet.getAAList(Location.Anywhere);
        int[] aaIndices = getAAIndices(aaList);
        for (int i = 1; i < pmNode.getNominalMass(); i++)
            makeForwardEdges(new NominalMass(i), aaList, aaIndices, false);
    }

    private int[] getAAIndices(ArrayList<AminoAcid> aaList) {
        int[] aaIndices = new int[aaList.size()];
        for (int i = 0; i < aaIndices.length; i++)
            aaIndices[i] = aaSet.getIndex(aaList.get(i));
        return aaIndices;
    }

    private void setBackwardEdgesFromSink() {
        Location location;
        if (direction) {
            if (!this.useProtCTerm)
                location = Location.C_Term;
            else
                location = Location.Protein_C_Term;
        } else {
            if (!this.useProtNTerm)
                location = Location.N_Term;
            else
                location = Location.Protein_N_Term;
        }

        ArrayList<AminoAcid> aaList = aaSet.getAAList(location);

//		if(enzymaticCleavageOnly && direction != enzyme.isNTerm())
//			aaList = aaSet.getEnzymeAAList();
//		else
//			aaList = aaSet.getAAList(location);

        int peptideNominalMass = pmNode.getNominalMass();
        ArrayList<DeNovoGraph.Edge<NominalMass>> edges = new ArrayList<DeNovoGraph.Edge<NominalMass>>();
        for (AminoAcid aa : aaList) {
            NominalMass prevNode = new NominalMass(peptideNominalMass - aa.getNominalMass());
            if (getEdges(prevNode) != null) {
                DeNovoGraph.Edge<NominalMass> edge = new DeNovoGraph.Edge<NominalMass>(prevNode, aa.getProbability(), aaSet.getIndex(aa), aa.getMass());
                edges.add(edge);
                if (enzyme != null && direction != enzyme.isNTerm()) {
                    if (enzyme.isCleavable(aa))
                        edge.setCleavageScore(aaSet.getPeptideCleavageCredit());
                    else
                        edge.setCleavageScore(aaSet.getPeptideCleavagePenalty());
                }
                if (aa.isModified())
                    edge.setErrorScore(MODIFIED_EDGE_PENALTY);
            }
        }
        edgeMap[getNodeIndex(pmNode)] = edges;
    }

    private void makeForwardEdges(NominalMass curNode, ArrayList<AminoAcid> aaList, int[] aaIndices, boolean addCleavageScore) {
        if (getEdges(curNode) == null)
            return;
        int curNominalMass = curNode.getNominalMass();
        for (int aaIndex = 0; aaIndex < aaList.size(); aaIndex++) {
            AminoAcid aa = aaList.get(aaIndex);
            int nextNodeNominalMass = curNominalMass + aa.getNominalMass();
            if (nextNodeNominalMass >= pmNode.getNominalMass())
                continue;
            NominalMass nextNode = new NominalMass(nextNodeNominalMass);
            int nextNodeIndex = getNodeIndex(nextNode);
            ArrayList<DeNovoGraph.Edge<NominalMass>> edges = edgeMap[nextNodeIndex];
            if (edges == null) {
                edges = new ArrayList<DeNovoGraph.Edge<NominalMass>>();
                edgeMap[nextNodeIndex] = edges;
            }

            DeNovoGraph.Edge<NominalMass> edge = new DeNovoGraph.Edge<NominalMass>(
                    curNode,
                    aa.getProbability(),
                    aaIndices[aaIndex],
                    aa.getMass());
//			if(curNode.getNominalMass() == 57 && nextNode.getNominalMass() == 114)
//				System.out.println("Debug");
            int errorScore = scoredSpec.getEdgeScore(nextNode, curNode, aa.getMass());
//			if(aa.isModified())
//				errorScore += MODIFIED_EDGE_PENALTY;
            if (errorScore < -100 || errorScore > 100) {
                System.err.println("Warning, invalid ErrorScore: " + errorScore);

                // Could abort the search
                // System.exit(-1);

                // Instead, use a score of -4
                errorScore = -4;
            }
            edge.setErrorScore(errorScore);
            if (addCleavageScore) {
                if (enzyme.isCleavable(aa))
                    edge.setCleavageScore(aaSet.getPeptideCleavageCredit());
                else
                    edge.setCleavageScore(aaSet.getPeptideCleavagePenalty());
            }

            edges.add(edge);
        }
    }

//	private void computeEdgeScores()
//	{
//		for(NominalMass curNode : intermediateNodes)
//		{
//			ArrayList<DeNovoGraph.Edge<NominalMass>> edges = edgeMap.get(curNode);
//			for(DeNovoGraph.Edge<NominalMass> edge : edges)
//			{
//				NominalMass prevNode = edge.getPrevNode();
//				int errorScore = scoredSpec.getEdgeScore(curNode, prevNode, edge.getEdgeMass());
//				assert(errorScore == edge.getErrorScore());
//				edge.setErrorScore(errorScore);
//			}
//			edgeMap.put(curNode, edges);
//		}
//	}	
}
