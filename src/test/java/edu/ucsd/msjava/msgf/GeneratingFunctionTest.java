package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msutil.*;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.ui.MSGFPlus;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.File;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.security.MessageDigest;
import java.util.HashMap;

import static org.junit.Assert.*;

/**
 * The amino acid graph and the generating function computed on it, against digests recorded before their tables
 * were changed from hash maps to arrays (commit f0bc79e). A digest covers the nodes in order, their scores and edges,
 * and every probability and count of the score distribution.
 */
public class GeneratingFunctionTest {

    @Rule
    public TemporaryFolder temporaryFolder = new TemporaryFolder();

    // Deterministic node and edge scores, negative and positive
    private static class TestSpectrum implements ScoredSpectrum<NominalMass> {
        private final boolean mainIonDirection;

        TestSpectrum(boolean mainIonDirection) {
            this.mainIonDirection = mainIonDirection;
        }

        public int getNodeScore(NominalMass prefixMass, NominalMass suffixMass) {
            return Math.floorMod(prefixMass.getNominalMass() * 7 + suffixMass.getNominalMass(), 11) - 5;
        }

        public float getNodeScore(NominalMass node, boolean isPrefix) {
            return 0;
        }

        public int getEdgeScore(NominalMass curNode, NominalMass prevNode, float theoMass) {
            return Math.floorMod(curNode.getNominalMass() + prevNode.getNominalMass(), 5) - 2;
        }

        public boolean getMainIonDirection() {
            return mainIonDirection;
        }

        public Peak getPrecursorPeak() {
            return new Peak(500, 1, 2);
        }

        public ActivationMethod[] getActivationMethodArr() {
            return new ActivationMethod[]{ActivationMethod.CID};
        }

        public int[] getScanNumArr() {
            return new int[]{1};
        }
    }

    private static final int NO_THRESHOLD = Integer.MIN_VALUE;

    private static DeNovoGraph<NominalMass> getGraph(AminoAcidSet aaSet, boolean generic, boolean direction, int nominalPeptideMass) {
        if (generic)
            return new GenericDeNovoGraph<NominalMass>(new NominalMassFactory(aaSet, direction ? Enzyme.LysN : Enzyme.TRYPSIN, 20),
                    nominalPeptideMass + (float) Composition.H2O, new Tolerance(1.5f, false), null, new TestSpectrum(direction));
        return new FlexAminoAcidGraph(aaSet, nominalPeptideMass, null, new TestSpectrum(direction), true, true);
    }

    private static void update(MessageDigest digest, Object value) {
        digest.update((value + "\n").getBytes(StandardCharsets.UTF_8));
    }

    // Both ion directions, three peptide masses, a table capacity that evicts at every node, one that evicts at some
    // and one that never does, and no score threshold, a low one and one that prunes
    private static String getDigest(AminoAcidSet aaSet, boolean generic) throws Exception {
        return getDigest(aaSet, generic, new int[]{300, 650, 1200});
    }

    private static String getDigest(AminoAcidSet aaSet, boolean generic, int[] nominalPeptideMasses) throws Exception {
        MessageDigest digest = MessageDigest.getInstance("SHA-256");
        for (boolean direction : new boolean[]{false, true}) {
            for (int nominalPeptideMass : nominalPeptideMasses) {
                for (int capacity : new int[]{1, 100, 10000}) {
                    for (int threshold : new int[]{NO_THRESHOLD, 0, 30}) {
                        DeNovoGraph<NominalMass> graph = getGraph(aaSet, generic, direction, nominalPeptideMass);
                        for (NominalMass node : graph.getIntermediateNodeList()) {
                            update(digest, node.getNominalMass());
                            update(digest, graph.getNodeScore(node));
                            for (DeNovoGraph.Edge<NominalMass> edge : graph.getEdges(node))
                                update(digest, edge.getPrevNode().getNominalMass() + "/" + edge.getEdgeScore() + "/" + edge.getEdgeIndex() + "/" + edge.getEdgeProbability());
                        }
                        GeneratingFunction<NominalMass> gf = new GeneratingFunction<NominalMass>(graph).enzyme(null).gfTableCapacity(capacity);
                        if (threshold != NO_THRESHOLD)
                            gf.setUpScoreThreshold(threshold);
                        boolean isComputed = gf.computeGeneratingFunction();
                        update(digest, isComputed);
                        if (isComputed) {
                            ScoreDist dist = gf.getScoreDist();
                            update(digest, dist.getMinScore());
                            update(digest, dist.getMaxScore());
                            for (int score = dist.getMinScore(); score < dist.getMaxScore(); score++)
                                update(digest, Double.doubleToLongBits(dist.getProbability(score)) + "/" + Double.doubleToLongBits(dist.getNumberRecs(score)));
                        }
                    }
                }
            }
        }
        StringBuilder hex = new StringBuilder();
        for (byte b : digest.digest())
            hex.append(String.format("%02x", b));
        return hex.toString();
    }

    // A modification at the protein N-terminus, one at the peptide C-terminus that removes its residue, and one
    // heavier than its residue, which gives the graph nodes of negative mass
    private static AminoAcidSet getAASetWithTerminalAndNegativeMods() throws Exception {
        File modFile = File.createTempFile("GeneratingFunctionTestMods", ".txt");
        modFile.deleteOnExit();
        PrintWriter writer = new PrintWriter(modFile, "UTF-8");
        writer.println("NumMods=3");
        writer.println("57.021464,C,fix,any,Carbamidomethyl");
        writer.println("15.994915,M,opt,any,Oxidation");
        writer.println("42.010565,*,opt,Prot-N-term,Acetyl");
        writer.println("-128.094963,K,opt,C-term,Lys-Loss");
        writer.println("-190.164215,K,opt,any,NegativeMass");
        writer.close();
        ParamManager paramManager = new ParamManager("MS-GF+", "test", "test", "test");
        paramManager.addMSGFPlusParams();
        return AminoAcidSet.getAminoAcidSetFromModFile(modFile.getPath(), paramManager);
    }

    @Test
    public void negativeSinkMasses() throws Exception {
        int[] masses = {-100, -62, -61, -1, 0, 1};
        assertEquals("57e84754ad8c893e9a3c64340c9c0340505718440ce5e210e0dc205473c24458", getDigest(AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys(), false, masses));
        assertEquals("328e0eb866481ff65b96eeb0767c8a6a38e1ad50566a281b44fd1f9ea279b835", getDigest(getAASetWithTerminalAndNegativeMods(), false, masses));
    }

    @Test
    public void precursorWindowCrossingZeroCompletes() throws Exception {
        File spectra = temporaryFolder.newFile("wide.mgf");
        try (PrintWriter writer = new PrintWriter(spectra, "UTF-8")) {
            writer.println("BEGIN IONS\nTITLE=wide\nPEPMASS=100\nCHARGE=2+\nSCANS=1");
            for (int mass : new int[]{40, 51, 63, 78, 87, 95, 105, 115, 122, 134, 142, 155, 162, 177, 188})
                writer.println(mass + " 50");
            writer.println("END IONS");
        }
        File database = temporaryFolder.newFile("wide.fasta");
        try (PrintWriter writer = new PrintWriter(database, "UTF-8")) {
            writer.println(">protein\nAAAAKR");
        }
        File result = new File(temporaryFolder.getRoot(), "wide.mzid");
        ParamManager params = new ParamManager("MS-GF+", "test", "test", "test");
        params.addMSGFPlusParams();
        assertNull(params.parseParams(new String[]{"-s", spectra.getPath(), "-d", database.getPath(),
                "-o", result.getPath(), "-t", "200Da", "-minLength", "3", "-maxLength", "6",
                "-e", "0", "-ntt", "0", "-thread", "1", "-tasks", "1"}));
        assertNull(MSGFPlus.runMSGFPlus(params));
        assertTrue(new String(Files.readAllBytes(result.toPath()), StandardCharsets.UTF_8).contains("<SpectrumIdentificationItem"));
    }

    @Test
    public void aminoAcidGraph() throws Exception {
        assertEquals("80bb1f6b1fe5c86c9d267d03ee8831e15131bc75a3b52d44ccf39863ac856e97", getDigest(AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys(), false));
    }

    @Test
    public void aminoAcidGraphWithTerminalAndNegativeMods() throws Exception {
        assertEquals("f1952729902519ffd937eedbccedc760196a5eebcd4fd3afbff07ff927f17e60", getDigest(getAASetWithTerminalAndNegativeMods(), false));
    }

    // GenericDeNovoGraph does not override the node index that the generating function's tables use
    @Test
    public void genericGraph() throws Exception {
        assertEquals("82548a5a3c5b497e1422681a0b3a0f69cd7b6ce82f6da129aa28a8e9cdc1580f", getDigest(AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys(), true));
    }

    // ProfileGF asks for the forward table at every node, so the table must not be rebuilt at every call;
    // computing the generating function again must replace it
    @Test
    public void forwardTableIsKeptUntilTheNextComputation() {
        AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
        // ProfileGF fails on negative node scores and on edge scores
        DeNovoGraph<NominalMass> graph = new FlexAminoAcidGraph(aaSet, 650, null, new TestSpectrum(true) {
            public int getNodeScore(NominalMass prefixMass, NominalMass suffixMass) {
                return 0;
            }

            public int getEdgeScore(NominalMass curNode, NominalMass prevNode, float theoMass) {
                return 0;
            }
        }, true, true);
        GeneratingFunction<NominalMass> gf = new GeneratingFunction<NominalMass>(graph).enzyme(null);
        assertTrue(gf.computeGeneratingFunction());

        HashMap<NominalMass, ScoreDist> fwdTable = gf.getFwdTable();
        assertSame(fwdTable, gf.getFwdTable());
        assertNotNull(fwdTable.get(graph.getSource()));
        for (NominalMass sink : graph.getSinkList())
            assertNotNull(fwdTable.get(sink));
        int numNodesWithDist = 0;
        for (NominalMass node : graph.getIntermediateNodeList()) {
            if (fwdTable.get(node) != null)
                numNodesWithDist++;
        }
        assertTrue(numNodesWithDist > 0);

        float sumProbabilities = new ProfileGF<NominalMass>(gf).computeProfile(0).getSpectralProfile().getSumProbabilities();
        assertTrue(sumProbabilities > 0);
        assertSame(fwdTable, gf.getFwdTable());

        assertTrue(gf.computeGeneratingFunction());
        assertNotSame(fwdTable, gf.getFwdTable());
        assertEquals(sumProbabilities, new ProfileGF<NominalMass>(gf).computeProfile(0).getSpectralProfile().getSumProbabilities(), 0);
    }
}
