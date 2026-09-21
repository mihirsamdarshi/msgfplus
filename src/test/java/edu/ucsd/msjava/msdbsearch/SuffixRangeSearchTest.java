package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msscorer.NewScorerFactory;
import edu.ucsd.msjava.msscorer.SimpleDBSearchScorer;
import edu.ucsd.msjava.msutil.*;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.PriorityQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import static org.junit.Assert.*;

/**
 * A database search split into suffix ranges must leave every spectrum's match queue as a search of the whole
 * suffix array leaves it, including which of the matches tied at the worst kept score stay and in which order.
 */
public class SuffixRangeSearchTest {

    private static CompactSuffixArray sa;

    @BeforeClass
    public static void buildSuffixArray() throws Exception {
        File dbFile = new File(SuffixRangeSearchTest.class.getClassLoader().getResource("Tryp_Pig_Bov.fasta").toURI());
        sa = new CompactSuffixArray(new CompactFastaSequence(dbFile.getPath()), 12);
    }

    // Scores from a small set of values, so that most queues hold matches tied at the worst kept score
    private static class TiedScorer implements SimpleDBSearchScorer<NominalMass> {
        public int getScore(double[] prefixMassArr, int[] nominalPrefixMassArr, int fromIndex, int toIndex, int numMods) {
            return (nominalPrefixMassArr[Math.min(3, toIndex - 1)] + 3 * toIndex + numMods) % 13;
        }

        public int getNodeScore(NominalMass prefixMass, NominalMass suffixMass) {
            return 0;
        }

        public float getNodeScore(NominalMass node, boolean isPrefix) {
            return 0;
        }

        public int getEdgeScore(NominalMass curNode, NominalMass prevNode, float theoMass) {
            return 0;
        }

        public boolean getMainIonDirection() {
            return true;
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

    private static final List<SpecKey> SPEC_KEYS = Arrays.asList(new SpecKey(1, 2), new SpecKey(2, 3));

    // Two spectra near 1000 Da with a 300 Da tolerance, so that most peptides of the database are candidates
    private static DBScanner newScanner(Enzyme enzyme, int numPeptidesPerSpec) {
        ScoredSpectraMap specScanner = new ScoredSpectraMap(null, SPEC_KEYS, new Tolerance(300, false), new Tolerance(300, false), 0, 0,
                new NewScorerFactory.SpecDataType(ActivationMethod.CID, InstrumentType.LOW_RESOLUTION_LTQ, enzyme), false, false);
        for (int i = 0; i < SPEC_KEYS.size(); i++) {
            specScanner.getPepMassSpecKeyMap().put(1000.0 + i, SPEC_KEYS.get(i));
            specScanner.getSpecKeyScorerMap().put(SPEC_KEYS.get(i), new TiedScorer());
        }
        AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
        DBScanner scanner = new DBScanner(specScanner, sa, enzyme, aaSet, numPeptidesPerSpec, 3, 12, 128, 0, false, 2);
        scanner.setPrintStream(new PrintStream(new ByteArrayOutputStream()));
        return scanner;
    }

    // Every match of every queue, in the queue's iteration order
    private static String getMatches(DBScanner scanner) {
        StringBuilder matches = new StringBuilder();
        for (SpecKey specKey : SPEC_KEYS) {
            matches.append(specKey.getSpecKeyString()).append(':');
            PriorityQueue<DatabaseMatch> matchQueue = scanner.getSpecKeyDBMatchMap().get(specKey);
            if (matchQueue != null) {
                for (DatabaseMatch m : matchQueue)
                    matches.append(m.getIndex()).append('/').append(m.getLength()).append('/').append(m.getScore()).append('/')
                            .append(m.getPepSeq()).append('/').append(m.getIndices()).append('/').append(m.isProteinNTerm()).append('/')
                            .append(m.isProteinCTerm()).append('/').append(m.isNTermMetCleaved()).append(';');
            }
            matches.append('\n');
        }
        return matches.toString();
    }

    private static String searchWholeArray(Enzyme enzyme, int nnet, int numPeptidesPerSpec) {
        DBScanner scanner = newScanner(enzyme, numPeptidesPerSpec);
        scanner.dbSearch(nnet, 0, sa.getSize(), false);
        return getMatches(scanner);
    }

    private static String searchRanges(Enzyme enzyme, final int nnet, int numPeptidesPerSpec, int numRanges) throws Exception {
        final DBScanner scanner = newScanner(enzyme, numPeptidesPerSpec);
        final int[] boundaries = scanner.getSuffixRangeBoundaries(numRanges, nnet);
        assertEquals(0, boundaries[0]);
        assertEquals(sa.getSize(), boundaries[numRanges]);
        scanner.startRangeSearch(numRanges, 8);
        ExecutorService pool = Executors.newFixedThreadPool(4);
        List<Future<?>> futures = new ArrayList<Future<?>>();
        for (int i = 0; i < numRanges; i++) {
            assertTrue(boundaries[i] <= boundaries[i + 1]);
            final int rangeIndex = i;
            futures.add(pool.submit(new Runnable() {
                public void run() {
                    scanner.dbSearch(nnet, boundaries[rangeIndex], boundaries[rangeIndex + 1], false, rangeIndex);
                }
            }));
        }
        pool.shutdown();
        for (Future<?> future : futures)
            future.get(60, TimeUnit.SECONDS);
        scanner.mergeRangeDBMatches();
        return getMatches(scanner);
    }

    @Test
    public void rangesGiveTheMatchesOfTheWholeArray() throws Exception {
        for (Enzyme enzyme : new Enzyme[]{Enzyme.UnspecificCleavage, Enzyme.TRYPSIN, Enzyme.LysN}) {
            for (int nnet = 0; nnet <= 2; nnet++) {
                if (enzyme == Enzyme.UnspecificCleavage && nnet != 2)
                    continue;
                for (int numPeptidesPerSpec : new int[]{1, 3}) {
                    String expected = searchWholeArray(enzyme, nnet, numPeptidesPerSpec);
                    assertTrue("no matches to compare", expected.contains(";"));
                    // one range, a few, more ranges than the admission window, and more ranges than suffixes
                    for (int numRanges : new int[]{1, 4, 31, sa.getSize() + 3}) {
                        String message = enzyme.getName() + " nnet=" + nnet + " numPeptidesPerSpec=" + numPeptidesPerSpec + " numRanges=" + numRanges;
                        assertEquals(message, expected, searchRanges(enzyme, nnet, numPeptidesPerSpec, numRanges));
                    }
                }
            }
        }
    }

    private static Thread startRange(final DBScanner scanner, final int[] boundaries, final int rangeIndex) {
        Thread thread = new Thread(new Runnable() {
            public void run() {
                scanner.dbSearch(0, boundaries[rangeIndex], boundaries[rangeIndex + 1], false, rangeIndex);
            }
        });
        thread.start();
        return thread;
    }

    private static void waitUntilWaiting(Thread thread) throws InterruptedException {
        long deadline = System.currentTimeMillis() + 30000;
        while (thread.getState() != Thread.State.WAITING) {
            assertTrue("the range did not wait for the ranges before it", thread.isAlive() && System.currentTimeMillis() < deadline);
            Thread.sleep(10);
        }
    }

    @Test
    public void rangeWaitsUntilItIsWithinTheWindowOfTheMergedRanges() throws Exception {
        String expected = searchWholeArray(Enzyme.TRYPSIN, 0, 3);

        DBScanner scanner = newScanner(Enzyme.TRYPSIN, 3);
        int numRanges = 4;
        int[] boundaries = scanner.getSuffixRangeBoundaries(numRanges, 0);
        scanner.startRangeSearch(numRanges, 2);

        // nothing is merged yet, so ranges 0 and 1 may run and range 2 has to wait
        Thread range2 = startRange(scanner, boundaries, 2);
        waitUntilWaiting(range2);
        // range 1 finishes but cannot be merged before range 0, so range 2 still waits
        scanner.dbSearch(0, boundaries[1], boundaries[2], false, 1);
        Thread.sleep(100);
        assertEquals(Thread.State.WAITING, range2.getState());
        // range 0 finishes; ranges 0 and 1 are merged and range 2 runs
        scanner.dbSearch(0, boundaries[0], boundaries[1], false, 0);
        range2.join(30000);
        assertFalse(range2.isAlive());
        scanner.dbSearch(0, boundaries[3], boundaries[4], false, 3);
        scanner.mergeRangeDBMatches();
        assertEquals(expected, getMatches(scanner));
    }

    @Test
    public void mergeFailsWhenARangeWasInterrupted() throws Exception {
        DBScanner scanner = newScanner(Enzyme.TRYPSIN, 3);
        int numRanges = 2;
        int[] boundaries = scanner.getSuffixRangeBoundaries(numRanges, 0);
        scanner.startRangeSearch(numRanges, 1);

        Thread range1 = startRange(scanner, boundaries, 1);
        waitUntilWaiting(range1);
        range1.interrupt();
        range1.join(30000);
        assertFalse(range1.isAlive());

        scanner.dbSearch(0, boundaries[0], boundaries[1], false, 0);
        try {
            scanner.mergeRangeDBMatches();
            fail("merged although range 1 was not searched");
        } catch (IllegalStateException e) {
            // expected
        }
    }

    @Test
    public void skipFullySkipsPastTheIntRange() throws IOException {
        // the offset of suffix 536,870,910 in the index file: 2 ints of header and 4 bytes per suffix, more than 2^31
        long offset = CompactSuffixArray.INT_BYTE_SIZE * 2 + (long) CompactSuffixArray.INT_BYTE_SIZE * 536870910;
        final long[] numBytesLeft = {offset};
        InputStream in = new InputStream() {
            // skips less than asked, as a BufferedInputStream does
            public long skip(long n) {
                long skipped = Math.min(Math.min(n, 1000000000L), numBytesLeft[0]);
                numBytesLeft[0] -= skipped;
                return skipped;
            }

            public int read() {
                if (numBytesLeft[0] == 0)
                    return -1;
                numBytesLeft[0]--;
                return 0;
            }
        };
        DBScanner.skipFully(in, offset);
        assertEquals(0, numBytesLeft[0]);
    }

    @Test
    public void skipFullyReadsWhenSkipMakesNoProgressAndFailsAtTheEnd() throws IOException {
        InputStream in = new ByteArrayInputStream(new byte[7]) {
            public synchronized long skip(long n) {
                return 0;
            }
        };
        DBScanner.skipFully(in, 7);
        assertEquals(-1, in.read());
        try {
            DBScanner.skipFully(in, 1);
            fail("skipped past the end of the stream");
        } catch (EOFException e) {
            // expected
        }
    }
}
