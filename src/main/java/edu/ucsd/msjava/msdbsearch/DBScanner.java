package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.misc.ProgressData;
import edu.ucsd.msjava.msgf.*;
import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.msscorer.SimpleDBSearchScorer;
import edu.ucsd.msjava.msutil.*;
import edu.ucsd.msjava.msutil.Modification.Location;
import edu.ucsd.msjava.parser.BufferedLineReader;
import edu.ucsd.msjava.sequences.Constants;

import java.io.*;
import java.util.*;
import java.util.Map.Entry;

public class DBScanner {

    private int minPeptideLength;
    private int maxPeptideLength;
    private int maxMissedCleavages;

    /**
     * Number of isoforms to consider per peptide.
     * NUM_VARIANTS_PER_PEPTIDE is 128 in Constants.java
     */
    private int maxNumVariantsPerPeptide;

    private AminoAcidSet aaSet;
    private double[] aaMass;
    private int[] intAAMass;

    private Enzyme enzyme;
    private int numPeptidesPerSpec;

    private final CompactSuffixArray sa;
    private final int size;
    // to scan the database partially
    // Input spectra
    private final ScoredSpectraMap specScanner;

    private int minDeNovoScore;
    private boolean ignoreNTermMetCleavage;

    // DB search results
    private Map<SpecKey, PriorityQueue<DatabaseMatch>> specKeyDBMatchMap;
    private Map<Integer, PriorityQueue<DatabaseMatch>> specIndexDBMatchMap;

    // When the database is searched by suffix range: the matches of the ranges that have finished but cannot be
    // merged yet because an earlier range is still running, per spectrum key, in the order the walk created them;
    // the number of ranges merged so far; and the merged queues (see addRangeDBMatches)
    private List<Map<SpecKey, List<DatabaseMatch>>> rangeMatchLists;
    private int numMergedRanges;
    private int maxRangesAhead;
    private Map<SpecKey, PriorityQueue<DatabaseMatch>> mergedRangeDBMatchMap;

    // Unsynchronized copies of the spectrum maps for the search, which only reads them; the synchronized
    // wrappers in ScoredSpectraMap would serialize suffix ranges searched in parallel on the same scanner.
    // The peptide mass map is kept as sorted arrays, so that a mass window is found by binary search
    // (the masses are positive and never NaN, so primitive comparison orders them as Double.compareTo did).
    private double[] searchPepMassKeys;
    private List<SpecKey> searchPepMassSpecKeyList;
    private Map<SpecKey, SimpleDBSearchScorer<NominalMass>> searchSpecKeyScorerMap;

    private ProgressData progress;
    private PrintStream output;

    // For output
    private String threadName = "";

    public DBScanner(
            ScoredSpectraMap specScanner,
            CompactSuffixArray sa,
            Enzyme enzyme,
            AminoAcidSet aaSet,
            int numPeptidesPerSpec,
            int minPeptideLength,
            int maxPeptideLength,
            int maxNumVariantsPerPeptide,
            int minDeNovoScore,
            boolean ignoreNTermMetCleavage,
            int maxMissedCleavages
    ) {
        this.specScanner = specScanner;
        this.sa = sa;
        this.size = sa.getSize();
        this.aaSet = aaSet;
        this.enzyme = enzyme;
        this.numPeptidesPerSpec = numPeptidesPerSpec;
        this.minPeptideLength = minPeptideLength;
        this.maxPeptideLength = maxPeptideLength;
        this.maxMissedCleavages = maxMissedCleavages;
        this.maxNumVariantsPerPeptide = maxNumVariantsPerPeptide;
        this.minDeNovoScore = minDeNovoScore;
        this.ignoreNTermMetCleavage = ignoreNTermMetCleavage;

        // Initialize mass arrays for a faster search
        aaMass = new double[aaSet.getMaxResidue()];
        intAAMass = new int[aaSet.getMaxResidue()];
        for (int i = 0; i < aaMass.length; i++) {
            aaMass[i] = -1;
            intAAMass[i] = -1;
        }
        for (AminoAcid aa : aaSet.getAllAminoAcidArr()) {
            aaMass[aa.getResidue()] = aa.getAccurateMass();
            intAAMass[aa.getResidue()] = aa.getNominalMass();
        }

        specKeyDBMatchMap = Collections.synchronizedMap(new HashMap<SpecKey, PriorityQueue<DatabaseMatch>>());
        specIndexDBMatchMap = Collections.synchronizedMap(new HashMap<Integer, PriorityQueue<DatabaseMatch>>());

        progress = null;
        output = System.out;
    }

    // builder
    public DBScanner maxPeptideLength(int maxPeptideLength) {
        this.maxPeptideLength = maxPeptideLength;
        return this;
    }

    // builder
    public DBScanner minPeptideLength(int minPeptideLength) {
        if (minPeptideLength > 1)
            this.minPeptideLength = minPeptideLength;
        else
            minPeptideLength = 1;
        return this;
    }

    public DBScanner setThreadName(String threadName) {
        this.threadName = threadName;
        return this;
    }

    public synchronized void addDBMatches(Map<SpecKey, PriorityQueue<DatabaseMatch>> map) {
        if (map == null)
            return;
        Iterator<Entry<SpecKey, PriorityQueue<DatabaseMatch>>> itr = map.entrySet().iterator();
        while (itr.hasNext()) {
            Entry<SpecKey, PriorityQueue<DatabaseMatch>> entry = itr.next();
            SpecKey specKey = entry.getKey();
            PriorityQueue<DatabaseMatch> queue = entry.getValue();

            PriorityQueue<DatabaseMatch> existingQueue = specKeyDBMatchMap.get(entry.getKey());
            if (existingQueue == null) {
                existingQueue = new PriorityQueue<DatabaseMatch>();
                specKeyDBMatchMap.put(specKey, existingQueue);
            }
            existingQueue.addAll(queue);
        }
    }

    /**
     * Add a match to a per-spectrum queue, keeping the best numPeptidesPerSpec scores and every match tied with the worst kept score.
     * Matches scoring below the worst kept score are dropped.
     */
    private void addMatch(PriorityQueue<DatabaseMatch> matchQueue, DatabaseMatch match) {
        int score = match.getScore();
        if (matchQueue.size() < this.numPeptidesPerSpec || score == matchQueue.peek().getScore()) {
            matchQueue.add(match);
        } else {
            int worstScore = matchQueue.peek().getScore();
            if (score > worstScore) {
                List<DatabaseMatch> removed = new ArrayList<DatabaseMatch>();
                while (!matchQueue.isEmpty() && matchQueue.peek().getScore() == worstScore) {
                    removed.add(matchQueue.poll());
                }
                matchQueue.add(match);

                if (matchQueue.size() < this.numPeptidesPerSpec) {
                    for (DatabaseMatch m : removed)
                        matchQueue.add(m);
                }
            }
        }
    }

    private synchronized void copySpectrumMaps() {
        if (searchPepMassKeys == null) {
            SortedMap<Double, SpecKey> pepMassSpecKeyMap = specScanner.getPepMassSpecKeyMap();
            synchronized (pepMassSpecKeyMap) {
                double[] keys = new double[pepMassSpecKeyMap.size()];
                SpecKey[] specKeys = new SpecKey[keys.length];
                int i = 0;
                for (Entry<Double, SpecKey> entry : pepMassSpecKeyMap.entrySet()) {
                    keys[i] = entry.getKey();
                    specKeys[i] = entry.getValue();
                    i++;
                }
                searchPepMassSpecKeyList = Arrays.asList(specKeys);
                searchPepMassKeys = keys;
            }
            searchSpecKeyScorerMap = new HashMap<SpecKey, SimpleDBSearchScorer<NominalMass>>(specScanner.getSpecKeyScorerMap());
        }
    }

    /**
     * @return the index of the first key that is not less than mass, or keys.length if there is none
     */
    private static int firstKeyNotLessThan(double[] keys, double mass) {
        int low = 0;
        int high = keys.length;
        while (low < high) {
            int mid = (low + high) >>> 1;
            if (keys[mid] < mass)
                low = mid + 1;
            else
                high = mid;
        }
        return low;
    }

    /**
     * Prepare for searching the database as numRanges suffix ranges, each with dbSearch(..., rangeIndex),
     * followed by mergeRangeDBMatches. At most maxRangesAhead ranges beyond the merged ones are searched or
     * held for merging at a time, which bounds the matches kept aside.
     */
    public void startRangeSearch(int numRanges, int maxRangesAhead) {
        rangeMatchLists = new ArrayList<Map<SpecKey, List<DatabaseMatch>>>(Collections.<Map<SpecKey, List<DatabaseMatch>>>nCopies(numRanges, null));
        numMergedRanges = 0;
        this.maxRangesAhead = maxRangesAhead;
        mergedRangeDBMatchMap = new HashMap<SpecKey, PriorityQueue<DatabaseMatch>>();
    }

    /**
     * Ranges start in suffix-array order, so the first unmerged range is always being searched and this cannot deadlock.
     */
    private synchronized void waitUntilRangeMayStart(int rangeIndex) throws InterruptedException {
        while (rangeIndex >= numMergedRanges + maxRangesAhead)
            wait();
    }

    /**
     * Take the matches of a finished suffix range. The ranges are merged in suffix-array order, each as soon as
     * every range before it has finished, so the matches of only a few ranges are held at any time.
     * The matches are replayed in the order a search of the whole array would have created them, through the
     * same queue rule, so the merged queues end up identical to a whole-array search.
     */
    private synchronized void addRangeDBMatches(int rangeIndex, Map<SpecKey, List<DatabaseMatch>> matchListMap) {
        rangeMatchLists.set(rangeIndex, matchListMap);
        while (numMergedRanges < rangeMatchLists.size() && rangeMatchLists.get(numMergedRanges) != null) {
            for (Entry<SpecKey, List<DatabaseMatch>> entry : rangeMatchLists.get(numMergedRanges).entrySet()) {
                SpecKey specKey = entry.getKey();
                PriorityQueue<DatabaseMatch> matchQueue = mergedRangeDBMatchMap.get(specKey);
                if (matchQueue == null) {
                    matchQueue = new PriorityQueue<DatabaseMatch>();
                    mergedRangeDBMatchMap.put(specKey, matchQueue);
                }
                for (DatabaseMatch match : entry.getValue())
                    addMatch(matchQueue, match);
            }
            rangeMatchLists.set(numMergedRanges, null);
            numMergedRanges++;
        }
        notifyAll();
    }

    static void skipFully(InputStream in, long numBytes) throws IOException {
        while (numBytes > 0) {
            long skipped = in.skip(numBytes);
            if (skipped <= 0) {
                if (in.read() < 0)
                    throw new EOFException();
                skipped = 1;
            }
            numBytes -= skipped;
        }
    }

    /**
     * Split the suffix array into numRanges ranges of similar size.
     * Every range starts at a suffix that shares at most minPeptideLength residues with the suffix before it,
     * so no peptide is shared across a boundary and searching the ranges separately visits exactly the peptides
     * a search of the whole array visits.
     *
     * @return numRanges+1 boundaries; range i is [boundaries[i], boundaries[i+1])
     */
    public int[] getSuffixRangeBoundaries(int numRanges) throws IOException {
        int[] boundaries = new int[numRanges + 1];
        boundaries[numRanges] = size;
        try (DataInputStream nlcps = new DataInputStream(new BufferedInputStream(new FileInputStream(sa.getNeighboringLcpFile())))) {
            // skip size and id
            skipFully(nlcps, CompactSuffixArray.INT_BYTE_SIZE * 2);
            int position = 0;    // the next suffix whose neighboring LCP will be read from the file
            for (int i = 1; i < numRanges; i++) {
                int target = (int) ((long) size * i / numRanges);
                if (target < position)
                    target = position;
                skipFully(nlcps, target - position);
                position = target;
                // Move on to the first suffix that shares at most minPeptideLength residues with its predecessor: dbSearch
                // resets the LCP at the start of a range, so a candidate of minPeptideLength or more (plus its preceding
                // residue) can never span the boundary
                while (position < size) {
                    int lcp = nlcps.readByte();
                    position++;
                    if (lcp <= minPeptideLength) {
                        break;
                    }
                }
                // position is one past the suffix just read; if the file ended before such a suffix, the range is empty
                boundaries[i] = position < size ? position - 1 : size;
            }
        }
        return boundaries;
    }

    /**
     * Add the merged matches of the suffix ranges searched with dbSearch(..., rangeIndex) to specKeyDBMatchMap,
     * once every range has finished.
     */
    public void mergeRangeDBMatches() {
        if (numMergedRanges < rangeMatchLists.size())
            throw new IllegalStateException("Only " + numMergedRanges + " of " + rangeMatchLists.size() + " suffix ranges finished");
        rangeMatchLists = null;
        this.addDBMatches(mergedRangeDBMatchMap);
        mergedRangeDBMatchMap = null;
    }

    public Map<SpecKey, PriorityQueue<DatabaseMatch>> getSpecKeyDBMatchMap() {
        return specKeyDBMatchMap;
    }

    public Map<Integer, PriorityQueue<DatabaseMatch>> getSpecIndexDBMatchMap() {
        return specIndexDBMatchMap;
    }

    public void setProgressObj(ProgressData progObj) {
        progress = progObj;
    }

    public ProgressData getProgressObj() {
        return progress;
    }

    public void setPrintStream(PrintStream out) {
        if (out == null) {
            output = System.out;
        } else {
            output = out;
        }
    }

    public PrintStream getPrintStream() {
        return output;
    }

    public void dbSearchCTermEnzymeNoMod(int numberOfAllowableNonEnzymaticTermini, boolean verbose) {
        dbSearch(numberOfAllowableNonEnzymaticTermini, 0, size, verbose);
    }

    public void dbSearchCTermEnzyme(int numberOfAllowableNonEnzymaticTermini, boolean verbose) {
        dbSearch(numberOfAllowableNonEnzymaticTermini, 0, size, verbose);
    }

    public void dbSearchNTermEnzyme(int numberOfAllowableNonEnzymaticTermini, boolean verbose) {
        dbSearch(numberOfAllowableNonEnzymaticTermini, 0, size, verbose);
    }

    public void dbSearchNoEnzyme(boolean verbose) {
        dbSearch(2, 0, size, verbose);
    }

    public void dbSearch(int numberOfAllowableNonEnzymaticTermini) {
        dbSearch(numberOfAllowableNonEnzymaticTermini, 0, size, true);
    }

    public void dbSearch(int numberOfAllowableNonEnzymaticTermini, int fromIndex, int toIndex, boolean verbose) {
        dbSearch(numberOfAllowableNonEnzymaticTermini, fromIndex, toIndex, verbose, -1);
    }

    /**
     * Search the suffixes [fromIndex, toIndex) of the database.
     *
     * @param rangeIndex when the search is split into suffix ranges (see startRangeSearch), the index of this range;
     *                   the matches are then kept aside for mergeRangeDBMatches instead of being added to the results.
     */
    public void dbSearch(int numberOfAllowableNonEnzymaticTermini, int fromIndex, int toIndex, boolean verbose, int rangeIndex) {
        if (progress == null) {
            progress = new ProgressData();
        }
        if (rangeIndex >= 0) {
            try {
                waitUntilRangeMayStart(rangeIndex);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                return;
            }
        }

        Map<SpecKey, PriorityQueue<DatabaseMatch>> curSpecKeyDBMatchMap = new HashMap<SpecKey, PriorityQueue<DatabaseMatch>>();
        // in key order of first match, so that the merged map has the same iteration order as a whole-array search
        Map<SpecKey, List<DatabaseMatch>> curMatchListMap = rangeIndex >= 0 ? new LinkedHashMap<SpecKey, List<DatabaseMatch>>() : null;

        CandidatePeptideGrid candidatePepGrid;
        if (enzyme != null && !ignoreNTermMetCleavage)
            candidatePepGrid = new CandidatePeptideGridConsideringMetCleavage(aaSet, enzyme, maxPeptideLength, maxNumVariantsPerPeptide, maxMissedCleavages);
        else
            candidatePepGrid = new CandidatePeptideGrid(aaSet, enzyme, maxPeptideLength, maxNumVariantsPerPeptide, maxMissedCleavages);

        int peptideLengthIndex = Integer.MAX_VALUE - 1000;

        boolean enzymaticSearch;
        enzymaticSearch = numberOfAllowableNonEnzymaticTermini != 2;

        int neighboringAACleavageCredit = aaSet.getNeighboringAACleavageCredit();
        int neighboringAACleavagePenalty = aaSet.getNeighboringAACleavagePenalty();
        int peptideCleavageCredit = aaSet.getPeptideCleavageCredit();
        int peptideCleavagePenalty = aaSet.getPeptideCleavagePenalty();

        boolean containsCTermMod = aaSet.containsCTermModification();

        copySpectrumMaps();
        double[] pepMassKeys = this.searchPepMassKeys;
        List<SpecKey> pepMassSpecKeyList = this.searchPepMassSpecKeyList;
        Map<SpecKey, SimpleDBSearchScorer<NominalMass>> specKeyScorerMap = this.searchSpecKeyScorerMap;

        try {
            DataInputStream indices = new DataInputStream(new BufferedInputStream(new FileInputStream(sa.getIndexFile())));

            // skip size and id
            skipFully(indices, CompactSuffixArray.INT_BYTE_SIZE * 2 + (long) CompactSuffixArray.INT_BYTE_SIZE * fromIndex);

            DataInputStream nlcps = new DataInputStream(new BufferedInputStream(new FileInputStream(sa.getNeighboringLcpFile())));

            // skip size
            skipFully(nlcps, CompactSuffixArray.INT_BYTE_SIZE * 2 + (long) fromIndex);
            CompactFastaSequence sequence = sa.getSequence();

            boolean isProteinNTerm = true;
            int nTermCleavageScore = 0;

            boolean isExtensionAtTheSameIndex;

            // number of non-enzymatic termini
            int numNonEnzTermini = 0;

            int numIndices = toIndex - fromIndex;

            class MatchList extends ArrayList<DatabaseMatch> {
                private static final long serialVersionUID = 1L;
            }
            MatchList[] prevMatchList = new MatchList[maxPeptideLength + 2];

            for (int bufferIndex = 0; bufferIndex < numIndices; bufferIndex++) {
                // Print out the progress
                if (verbose && bufferIndex % 2000000 == 0) {
                    output.print(threadName + ": Database search progress... ");
                    output.format("%.1f%% complete\n", bufferIndex / (float) numIndices * 100);
                }
                progress.report(bufferIndex, numIndices);
                isExtensionAtTheSameIndex = false;
                int index = indices.readInt();
                int lcp = nlcps.readByte();
                if (bufferIndex == 0)
                    lcp = 0;

                // For debugging
//				System.out.println(index+": " +sequence.getSubsequence(index, sequence.getSize()));
//				if(index == 4)
//					System.out.println("Debug");
                // skip redundant peptides

                if (Thread.currentThread().isInterrupted()) {
                    return;
                }

                // lcp: shared prefix length
                for (int peptideLength = minPeptideLength; peptideLength < prevMatchList.length; peptideLength++) {
                    if (Thread.currentThread().isInterrupted()) {
                        return;
                    }

                    if (lcp >= peptideLength + 2)    // peptide, N-term, C-term are shared
                    {
                        if (prevMatchList[peptideLength] != null) {
                            for (DatabaseMatch m : prevMatchList[peptideLength]) {
                                m.addIndex(index);
                            }
                        }
                    } else if (lcp == peptideLength + 1) {
                        if (prevMatchList[peptideLength] != null) {
                            for (DatabaseMatch m : prevMatchList[peptideLength]) {
                                if (Thread.currentThread().isInterrupted()) {
                                    return;
                                }

                                if (!m.isProteinCTerm() || enzyme == null || enzyme.isNTerm() || numberOfAllowableNonEnzymaticTermini == 2) {
                                    m.addIndex(index);
                                    continue;
                                }

                                char pre = sequence.getCharAt(index);
                                if (numberOfAllowableNonEnzymaticTermini == 1 && enzyme.isCleavable(pre)) {
                                    m.addIndex(index);
                                    continue;
                                }

                                // C-term should be enzymatic
                                char cTermResidue = sequence.getCharAt(index + peptideLength);
                                if (enzyme.isCleavable(cTermResidue)) {
                                    m.addIndex(index);
                                    continue;
                                }

                                // post should be protein c term
                                char post = sequence.getCharAt(index + peptideLength + 1);
                                if (post == Constants.TERMINATOR_CHAR) {
                                    m.addIndex(index);
                                }
                            }
                        }
                    } else
                        prevMatchList[peptideLength] = null;
                }

                if (lcp >= peptideLengthIndex + 2 ||
                        lcp == peptideLengthIndex + 1 && (enzyme == null || enzyme.isCTerm())) {
                    continue;
                }
                else if (lcp == 0)    // preceding aa is changed
                {
                    char precedingAA = sequence.getCharAt(index);
                    isProteinNTerm = precedingAA == Constants.TERMINATOR_CHAR;

                    // determine neighboring N-term score
                    if (enzyme == null || enzyme.isNTerm()) {
                        nTermCleavageScore = 0;
                    } else if (enzyme.isCTerm()) {
                        if (isProteinNTerm || enzyme.isCleavable(precedingAA))// || precedingAA == Constants.INVALID_CHAR)
                        {
                            nTermCleavageScore = neighboringAACleavageCredit;
                            if (enzymaticSearch)
                                numNonEnzTermini = 0;
                        } else {
                            nTermCleavageScore = neighboringAACleavagePenalty;
                            if (enzymaticSearch) {
                                numNonEnzTermini = 1;
                                if (numNonEnzTermini > numberOfAllowableNonEnzymaticTermini) {
                                    peptideLengthIndex = 0;
                                    continue;
                                }
                            }
                        }
                    }
                }    // end lcp=0

                if (lcp == 0)
                    peptideLengthIndex = 1;
                    //else if(lcp < peptideLengthIndex + 1)
                else {
                    if (enzyme != null && enzyme.isNTerm()) {
                        if (lcp > 1)
                            peptideLengthIndex = lcp - 1;
                        else
                            peptideLengthIndex = 1;
                    } else {
                        peptideLengthIndex = lcp;
                    }
                }

                for (; peptideLengthIndex <= maxPeptideLength && index + peptideLengthIndex < size - 1; peptideLengthIndex++)    // ith character of a peptide
                {
                    if (Thread.currentThread().isInterrupted()) {
                        return;
                    }

                    char residue = sequence.getCharAt(index + peptideLengthIndex);
                    boolean isProteinCTerm = false;
                    if (peptideLengthIndex == 1)    // N-term residue
                    {
                        if (enzyme != null && enzyme.isNTerm()) {
                            if (isProteinNTerm || enzyme.isCleavable(residue))    // || sequence.getCharAt(index) == Constants.INVALID_CHAR)
                            {
                                nTermCleavageScore = peptideCleavageCredit;
                                if (enzymaticSearch)
                                    numNonEnzTermini = 0;
                            } else {
                                nTermCleavageScore = peptideCleavagePenalty;
                                if (enzymaticSearch) {
                                    numNonEnzTermini = 1;
                                    if (numNonEnzTermini > numberOfAllowableNonEnzymaticTermini)
                                        break;
                                }
                            }
                        }

                        if (isProteinNTerm) {
                            if (candidatePepGrid.addProtNTermResidue(residue) == false)
                                break;
                        } else {
                            if (candidatePepGrid.addNTermResidue(residue) == false)
                                break;
                        }
                    } else {
                        if (!containsCTermMod) {
                            if (candidatePepGrid.addResidue(peptideLengthIndex, residue) == false)
                                break;
                        } else {
                            if (peptideLengthIndex < minPeptideLength) {
                                if (candidatePepGrid.addResidue(peptideLengthIndex, residue) == false)
                                    break;
                                else
                                    continue;
                            } else {
                                if (isExtensionAtTheSameIndex && peptideLengthIndex > minPeptideLength)
                                    candidatePepGrid.addResidue(peptideLengthIndex - 1, sequence.getCharAt(index + peptideLengthIndex - 1));
                                boolean success;
                                if (isProteinCTerm = (sequence.getCharAt(index + peptideLengthIndex + 1) == Constants.TERMINATOR_CHAR))    // protein C-term
                                    success = candidatePepGrid.addProtCTermResidue(peptideLengthIndex, residue);
                                else    // peptide C-term
                                    success = candidatePepGrid.addCTermResidue(peptideLengthIndex, residue);
                                if (!success)
                                    break;
                            }
                        }
                    }

                    if (peptideLengthIndex < minPeptideLength)
                        continue;

//					System.out.println(sequence.getSubsequence(index+1, index+i+1));
//					if(sequence.getSubsequence(index+1, index+i+1).equalsIgnoreCase("KYPCRYCEK"))
//					{
//						System.out.println("DebugSequence: " + sequence.getSubsequence(index, index+i+1));
//					}

                    int cTermCleavageScore = 0;
                    if (enzyme != null) {
                        char cTermNeighboringResidue = sequence.getCharAt(index + peptideLengthIndex + 1);
                        isProteinCTerm = (cTermNeighboringResidue == Constants.TERMINATOR_CHAR);
                        if (enzyme.isCTerm()) {
//							if(isProteinCTerm || enzyme.isCleavable(residue)) // || cTermNeighboringResidue == Constants.INVALID_CHAR)
                            if (enzyme.isCleavable(residue)) // || cTermNeighboringResidue == Constants.INVALID_CHAR)	// changed by Sangtae to avoid SpecProb=0
                                cTermCleavageScore = peptideCleavageCredit;
                            else {
                                cTermCleavageScore = peptideCleavagePenalty;
                                if (!isProteinCTerm && numNonEnzTermini + 1 > numberOfAllowableNonEnzymaticTermini) {
                                    isExtensionAtTheSameIndex = true;
                                    continue;
                                }
                            }
                        } else if (enzyme.isNTerm()) {
                            if (isProteinCTerm || enzyme.isCleavable(cTermNeighboringResidue)) // || cTermNeighboringResidue == Constants.INVALID_CHAR)
                                cTermCleavageScore = neighboringAACleavageCredit;
                            else {
                                cTermCleavageScore = neighboringAACleavagePenalty;
                                if (numNonEnzTermini + 1 > numberOfAllowableNonEnzymaticTermini) {
                                    isExtensionAtTheSameIndex = true;
                                    continue;
                                }
                            }
                        }
                    }

                    int cleavageScore = nTermCleavageScore + cTermCleavageScore;

                    for (int j = 0; j < candidatePepGrid.size(); j++) {
                        if (Thread.currentThread().isInterrupted()) {
                            return;
                        }

                        /*
                         * Check for edge case where peptides derived from the
                         * start of a protein sequence containing an N-terminus
                         * methionine may have more missed cleavages than the
                         * peptides derived from removing the methionine when
                         * digesting with N-term enzymes.
                         *
                         * E.g., a grid that considers methionine cleavage on
                         * protein sequence 'MDT' will return peptides
                         * ['MDT','DT']. If we are using AspN as the enzyme
                         * the MDT peptide has one missed cleavage and the DT
                         * peptide has zero. We want to skip the peptides that
                         * are over the maximum number of missed cleavages.
                         *
                         */
                        if (candidatePepGrid.gridIsOverMaxMissedCleavages(j))
                            continue;

                        float theoPeptideMass = candidatePepGrid.getPeptideMass(j);
//						/// Debug
//						System.out.println("PepStr: " + candidatePepGrid.getPeptideSeq(j) + " GridSize:" + candidatePepGrid.size());
//						///
                        int nominalPeptideMass = candidatePepGrid.getNominalPeptideMass(j);
                        float tolDaLeft = specScanner.getLeftPrecursorMassTolerance().getToleranceAsDa(theoPeptideMass);
                        float tolDaRight = specScanner.getRightPrecursorMassTolerance().getToleranceAsDa(theoPeptideMass);

                        double leftThr = (double) (theoPeptideMass - tolDaLeft);
                        double rightThr = (double) (theoPeptideMass + tolDaRight);

//						float tolDaLeft = specScanner.getLeftPrecursorMassTolerance().getToleranceAsDa(peptideMass);
//						float tolDaRight = specScanner.getRightPrecursorMassTolerance().getToleranceAsDa(peptideMass);
//						int maxPeptideMassIndex, minPeptideMassIndex;
//						
//						maxPeptideMassIndex = maxNominalPeptideMass + Math.round(tolDaLeft-0.4999f);
//						minPeptideMassIndex = minNominalPeptideMass - Math.round(tolDaRight-0.4999f);

                        if (leftThr < 1 || rightThr < 1) {
                            // Either or both of the thresholds is less than 1 (and probably negative)
                            // This can happen when a dynamic mod with a large negative mass is defined and is applied to a small peptide

                            // For example:
                            //  DynamicMod=304.207146,  *,  opt, N-term,    TMTpro                # 16-plex TMT
                            //  DynamicMod=304.207146,  K,  opt, any,       TMTpro                # 16-plex TMT
                            //  DynamicMod=-190.164215, K,  opt, any,       UbNoTMT16             # Residue tagged by MS-GF+ with TMT16, but is actually ubiquitinated and does not have TMT (+114.042931 - 304.207146)
                            continue;
                        }

                        Collection<SpecKey> matchedSpecKeyList = pepMassSpecKeyList.subList(firstKeyNotLessThan(pepMassKeys, leftThr), firstKeyNotLessThan(pepMassKeys, rightThr));
                        if (matchedSpecKeyList.size() > 0) {
                            //////
//							System.out.println("\tMatch: " + sequence.getCharAt(index)+"."+sequence.getSubsequence(index+1, index+i+1)+"."+sequence.getCharAt(index+i+1));
                            ///////
                            boolean isNTermMetCleaved = candidatePepGrid.isNTermMetCleaved(j);
//							int pepLength = i;
                            int pepLength;
                            if (!isNTermMetCleaved)
                                pepLength = peptideLengthIndex;
                            else
                                pepLength = peptideLengthIndex - 1;

                            if (pepLength < minPeptideLength)
                                continue;

                            for (SpecKey specKey : matchedSpecKeyList) {
                                if (Thread.currentThread().isInterrupted()) {
                                    return;
                                }

//								Tolerance specSpecificTol;
//								if((specSpecificTol = specScanner.getSpectrumSpecificPrecursorTolerance(specKey)) != null)
//								{
//								}

                                SimpleDBSearchScorer<NominalMass> scorer = specKeyScorerMap.get(specKey);
//								if(sequence.getSubsequence(index, index+i+1).equalsIgnoreCase("SRDTAIKT"))
//									System.out.println("Debug");
                                int score = cleavageScore + scorer.getScore(candidatePepGrid.getPRMGrid(j), candidatePepGrid.getNominalPRMGrid(j), 1, pepLength + 1, candidatePepGrid.getNumMods(j));
                                PriorityQueue<DatabaseMatch> prevMatchQueue = curSpecKeyDBMatchMap.get(specKey);
                                if (prevMatchQueue == null) {
                                    prevMatchQueue = new PriorityQueue<DatabaseMatch>();
                                    curSpecKeyDBMatchMap.put(specKey, prevMatchQueue);
                                }

                                if (prevMatchQueue.size() < this.numPeptidesPerSpec || score >= prevMatchQueue.peek().getScore()) {
                                    DatabaseMatch dbMatch = new DatabaseMatch(index, (byte) (pepLength + 2), score, theoPeptideMass, nominalPeptideMass, specKey.getCharge(), candidatePepGrid.getPeptideSeq(j), scorer.getActivationMethodArr()).setProteinNTerm(isProteinNTerm).setProteinCTerm(isProteinCTerm);
                                    dbMatch.setNTermMetCleaved(isNTermMetCleaved);
                                    addMatch(prevMatchQueue, dbMatch);
                                    if (curMatchListMap != null) {
                                        List<DatabaseMatch> matchList = curMatchListMap.get(specKey);
                                        if (matchList == null) {
                                            matchList = new ArrayList<DatabaseMatch>();
                                            curMatchListMap.put(specKey, matchList);
                                        }
                                        matchList.add(dbMatch);
                                    }
                                    if (prevMatchList[peptideLengthIndex] == null)
                                        prevMatchList[peptideLengthIndex] = new MatchList();
                                    prevMatchList[peptideLengthIndex].add(dbMatch);
                                }
                            }
                        }
                    }
                    isExtensionAtTheSameIndex = true;
                }
            }
            if (rangeIndex >= 0)
                addRangeDBMatches(rangeIndex, curMatchListMap);
            else
                this.addDBMatches(curSpecKeyDBMatchMap);
            indices.close();
            nlcps.close();
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }

    public void computeSpecEValue(boolean storeScoreDist) {
        computeSpecEValue(storeScoreDist, 0, specScanner.getSpecKeyList().size());
    }

    public void computeSpecEValue(boolean storeScoreDist, int fromIndex, int toIndex) {
        if (progress == null) {
            progress = new ProgressData();
        }
        List<SpecKey> specKeyList = specScanner.getSpecKeyList().subList(fromIndex, toIndex);

        int numSpecs = toIndex - fromIndex;
        int numProcessedSpecs = 0;
        for (SpecKey specKey : specKeyList) {
            if (Thread.currentThread().isInterrupted()) {
                return;
            }
            numProcessedSpecs++;
            if (numProcessedSpecs % 1000 == 0) {
                output.print(threadName + ": Computing spectral E-values... ");
                output.format("%.1f%% complete\n", numProcessedSpecs / (float) numSpecs * 100);
            }
            progress.report(numProcessedSpecs, numSpecs);

            PriorityQueue<DatabaseMatch> matchQueue = specKeyDBMatchMap.get(specKey);
            if (matchQueue == null)
                continue;

            int specIndex = specKey.getSpecIndex();

            boolean useProtNTerm = false;
            boolean useProtCTerm = false;
            int minScore = Integer.MAX_VALUE;
            for (DatabaseMatch m : matchQueue) {
                if (m.isProteinNTerm())
                    useProtNTerm = true;
                if (m.isProteinCTerm())
                    useProtCTerm = true;
                if (m.getScore() < minScore)
                    minScore = m.getScore();
            }

            GeneratingFunctionGroup<NominalMass> gf = new GeneratingFunctionGroup<NominalMass>();
            SimpleDBSearchScorer<NominalMass> scoredSpec = specScanner.getSpecKeyScorerMap().get(specKey);
            float peptideMass = scoredSpec.getPrecursorPeak().getMass() - (float) Composition.H2O;
            int nominalPeptideMass = NominalMass.toNominalMass(peptideMass);
            int minNominalPeptideMass = nominalPeptideMass - specScanner.getMaxIsotopeError();
            int maxNominalPeptideMass = nominalPeptideMass - specScanner.getMinIsotopeError();

            // Debug
//			System.out.println("***Debug " + threadName + ": " + scoredSpec.getPrecursorPeak().getMz() + "," + scoredSpec.getPrecursorPeak().getCharge());
            ///

            float tolDaLeft = specScanner.getLeftPrecursorMassTolerance().getToleranceAsDa(peptideMass);
            float tolDaRight = specScanner.getRightPrecursorMassTolerance().getToleranceAsDa(peptideMass);
            int maxPeptideMassIndex, minPeptideMassIndex;

            maxPeptideMassIndex = maxNominalPeptideMass + Math.round(tolDaLeft - 0.4999f);
            minPeptideMassIndex = minNominalPeptideMass - Math.round(tolDaRight - 0.4999f);

            for (int peptideMassIndex = minPeptideMassIndex; peptideMassIndex <= maxPeptideMassIndex; peptideMassIndex++) {
                DeNovoGraph<NominalMass> graph = new FlexAminoAcidGraph(
                        aaSet,
                        peptideMassIndex,
                        enzyme,
                        scoredSpec,
                        useProtNTerm,
                        useProtCTerm
                );

                GeneratingFunction<NominalMass> gfi = new GeneratingFunction<NominalMass>(graph)
                        .doNotBacktrack()
                        .doNotCalcNumber();
                gfi.setUpScoreThreshold(minScore);
                gf.registerGF(graph.getPMNode(), gfi);
            }

            boolean isGFComputed = gf.computeGeneratingFunction();

            for (DatabaseMatch match : matchQueue) {
                if (!isGFComputed || match.getNominalPeptideMass() < minPeptideMassIndex || match.getNominalPeptideMass() > maxPeptideMassIndex) {
                    match.setDeNovoScore(Integer.MIN_VALUE);
                    match.setSpecProb(1);
                } else {
                    match.setDeNovoScore(gf.getMaxScore() - 1);
                    int score = match.getScore();
                    double specProb = gf.getSpectralProbability(score);
                    assert (specProb > 0) : specIndex + ": " + match.getDeNovoScore() + " " + match.getScore() + " " + specProb;
                    match.setSpecProb(specProb);
                    if (storeScoreDist)
                        match.setScoreDist(gf.getScoreDist());
                }
            }
        }
    }

    public synchronized void generateSpecIndexDBMatchMap() {
        Iterator<Entry<SpecKey, PriorityQueue<DatabaseMatch>>> itr = specKeyDBMatchMap.entrySet().iterator();
        int numPeptidesPerSpec = this.numPeptidesPerSpec;

        while (itr.hasNext()) {
            Entry<SpecKey, PriorityQueue<DatabaseMatch>> entry = itr.next();
            SpecKey specKey = entry.getKey();
            PriorityQueue<DatabaseMatch> matchQueue = entry.getValue();
            if (matchQueue == null || matchQueue.size() == 0)
                continue;
            else {
                Map<String, DatabaseMatch> pepSeqMap = new HashMap<String, DatabaseMatch>();
                for (DatabaseMatch m : matchQueue) {
                    String pepSeq = m.getPepSeq();
//					int index = m.getIndex();
//					char pre = sa.getSequence().getCharAt(index);
//					char post;
//					if(m.isNTermMetCleaved())
//						post = sa.getSequence().getCharAt(index+m.getLength());
//					else
//						post = sa.getSequence().getCharAt(index+m.getLength()-1);
//					String key = pre+pepSeq+post;

                    String key = pepSeq + m.getScore();
                    DatabaseMatch existingMatch = pepSeqMap.get(key);
                    if (existingMatch == null)
                        pepSeqMap.put(key, m);
                    else {
                        for (int index : m.getIndices())
                            existingMatch.addIndex(index);
                    }
                }
                matchQueue = new PriorityQueue<DatabaseMatch>(pepSeqMap.values());
                pepSeqMap = null;
            }


            int specIndex = specKey.getSpecIndex();
            PriorityQueue<DatabaseMatch> existingQueue = specIndexDBMatchMap.get(specIndex);
            if (existingQueue == null) {
                existingQueue = new PriorityQueue<DatabaseMatch>(numPeptidesPerSpec, new DatabaseMatch.SpecProbComparator());
                specIndexDBMatchMap.put(specIndex, existingQueue);
            }

            for (DatabaseMatch match : matchQueue) {
                double curEValue = match.getSpecEValue();
                if (existingQueue.size() < numPeptidesPerSpec || curEValue == existingQueue.peek().getSpecEValue()) {
                    existingQueue.add(match);
                } else {
                    double prevEValue = existingQueue.peek().getSpecEValue();
                    if (curEValue < prevEValue) {
                        while (!existingQueue.isEmpty() && existingQueue.peek().getSpecEValue() == prevEValue)
                            existingQueue.poll();
                        existingQueue.add(match);
                    }
                }
            }
        }
    }

    public synchronized void addResultsToList(List<MSGFPlusMatch> resultList) {
        Iterator<Entry<Integer, PriorityQueue<DatabaseMatch>>> itr = specIndexDBMatchMap.entrySet().iterator();
        while (itr.hasNext()) {
            Entry<Integer, PriorityQueue<DatabaseMatch>> entry = itr.next();
            resultList.add(new MSGFPlusMatch(entry.getKey(), entry.getValue()));
        }
    }

    public void addAdditionalFeatures() {
        Iterator<Entry<Integer, PriorityQueue<DatabaseMatch>>> itr = specIndexDBMatchMap.entrySet().iterator();
        while (itr.hasNext()) {
            Entry<Integer, PriorityQueue<DatabaseMatch>> entry = itr.next();
            int specIndex = entry.getKey();

            PriorityQueue<DatabaseMatch> matchQueue = entry.getValue();
            if (matchQueue == null || matchQueue.size() == 0)
                continue;

            Spectrum spec = specScanner.getSpectraAccessor().getSpectrumBySpecIndex(specIndex);
            for (DatabaseMatch match : matchQueue) {
                NewRankScorer scorer = specScanner.getRankScorer(new SpecKey(specIndex, match.getCharge()));
                if (scorer == null)
                    continue;

                spec.setCharge(match.getCharge());
                PSMFeatureFinder addFeatures = new PSMFeatureFinder(spec, aaSet.getPeptide(match.getPepSeq()), scorer);
                for (Pair<String, String> feature : addFeatures.getAllFeatures())
                    match.addAdditionalFeature(feature.getFirst(), feature.getSecond());
            }
        }
    }

    // for MS-GFDB
    public synchronized void addDBSearchResults(List<MSGFDBResultGenerator.DBMatch> gen, String specFileName, boolean replicateMergedResults) {
        Map<Integer, PriorityQueue<DatabaseMatch>> specIndexDBMatchMap = new HashMap<Integer, PriorityQueue<DatabaseMatch>>();

        Iterator<Entry<SpecKey, PriorityQueue<DatabaseMatch>>> itr = specKeyDBMatchMap.entrySet().iterator();
        while (itr.hasNext()) {
            Entry<SpecKey, PriorityQueue<DatabaseMatch>> entry = itr.next();
            SpecKey specKey = entry.getKey();
            PriorityQueue<DatabaseMatch> matchQueue = entry.getValue();
            if (matchQueue == null || matchQueue.size() == 0)
                continue;

            int specIndex = specKey.getSpecIndex();
            PriorityQueue<DatabaseMatch> existingQueue = specIndexDBMatchMap.get(specIndex);
            if (existingQueue == null) {
                existingQueue = new PriorityQueue<DatabaseMatch>(this.numPeptidesPerSpec, new DatabaseMatch.SpecProbComparator());
                specIndexDBMatchMap.put(specIndex, existingQueue);
            }

            for (DatabaseMatch match : matchQueue) {
                if (existingQueue.size() < this.numPeptidesPerSpec) {
                    existingQueue.add(match);
                } else if (existingQueue.size() >= this.numPeptidesPerSpec) {
                    if (match.getSpecEValue() < existingQueue.peek().getSpecEValue()) {
                        existingQueue.poll();
                        existingQueue.add(match);
                    }
                }
            }
        }

        Iterator<Entry<Integer, PriorityQueue<DatabaseMatch>>> itr2 = specIndexDBMatchMap.entrySet().iterator();
        while (itr2.hasNext()) {
            Entry<Integer, PriorityQueue<DatabaseMatch>> entry = itr2.next();
            int specIndex = entry.getKey();
            PriorityQueue<DatabaseMatch> matchQueue = entry.getValue();
            if (matchQueue == null)
                continue;

            ArrayList<DatabaseMatch> matchList = new ArrayList<DatabaseMatch>(matchQueue);
            if (matchList.size() == 0)
                continue;

            for (int i = matchList.size() - 1; i >= 0; --i) {
                DatabaseMatch match = matchList.get(i);

                if (match.getDeNovoScore() < minDeNovoScore)
                    continue;

                int index = match.getIndex();
                int length = match.getLength();
                int charge = match.getCharge();

                String peptideStr = match.getPepSeq();
                if (peptideStr == null)
                    peptideStr = sa.getSequence().getSubsequence(index + 1, index + length - 1);
                Peptide pep = aaSet.getPeptide(peptideStr);
                String annotationStr = sa.getSequence().getCharAt(index) + "." + pep + "." + sa.getSequence().getCharAt(index + length - 1);
                SimpleDBSearchScorer<NominalMass> scorer = specScanner.getSpecKeyScorerMap().get(new SpecKey(specIndex, charge));
                ArrayList<Integer> specIndexList = specScanner.getSpecKey(specIndex, charge).getSpecIndexList();
                if (specIndexList == null) {
                    specIndexList = new ArrayList<Integer>();
                    specIndexList.add(specIndex);
                }

                float expMass = scorer.getPrecursorPeak().getMass();
//				float theoMass = pep.getParentMass();
                float peptideMass = match.getPeptideMass();
                float pmError = Float.MAX_VALUE;
                float theoMass = peptideMass + (float) Composition.H2O;

                for (int delta = specScanner.getMinIsotopeError(); delta <= specScanner.getMaxIsotopeError(); delta++) {
                    float error = expMass - theoMass - (float) (Composition.ISOTOPE) * delta;
                    if (Math.abs(error) < Math.abs(pmError)) {
                        pmError = error;
                    }
                }
//				if(pmError > )
                if (specScanner.getRightPrecursorMassTolerance().isTolerancePPM())
                    pmError = pmError / theoMass * 1e6f;

                String protein = sa.getAnnotation(index + 1);

                int score = match.getScore();
                double specProb = match.getSpecEValue();
                int numPeptides = sa.getNumDistinctPeptides(peptideStr.length() + 1);
                double pValue = MSGFDBResultGenerator.DBMatch.getPValue(specProb, numPeptides);
                String specProbStr;
                if (specProb < Float.MIN_NORMAL)
                    specProbStr = String.valueOf(specProb);
                else
                    specProbStr = String.valueOf((float) specProb);
                String pValueStr;
                if (specProb < Float.MIN_NORMAL)
                    pValueStr = String.valueOf(pValue);
                else
                    pValueStr = String.valueOf((float) pValue);

                if (!replicateMergedResults) {
                    StringBuffer specIndexStrBuf = new StringBuffer();
                    StringBuffer scanNumStrBuf = new StringBuffer();
                    StringBuffer actMethodStrBuf = new StringBuffer();
                    specIndexStrBuf.append(specIndexList.get(0));
                    actMethodStrBuf.append(scorer.getActivationMethodArr()[0]);
                    scanNumStrBuf.append(scorer.getScanNumArr()[0]);
                    for (int j = 1; j < scorer.getActivationMethodArr().length; j++) {
                        specIndexStrBuf.append("/" + specIndexList.get(j));
                        scanNumStrBuf.append("/" + scorer.getScanNumArr()[j]);
                        actMethodStrBuf.append("/" + scorer.getActivationMethodArr()[j]);
                    }

                    String resultStr =
                            specFileName + "\t"
                                    + specIndexStrBuf.toString() + "\t"
                                    + scanNumStrBuf.toString() + "\t"
                                    + actMethodStrBuf.toString() + "\t"
                                    + scorer.getPrecursorPeak().getMz() + "\t"
                                    + pmError + "\t"
                                    + match.getCharge() + "\t"
                                    + annotationStr + "\t"
                                    + protein + "\t"
                                    + match.getDeNovoScore() + "\t"
                                    + score + "\t"
                                    + specProbStr + "\t"
                                    + pValueStr;
                    MSGFDBResultGenerator.DBMatch dbMatch = new MSGFDBResultGenerator.DBMatch(specProb, numPeptides, resultStr, match.getScoreDist());
                    gen.add(dbMatch);
                } else {
                    for (int j = 0; j < scorer.getActivationMethodArr().length; j++) {
                        String resultStr =
                                specFileName + "\t"
                                        + specIndexList.get(j) + "\t"
                                        + scorer.getScanNumArr()[j] + "\t"
                                        + scorer.getActivationMethodArr()[j] + "\t"
                                        + scorer.getPrecursorPeak().getMz() + "\t"
                                        + pmError + "\t"
                                        + match.getCharge() + "\t"
                                        + annotationStr + "\t"
                                        + protein + "\t"
                                        + match.getDeNovoScore() + "\t"
                                        + score + "\t"
                                        + specProbStr + "\t"
                                        + pValueStr;
                        MSGFDBResultGenerator.DBMatch dbMatch = new MSGFDBResultGenerator.DBMatch(specProb, numPeptides, resultStr, match.getScoreDist());
                        gen.add(dbMatch);
                    }
                }
            }
        }
    }

    public static void setAminoAcidProbabilities(String databaseFileName, AminoAcidSet aaSet) {
        BufferedLineReader in = null;
        try {
            in = new BufferedLineReader(databaseFileName);
        } catch (IOException e) {
            e.printStackTrace();
        }

        long[] aaCount = new long[128];
        String s;
        while ((s = in.readLine()) != null) {
            if (s.startsWith(">"))    // annotation
                continue;
            for (int i = 0; i < s.length(); i++) {
                char residue = s.charAt(i);
                //if(aaSet.getAminoAcid(residue) != null)
                if (Character.isLetter(residue))
                    aaCount[residue]++;
            }
        }
        long totalAACount = 0;
        for (AminoAcid aa : aaSet.getAAList(Location.Anywhere))
            if (!aa.isModified())
                totalAACount += aaCount[aa.getResidue()];

        boolean success = true;
        for (AminoAcid aa : aaSet.getAllAminoAcidArr()) {
            long count = aaCount[aa.getUnmodResidue()];
            if (count == 0 && AminoAcid.isStdAminoAcid(aa.getUnmodResidue())) {
                success = false;
                break;
            }
            aa.setProbability(count / (float) totalAACount);
        }
        for (int i = 0; i < 128; i++) {
            if (!aaSet.contains((char) i) && aaCount[i] > 0) {
                System.out.println("Warning: Sequence database contains " +
                        aaCount[i] + " counts of letter '" + (char) i +
                        "', which does not correspond to an amino acid.");
            }
        }

        if (!success) {
            System.out.println("Warning: database does not contain all standard amino acids. " +
                    "Probability 0.05 will be used for all amino acids.");
            for (AminoAcid aa : aaSet.getAllAminoAcidArr())
                aa.setProbability(0.05f);
        }
    }
}
