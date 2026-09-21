package edu.ucsd.msjava.ui;

import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msutil.SpecKey;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;

/**
 * The spectra searched together are limited by the memory of their scorers, which grows with the peptide mass:
 * about 13 KB for a key at 1,000 Da and 73 KB at 6,000 Da.
 */
public class SpectrumBatchTest {

    private static final Tolerance TOLERANCE = new Tolerance(20, true);
    private static final int MIN_ISOTOPE_ERROR = -1;

    // numSpectra spectra of numCharges charge states each, all of about peptideMass Da
    private static List<SpecKey> getSpecKeyList(int numSpectra, int numCharges, float peptideMass) {
        List<SpecKey> specKeyList = new ArrayList<SpecKey>();
        for (int specIndex = 1; specIndex <= numSpectra; specIndex++) {
            for (int charge = 2; charge < 2 + numCharges; charge++)
                specKeyList.add(new SpecKey(specIndex, charge, (peptideMass + 18.0106f) / charge + 1.00727f));
        }
        return specKeyList;
    }

    @Test
    public void smallHeapsStillHoldManySpectra() {
        List<SpecKey> specKeyList = getSpecKeyList(500, 1, 6000);
        long scorerMemory = MSGFPlus.getScorerMemory(256L << 20, 100000);
        assertEquals(500, MSGFPlus.getBatchEndIndex(specKeyList, 0, 500, scorerMemory, TOLERANCE, MIN_ISOTOPE_ERROR));
    }

    @Test
    public void reserveRemainsCappedOnLargeHeaps() {
        assertEquals((4L << 30) - (1L << 27) - 50000, MSGFPlus.getScorerMemory(8L << 30, 100000));
        assertEquals(0, MSGFPlus.getScorerMemory(256L << 20, 256L << 20));
    }

    @Test
    public void allSpectraWhenTheyFit() {
        List<SpecKey> specKeyList = getSpecKeyList(10, 1, 1000);
        assertEquals(10, MSGFPlus.getBatchEndIndex(specKeyList, 0, 10, 1L << 30, TOLERANCE, MIN_ISOTOPE_ERROR));
        assertEquals(10, MSGFPlus.getBatchEndIndex(specKeyList, 4, 10, 1L << 30, TOLERANCE, MIN_ISOTOPE_ERROR));
    }

    @Test
    public void asManySpectraAsFitInTheMemory() {
        List<SpecKey> specKeyList = getSpecKeyList(10, 1, 1000);
        assertEquals(7, MSGFPlus.getBatchEndIndex(specKeyList, 0, 10, 100000, TOLERANCE, MIN_ISOTOPE_ERROR));
        assertEquals(10, MSGFPlus.getBatchEndIndex(specKeyList, 7, 10, 100000, TOLERANCE, MIN_ISOTOPE_ERROR));
    }

    @Test
    public void fewerSpectraWhenTheyAreHeavier() {
        List<SpecKey> specKeyList = getSpecKeyList(10, 1, 6000);
        assertEquals(1, MSGFPlus.getBatchEndIndex(specKeyList, 0, 10, 100000, TOLERANCE, MIN_ISOTOPE_ERROR));
        assertEquals(4, MSGFPlus.getBatchEndIndex(specKeyList, 0, 10, 300000, TOLERANCE, MIN_ISOTOPE_ERROR));
    }

    @Test
    public void atLeastOneSpectrum() {
        List<SpecKey> specKeyList = getSpecKeyList(10, 1, 6000);
        assertEquals(1, MSGFPlus.getBatchEndIndex(specKeyList, 0, 10, 0, TOLERANCE, MIN_ISOTOPE_ERROR));
        assertEquals(1, MSGFPlus.getBatchEndIndex(specKeyList, 0, 10, -1L << 30, TOLERANCE, MIN_ISOTOPE_ERROR));
        assertEquals(6, MSGFPlus.getBatchEndIndex(specKeyList, 5, 10, 0, TOLERANCE, MIN_ISOTOPE_ERROR));
    }

    @Test
    public void noMoreThanTheGivenNumberOfKeys() {
        List<SpecKey> specKeyList = getSpecKeyList(10, 1, 1000);
        assertEquals(3, MSGFPlus.getBatchEndIndex(specKeyList, 0, 3, 1L << 30, TOLERANCE, MIN_ISOTOPE_ERROR));
        assertEquals(6, MSGFPlus.getBatchEndIndex(specKeyList, 3, 3, 1L << 30, TOLERANCE, MIN_ISOTOPE_ERROR));
    }

    @Test
    public void chargeStatesOfASpectrumStayTogether() {
        // keys 0-2 belong to the first spectrum, 3-5 to the second
        List<SpecKey> specKeyList = getSpecKeyList(4, 3, 1000);
        assertEquals(3, MSGFPlus.getBatchEndIndex(specKeyList, 0, 1, 1L << 30, TOLERANCE, MIN_ISOTOPE_ERROR));
        assertEquals(6, MSGFPlus.getBatchEndIndex(specKeyList, 0, 4, 1L << 30, TOLERANCE, MIN_ISOTOPE_ERROR));
        assertEquals(3, MSGFPlus.getBatchEndIndex(specKeyList, 0, 12, 0, TOLERANCE, MIN_ISOTOPE_ERROR));
    }
}
