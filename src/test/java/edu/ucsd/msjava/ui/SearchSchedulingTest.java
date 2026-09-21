package edu.ucsd.msjava.ui;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.File;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static org.junit.Assert.*;

public class SearchSchedulingTest {
    @Rule
    public TemporaryFolder temporaryFolder = new TemporaryFolder();

    // Run a real search in a separate JVM so the memory-budget calculation sees the requested heap limit.
    private String search(String heap, int tasks) throws Exception {
        File spectra = temporaryFolder.newFile("spectra.mgf");
        try (PrintWriter writer = new PrintWriter(spectra, "UTF-8")) {
            for (int i = 1; i <= 500; i++) {
                writer.println("BEGIN IONS\nTITLE=spectrum" + i + "\nPEPMASS=500\nCHARGE=2+\nSCANS=" + i);
                for (int mass = 100; mass <= 400; mass += 20)
                    writer.println(mass + " 50");
                writer.println("END IONS");
            }
        }
        File database = temporaryFolder.newFile("database.fasta");
        try (PrintWriter writer = new PrintWriter(database, "UTF-8")) {
            writer.println(">protein\nAAAAKR");
        }
        File result = new File(temporaryFolder.getRoot(), "result.mzid");
        File logFile = temporaryFolder.newFile("search.log");
        List<String> command = Arrays.asList(
                new File(System.getProperty("java.home"), "bin/java").getPath(), "-Xmx" + heap,
                "-cp", System.getProperty("java.class.path"), MSGFPlus.class.getName(),
                "-s", spectra.getPath(), "-d", database.getPath(), "-o", result.getPath(),
                "-thread", "1", "-tasks", Integer.toString(tasks));
        Process process = new ProcessBuilder(command).redirectErrorStream(true).redirectOutput(logFile).start();
        try {
            assertTrue("search did not finish", process.waitFor(60, TimeUnit.SECONDS));
            String log = new String(Files.readAllBytes(logFile.toPath()), StandardCharsets.UTF_8);
            assertEquals(log, 0, process.exitValue());
            assertTrue(log, log.contains("MS-GF+ complete"));
            assertTrue(log, result.isFile());
            return log;
        } finally {
            process.destroyForcibly();
        }
    }

    @Test
    public void smallHeapDoesNotForceOneSpectrumPerBatch() throws Exception {
        String log = search("256m", 1);
        Matcher batches = Pattern.compile("(?m)^Spectrum [0-9]+-[0-9]+ ").matcher(log);
        assertTrue(log, batches.find());
        assertFalse(log, batches.find());
        assertTrue(log, log.contains("Spectrum 0-499 (total: 500)"));
    }

    @Test
    public void explicitTasksLimitBatchSizeWithoutMultiplyingJobs() throws Exception {
        String log = search("512m", 64);
        assertTrue(log, log.contains("Spectrum 0-6 (total: 500)"));
        assertTrue(log, log.contains("Spectrum 497-499 (total: 500)"));
        Matcher progress = Pattern.compile("(?m)^(Preprocessing|Search|E-value) progress: ([0-9]+) / ([0-9]+) tasks,").matcher(log);
        int reports = 0;
        while (progress.find()) {
            // In-flight executor counts are approximate; inspect completed phase reports.
            if (!progress.group(2).equals(progress.group(3)))
                continue;
            // The database has only eight suffixes, so even the search must submit fewer than twelve jobs.
            int limit = progress.group(1).equals("Search") ? 8 : 3;
            assertTrue(progress.group(), Integer.parseInt(progress.group(3)) <= limit);
            reports++;
        }
        assertTrue("no phase progress was reported", reports > 0);
    }
}
