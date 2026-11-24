package pgl.app.popdep;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import pgl.AppAbstract;
import pgl.PGLAPPEntrance;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Population depth profiling with one sequential BAM scan per taxon across all chromosomes.
 *
 * <p>Unlike {@link PopDepFull}, which re-reads each taxon's BAM once per chromosome,
 * PopDepCrossChr runs a single {@code samtools depth} call per taxon over the full
 * alignment file and accumulates mean/SD statistics for every chromosome in one pass.
 * On HDD this reduces total disk reads by roughly the number of chromosomes (e.g. ~43×
 * for a typical chicken genome).</p>
 *
 * <p>Output format per chromosome matches PopDep / PopDepFull:
 * {@code Position\tDepth_Mean\tDepth_SD}, written as {@code {outdir}/{chr}.popdep.txt.gz}.</p>
 *
 * <p>Chromosome length file format ({@code -b}): {@code Chr\\tLength\\tnTaxa}. The third column
 * is the per-chromosome taxon count used as the denominator for mean/SD (matching PopDepFull's
 * tb.A / tb.B / tb.D / tb.ALL row counts). If omitted, defaults to {@code taxa.length} for all
 * chromosomes.</p>
 *
 * <p>Checkpoint / resume ({@code -k}): when a checkpoint directory is given, taxa are processed in
 * batches; after each batch the shared {@code sum}/{@code sumSq} accumulators and the list of
 * completed taxa are written atomically to {@code checkpoint.bin}. On restart with the same
 * {@code -k} directory, the accumulators are reloaded and already-completed taxa are skipped, so a
 * crash only loses the in-progress batch. The raw {@code samtools depth} stream is never stored.</p>
 */
public class PopDepCrossChr extends AppAbstract {
    String taxaBamFileS = null;
    String chrLengthFileS = null;
    String samPath = null;
    /** Minimum samtools mapping quality (MAPQ); passed to depth as {@code -q}. */
    int minMapq = 0;
    int threadNum = 16;
    String outDirS = null;
    String[] taxa = null;
    HashMap<String, String[]> taxaBamPathsMap = new HashMap<>();
    /** Per-taxon genome-wide mean depth (column 2 of taxaBamMap), used to normalize to relative depth. */
    HashMap<String, Double> taxaCoverageMap = new HashMap<>();
    /** Coverage aligned to the final {@link #taxa} order (after orderTaxaForIO). */
    double[] taxaCoverage = null;

    /**
     * Taxa submission order: {@code interleave} (round-robin by BAM directory, default) spreads
     * concurrent samtools across physical disks; {@code shuffle} randomizes (fixed seed);
     * {@code sorted} keeps lexical order.
     */
    String taxaOrder = "interleave";

    /** Checkpoint directory; null disables checkpointing (single continuous pass). */
    String checkpointDir = null;
    /** Number of taxa completed between checkpoints (batch size when checkpointing). */
    int checkpointIntervalTaxa = 500;

    /** Lock stripes per chromosome to reduce contention when many taxa run in parallel. */
    private static final int LOCK_SEGMENTS = 256;
    /** Reader buffer for samtools stdout; large buffer reduces pipe read syscalls. */
    private static final int READ_BUFFER_BYTES = 1 << 20;
    /** Checkpoint file name inside the checkpoint directory. */
    private static final String CHECKPOINT_NAME = "checkpoint.bin";
    /** Bulk I/O buffer for reading/writing accumulator arrays (must be a multiple of 8). */
    private static final int CKPT_BUFFER_BYTES = 8 << 20;
    private static final int CKPT_MAGIC = 0x50445043; // "PDPC"
    private static final int CKPT_VERSION = 2; // v2 adds relSum/relSumSq arrays

    short[] chromosomes = null;
    int[] chrLengths = null;
    /** Per-chromosome taxon count for mean/SD denominator (from length file column 3). */
    int[] nTaxaPerChr = null;
    HashMap<String, Integer> refToChrIndex = new HashMap<>();
    Object[][] segmentLocks = null;

    /** Shared accumulators, allocated once before the worker pool starts. */
    double[][] sum = null;
    double[][] sumSq = null;
    /** Relative-depth accumulators (depth normalized by per-taxon coverage). */
    double[][] relSum = null;
    double[][] relSumSq = null;

    /** Chunk size for progress messages while writing output. */
    int progressWindowSize = 500_000;

    public PopDepCrossChr(String[] args) {
        this.creatAppOptions();
        this.retrieveAppParameters(args);
        this.profileDepth();
    }

    public void profileDepth() {
        int nChr = chromosomes.length;
        sum = new double[nChr][];
        sumSq = new double[nChr][];
        relSum = new double[nChr][];
        relSumSq = new double[nChr][];
        segmentLocks = new Object[nChr][LOCK_SEGMENTS];
        for (int i = 0; i < nChr; i++) {
            sum[i] = new double[chrLengths[i]];
            sumSq[i] = new double[chrLengths[i]];
            relSum[i] = new double[chrLengths[i]];
            relSumSq[i] = new double[chrLengths[i]];
            for (int s = 0; s < LOCK_SEGMENTS; s++) {
                segmentLocks[i][s] = new Object();
            }
        }

        String[] commands = buildDepthCommands();
        System.out.println("Example depth command: " + commands[0]);

        LinkedHashSet<String> doneSet = new LinkedHashSet<>();
        if (checkpointDir != null) {
            File ckptFile = new File(checkpointDir, CHECKPOINT_NAME);
            if (ckptFile.exists()) {
                try {
                    loadCheckpoint(ckptFile, doneSet);
                    System.out.println("Resumed from checkpoint " + ckptFile.getAbsolutePath()
                            + ": " + doneSet.size() + " taxa already accumulated.");
                }
                catch (Exception e) {
                    throw new RuntimeException("Failed to load checkpoint " + ckptFile.getAbsolutePath()
                            + " (delete it to restart from scratch): " + e.getMessage(), e);
                }
            }
            else {
                System.out.println("No existing checkpoint; starting fresh. Checkpoint dir: " + checkpointDir);
            }
        }

        try {
            if (checkpointDir == null) {
                runContinuous(commands);
            }
            else {
                runBatchedWithCheckpoint(commands, doneSet);
            }
            writeAllChromosomeOutputs(sum, sumSq);
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("PopDepCrossChr is finished.");
    }

    /** Single continuous pass over all taxa (no checkpointing). */
    private void runContinuous(String[] commands) throws Exception {
        int total = taxa.length;
        AtomicInteger finished = new AtomicInteger(0);
        AtomicInteger errors = new AtomicInteger(0);
        List<Future<?>> futures = new ArrayList<>();
        ExecutorService pool = Executors.newFixedThreadPool(threadNum);
        for (int j = 0; j < total; j++) {
            final int taxonIndex = j;
            futures.add(pool.submit(() -> {
                try {
                    processTaxon(commands[taxonIndex], taxonIndex, finished, total);
                }
                catch (RuntimeException ex) {
                    errors.incrementAndGet();
                    throw ex;
                }
            }));
        }
        pool.shutdown();
        try {
            for (Future<?> future : futures) {
                future.get();
            }
            if (!pool.awaitTermination(1, TimeUnit.HOURS)) {
                pool.shutdownNow();
            }
        }
        finally {
            pool.shutdownNow();
        }
        if (errors.get() > 0) {
            throw new RuntimeException("PopDepCrossChr failed with " + errors.get() + " taxon error(s)");
        }
    }

    /**
     * Batched pass with a checkpoint after every {@link #checkpointIntervalTaxa} completed taxa.
     * The batch barrier guarantees a consistent {@code sum}/{@code sumSq} snapshot (no taxon is
     * mid-application) before each checkpoint is written.
     */
    private void runBatchedWithCheckpoint(String[] commands, LinkedHashSet<String> doneSet) throws Exception {
        int total = taxa.length;
        List<Integer> todo = new ArrayList<>();
        for (int j = 0; j < total; j++) {
            if (!doneSet.contains(taxa[j])) {
                todo.add(j);
            }
        }
        AtomicInteger finished = new AtomicInteger(doneSet.size());
        System.out.println("Checkpointing every " + checkpointIntervalTaxa + " taxa; "
                + todo.size() + " to process, " + doneSet.size() + " already done.");

        ExecutorService pool = Executors.newFixedThreadPool(threadNum);
        try {
            int start = 0;
            while (start < todo.size()) {
                int end = Math.min(start + checkpointIntervalTaxa, todo.size());
                List<Future<?>> futures = new ArrayList<>();
                AtomicInteger errors = new AtomicInteger(0);
                for (int bi = start; bi < end; bi++) {
                    final int taxonIndex = todo.get(bi);
                    futures.add(pool.submit(() -> {
                        try {
                            processTaxon(commands[taxonIndex], taxonIndex, finished, total);
                        }
                        catch (RuntimeException ex) {
                            errors.incrementAndGet();
                            throw ex;
                        }
                    }));
                }
                for (Future<?> future : futures) {
                    future.get();
                }
                if (errors.get() > 0) {
                    throw new RuntimeException("Batch failed; last checkpoint preserved for resume.");
                }
                for (int bi = start; bi < end; bi++) {
                    doneSet.add(taxa[todo.get(bi)]);
                }
                long t0 = System.currentTimeMillis();
                writeCheckpoint(doneSet);
                System.out.println("Checkpoint written after " + doneSet.size() + " / " + total
                        + " taxa (" + (System.currentTimeMillis() - t0) + " ms).");
                start = end;
            }
        }
        finally {
            pool.shutdownNow();
        }
    }

    private void processTaxon(String command, int taxonIndex, AtomicInteger finished, int total) {
        try {
            long linesRead = streamDepthAndAccumulate(command, taxaCoverage[taxonIndex]);
            int done = finished.incrementAndGet();
            System.out.println("Finished taxa " + done + " / " + total + " ("
                    + taxa[taxonIndex] + ", depth lines=" + linesRead + ")");
        }
        catch (Exception e) {
            throw new RuntimeException("Depth scan failed for taxon " + taxa[taxonIndex], e);
        }
    }

    /**
     * Atomically persist the current accumulators plus the set of completed taxa. The file is
     * written to a temporary path and renamed, so an interrupted write never corrupts the
     * resumable checkpoint. Arrays are written via a bulk {@link ByteBuffer} (not per-double) to
     * stay disk-bandwidth bound rather than CPU bound.
     */
    private void writeCheckpoint(LinkedHashSet<String> doneSet) throws Exception {
        File dir = new File(checkpointDir);
        if (!dir.exists() && !dir.mkdirs()) {
            throw new IOException("Failed to create checkpoint directory: " + checkpointDir);
        }
        Path finalPath = new File(dir, CHECKPOINT_NAME).toPath();
        Path tmpPath = new File(dir, CHECKPOINT_NAME + ".tmp").toPath();

        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        try (DataOutputStream dos = new DataOutputStream(baos)) {
            dos.writeInt(CKPT_MAGIC);
            dos.writeInt(CKPT_VERSION);
            dos.writeInt(chromosomes.length);
            for (int c = 0; c < chromosomes.length; c++) {
                dos.writeShort(chromosomes[c]);
                dos.writeInt(chrLengths[c]);
                dos.writeInt(nTaxaPerChr[c]);
            }
            dos.writeInt(doneSet.size());
            for (String s : doneSet) {
                dos.writeUTF(s);
            }
        }
        byte[] header = baos.toByteArray();

        try (FileChannel ch = FileChannel.open(tmpPath, StandardOpenOption.CREATE,
                StandardOpenOption.WRITE, StandardOpenOption.TRUNCATE_EXISTING)) {
            ByteBuffer hb = ByteBuffer.allocate(4 + header.length);
            hb.putInt(header.length);
            hb.put(header);
            hb.flip();
            while (hb.hasRemaining()) {
                ch.write(hb);
            }
            ByteBuffer buf = ByteBuffer.allocateDirect(CKPT_BUFFER_BYTES);
            for (int c = 0; c < chromosomes.length; c++) {
                writeDoubleArray(ch, buf, sum[c]);
            }
            for (int c = 0; c < chromosomes.length; c++) {
                writeDoubleArray(ch, buf, sumSq[c]);
            }
            for (int c = 0; c < chromosomes.length; c++) {
                writeDoubleArray(ch, buf, relSum[c]);
            }
            for (int c = 0; c < chromosomes.length; c++) {
                writeDoubleArray(ch, buf, relSumSq[c]);
            }
            ch.force(true);
        }
        Files.move(tmpPath, finalPath, StandardCopyOption.ATOMIC_MOVE, StandardCopyOption.REPLACE_EXISTING);
    }

    /** Load accumulators and completed-taxa set; validates the genome signature before trusting it. */
    private void loadCheckpoint(File ckptFile, LinkedHashSet<String> doneSetOut) throws Exception {
        try (FileChannel ch = FileChannel.open(ckptFile.toPath(), StandardOpenOption.READ)) {
            ByteBuffer lenBuf = ByteBuffer.allocate(4);
            readFully(ch, lenBuf);
            lenBuf.flip();
            int headerLen = lenBuf.getInt();
            if (headerLen <= 0 || headerLen > (64 << 20)) {
                throw new IOException("Invalid checkpoint header length: " + headerLen);
            }
            ByteBuffer hb = ByteBuffer.allocate(headerLen);
            readFully(ch, hb);
            hb.flip();
            byte[] headerBytes = new byte[headerLen];
            hb.get(headerBytes);

            try (DataInputStream dis = new DataInputStream(new ByteArrayInputStream(headerBytes))) {
                if (dis.readInt() != CKPT_MAGIC) {
                    throw new IOException("Not a PopDepCrossChr checkpoint (bad magic).");
                }
                if (dis.readInt() != CKPT_VERSION) {
                    throw new IOException("Unsupported checkpoint version.");
                }
                int nChr = dis.readInt();
                if (nChr != chromosomes.length) {
                    throw new IOException("Checkpoint chromosome count " + nChr
                            + " != current " + chromosomes.length);
                }
                for (int c = 0; c < nChr; c++) {
                    short chrom = dis.readShort();
                    int len = dis.readInt();
                    int nt = dis.readInt();
                    if (chrom != chromosomes[c] || len != chrLengths[c] || nt != nTaxaPerChr[c]) {
                        throw new IOException("Checkpoint chromosome metadata mismatch at index " + c
                                + " (chr/len/nTaxa). Length file or taxa set changed?");
                    }
                }
                int doneCount = dis.readInt();
                for (int k = 0; k < doneCount; k++) {
                    doneSetOut.add(dis.readUTF());
                }
            }

            ByteBuffer buf = ByteBuffer.allocateDirect(CKPT_BUFFER_BYTES);
            for (int c = 0; c < chromosomes.length; c++) {
                readDoubleArray(ch, buf, sum[c]);
            }
            for (int c = 0; c < chromosomes.length; c++) {
                readDoubleArray(ch, buf, sumSq[c]);
            }
            for (int c = 0; c < chromosomes.length; c++) {
                readDoubleArray(ch, buf, relSum[c]);
            }
            for (int c = 0; c < chromosomes.length; c++) {
                readDoubleArray(ch, buf, relSumSq[c]);
            }
        }
    }

    private static void writeDoubleArray(FileChannel ch, ByteBuffer buf, double[] arr) throws IOException {
        int chunk = buf.capacity() / 8;
        int off = 0;
        while (off < arr.length) {
            int m = Math.min(chunk, arr.length - off);
            buf.clear();
            DoubleBuffer db = buf.asDoubleBuffer();
            db.put(arr, off, m);
            buf.position(0).limit(m * 8);
            while (buf.hasRemaining()) {
                ch.write(buf);
            }
            off += m;
        }
    }

    private static void readDoubleArray(FileChannel ch, ByteBuffer buf, double[] arr) throws IOException {
        int chunk = buf.capacity() / 8;
        int off = 0;
        while (off < arr.length) {
            int m = Math.min(chunk, arr.length - off);
            buf.clear();
            buf.limit(m * 8);
            while (buf.hasRemaining()) {
                if (ch.read(buf) < 0) {
                    throw new EOFException("Checkpoint truncated while reading accumulator arrays.");
                }
            }
            buf.flip();
            DoubleBuffer db = buf.asDoubleBuffer();
            db.get(arr, off, m);
            off += m;
        }
    }

    private static void readFully(FileChannel ch, ByteBuffer b) throws IOException {
        while (b.hasRemaining()) {
            if (ch.read(b) < 0) {
                throw new EOFException("Checkpoint truncated.");
            }
        }
    }

    /** Map position to lock stripe; use long arithmetic to avoid int overflow on large chromosomes. */
    private static int segmentIndex(int pos, int chrLength) {
        return (int) Math.min(LOCK_SEGMENTS - 1, (long) (pos - 1) * LOCK_SEGMENTS / chrLength);
    }

    /** Write one {@code {chr}.popdep.txt.gz} per chromosome in parallel (each file is independent). */
    private void writeAllChromosomeOutputs(double[][] sum, double[][] sumSq) throws Exception {
        File outDir = new File(outDirS);
        if (!outDir.exists() && !outDir.mkdirs()) {
            throw new RuntimeException("Failed to create output directory: " + outDirS);
        }

        int nChr = chromosomes.length;
        int writeThreads = Math.min(threadNum, nChr);
        System.out.println("Writing " + nChr + " chromosome output file(s) with " + writeThreads + " thread(s).");

        AtomicInteger errors = new AtomicInteger(0);
        List<Future<?>> futures = new ArrayList<>();
        ExecutorService pool = Executors.newFixedThreadPool(writeThreads);
        try {
            for (int c = 0; c < nChr; c++) {
                final int chrIndex = c;
                futures.add(pool.submit(() -> {
                    try {
                        writeChromosomeOutput(chrIndex, sum, sumSq, outDir);
                    }
                    catch (RuntimeException ex) {
                        errors.incrementAndGet();
                        throw ex;
                    }
                    catch (Exception ex) {
                        errors.incrementAndGet();
                        throw new RuntimeException(ex);
                    }
                }));
            }
            for (Future<?> future : futures) {
                future.get();
            }
            pool.shutdown();
            if (!pool.awaitTermination(1, TimeUnit.HOURS)) {
                pool.shutdownNow();
            }
        }
        finally {
            pool.shutdownNow();
        }
        if (errors.get() > 0) {
            throw new RuntimeException("PopDepCrossChr failed writing " + errors.get() + " chromosome output(s)");
        }
    }

    private void writeChromosomeOutput(int c, double[][] sum, double[][] sumSq, File outDir) throws Exception {
        short chr = chromosomes[c];
        int len = chrLengths[c];
        int nTaxa = nTaxaPerChr[c];
        String outFile = new File(outDir, chr + ".popdep.txt.gz").getPath();
        StringBuilder sb = new StringBuilder();
        BufferedWriter bw = IOUtils.getTextGzipWriter(outFile);
        bw.write("Position\tDepth_Mean\tDepth_SD\tRelativeDepth_Mean\tRelativeDepth_SD");
        bw.newLine();

        int[][] progressWindows = PArrayUtils.getSubsetsIndicesBySubsetSize(len, progressWindowSize);
        for (int i = 0; i < progressWindows.length; i++) {
            int start = progressWindows[i][0];
            int end = progressWindows[i][1];
            for (int pos = start; pos < end; pos++) {
                sb.setLength(0);
                double mean = sum[c][pos] / nTaxa;
                double var = (sumSq[c][pos] - sum[c][pos] * sum[c][pos] / nTaxa) / (nTaxa - 1);
                double sd = var > 0 ? Math.sqrt(var) : 0;
                double relMean = relSum[c][pos] / nTaxa;
                double relVar = (relSumSq[c][pos] - relSum[c][pos] * relSum[c][pos] / nTaxa) / (nTaxa - 1);
                double relSd = relVar > 0 ? Math.sqrt(relVar) : 0;
                sb.append(pos + 1).append("\t").append((float) mean).append("\t").append((float) sd)
                        .append("\t").append((float) relMean).append("\t").append((float) relSd);
                bw.write(sb.toString());
                bw.newLine();
            }
            System.out.println("Wrote positions 1-" + end + " for chromosome " + chr
                    + " (nTaxa=" + nTaxa + ") -> " + outFile);
        }
        bw.flush();
        bw.close();
    }

    /**
     * Run one samtools depth process and accumulate directly into {@link #sum}/{@link #sumSq}.
     *
     * <p>Hot path is hand-written: each line is parsed with a single character scan (no
     * {@code fastSplit} List allocation, no {@code Integer.parseInt}), and the reference name is
     * cached across consecutive lines (depth output keeps one chromosome contiguous), so the
     * {@code refToChrIndex} lookup happens only when the chromosome changes. This removes the
     * per-stream single-thread parsing bottleneck that previously starved samtools (pipe_write).</p>
     *
     * @param coverage the taxon's genome-wide mean depth, used to normalize to relative depth
     * @return number of depth lines read
     */
    private long streamDepthAndAccumulate(String command, double coverage) throws Exception {
        List<String> cmd = parseCommand(command);
        ProcessBuilder pb = new ProcessBuilder(cmd);
        pb.redirectErrorStream(false);
        Process p = pb.start();
        try {
            p.getOutputStream().close();
        } catch (Exception ignored) {}

        Thread errDrainer = new Thread(() -> {
            try (BufferedReader ebr = new BufferedReader(new InputStreamReader(p.getErrorStream()))) {
                String errLine;
                while ((errLine = ebr.readLine()) != null) {
                    System.err.println("[samtools std err]: " + errLine);
                }
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        }, "samtools-stderr-drainer");
        errDrainer.setDaemon(true);
        errDrainer.start();

        long linesRead = 0;
        // Cache of the most recent reference name region -> chromosome index.
        String cachedRef = null;
        int cachedChrIdx = -1;
        int cachedLen = 0;

        try (BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()), READ_BUFFER_BYTES)) {
            String line;
            while ((line = br.readLine()) != null) {
                linesRead++;
                int n = line.length();
                if (n == 0) {
                    continue;
                }

                // Field 1: reference name region [0, refEnd)
                int i = 0;
                while (i < n && line.charAt(i) != '\t') {
                    i++;
                }
                int refEnd = i;
                if (i >= n) {
                    continue;
                }

                // Resolve chromosome index, reusing cache while the reference name is unchanged.
                int chrIdx;
                int len;
                if (regionEquals(line, 0, refEnd, cachedRef)) {
                    chrIdx = cachedChrIdx;
                    len = cachedLen;
                }
                else {
                    String ref = line.substring(0, refEnd);
                    Integer idxObj = refToChrIndex.get(ref);
                    cachedRef = ref;
                    cachedChrIdx = (idxObj == null) ? -1 : idxObj;
                    cachedLen = (idxObj == null) ? 0 : chrLengths[cachedChrIdx];
                    chrIdx = cachedChrIdx;
                    len = cachedLen;
                }

                // Field 2: position (1-based)
                i++; // skip tab after reference name
                int pos = 0;
                while (i < n) {
                    char c = line.charAt(i);
                    if (c == '\t') {
                        break;
                    }
                    pos = pos * 10 + (c - '0');
                    i++;
                }

                // Fields 3..: per-BAM depths, summed
                int depth = 0;
                while (i < n) {
                    i++; // skip tab
                    int v = 0;
                    while (i < n) {
                        char c = line.charAt(i);
                        if (c == '\t') {
                            break;
                        }
                        v = v * 10 + (c - '0');
                        i++;
                    }
                    depth += v;
                }

                if (chrIdx < 0 || pos < 1 || pos > len) {
                    continue;
                }
                int seg = segmentIndex(pos, len);
                double rv = depth / coverage;
                synchronized (segmentLocks[chrIdx][seg]) {
                    int idx = pos - 1;
                    sum[chrIdx][idx] += depth;
                    sumSq[chrIdx][idx] += (double) depth * depth;
                    relSum[chrIdx][idx] += rv;
                    relSumSq[chrIdx][idx] += rv * rv;
                }
            }
        }

        errDrainer.join();
        int exit = p.waitFor();
        if (exit != 0) {
            throw new RuntimeException("Command failed (exit=" + exit + "): " + command);
        }
        return linesRead;
    }

    /** True if {@code s} equals the substring {@code line[start, end)} without allocating. */
    private static boolean regionEquals(String line, int start, int end, String s) {
        if (s == null) {
            return false;
        }
        int len = end - start;
        if (s.length() != len) {
            return false;
        }
        for (int k = 0; k < len; k++) {
            if (line.charAt(start + k) != s.charAt(k)) {
                return false;
            }
        }
        return true;
    }

    private List<String> parseCommand(String command) {
        List<String> parts = new ArrayList<>();
        StringBuilder cur = new StringBuilder();
        boolean inQuote = false;
        for (int i = 0; i < command.length(); i++) {
            char c = command.charAt(i);
            if (c == '"') {
                inQuote = !inQuote;
            } else if (Character.isWhitespace(c) && !inQuote) {
                if (cur.length() > 0) {
                    parts.add(cur.toString());
                    cur.setLength(0);
                }
            } else {
                cur.append(c);
            }
        }
        if (cur.length() > 0) {
            parts.add(cur.toString());
        }
        return parts;
    }

    private String[] buildDepthCommands() {
        String[] commands = new String[taxa.length];
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < commands.length; i++) {
            sb.setLength(0);
            sb.append(this.samPath);
            sb.append(" depth -Q 20 -q ").append(this.minMapq);
            String[] paths = this.taxaBamPathsMap.get(taxa[i]);
            for (String path : paths) {
                sb.append(' ').append(path);
            }
            commands[i] = sb.toString();
        }
        return commands;
    }

    private void loadChromosomeLengths() throws Exception {
        BufferedReader br = IOUtils.getTextReader(chrLengthFileS);
        List<Short> chrList = new ArrayList<>();
        List<Integer> lenList = new ArrayList<>();
        List<Integer> nTaxaList = new ArrayList<>();
        boolean hasExplicitNTaxa = false;
        String temp;
        while ((temp = br.readLine()) != null) {
            temp = temp.trim();
            if (temp.isEmpty() || temp.startsWith("#")) {
                continue;
            }
            List<String> l = PStringUtils.fastSplit(temp);
            if (l.size() < 2) {
                continue;
            }
            String chrToken = l.get(0);
            String lenToken = l.get(1);
            if (!isInteger(lenToken)) {
                continue;
            }
            chrList.add(Short.parseShort(chrToken));
            lenList.add(Integer.parseInt(lenToken));
            if (l.size() >= 3 && isInteger(l.get(2))) {
                nTaxaList.add(Integer.parseInt(l.get(2)));
                hasExplicitNTaxa = true;
            }
            else {
                nTaxaList.add(-1);
            }
        }
        br.close();

        if (chrList.isEmpty()) {
            throw new IllegalArgumentException("No chromosome lengths loaded from " + chrLengthFileS);
        }

        int n = chrList.size();
        chromosomes = new short[n];
        chrLengths = new int[n];
        nTaxaPerChr = new int[n];
        refToChrIndex.clear();
        for (int i = 0; i < n; i++) {
            chromosomes[i] = chrList.get(i);
            chrLengths[i] = lenList.get(i);
            nTaxaPerChr[i] = nTaxaList.get(i);
            refToChrIndex.put(String.valueOf(chromosomes[i]), i);
        }
        if (!hasExplicitNTaxa) {
            System.out.println("Chromosome length file has no nTaxa column; "
                    + "will use taxa map size as denominator for all chromosomes.");
        }
    }

    /**
     * Decide the order in which taxa (and thus samtools processes) are submitted to the pool.
     *
     * <p>The default {@code interleave} round-robins taxa by their BAM parent directory. Because
     * same-prefix BAMs tend to live on the same physical disk (linear LVM), lexical sorting packs
     * all concurrently running samtools onto one spindle, leaving the other disks idle and starving
     * each process at the pipe. Round-robin spreads the first {@code threadNum} processes across
     * directories/disks so aggregate read bandwidth scales with the number of spindles.</p>
     *
     * <p>Result correctness and reproducibility are unaffected: depth values are integers and the
     * accumulation in {@code double} is exact, so any submission order yields identical sum/sumSq.</p>
     */
    private void orderTaxaForIO() {
        Arrays.sort(taxa); // deterministic baseline before reordering
        if ("sorted".equalsIgnoreCase(taxaOrder)) {
            System.out.println("Taxa order: sorted");
            return;
        }
        if ("shuffle".equalsIgnoreCase(taxaOrder)) {
            Random rng = new Random(42);
            for (int i = taxa.length - 1; i > 0; i--) {
                int j = rng.nextInt(i + 1);
                String t = taxa[i];
                taxa[i] = taxa[j];
                taxa[j] = t;
            }
            System.out.println("Taxa order: shuffle (seed=42)");
            return;
        }
        // interleave (default): round-robin by BAM parent directory
        LinkedHashMap<String, List<String>> byDir = new LinkedHashMap<>();
        for (String t : taxa) {
            String[] bams = taxaBamPathsMap.get(t);
            String dir = (bams != null && bams.length > 0) ? parentDirOf(bams[0]) : "";
            byDir.computeIfAbsent(dir, k -> new ArrayList<>()).add(t);
        }
        List<List<String>> groups = new ArrayList<>(byDir.values());
        String[] result = new String[taxa.length];
        int idx = 0;
        int round = 0;
        boolean placed = true;
        while (placed) {
            placed = false;
            for (List<String> g : groups) {
                if (round < g.size()) {
                    result[idx++] = g.get(round);
                    placed = true;
                }
            }
            round++;
        }
        taxa = result;
        System.out.println("Taxa order: interleave by BAM directory (" + groups.size()
                + " directories, spreads concurrent samtools across disks)");
    }

    private static String parentDirOf(String path) {
        String parent = new File(path).getParent();
        return parent == null ? "" : parent;
    }

    private void finalizeNTaxaPerChr() {
        int defaultNTaxa = taxa.length;
        for (int i = 0; i < nTaxaPerChr.length; i++) {
            if (nTaxaPerChr[i] < 0) {
                nTaxaPerChr[i] = defaultNTaxa;
            }
            if (nTaxaPerChr[i] < 2) {
                throw new IllegalArgumentException("nTaxa for chromosome " + chromosomes[i]
                        + " must be >= 2 for sample SD, got " + nTaxaPerChr[i]);
            }
        }
    }

    private static boolean isInteger(String s) {
        if (s == null || s.isEmpty()) {
            return false;
        }
        for (int i = 0; i < s.length(); i++) {
            if (!Character.isDigit(s.charAt(i))) {
                return false;
            }
        }
        return true;
    }

    @Override
    public void creatAppOptions() {
        options.addOption("app", true, "App name.");
        options.addOption("a", true, "The taxaBamMap file: Taxa\\tCoverage\\tBam1[\\tBam2...] per line " +
                "(column 2 is the genome-wide mean depth of the taxon, used for relative depth). " +
                "Bam files must have a .bai index in the same folder.");
        options.addOption("b", true, "Chromosome length file: Chr\\tLength\\tnTaxa per line. " +
                "nTaxa is the taxon count for mean/SD on that chromosome (PopDepFull tb row count). " +
                "If nTaxa omitted, uses taxa map size for all chromosomes.");
        options.addOption("d", true, "The path of samtools.");
        options.addOption("e", true, "Number of concurrent taxa (samtools depth processes).");
        options.addOption("f", true, "Output directory; writes {chr}.popdep.txt.gz for each chromosome.");
        options.addOption("g", true, "Minimum samtools mapping quality (MAPQ, depth -q). Default 0.");
        options.addOption("o", true, "Taxa submission order: interleave (default, round-robin by BAM " +
                "directory to spread samtools across disks), shuffle (fixed seed), or sorted.");
        options.addOption("k", true, "Checkpoint directory (enables resume). Writes checkpoint.bin " +
                "after each batch; restart with the same dir to skip completed taxa. Prefer a disk " +
                "separate from the BAMs. If omitted, runs a single pass with no checkpointing.");
        options.addOption("ci", true, "Taxa per checkpoint batch (default 500). Larger = less " +
                "checkpoint I/O overhead but more work lost on crash. Only used with -k.");
    }

    @Override
    public void retrieveAppParameters(String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            this.taxaBamFileS = line.getOptionValue("a");
            this.chrLengthFileS = line.getOptionValue("b");
            this.samPath = line.getOptionValue("d");
            this.threadNum = Integer.parseInt(line.getOptionValue("e"));
            this.outDirS = line.getOptionValue("f");
            if (line.hasOption("g")) {
                this.minMapq = Integer.parseInt(line.getOptionValue("g"));
                if (this.minMapq < 0) {
                    throw new IllegalArgumentException("-g must be >= 0");
                }
            }
            if (line.hasOption("o")) {
                this.taxaOrder = line.getOptionValue("o");
            }
            if (line.hasOption("k")) {
                this.checkpointDir = line.getOptionValue("k");
            }
            if (line.hasOption("ci")) {
                this.checkpointIntervalTaxa = Integer.parseInt(line.getOptionValue("ci"));
                if (this.checkpointIntervalTaxa < 1) {
                    throw new IllegalArgumentException("-ci must be >= 1");
                }
            }

            loadChromosomeLengths();

            BufferedReader br = IOUtils.getTextReader(this.taxaBamFileS);
            br.readLine();
            List<String> l;
            String temp;
            while ((temp = br.readLine()) != null) {
                if (temp.trim().isEmpty()) continue;
                l = PStringUtils.fastSplit(temp);
                String[] bams = new String[l.size() - 2];
                for (int i = 0; i < bams.length; i++) {
                    bams[i] = l.get(i + 2);
                }
                this.taxaBamPathsMap.put(l.get(0), bams);
                this.taxaCoverageMap.put(l.get(0), Double.parseDouble(l.get(1)));
            }
            Set<String> tSet = this.taxaBamPathsMap.keySet();
            this.taxa = tSet.toArray(new String[tSet.size()]);
            br.close();

            orderTaxaForIO();
            this.taxaCoverage = new double[taxa.length];
            for (int i = 0; i < taxa.length; i++) {
                this.taxaCoverage[i] = this.taxaCoverageMap.get(taxa[i]);
            }
            finalizeNTaxaPerChr();

            long totalBp = 0;
            for (int len : chrLengths) {
                totalBp += len;
            }
            System.out.println("PopDepCrossChr: " + chromosomes.length + " chromosomes, "
                    + totalBp + " bp total, " + taxa.length + " taxa in map, threadNum=" + threadNum
                    + ", minMapq=" + minMapq
                    + (checkpointDir == null ? ", checkpoint=off"
                        : ", checkpoint=" + checkpointDir + " (every " + checkpointIntervalTaxa + " taxa)"));
            for (int i = 0; i < chromosomes.length; i++) {
                System.out.println("  chr " + chromosomes[i] + ": length=" + chrLengths[i]
                        + ", nTaxa=" + nTaxaPerChr[i]);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.out.println("\nThere are input errors in the command line. Program stops.");
            this.printInstructionAndUsage();
            System.exit(0);
        }
    }

    @Override
    public void printInstructionAndUsage() {
        System.out.println(PGLAPPEntrance.getTIGERIntroduction());
        System.out.println("Below are the commands of PopDepCrossChr.");
        this.printUsage();
    }
}
