package pgl.infra.dna.genot.summa;

import pgl.infra.dna.genot.GenotypeRows;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.GenotypeTable;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;

import java.io.BufferedWriter;
import java.util.List;

public class SumTaxaDivergence {
    float[][] dxyMatrix = null;
    String[] taxa;

    public SumTaxaDivergence (GenotypeTable gt) {
        long start = System.nanoTime();
        if (gt instanceof GenotypeGrid) {
            dxyMatrix = gt.getIBSDistanceMatrix();
        }
        else if (gt instanceof GenotypeRows) {
            dxyMatrix = ((GenotypeRows) gt).getDxyMatrixFast10K();
        }
        else {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
        StringBuilder sb = new StringBuilder("IBS distance matrix calculation takes ");
        sb.append(Benchmark.getTimeSpanSeconds(start)).append(" seconds");
        System.out.println(sb.toString());
        taxa = gt.getTaxaNames();
    }

    private String getDxyMatrixHeader() {
        StringBuilder sb = new StringBuilder("Dxy");
        for (int i = 0; i < taxa.length; i++) {
            sb.append("\t").append(taxa[i]);
        }
        return sb.toString();
    }

    private String[] getDxyMatrixLine() {
        String[] content = new String[taxa.length];
        List<Integer> l = PArrayUtils.getIndexList(taxa.length);
        l.parallelStream().forEach(i ->{
            StringBuilder sb = new StringBuilder();
            sb.append(taxa[i]);
            for (int j = 0; j < taxa.length; j++) {
                sb.append("\t").append(String.format("%.3f", dxyMatrix[i][j]));
            }
            content[i] = sb.toString();
        });
        return content;
    }

    public void writeDxyMatrix(String outfileS, IOFileFormat format) {
        try {
            BufferedWriter bw = null;
            if (format == IOFileFormat.Text) {
                bw = IOUtils.getTextWriter(outfileS);
            }
            else if (format == IOFileFormat.TextGzip) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }
            else {
                throw new UnsupportedOperationException("Not supported yet.");
            }
            bw.write(this.getDxyMatrixHeader());
            bw.newLine();
            String[] content = this.getDxyMatrixLine();
            for (int i = 0; i < content.length; i++) {
                bw.write(content[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(0);
        }
    }
}
