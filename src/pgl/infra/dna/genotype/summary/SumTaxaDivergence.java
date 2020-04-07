package pgl.infra.dna.genotype.summary;

import pgl.infra.dna.genotype.GenotypeTable;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.List;

public class SumTaxaDivergence {
    float[][] dxyMatrix = null;
    String[] taxa;

    public SumTaxaDivergence (GenotypeTable gt) {
        dxyMatrix = gt.getDxyMatrix();
        taxa = gt.getTaxaNames();
    }

    private String getSumHeader () {
        StringBuilder sb = new StringBuilder("SumTaxaDivergence");
        for (int i = 0; i < taxa.length; i++) {
            sb.append("\t").append(taxa[i]);
        }
        return sb.toString();
    }

    public String[] getSumContent () {
        String[] content = new String[taxa.length];
        List<Integer> l = PArrayUtils.getIndexList(taxa.length);
        l.parallelStream().forEach(i ->{
            StringBuilder sb = new StringBuilder();
            sb.append(taxa[i]);
            for (int j = 0; j < taxa.length; j++) {
                sb.append("\t").append((float)dxyMatrix[i][j]);
            }
            content[i] = sb.toString();
        });
        return content;
    }

    public void writeSummary (String outfileS, IOFileFormat format) {
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
            bw.write(this.getSumHeader());
            bw.newLine();
            String[] content = this.getSumContent();
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
