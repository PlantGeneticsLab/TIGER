package pgl.infra.dna.genot;

import pgl.app.grt.GRTVCFUtils;
import pgl.infra.utils.IOUtils;

import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.nio.ByteBuffer;

/**
 * Data output utilities for genotype object implemented {@link GenotypeTable}
 *
 * @author feilu
 */
public class GenotypeExport {

    public static void output(GenotypeTable gt, String outfileS, GenoIOFormat format) {
        if (format == GenoIOFormat.VCF) {
            toVCF(gt, outfileS,false);
        }
        else if (format == GenoIOFormat.VCF_GZ) {
            toVCF(gt, outfileS,true);
        }
        else if (format == GenoIOFormat.Binary) {
            toBinary(gt, outfileS, false);
        }
        else if (format == GenoIOFormat.Binary_GZ) {
            toBinary(gt, outfileS, true);
        }
        else if (format == GenoIOFormat.HDF5) {

        }
        System.out.println("Genotype table exported to "+ outfileS);
    }

    private static void toBinary (GenotypeTable gt, String outfileS, boolean ifGZ) {
        try {
            DataOutputStream dos = null;
            if (ifGZ) {
                dos = IOUtils.getBinaryGzipWriter(outfileS);
            }
            else  {
                dos = IOUtils.getBinaryWriter(outfileS);
            }
            int siteNumber = gt.getSiteNumber();
            int taxaNumber = gt.getTaxaNumber();
            dos.writeInt(siteNumber);
            dos.writeInt(taxaNumber);
            for (int i = 0; i < taxaNumber; i++) {
                dos.writeUTF(gt.getTaxonName(i));
            }
            int byteSize = getByteSizeOfSiteInBinary(gt.getTaxaNumber());
            ByteBuffer bb = ByteBuffer.allocate(byteSize);
            for (int i = 0; i < siteNumber; i++) {
                bb = gt.getBinaryOutput(i, bb);
                bb.flip();
                dos.write(bb.array());
                bb.clear();
            }
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Return the number of byte for a {@link pgl.infra.dna.snp.BiSNP}
     * @return
     */
    public static int getByteSizeOfSNPInBinary () {
        //short chr, int pos, genotype (ref+alt), byte feature *2
        return 9;
    }

    /**
     * Return the number of byte for genotype of a site (chr, pos, ref+alt, alternative at phase 1 and phase 2, and missing)
     * @param taxaNumber
     * @return
     */
    public static int getByteSizeOfSiteInBinary (int taxaNumber) {
        //short chr, int pos, genotype (ref+alt), byte feature *2, BitSet phase1, BitSet phase2, BitSet missing
        int n = getByteSizeOfSNPInBinary();
        if (taxaNumber%64 == 0) {
            n += (taxaNumber/64)*8*3;
        }
        else  {
            n += (taxaNumber/64+1)*8*3;
        }
        return n;
    }

    private static void toVCF (GenotypeTable gt, String outfileS, boolean ifGZ) {
        BufferedWriter bw = null;
        if (ifGZ) bw = IOUtils.getTextGzipWriter(outfileS);
        else bw = IOUtils.getTextWriter(outfileS);
        try {
            bw.write(VCFUtils.getVCFAnnotation());
            bw.write(VCFUtils.getVCFHeader(gt.getTaxaNames()));
            bw.newLine();
            for (int i = 0; i < gt.getSiteNumber(); i++) {
                bw.write(gt.getUnphasedVCFRecord(i));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
