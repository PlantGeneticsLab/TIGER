package pgl.infra.germ;

public class GermRecord implements Comparable<GermRecord>{
    public String gid = null;
    public String accession = null;
    public String chineseName = null;
    public String englishName = null;
    public String taxonName = null;
    byte ifContamination = Byte.MIN_VALUE;
    byte ifSeedAvailable = Byte.MIN_VALUE;

    public GermRecord (String inputline) {
        String[] temp = inputline.split("\t");
        gid = temp[0];
        accession = temp[1];
        chineseName = temp[2];
        englishName = temp[3];
        taxonName = temp[4];
        ifContamination = Byte.parseByte(temp[5]);
        ifSeedAvailable = Byte.parseByte(temp[6]);
    }


    @Override
    public int compareTo(GermRecord o) {
        return this.gid.compareTo(o.gid);
    }
}
