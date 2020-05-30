package pgl.infra.germ;

import pgl.infra.utils.IOUtils;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class GermBank {
    String header = "GID\tAccession\tChineseName\tEnglishName\tTaxaName\tContamination\tSeedStorage";
    List<GermRecord> recordList = null;

    public GermBank (String infileS) {
        this.readIn(infileS);
    }

    private void readIn (String infileS) {
        recordList = new ArrayList<>();
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = br.readLine();
            while ((temp = br.readLine()) != null) {
                GermRecord gr = new GermRecord(temp);
                recordList.add(gr);
            }
            br.close();
            System.out.println(this.getSummary());
            this.sortByGID();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    public String getSummary () {
        StringBuilder sb = new StringBuilder();
        sb.append("There are ").append(this.getGermplasmNumber()).append(" germplasm records in the current germplasm bank\n");

        return sb.toString();
    }

    public int getGermplasmNumber () {
        return this.recordList.size();
    }

    public void sortByGID () {
        Collections.sort(this.recordList);
    }
}
