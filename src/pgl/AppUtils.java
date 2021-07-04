package pgl;

import pgl.infra.utils.Dyad;
import pgl.infra.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

/**
 * Utilities of TIGER apps
 * @author feilu
 */
public class AppUtils {
    /**
     * Return the list of steps and the list of parameters
     * @param parameterFileS
     * @return
     */
    public static Dyad<List<String>, List<String>> getParameterList (String parameterFileS) {
        List<String> pLineList = new ArrayList<>();
        List<String> sLineList = new ArrayList<>();
        Dyad<List<String>, List<String>> d = new Dyad<>(pLineList, sLineList);
        try {
            BufferedReader br = IOUtils.getTextReader(parameterFileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) continue;
                if (temp.startsWith("@")) {
                    sLineList.add(temp);
                    continue;
                }
                if (temp.isEmpty()) continue;
                pLineList.add(temp);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        return d;
    }

    public static void creatParameterFile (String sampleParameterFileS, List<String> sLineList, List<String> pLineList, String outfileS) {
        try {
            int sCnt = 0;
            int pCnt = 0;
            BufferedReader br = IOUtils.getTextReader(sampleParameterFileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    bw.write(temp);
                }
                else if (temp.startsWith("@")) {
                    bw.write(sLineList.get(sCnt));
                    sCnt++;
                }
                else if (temp.length() == 0) {

                }
                else if (temp.startsWith("\\s+")) {
                    bw.write(temp);
                }
                else {
                    bw.write(pLineList.get(pCnt));
                    pCnt++;
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
