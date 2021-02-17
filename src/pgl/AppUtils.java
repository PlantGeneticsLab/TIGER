package pgl;

import pgl.infra.utils.Dyad;
import pgl.infra.utils.IOUtils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

/**
 * Utilities of TIGER apps
 * @author feilu
 */
public class AppUtils {

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
}
