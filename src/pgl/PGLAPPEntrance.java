/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl;

import org.apache.commons.cli.Options;
import pgl.app.fastCall.FastCall;
import pgl.app.fastCall2.FastCall2;
import pgl.app.fastCall3.FastCall3;
import pgl.app.hapScanner.HapScanner;
import pgl.app.popdep.PopDep;

/**
 * This provides the interface between users and TIGER apps
 * @author feilu
 */
public class PGLAPPEntrance {
    Options options = new Options();
    static String introduction = createIntroduction();
    String app = null;
    String parameterPath = null;

    
    public PGLAPPEntrance (String[] args) {
        this.selectApp(args);
        this.runApp(args);
    }

    /**
     * This method is used to select the app according to the command line
     * parameter -app.
     * @param args the command line parameters
     */
    private void selectApp (String[] args) {
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-app")) {
                app = args[i+1];
                break;
            }
        }
        if (app == null) {
            System.out.println("App does not exist. The program stops.");
            System.out.println(introduction);
            System.exit(0);
        }

    }

    /**
     * This method runs the specified app with given parameters
     * @param args the parameters for the app
     */
    private void runApp(String[] args) {
        options = new Options();
        if (app.equals(AppNames.FastCall.getName())) {
            new FastCall (args);
        }
        else if (app.equals(AppNames.FastCall2.getName())) {
            new FastCall2(args);
        }
        else if (app.equals(AppNames.FastCall3.getName())) {
            new FastCall3(args);
        }
        else if (app.equals(AppNames.PopDep.getName())) {
            new PopDep(this.parameterPath);
        }
        else if (app.equals(AppNames.HapScanner.getName())) {
            new HapScanner(this.parameterPath);
        }
        else {
            System.out.println("App does not exist. Programs stops.");
            System.out.println(introduction);
            System.exit(0);
        }
    }

    /**
     * Create an introduction string for the main method.
     * @return an introduction string
     */
    private static String createIntroduction() {
        StringBuilder sb = new StringBuilder();
        sb.append("\nToolkits Integrated for Genetic and Evolutionary Research (TIGER) is designed to simplify its usage.\n");
        sb.append("It uses two sets of options to run its apps. \"-app\" is used to select an app in TIGER. The other options are used to set parameters of a specific app.\n");
        sb.append("e.g. The command line usage of the app FastCall is: ");
        sb.append("java -Xmx100g -jar TIGER.jar -app FastCall -a chr001.fa -b taxaBamMap.txt -c 30 -d 20 -e 2 -f 0.2 -g 3 -h 0.8 -i 0.4 -j 0.2 -k 1 -l 32 -m /ing -n /usr/local/bin/samtools > log.txt &\n");
        sb.append("\nAvailable apps in TIGER include,\n");
        for (int i = 0; i < AppNames.values().length; i++) {
            sb.append(AppNames.values()[i].getName()).append("\n");
        }
        sb.append("\nPlease visit https://github.com/PlantGeneticsLab/TIGER for details.\n");
        return sb.toString();
    }

    /**
     * Return the introduction string for the main method.
     * @return the introduction string
     */
    public static String getTIGERIntroduction() {
        return introduction;
    }

    /**
     * This is the main method of TIGER. It creates an instance of {@link PGLAPPEntrance} to run TIGER.
     * @param args the command line parameters
     */
    public static void main (String[] args) {
        new PGLAPPEntrance(args);
    }
}
