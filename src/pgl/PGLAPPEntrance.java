/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl;

import java.io.File;
import java.util.Arrays;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import pgl.app.fastCall.FastCall;
import pgl.app.popdep.PopDep;
import pgl.app.speedCall.SpeedCall;
import pgl.infra.utils.CLIInterface;

/**
 *
 * @author feilu
 */
public class PGLAPPEntrance implements CLIInterface {
    Options options = new Options();
    HelpFormatter optionFormat = new HelpFormatter();
    String introduction = this.createIntroduction();
    String app = null;
    String parameterPath = null;

    
    public PGLAPPEntrance (String[] args) {
        this.createOptions();
        this.retrieveParameters (args);
    }
    
    @Override
    public void createOptions() {
        options = new Options();
        options.addOption("a", true, "App. e.g. -a FastCall");
        options.addOption("p", true, "Parameter file path of an app. e.g. parameter_fastcall.txt");
    }

    @Override
    public void retrieveParameters(String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            app = line.getOptionValue("a");
            parameterPath = line.getOptionValue("p");
        }
        catch(Exception e) {
            e.printStackTrace();
            System.exit(0);
        }
        if (app == null) {
            System.out.println("App does not exist");
            this.printIntroductionAndUsage();
            System.exit(0);
        }
        if (app.equals(AppNames.FastCall.getName())) {
            new FastCall (this.parameterPath);
        }
        else if (app.equals(AppNames.PopDep.getName())) {
            new PopDep(this.parameterPath);
        }
        else if (app.equals(AppNames.HapScanner.getName())) {

        }
        else {
            System.out.println("App does not exist");
            this.printIntroductionAndUsage();
            System.exit(0);
        }
        if (this.parameterPath == null) {
            System.out.println("Parametar file does not exist");
            this.printIntroductionAndUsage();
            System.exit(0);
        }
        File f = new File (this.parameterPath);
        if (!f.exists()) {
            System.out.println("Parametar file does not exist");
            this.printIntroductionAndUsage();
            System.exit(0);
        }
    }

    @Override
    public void printIntroductionAndUsage() {
        System.out.println("Incorrect options input. Program stops.");
        System.out.println(introduction);
        optionFormat.printHelp("TIGER.jar", options );
    }

    @Override
    public String createIntroduction() {
        StringBuilder sb = new StringBuilder();
        sb.append("\nToolkits Integrated for Genetic and Evolutionary Research (TIGER) is designed to simplify its usage.\n");
        sb.append("It uses two options to run its apps. \"-a\" is used to select an app. \"-p\" is used to set parameters of an app.\n");
        sb.append("e.g. The command line usage of the app FastCall is: ");
        sb.append("java -Xmx100g -jar TIGER.jar -a FastCall -p parameter_fastcall.txt > log.txt &\n");
        sb.append("\nAvailable apps in TIGER include,\n");
        for (int i = 0; i < AppNames.values().length; i++) {
            sb.append(AppNames.values()[i].getName()).append("\n");
        }
        sb.append("\nPlease visit https://github.com/PlantGeneticsLab/TIGER for details.\n");
        return sb.toString();
    }
    
    public static void main (String[] args) {
        new PGLAPPEntrance(args);
    }

}
