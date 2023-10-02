package pgl.app.fastCall2;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import pgl.AppAbstract;
import pgl.PGLAPPEntrance;
import pgl.infra.table.RowTable;

import java.util.Arrays;

class CustomizeVariationLibrary extends AppAbstract {
    // the genetic variation file in binary format
    String inputLibFileS = null;
    String positionFileS = null;
    // the genetic variation file in text format
    String outputLibFileS = null;

    public CustomizeVariationLibrary(String[] args) {
        this.creatAppOptions();
        this.retrieveAppParameters(args);
        this.customizeLibrary();
    }

    @Override
    public void creatAppOptions() {
        options.addOption("app", true, "App name.");
        options.addOption("mod", true, "Module name of FastCall 2.");
        options.addOption("a", true, "The input genetic variation library file in binary format.");
        options.addOption("b", true, "The user provided file with custom positions.");
        options.addOption("c", true, "The custom genetic variation library file.");
    }

    @Override
    public void retrieveAppParameters(String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            this.inputLibFileS = line.getOptionValue("a");
            this.positionFileS = line.getOptionValue("b");
            this.outputLibFileS = line.getOptionValue("c");
        }
        catch(Exception e) {
            e.printStackTrace();
            System.out.println("\nThere are input errors in the command line. Program stops.");
            this.printInstructionAndUsage();
            System.exit(0);
        }
    }

    public void customizeLibrary() {
        VariationLibrary vl = new VariationLibrary(this.inputLibFileS);
        RowTable<String> table = new RowTable<String>(this.positionFileS);
        int[] customPositions = table.getColumnAsIntArray(0);
        Arrays.sort(customPositions);
        vl.writeBinaryFileS(this.outputLibFileS, customPositions);
    }

    @Override
    public void printInstructionAndUsage() {
        System.out.println(PGLAPPEntrance.getTIGERIntroduction());
        System.out.println("Below are the commands of module \"clib\" in FastCall 2.");
        this.printUsage();
    }
}
