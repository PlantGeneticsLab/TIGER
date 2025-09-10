package pgl.app.fastCall3;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import pgl.AppAbstract;
import pgl.PGLAPPEntrance;
import pgl.infra.table.RowTable;

import java.util.Arrays;

/**
 * A utility class for creating customized variation libraries from an existing library
 * based on user-specified genomic positions.
 * <p>
 * This class extends {@link AppAbstract} to provide functionality for filtering
 * and extracting specific variations from a pre-built variation library. It allows
 * users to create smaller, focused variation libraries containing only positions
 * of interest, which can be useful for targeted analysis or reducing data size.
 *
 * <p>Key features include:
 * <ul>
 *   <li>Reading from a binary variation library file</li>
 *   <li>Filtering variations based on user-provided positions</li>
 *   <li>Generating a new, smaller variation library file</li>
 *   <li>Command-line interface for easy integration into pipelines</li>
 * </ul>
 *
 * @author Fei Lu
 * @version 3.0
 * @since 1.0
 */
class CustomizeVariationLibraryF3 extends AppAbstract {
    // the genetic variation file in binary format
    String inputLibFileS = null;
    String positionFileS = null;
    // the genetic variation file in text format
    String outputLibFileS = null;

    public CustomizeVariationLibraryF3(String[] args) {
        this.creatAppOptions();
        this.retrieveAppParameters(args);
        this.customizeLibrary();
    }

    @Override
    public void creatAppOptions() {
        options.addOption("app", true, "App name.");
        options.addOption("mod", true, "Module name of FastCall 3.");
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

    /**
     * A method to customize a genetic variation library based on user-provided positions.
     * <p>
     * This method reads a pre-built variation library from a binary file, reads custom
     * positions from a text file, filters the library and generates a new variation
     * library file with only the custom positions.
     *
     * @param inputLibFileS
     *            the input genetic variation library file in binary format
     * @param positionFileS
     *            the user provided file with custom positions
     * @param outputLibFileS
     *            the custom genetic variation library file
     */
    public void customizeLibrary() {
        VariationLibraryF3 vl = new VariationLibraryF3(this.inputLibFileS);
        RowTable<String> table = new RowTable<String>(this.positionFileS);
        int[] customPositions = table.getColumnAsIntArray(0);
        Arrays.sort(customPositions);
        vl.writeBinaryFileS(this.outputLibFileS, customPositions);
    }

    @Override
    public void printInstructionAndUsage() {
        System.out.println(PGLAPPEntrance.getTIGERIntroduction());
        System.out.println("Below are the commands of module \"clib\" in FastCall 3.");
        this.printUsage();
    }
}
