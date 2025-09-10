package pgl.app.fastCall3;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import pgl.AppAbstract;
import pgl.PGLAPPEntrance;

/**
 * A utility class for viewing and converting genetic variation library files between binary and text formats.
 * 
 * <p>This class extends {@link AppAbstract} and provides functionality to:
 * <ul>
 *   <li>Convert binary variation library files to human-readable text format</li>
 *   <li>Handle command-line arguments for input/output file specification</li>
 *   <li>Integrate with the FastCall3 pipeline for genetic variation analysis</li>
 * </ul>
 *
 * <p>Usage example:
 * <pre>
 * ViewVariationLibraryF3 -i input.bin -o output.txt
 * </pre>
 *
 * <p>The text output format is structured as follows:
 * <ul>
 *   <li>Each line represents a single genetic variant</li>
 *   <li>Fields are tab-separated</li>
 *   <li>Format: Chromosome\tPosition\tReference_Allele\tAlternative_Alleles</li>
 * </ul>
 *
 * @author Fei Lu
 * @version 3.0
 * @since 1.0
 * @see AppAbstract
 * @see VariationLibraryF3
 * @see FastCall3
 */
class ViewVariationLibraryF3 extends AppAbstract {
    // the genetic variation file in binary format
    String binaryLibFileS = null;
    // the genetic variation file in text format
    String textLibFileS = null;

    public ViewVariationLibraryF3(String[] args) {
        this.creatAppOptions();
        this.retrieveAppParameters(args);
        this.convertLibrary();
    }

    @Override
    public void creatAppOptions() {
        options.addOption("app", true, "App name.");
        options.addOption("mod", true, "Module name of FastCall 3.");
        options.addOption("a", true, "The input genetic variation library file in binary format.");
        options.addOption("b", true, "The output genetic variation library file in text format");
    }

    @Override
    public void retrieveAppParameters(String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            this.binaryLibFileS = line.getOptionValue("a");
            this.textLibFileS = line.getOptionValue("b");
        }
        catch(Exception e) {
            e.printStackTrace();
            System.out.println("\nThere are input errors in the command line. Program stops.");
            this.printInstructionAndUsage();
            System.exit(0);
        }
    }

    /**
     * Converts the genetic variation library file from binary format to text format.
     * <p>
     * The method reads the genetic variation library file in binary format and
     * writes it to a text format file.
     */
    public void convertLibrary () {
        VariationLibraryF3 vl = new VariationLibraryF3(this.binaryLibFileS);
        vl.writeTextFileS(this.textLibFileS);
    }

    @Override
    public void printInstructionAndUsage() {
        System.out.println(PGLAPPEntrance.getTIGERIntroduction());
        System.out.println("Below are the commands of module \"vlib\" in FastCall 3.");
        this.printUsage();
    }

}
