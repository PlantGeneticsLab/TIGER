package pgl.app.fastCall2;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import pgl.AppAbstract;
import pgl.PGLAPPEntrance;
import pgl.PGLConstraints;
import pgl.infra.dna.FastaBit;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.Dyad;
import pgl.infra.utils.IOUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.*;
import java.util.concurrent.atomic.LongAdder;

class ViewVariationLibrary extends AppAbstract {
    // the genetic variation file in binary format
    String binaryLibFileS = null;
    // the genetic variation file in text format
    String textLibFileS = null;

    public ViewVariationLibrary(String[] args) {
        this.creatAppOptions();
        this.retrieveAppParameters(args);
        this.convertLibrary();
    }

    @Override
    public void creatAppOptions() {
        options.addOption("app", true, "App name.");
        options.addOption("mod", true, "Module name of FastCall 2.");
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

    public void convertLibrary () {
        VariationLibrary vl = new VariationLibrary(this.binaryLibFileS);
        vl.writeTextFileS(this.textLibFileS);
    }

    @Override
    public void printInstructionAndUsage() {
        System.out.println(PGLAPPEntrance.getTIGERIntroduction());
        System.out.println("Below are the commands of module \"vlib\" in FastCall 2.");
        this.printUsage();
    }

}
