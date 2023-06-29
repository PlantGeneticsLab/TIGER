package pgl;

import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

public abstract class AppAbstract implements AppInterface{
    protected Options options = new Options();
    protected HelpFormatter optionFormat = new HelpFormatter();

    public void printUsage() {
        optionFormat.printHelp(" ", options );
    }
}
