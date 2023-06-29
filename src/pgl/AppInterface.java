package pgl;

import org.apache.commons.cli.Options;

public interface AppInterface {

    public void creatAppOptions ();

    public void retrieveAppParameters (String[] args);

    public void printUsage ();

    public void printInstructionAndUsage ();
}
