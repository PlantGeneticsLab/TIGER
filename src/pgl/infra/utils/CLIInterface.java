package pgl.infra.utils;

/**
 * This interface is a template of building command line software from PlantGenetics project
 * @author feilu
 */
public interface CLIInterface {
    /**
     * Build options for a CLI
     * */
    void createOptions ();

    /**
     * Retrieve parameters from command line
     * @param args
     */
    void retrieveParameters (String[] args);

    /**
     * Output software introduction and usage
     */
    void printIntroductionAndUsage ();

    /**
     * Create introduction of the software
     * @return 
     */
    String createIntroduction ();
}
