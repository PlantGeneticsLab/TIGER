package pgl.infra.dna;

/**
 * Fasta record with 3 bits storage.
 * @author feilu
 */
public class FastaRecordBit extends Sequence3Bit implements FastaRecordInterface {
    /**
     * The description of Fasta record
     */
    String description;
    /**
     * The ID of Fasta record
     */
    int id;

    /**
     * Constructs an object with DNA sequence {@link String}.
     * @param description
     * @param seq
     * @param id
     */
    public FastaRecordBit(String description, String seq, int id) {
        super(seq);
        this.description = description;
        this.id = id;
    }

    /**
     * Constructs an object with {@link Sequence3Bit}.
     * @param description
     * @param s3
     * @param id
     */
    public FastaRecordBit (String description, Sequence3Bit s3, int id) {
        super.sequenceLength = s3.sequenceLength;
        super.seqS = s3.seqS;
        this.description = description;
        this.id = id;
    }
    
    @Override
    public String getDescription() {
        return description;
    }

    @Override
    public int getID() {
        return id;
    }

    @Override
    public void setDescription(String description) {
        this.description = description;
    }

    @Override
    public void setID(int id) {
        this.id = id;
    }
}
